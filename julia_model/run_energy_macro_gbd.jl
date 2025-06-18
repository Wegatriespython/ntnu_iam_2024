# run_energy_macro_gbd.jl - Generalised Benders reformulation
# Energy systems x macro model (MESSAGE + MACRO)
# Conventions used below:
#   • All monetary variables in the **macro** block are in trillion US$.
#   • The energy sub-problem (MESSAGE) returns TOTAL_COST in **billion US$**.
#     We convert to trillion by dividing by 1000 exactly once when
#     constructing every Benders cut.
#   • No other ad-hoc scaling factors appear.

using JuMP, Ipopt, LinearAlgebra
println("Starting Energy-Macro GBD decomposition...")

# include shared definitions
include("shared.jl");
println("✓ loaded shared definitions");
include("energy_model_world.jl");
println("✓ loaded energy model");
include("macro_data_load.jl");
println("✓ loaded macro data");
include("macro_core.jl");
println("✓ loaded macro core");
include("macro_presolve.jl");
println("✓ loaded macro preprocessing");

# 1 Energy sub-problem (fixed service demands S_bar)
function create_energy_subproblem(S_bar::Dict{Tuple{String,Int},Float64})
  model = Model(Ipopt.Optimizer)
  set_optimizer_attribute(model, "print_level", 0)
  set_optimizer_attribute(model, "tol", 1e-6)
  set_optimizer_attribute(model, "max_iter", 1000)

  # decision variables
  @variable(model, ACT[technology, year_all] >= 0)
  @variable(model, CAP_NEW[technology, year_all] >= 0)
  @variable(model, EMISS[year_all])
  @variable(model, CUM_EMISS)
  @variable(model, TOTAL_COST)                # **billion US$**
  @variable(model, COST_ANNUAL[year_all])
  @variable(model, demand_slack[sector, year_all] >= 0)

  # bookkeeping for sector-split mapping
  n_map = Dict(s => sum(haskey(map_energy_sector, (e, l, s)) for e in energy, l in level)
               for s in sector)

  # energy balance
  bal = Dict{Tuple{String,String,Int},ConstraintRef}()
  for e in energy, l in level, y in year_all
    sector_key = findfirst(s -> haskey(map_energy_sector, (e, l, s)), sector)
    if sector_key === nothing
      continue
    end
    if haskey(S_bar, (sector_key, y))
      D = S_bar[(sector_key, y)] / n_map[sector_key]
      bal[(e, l, y)] = @constraint(model,
        sum(ACT[t, y] * (get(output, (t, e, l), 0) - get(input, (t, e, l), 0)) for t in technology)
        +
        demand_slack[sector_key, y] == D)
    else
      bal[(e, l, y)] = @constraint(model,
        sum(ACT[t, y] * (get(output, (t, e, l), 0) - get(input, (t, e, l), 0)) for t in technology)
        >=
        demand[(e, l)] * (gdp[y] / gdp[2020])^beta)
    end
  end

  # capacity & technical constraints
  for t in technology, y in year_all
    if haskey(hours, t) && haskey(lifetime, t)
      idx = findfirst(==(y), year_all)
      @constraint(model,
        ACT[t, y] <= sum(CAP_NEW[t, year_all[i]] * hours[t]
                        for i in 1:idx if (idx - i + 1) * duration_period <= lifetime[t]))
    end
  end
  for t in technology, y in year_all[2:end]
    if !haskey(diffusion_up, t)
      continue
    end
    y_prev = year_all[findfirst(==(y), year_all)-1]
    @constraint(model, CAP_NEW[t, y] <= CAP_NEW[t, y_prev] * (1 + diffusion_up[t])^duration_period + get(startup, t, 0))
  end

  # share constraints
  for s in share, y in year_all
    lhs = [t for t in technology if haskey(tec_share, (s, t))]
    rhs = [t for t in technology if haskey(tec_share_rhs, (s, t))]
    if haskey(share_up, s)
      @constraint(model, sum(ACT[t, y] for t in lhs) <= share_up[s] * sum(ACT[t, y] for t in rhs))
    end
    if haskey(share_lo, s)
      @constraint(model, sum(ACT[t, y] for t in lhs) >= share_lo[s] * sum(ACT[t, y] for t in rhs))
    end
  end

  # emissions
  @constraint(model, [y in year_all], sum(ACT[t, y] * get(CO2_emission, t, 0) for t in technology) == EMISS[y])
  @constraint(model, sum(EMISS[y] * duration_period for y in year_all) == CUM_EMISS)

  # cost accounting
  year_idx = Dict(y => i for (i, y) in enumerate(year_all))
  lifetime_r = Dict{Tuple{String,Int},Vector{Int}}()
  for t in technology, y in year_all
    if !haskey(lifetime, t)
      continue
    end
    idx, lt = year_idx[y], lifetime[t]
    lifetime_r[(t, y)] = [y2 for y2 in year_all if year_idx[y2] <= idx && (idx - year_idx[y2] + 1) * duration_period <= lt]
  end
  disc = Dict(y => (1 - drate)^(duration_period * (year_idx[y] - 1)) for y in year_all)
  @constraint(model, [y in year_all],
    sum(ACT[t, y] * get(vom, (t, y), 0) for t in technology) +
    sum(sum(CAP_NEW[t, y2] * get(cost_capacity, (t, y2), 0) for y2 in lifetime_r[(t, y)]) for t in keys(lifetime))
    ==
    COST_ANNUAL[y])
  @constraint(model, sum(COST_ANNUAL[y] * duration_period * disc[y] for y in year_all) == TOTAL_COST)

  # historical calibration
  for ((t, y), v) in energy_calibration
    fix(ACT[t, y], v; force=true)
  end
  for ((t, y), v) in energy_calibration_lo
    set_lower_bound(ACT[t, y], v)
  end
  for t in technology
    if haskey(startup, t)
      fix(CAP_NEW[t, 2020], startup[t]; force=true)
    end
  end

  @objective(model, Min, TOTAL_COST + 1e6 * sum(demand_slack))
  return model, bal, demand_slack
end

# 2 Master problem (MACRO + theta)
function create_master_problem()
  m = Model(Ipopt.Optimizer)
  set_optimizer_attribute(m, "print_level", 3)
  set_optimizer_attribute(m, "tol", 1e-5)
  set_optimizer_attribute(m, "max_iter", 2000)

  # shared & macro vars
  @variable(m, S[sector, year_all] >= 0)
  @variable(m, PHYSENE[sector, year_all] >= 0)
  @variable(m, theta >= 0, upper_bound = 50)          # trillion US$

  @variable(m, K[year_all] >= 0)
  @variable(m, KN[year_all] >= 0)
  @variable(m, Y[year_all] >= 0)
  @variable(m, YN[year_all] >= 0)
  @variable(m, PRODENE[sector, year_all] >= 0)
  @variable(m, NEWENE[sector, year_all] >= 0)
  @variable(m, C[year_all] >= 0)
  @variable(m, I[year_all] >= 0)
  @variable(m, UTILITY)
  @variable(m, EC[year_all])

  # utility
  util = sum(udf[y] * log(C[y]) * duration_period for y in year_all if y ∉ (2020, 2080))
  if 2080 in year_all
    util += udf[2080] * log(C[2080]) * (duration_period + 1 / finite_time_corr[2080])
  end
  @constraint(m, UTILITY == util)

  # accounting
  @constraint(m, [y in year_all], Y[y] == C[y] + I[y] + EC[y])
  @constraint(m, [y in year_all; y != 2020], KN[y] == duration_period * I[y])
  @constraint(m, [y in year_all; y != 2020],
    YN[y] == (a * KN[y]^(rho * kpvs) * newlab[y]^(rho * (1 - kpvs)) +
              b * NEWENE["ELEC", y]^(rho * elvs) * NEWENE["NELE", y]^(rho * (1 - elvs)))^(1 / rho))
  @constraint(m, [y in year_all; y != 2020],
    K[y] == K[year_all[findfirst(==(y), year_all)-1]] * (1 - depr)^duration_period + KN[y])
  @constraint(m, [y in year_all; y != 2020],
    Y[y] == Y[year_all[findfirst(==(y), year_all)-1]] * (1 - depr)^duration_period + YN[y])
  @constraint(m, [s in sector, y in year_all; y != 2020],
    PRODENE[s, y] == PRODENE[s, year_all[findfirst(==(y), year_all)-1]] * (1 - depr)^duration_period + NEWENE[s, y])

  # link S <-> PRODENE
  @constraint(m, [s in sector, y in year_all; y != 2020], S[s, y] == PRODENE[s, y] * aeei_factor[(s, y)])
  @constraint(m, [s in sector, y in year_all], PHYSENE[s, y] == S[s, y])

  # energy cost expansion
  disc = Dict(y => (1 - drate)^(duration_period * (findfirst(==(y), year_all) - 1)) for y in year_all)
  totW = sum(cost_MESSAGE[y] for y in year_all if y != 2020)
  @constraint(m, [y in year_all; y != 2020],
    EC[y] == theta * (cost_MESSAGE[y] / totW) / (disc[y] * duration_period))

  # terminal investment
  if 2080 in year_all
    @constraint(m, I[2080] >= K[2080] * (grow[2080] + depr))
  end

  # objective
  @objective(m, Min, -UTILITY + theta)       # both trillion $

  # starts
  set_start_value(theta, 40.0)
  for y in year_all
    set_start_value(C[y], 0.8 * c0)
    set_start_value(I[y], 0.8 * i0)
    set_start_value(Y[y], 0.8 * y0)
  end
  set_start_value(UTILITY, 0.0)
  return m
end

# 3 Benders infrastructure
mutable struct BendersCut
  type::String                       # "optimality" | "feasibility"
  v::Float64                         # cost value (trn $)
  lambda::Dict{Tuple{String,Int},Float64} # duals (trn $ / PWh)
  S_hat::Dict{Tuple{String,Int},Float64} # incumbent S
end

function add_opt_cut!(M::Model, cut::BendersCut)
  @constraint(M, M[:theta] >= cut.v + sum(cut.lambda[key] * (M[:S][key...] - cut.S_hat[key]) for key in keys(cut.lambda)))
end

function set_gbd_master_bounds_and_initial_values!(M::Model)
  # Set bounds and initial values for master problem variables
  # K, Y bounds
  set_lower_bound(M[:K][2020], k0 * 0.8)
  set_upper_bound(M[:K][2020], k0 * 1.2)
  fix(M[:K][2020], k0; force=true)
  
  # Energy service bounds
  for s in sector, y in year_all
    set_lower_bound(M[:S][s, y], 0.01)
    set_upper_bound(M[:S][s, y], 1000.0)
  end
end

# 4 Main GBD loop
function solve_gbd(maxit::Int=20, tol=1e-4)
  println("\n" * "="^60 * "\nGBD algorithm\n" * "="^60)

  S_curr = Dict((s, y) => enestart[(s, y)] for s in sector for y in year_all)
  cuts = BendersCut[]
  LB, UB = -Inf, Inf
  k = 0

  while k < maxit
    k += 1
    println("\n" * "-"^40 * "\nITERATION ", k)

    # sub-problem
    sp, bal, _ = create_energy_subproblem(S_curr)
    optimize!(sp)
    if termination_status(sp) ∉ (MOI.OPTIMAL, MOI.LOCALLY_SOLVED)
      error("energy sub-problem failed: " * string(termination_status(sp)))
    end
    v = objective_value(sp) / 1000.0             # billion → trillion
    lambda = Dict{Tuple{String,Int},Float64}((s, y) => 0.0 for s in sector, y in year_all)
    for ((e, l, y), con) in bal
      s_idx = findfirst(t -> haskey(map_energy_sector, (e, l, t)), sector)
      if s_idx !== nothing
        s = sector[s_idx]
        lambda[(s, y)] += -dual(con) / 1000.0           # trillion
      end
    end
    push!(cuts, BendersCut("opt", v, lambda, copy(S_curr)))
    UB = min(UB, v)
    println("sub-problem cost v=$(round(v, digits=3)) trn USD")

    # master
    M = create_master_problem()
    set_gbd_master_bounds_and_initial_values!(M)
    foreach(cut -> add_opt_cut!(M, cut), cuts)
    optimize!(M)
    if termination_status(M) ∉ (MOI.OPTIMAL, MOI.LOCALLY_SOLVED)
      error("master infeasible: " * string(termination_status(M)))
    end
    theta = value(M[:theta])
    util = value(M[:UTILITY])
    obj = objective_value(M)
    LB = max(LB, obj)

    gap = abs((theta - util) - v)
    rel = gap / max(abs(theta - util), abs(v), 1.0)
    println("   theta=$(round(theta, digits=3))  UT=$(round(util, digits=3))  gap=$(round(gap, digits=6)) trn USD  rel=$(round(rel*100, digits=4))%")
    if rel <= tol
      println("\n✓ converged in $k iterations")
      return true
    end

    # update S
    S_curr = Dict((s, y) => value(M[:S][s, y]) for s in sector for y in year_all)
  end
  println("⚠️  max iterations reached")
  false
end

# driver
if abspath(PROGRAM_FILE) == @__FILE__
  solve_gbd()
end