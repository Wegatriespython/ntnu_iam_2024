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

# Only include files if not already in module context
if !@isdefined(year_all)  # Check if shared.jl already loaded
    include("model_config.jl");
    println("✓ loaded model configuration");
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
else
    println("✓ using existing module definitions");
end

# 1 Energy sub-problem (fixed service demands S_bar)
function create_energy_subproblem(S_bar::Dict{Tuple{String,Int},Float64}, config::ModelConfig = default_config())
  model = Model(Ipopt.Optimizer)
  set_optimizer_attribute(model, "print_level", 0)
  set_optimizer_attribute(model, "tol", config.solver_tolerance)
  set_optimizer_attribute(model, "max_iter", config.solver_max_iter)

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
    # Check if this energy-level combination is valid
    if !haskey(energy_level, (e, l))
      continue
    end
    
    # Check if this maps to a sector (PHYSENE demand)
    # Fixed: Properly identify PHYSENE commodities
    if haskey(map_energy_sector, (e, l, "ELEC")) && haskey(S_bar, ("ELEC", y))
      # This is electricity useful energy - use ELEC demand
      bal[(e, l, y)] = @constraint(model,
        sum(ACT[t, y] * (get(output, (t, e, l), 0) - get(input, (t, e, l), 0)) for t in technology)
        >=
        S_bar[("ELEC", y)])
    elseif haskey(map_energy_sector, (e, l, "NELE")) && haskey(S_bar, ("NELE", y))
      # This is non-electric useful energy - use NELE demand  
      bal[(e, l, y)] = @constraint(model,
        sum(ACT[t, y] * (get(output, (t, e, l), 0) - get(input, (t, e, l), 0)) for t in technology)
        >=
        S_bar[("NELE", y)])
    elseif haskey(demand, (e, l))
      # This has a standard demand (non-PHYSENE commodities)
      bal[(e, l, y)] = @constraint(model,
        sum(ACT[t, y] * (get(output, (t, e, l), 0) - get(input, (t, e, l), 0)) for t in technology)
        >=
        demand[(e, l)] * (gdp[y] / gdp[2020])^beta)
    else
      # No demand - net production must be zero (no free disposal)
      bal[(e, l, y)] = @constraint(model,
        sum(ACT[t, y] * (get(output, (t, e, l), 0) - get(input, (t, e, l), 0)) for t in technology)
        ==
        0)
    end
  end

  # capacity & technical constraints
  for t in technology, y in year_all
    if haskey(hours, t) && haskey(lifetime, t)
      idx = findfirst(==(y), year_all)
      @constraint(model,
        ACT[t, y] <= sum(CAP_NEW[t, year_all[i]] * hours[t]
                        for i in 1:idx if (idx - i + 1) * period_length <= lifetime[t]))
    end
  end
  for t in technology, y in year_all[2:end]
    if !haskey(diffusion_up, t)
      continue
    end
    y_prev = year_all[findfirst(==(y), year_all)-1]
    @constraint(model, CAP_NEW[t, y] <= CAP_NEW[t, y_prev] * (1 + diffusion_up[t])^period_length + get(startup, t, 0))
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
  @constraint(model, sum(EMISS[y] * period_length for y in year_all) == CUM_EMISS)

  # cost accounting
  year_idx = Dict(y => i for (i, y) in enumerate(year_all))
  lifetime_r = Dict{Tuple{String,Int},Vector{Int}}()
  for t in technology, y in year_all
    if !haskey(lifetime, t)
      continue
    end
    idx, lt = year_idx[y], lifetime[t]
    lifetime_r[(t, y)] = [y2 for y2 in year_all if year_idx[y2] <= idx && (idx - year_idx[y2] + 1) * period_length <= lt]
  end
  disc = Dict(y => (1 - drate)^(period_length * (year_idx[y] - 1)) for y in year_all)
  @constraint(model, [y in year_all],
    sum(ACT[t, y] * get(vom, (t, y), 0) for t in technology) +
    sum(sum(CAP_NEW[t, y2] * get(cost_capacity, (t, y2), 0) for y2 in lifetime_r[(t, y)]) for t in keys(lifetime))
    ==
    COST_ANNUAL[y])
  @constraint(model, sum(COST_ANNUAL[y] * period_length * disc[y] for y in year_all) == TOTAL_COST)

  # historical calibration
  for ((t, y), v) in energy_calibration
    fix(ACT[t, y], v; force=true)
  end
  for ((t, y), v) in energy_calibration_lo
    set_lower_bound(ACT[t, y], v)
  end

  @objective(model, Min, TOTAL_COST + 1e6 * sum(demand_slack))
  return model, bal, demand_slack
end

# 2 Master problem (MACRO + theta)
function create_master_problem()
  m = Model(Ipopt.Optimizer)
  set_optimizer_attribute(m, "print_level", 0)
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

  # energy cost - simple linear relationship with theta
  # theta represents total discounted cost from energy subproblem
  # We distribute it across years proportionally to cost_MESSAGE
  totW = sum(cost_MESSAGE[y] for y in year_all if y != 2020)
  @constraint(m, [y in year_all; y != 2020],
    EC[y] == theta * (cost_MESSAGE[y] / totW))

  # terminal investment
  if 2080 in year_all
    @constraint(m, I[2080] >= K[2080] * (grow[2080] + depr))
  end

  # objective
  @objective(m, Min, -UTILITY + theta)       # both trillion $

  # starts
  # Calculate expected theta based on cost_MESSAGE
  expected_theta = sum(cost_MESSAGE[y] for y in year_all if y != 2020)
  set_start_value(theta, expected_theta)
  
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
  # Set bounds and initial values for master problem variables similar to integrated model
  
  # Set lower bounds on all variables to avoid singularities (like macro model)
  for y in year_all
    set_lower_bound(M[:K][y], lotol * k0)
    set_lower_bound(M[:Y][y], lotol * y0)
    set_lower_bound(M[:C][y], lotol * c0)
    set_lower_bound(M[:I][y], lotol * i0)
    
    if y != 2020
      set_lower_bound(M[:KN][y], lotol * i0 * duration_period)
      set_lower_bound(M[:YN][y], lotol * y0 * newlab[y])
    end
    
    for s in sector
      set_lower_bound(M[:PRODENE][s, y], lotol * enestart[(s, y)] / aeei_factor[(s, y)])
      if y != 2020
        set_lower_bound(M[:NEWENE][s, y], lotol * enestart[(s, y)] / aeei_factor[(s, y)])
      end
    end
  end
  
  # Fix base year values to historical values
  fix(M[:Y][2020], y0; force=true)
  fix(M[:K][2020], k0; force=true)
  fix(M[:C][2020], c0; force=true)
  fix(M[:I][2020], i0; force=true)
  fix(M[:EC][2020], y0 - i0 - c0; force=true)
  
  for s in sector
    fix(M[:PRODENE][s, 2020], demand_base[s] / aeei_factor[(s, 2020)]; force=true)
  end
  
  # Energy service bounds and starting values
  for s in sector, y in year_all
    # Set tighter bounds based on reasonable growth from base demands
    base_demand = demand_MESSAGE[(s, 2020)]
    year_factor = (y - 2020) / 10 + 1  # Allow for growth over time
    set_lower_bound(M[:S][s, y], 0.5 * base_demand)  # At least 50% of base
    set_upper_bound(M[:S][s, y], 5.0 * base_demand * year_factor)  # Max 5x base with growth
    # Set starting values for S based on demand_MESSAGE
    set_start_value(M[:S][s, y], demand_MESSAGE[(s, y)])
    set_start_value(M[:PHYSENE][s, y], demand_MESSAGE[(s, y)])
  end
  
  # Set starting values for other variables
  set_start_value(M[:UTILITY], 100.0)
  for y in year_all
    if y != 2020
      set_start_value(M[:KN][y], i0 * duration_period)
      set_start_value(M[:YN][y], y0 * 0.1)
      for s in sector
        set_start_value(M[:NEWENE][s, y], demand_MESSAGE[(s, y)] * 0.1)
      end
    end
  end
end

# 4 Main GBD loop
function solve_gbd(maxit::Int=20, tol=1e-4)
  println("\n" * "="^60 * "\nGBD algorithm\n" * "="^60)

  S_curr = Dict((s, y) => enestart[(s, y)] for s in sector for y in year_all)
  cuts = BendersCut[]
  LB, UB = -Inf, Inf
  k = 0
  
  println("Initial S_curr values (enestart):")
  for s in sector
    println("  $s: ", [round(S_curr[(s, y)], digits=1) for y in year_all])
  end

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
    
    # Check demand slack usage
    println("Demand slack usage:")
    total_slack = 0.0
    for s in sector, y in year_all
      slack_val = value(sp[:demand_slack][s, y])
      if slack_val > 0.001
        println("  demand_slack[($s, $y)] = $(round(slack_val, digits=3)) PWh")
        total_slack += slack_val
      end
    end
    if total_slack > 0.001
      println("  Total slack: $(round(total_slack, digits=3)) PWh")
      println("  Slack penalty in objective: $(round(1e6 * total_slack / 1000.0, digits=3)) trillion USD")
    end
    
    # Extract actual TOTAL_COST without slack penalty
    actual_cost = value(sp[:TOTAL_COST])
    println("  Actual energy system cost: $(round(actual_cost, digits=3)) billion USD")
    println("  Actual energy system cost: $(round(actual_cost / 1000.0, digits=3)) trillion USD")
    
    # Show cost breakdown by year
    println("  Annual costs (billion USD):")
    for y in year_all
      annual_cost = value(sp[:COST_ANNUAL][y])
      println("    $y: $(round(annual_cost, digits=1)) billion USD")
    end
    
    # Check key metrics
    println("\n  Key energy metrics for 2020:")
    # Just show a few key technologies
    for t in ["coal_ppl", "gas_ppl", "electricity_grid", "appliances"]
      if t in technology
        act = value(sp[:ACT][t, 2020])
        if act > 0.01
          println("    ACT[$t] = $(round(act, digits=2)) PWh")
        end
      end
    end
    
    lambda = Dict{Tuple{String,Int},Float64}((s, y) => 0.0 for s in sector, y in year_all)
    for ((e, l, y), con) in bal
      # Fixed: Directly map energy-level to sectors using map_energy_sector
      if haskey(map_energy_sector, (e, l, "ELEC"))
        # Units: shadow price billion USD/PWh -> trillion USD/PWh  
        lambda[("ELEC", y)] = dual(con) / 1000.0
      elseif haskey(map_energy_sector, (e, l, "NELE"))
        # Units: shadow price billion USD/PWh -> trillion USD/PWh
        lambda[("NELE", y)] = dual(con) / 1000.0
      end
    end
    push!(cuts, BendersCut("opt", v, lambda, copy(S_curr)))
    UB = min(UB, v)
    println("sub-problem cost v=$(round(v, digits=3)) trn USD")
    println("Shadow prices (lambda):")
    for s in sector, y in year_all
      if lambda[(s, y)] != 0.0
        println("  λ[($s, $y)] = $(round(lambda[(s, y)], digits=6)) trillion USD/PWh")
      end
    end

    # master
    M = create_master_problem()
    set_gbd_master_bounds_and_initial_values!(M)
    
    println("\nAdding $(length(cuts)) Benders cuts to master problem...")
    foreach(cut -> add_opt_cut!(M, cut), cuts)
    
    println("Master problem has $(length(cuts)) Benders cuts")
    println("Master problem variables: $(num_variables(M))")
    println("Master problem constraints: $(num_constraints(M; count_variable_in_set_constraints=false))")
    
    # Try to identify infeasibility issues
    if k == 1 && length(cuts) == 1
      println("\nFirst cut details:")
      cut = cuts[1]
      println("  Cut cost v = $(round(cut.v, digits=3)) trillion USD")
      println("  Cut S_hat values:")
      for s in sector, y in year_all
        if haskey(cut.S_hat, (s, y))
          println("    S_hat[($s, $y)] = $(round(cut.S_hat[(s, y)], digits=3)) PWh")
        end
      end
    end
    
    optimize!(M)
    status = termination_status(M)
    println("\nMaster termination status: $status")
    if status ∉ (MOI.OPTIMAL, MOI.LOCALLY_SOLVED)
      error("master infeasible: " * string(status))
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
    S_new = Dict((s, y) => value(M[:S][s, y]) for s in sector for y in year_all)
    println("\n  Master problem S values:")
    for s in sector
      println("    $s: ", [round(S_new[(s, y)], digits=1) for y in year_all])
    end
    
    # Check if S changed
    max_change = 0.0
    for s in sector, y in year_all
      change = abs(S_new[(s, y)] - S_curr[(s, y)])
      max_change = max(max_change, change)
    end
    println("  Maximum S change: $(round(max_change, digits=6))")
    
    S_curr = S_new
  end
  println("⚠️  max iterations reached")
  false
end

# driver
if abspath(PROGRAM_FILE) == @__FILE__
  solve_gbd()
end