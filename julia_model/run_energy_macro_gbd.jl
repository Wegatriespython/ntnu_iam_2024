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
  include("model_config.jl")
  println("✓ loaded model configuration")
  include("shared.jl")
  println("✓ loaded shared definitions")
  include("energy/energy_model_world.jl")
  println("✓ loaded energy model")
  include("energy_sub_problem.jl")
  println("✓ loaded energy subproblem")
  include("macro/macro_data_load.jl")
  println("✓ loaded macro data")
  include("macro/macro_core.jl")
  println("✓ loaded macro core")
  include("macro/macro_presolve.jl")
  println("✓ loaded macro preprocessing")
  include("run_energy_macro.jl")
  println("✓ loaded integrated model functions")
end

# 1 Energy sub-problem (fixed service demands S_bar) - now uses external file
# This function is a wrapper that calls the version in energy_sub_problem.jl

# 2 Master problem (MACRO + theta)
function create_master_problem()
  # Use the enhanced integrated model in GBD mode
  return create_integrated_model(gbd_mode=true)
end

# 3 Benders infrastructure
mutable struct BendersCut
  type::String                       # "optimality" | "feasibility"
  v::Float64                         # cost value (trn $)
  lambda::Dict{Tuple{String,Int},Float64} # duals (trn $ / PWh)
  S_hat::Dict{Tuple{String,Int},Float64} # incumbent S
end

function add_opt_cut!(M::Model, cut::BendersCut)
  # Filter out zero or near-zero shadow prices to avoid numerical issues
  active_keys = [key for key in keys(cut.lambda) if abs(cut.lambda[key]) > 1e-8]
  
  if length(active_keys) == 0
    println("  Warning: All shadow prices are zero or near-zero. Adding simple cut: theta >= $(cut.v)")
    @constraint(M, M[:theta] >= cut.v)
  else
    println("  Adding cut with $(length(active_keys)) active shadow prices")
    @constraint(M, M[:theta] >= cut.v + sum(cut.lambda[key] * (M[:PHYSENE][key...] - cut.S_hat[key]) for key in active_keys))
    
    # Debug: show the cut components
    println("  Cut components:")
    println("    Constant term: $(cut.v)")
    for key in active_keys
      if abs(cut.lambda[key]) > 1e-6  # Only show significant terms
        println("    λ[$(key)] * (PHYSENE[$(key)] - $(round(cut.S_hat[key], digits=1))) = $(round(cut.lambda[key], digits=6)) * (PHYSENE[$(key)] - $(round(cut.S_hat[key], digits=1)))")
      end
    end
  end
end

# This function is now replaced by the one in run_energy_macro.jl

# 4 Main GBD loop
function solve_gbd(maxit::Int=30, tol=1e-2, config::Union{ModelConfig,Nothing}=nothing)
  println("\n" * "="^60 * "\nGBD algorithm\n" * "="^60)

  # Use default config if not provided
  if config === nothing
    config = @isdefined(default_config) ? default_config() : ModelConfig()
  end

  # Use exact integrated solution PHYSENE values as starting point
  S_curr = Dict{Tuple{String,Int},Float64}()
  
  # Exact integrated solution values from integrated_debug_output.txt
  elec_values = [0.8, 21.9, 23.4, 26.1, 29.5, 29.8, 30.1]
  nele_values = [83.3, 130.9, 167.5, 201.9, 235.3, 240.4, 244.6]
  
  for (i, y) in enumerate(year_all)
    S_curr[("ELEC", y)] = elec_values[i]
    S_curr[("NELE", y)] = nele_values[i]
  end
  cuts = BendersCut[]
  LB, UB = -Inf, Inf
  util_curr = 0.0  # Track current utility for proper UB calculation
  k = 0

  println("Initial S_curr values (enestart):")
  for s in sector
    println("  $s: ", [round(S_curr[(s, y)], digits=1) for y in year_all])
  end

  while k < maxit
    k += 1
    println("\n" * "-"^40 * "\nITERATION ", k)

    # sub-problem
    sp, bal = create_energy_subproblem(S_curr; config=config)
    optimize!(sp)
    if termination_status(sp) ∉ (MOI.OPTIMAL, MOI.LOCALLY_SOLVED)
      error("energy sub-problem failed: " * string(termination_status(sp)))
    end
    v = objective_value(sp) / 1000.0             # billion → trillion

    # Check demand slack usage (commented out - no slack variables in current model)
    # println("Demand slack usage:")
    # total_slack = 0.0
    # for s in sector, y in year_all
    #   slack_val = value(sp[:demand_slack][s, y])
    #   if slack_val > 0.001
    #     println("  demand_slack[($s, $y)] = $(round(slack_val, digits=3)) PWh")
    #     total_slack += slack_val
    #   end
    # end
    # if total_slack > 0.001
    #   println("  Total slack: $(round(total_slack, digits=3)) PWh")
    #   println("  Slack penalty in objective: $(round(1e6 * total_slack / 1000.0, digits=3)) trillion USD")
    # end
    println("Energy subproblem solved successfully")

    # Extract actual TOTAL_COST without slack penalty
    actual_cost = value(sp[:TOTAL_COST])
    println("  Actual energy system cost: $(round(actual_cost, digits=3)) billion USD")
    println("  Actual energy system cost: $(round(actual_cost / 1000.0, digits=3)) trillion USD")
    
    # Store annual costs for proper EC constraint in master problem
    annual_costs = Dict{Int,Float64}()
    for y in year_all
      annual_costs[y] = value(sp[:COST_ANNUAL][y])
    end

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
      # Fixed: Sum all balance constraints that map to same sector (use += not =)
      if haskey(map_energy_sector, (e, l, "ELEC"))
        # Units: shadow price billion USD/PWh -> trillion USD/PWh  
        lambda[("ELEC", y)] += dual(con) / 1000.0
      elseif haskey(map_energy_sector, (e, l, "NELE"))
        # Units: shadow price billion USD/PWh -> trillion USD/PWh
        lambda[("NELE", y)] += dual(con) / 1000.0
      end
    end
    push!(cuts, BendersCut("opt", v, lambda, copy(S_curr)))
    
    # Store the most recent annual costs for the master problem
    global recent_annual_costs = annual_costs
    println("sub-problem cost v=$(round(v, digits=3)) trn USD")
    println("Shadow prices (lambda):")
    for s in sector, y in year_all
      if lambda[(s, y)] != 0.0
        println("  λ[($s, $y)] = $(round(lambda[(s, y)], digits=6)) trillion USD/PWh")
      end
    end

    # master
    M = create_master_problem()
    set_macro_bounds_and_initial_values!(M)  # Use existing function from macro model
    set_gbd_bounds_and_initial_values!(M)    # Use new function from run_energy_macro.jl

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

    # Debug: Check theta bounds and constraints before solving
    println("\n=== DEBUG: Master Problem Analysis Before Solving ===")
    theta_var = M[:theta]
    println("  Theta variable bounds: [$(lower_bound(theta_var)), $(upper_bound(theta_var))]")
    println("  Theta variable: $theta_var")
    
    # Check starting values of PHYSENE variables
    println("  PHYSENE starting values:")
    for s in sector
      start_values = []
      for y in year_all
        try
          start_val = start_value(M[:PHYSENE][s, y])
          push!(start_values, start_val !== nothing ? round(start_val, digits=1) : "none")
        catch
          push!(start_values, "error")
        end
      end
      println("    $s: $start_values")
    end
    
    # Show existing Benders cuts involving theta
    if length(cuts) > 0
      println("  Current Benders cuts requiring:")
      for (i, cut) in enumerate(cuts)
        println("    Cut $i: theta >= $(round(cut.v, digits=3)) + linear_terms")
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
    
    # Debug: Check theta constraint violations after solving
    println("\n=== DEBUG: Master Problem Solution Analysis ===")
    println("  Solved theta value: $(round(theta, digits=3)) trillion USD")
    println("  Energy subproblem cost: $(round(v, digits=3)) trillion USD")
    println("  Discrepancy (energy_cost - theta): $(round(v - theta, digits=3)) trillion USD")
    
    # Check if any Benders cuts are violated
    if length(cuts) > 0
      println("  Benders cut constraint checks:")
      for (i, cut) in enumerate(cuts)
        # Calculate the current cut value
        cut_rhs = cut.v
        println("    Cut $i details:")
        println("      Base cost v = $(round(cut.v, digits=3)) trillion USD")
        linear_term = 0.0
        for (s, y) in keys(cut.lambda)
          if haskey(cut.S_hat, (s, y)) && abs(cut.lambda[(s, y)]) > 1e-8
            physene_val = value(M[:PHYSENE][s, y])
            s_hat_val = cut.S_hat[(s, y)]
            lambda_val = cut.lambda[(s, y)]
            term_val = lambda_val * (physene_val - s_hat_val)
            linear_term += term_val
            cut_rhs += term_val
            println("      λ[($s,$y)] = $(round(lambda_val, digits=6)), PHYSENE = $(round(physene_val, digits=1)), S_hat = $(round(s_hat_val, digits=1)), term = $(round(term_val, digits=6))")
          end
        end
        println("      Total linear term = $(round(linear_term, digits=6))")
        println("      Cut RHS = $(round(cut.v, digits=3)) + $(round(linear_term, digits=6)) = $(round(cut_rhs, digits=3))")
        violation = theta - cut_rhs
        status_str = violation >= -1e-6 ? "✓ satisfied" : "❌ VIOLATED"
        println("      Constraint: theta($(round(theta, digits=3))) >= $(round(cut_rhs, digits=3)), violation = $(round(violation, digits=6)) - $status_str")
        
        # Special check for first iteration when PHYSENE should equal S_hat
        if k == 1
          println("      First iteration check - PHYSENE vs S_hat differences:")
          for (s, y) in keys(cut.S_hat)
            physene_val = value(M[:PHYSENE][s, y])
            s_hat_val = cut.S_hat[(s, y)]
            diff = abs(physene_val - s_hat_val)
            if diff > 1e-6
              println("        ⚠️ PHYSENE[($s,$y)] = $(round(physene_val, digits=3)) ≠ S_hat = $(round(s_hat_val, digits=3)), diff = $(round(diff, digits=6))")
            end
          end
        end
      end
    end
    
    # Check energy cost constraints (EC constraints) in master problem
    println("  Energy cost constraint status:")
    try
      if haskey(M.obj_dict, :EC)
        ec_constraints = M[:EC]
        println("    EC variables type: $(typeof(ec_constraints))")
        println("    EC variables indices: $(ec_constraints.axes)")
        
        # Use num_constraints to get total count
        total_constraints = num_constraints(M; count_variable_in_set_constraints=false)
        println("  Model reports $total_constraints total constraints")
        
        # Check all constraints in the model using list_of_constraint_types
        constraint_types = list_of_constraint_types(M)
        println("  Constraint types in model:")
        physene_constraint_count = 0
        
        for (F, S) in constraint_types
          constraints = all_constraints(M, F, S)
          println("    $(length(constraints)) constraints of type $F in $S")
          
          # Check first few constraints of each type for PHYSENE
          for (i, con_ref) in enumerate(constraints)
            if i <= 2  # Only check first 2 of each type
              try
                con_obj = constraint_object(con_ref)
                func_str = string(con_obj.func)
                if occursin("PHYSENE", func_str)
                  physene_constraint_count += 1
                  println("      PHYSENE constraint: $(con_obj.func) ∈ $(con_obj.set)")
                end
              catch
                # Skip if can't access
              end
            end
          end
        end
        println("  Total PHYSENE-related constraints found: $physene_constraint_count")
        
      else
        println("    No EC variables found in master problem")
        println("    Available keys: $(keys(M.obj_dict))")
      end
    catch e
      println("    EC constraint investigation failed: $e")
    end
    
    # Update current utility for proper UB calculation
    util_curr = util
    
    # Correct GBD gap calculation for Max(UTILITY) objective:
    # We want to maximize UTILITY subject to energy cost constraints
    # Upper Bound = Current utility from master problem 
    # Lower Bound = Best achievable utility (will improve with iterations)
    # When energy subproblem gives cost v, master finds utility util_curr with that constraint
    UB = max(UB, util_curr)  # Upper bound is best utility found so far
    LB = util_curr           # Lower bound is current iteration utility
    gap = abs(UB - LB)
    rel = gap / max(abs(UB), abs(LB), 1.0)
    println("   theta=$(round(theta, digits=3))  UT=$(round(util, digits=3))  obj=$(round(obj, digits=3))  UB=$(round(UB, digits=3))  LB=$(round(LB, digits=3))  gap=$(round(gap, digits=6)) trn USD  rel=$(round(rel*100, digits=4))%")
    if rel <= tol
      println("\n✓ converged in $k iterations")
      return true
    end

    # update S (use PHYSENE from the new model structure)
    S_new = Dict((s, y) => value(M[:PHYSENE][s, y]) for s in sector for y in year_all)
    println("\n  Master problem PHYSENE values:")
    for s in sector
      println("    $s: ", [round(S_new[(s, y)], digits=1) for y in year_all])
    end

    # Check if S changed
    max_change = 0.0
    for s in sector, y in year_all
      change = abs(S_new[(s, y)] - S_curr[(s, y)])
      max_change = max(max_change, change)
    end
    println("  Maximum PHYSENE change: $(round(max_change, digits=6))")
    
    # Check for stalling (no significant change in PHYSENE values)
    if max_change < 1e-6
      println("  ⚠️ Algorithm stalled - PHYSENE values not changing")
      println("  This suggests convergence to a suboptimal point or numerical issues")
      
      # Compare to integrated solution for debugging
      println("  Comparison to integrated solution:")
      integrated_elec = [0.8, 21.9, 23.4, 26.1, 29.5, 29.8, 30.1]
      integrated_nele = [83.3, 130.9, 167.5, 201.9, 235.3, 240.4, 244.6]
      
      for (i, y) in enumerate(year_all)
        gbd_elec = S_new[("ELEC", y)]
        gbd_nele = S_new[("NELE", y)]
        diff_elec = abs(gbd_elec - integrated_elec[i])
        diff_nele = abs(gbd_nele - integrated_nele[i])
        println("    $y: ELEC diff=$(round(diff_elec, digits=1)), NELE diff=$(round(diff_nele, digits=1))")
      end
    end

    S_curr = S_new
  end
  println("⚠️  max iterations reached")
  false
end

# driver
if abspath(PROGRAM_FILE) == @__FILE__
  solve_gbd()
end
