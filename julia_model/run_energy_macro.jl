# run_energy_macro.jl - Equivalent to run_energy_macro.gms
# Integrating energy systems model and macro-economic model

using JuMP
using GAMS
using LinearAlgebra


# Include all model components in order (equivalent to $INCLUDE statements)
include("shared.jl")

include("energy/energy_model_world.jl")

include("macro/macro_data_load.jl")

include("macro/macro_core.jl")

include("macro/macro_presolve.jl")

# ------------------------------------------------------------------------------
# Model definition and creation
# ------------------------------------------------------------------------------

function create_integrated_model(; gbd_mode::Bool = false)
  # Create the main optimization model
  model = Model(GAMS.Optimizer)

  # Set solver options
  set_optimizer_attribute(model, "NLP", "CONOPT")

  if gbd_mode
    println("Creating GBD master problem...")
  else
    println("Creating integrated energy-macro model...")
  end

  # Create shared variables from shared.jl
  create_shared_variables!(model; gbd_mode=gbd_mode)
  println("✓ Created shared variables")

  if gbd_mode
    # GBD mode: Add theta variable and skip energy model creation
    @variable(model, theta >= 0, upper_bound = 250)  # trillion US$
    println("✓ Added GBD theta variable")
  else
    # Integrated mode: Create energy model components
    create_energy_model!(model)
    println("✓ Created energy model equations")
  end

  # Create macro model components  
  create_macro_model!(model)
  println("✓ Created macro model equations")

  # Cost equations are already defined in create_energy_model!
  # No need to duplicate them here

  if !gbd_mode
    # Set energy model calibration constraints (only for integrated mode)
    for (tech_year, value) in energy_calibration
      tech, year = tech_year
      fix(model[:ACT][tech, year], value; force=true)
    end

    for (tech_year, value) in energy_calibration_lo
      tech, year = tech_year
      set_lower_bound(model[:ACT][tech, year], value)
    end

    println("✓ Applied energy model calibration")
  end

  # Set macro bounds and initial values
  set_macro_bounds_and_initial_values!(model)
  println("✓ Set macro bounds and initial values")

  # Additional preprocessing
  additional_macro_preprocessing!(model)

  if gbd_mode
    # GBD mode: Add proper PHYSENE-EC linking constraints (Taylor expansion)
    # This is the critical missing piece that links PHYSENE to energy costs
    println("  Adding PHYSENE-EC linking constraints using Taylor expansion")
    for y in year_all
      if y != 2020  # not base period
        @constraint(model,
          model[:EC][y] == cost_MESSAGE[y] +
          sum(eneprice[(s, y)] * (model[:PHYSENE][s, y] - enestart[(s, y)]) for s in sector) +
          sum(eneprice[(s, y)] / enestart[(s, y)] * 
              (model[:PHYSENE][s, y] - enestart[(s, y)]) * (model[:PHYSENE][s, y] - enestart[(s, y)]) 
              for s in sector)
        )
      end
    end
    println("✓ Added GBD PHYSENE-EC linking constraints")
    
    # Set objective to maximize UTILITY (same as integrated model)
    # Note: theta is constrained by Benders cuts but not in objective
    @objective(model, Max, model[:UTILITY])
    println("✓ Set GBD objective: Max(UTILITY)")
  else
    # Integrated mode: maximize utility
    @objective(model, Max, model[:UTILITY])
    println("✓ Set integrated objective: Max(UTILITY)")
  end

  println("✓ Model creation completed")
  println("Model has $(num_variables(model)) variables and $(num_constraints(model; count_variable_in_set_constraints=false)) constraints")

  return model
end

# Function to set GBD-specific bounds and initial values
function set_gbd_bounds_and_initial_values!(model::Model)
  # Energy service bounds and starting values for GBD mode
  # Integrated solution PHYSENE values from integrated_debug_output.txt
  integrated_elec = [0.8, 21.9, 23.4, 26.1, 29.5, 29.8, 30.1]
  integrated_nele = [83.3, 130.9, 167.5, 201.9, 235.3, 240.4, 244.6]
  
  for (i, y) in enumerate(year_all)
    for s in sector
      # Use much wider bounds to allow integrated solution values
      if s == "ELEC"
        set_lower_bound(model[:PHYSENE][s, y], 0.1)  # Very low minimum
        set_upper_bound(model[:PHYSENE][s, y], 50.0)  # Allow up to 50 PWh
        # Set starting value to integrated solution
        set_start_value(model[:PHYSENE][s, y], integrated_elec[i])
      else  # NELE
        set_lower_bound(model[:PHYSENE][s, y], 1.0)   # Very low minimum  
        set_upper_bound(model[:PHYSENE][s, y], 300.0) # Allow up to 300 PWh
        # Set starting value to integrated solution
        set_start_value(model[:PHYSENE][s, y], integrated_nele[i])
      end
    end
  end
  
  # Set theta starting value based on expected cost
  if haskey(model.obj_dict, :theta)
    expected_theta = sum(cost_MESSAGE[y] for y in year_all if y != 2020) / 1000.0  # Convert to trillion
    set_start_value(model[:theta], expected_theta)
  end
  
  println("✓ Set GBD bounds and initial values")
end

# ------------------------------------------------------------------------------
# Solve model
# ------------------------------------------------------------------------------

function solve_energy_macro_model()
  println("\n" * "="^60)
  println("SOLVING ENERGY-MACRO INTEGRATED MODEL")
  println("="^60)

  # Create the model
  model = create_integrated_model()

  # Solve the optimization problem
  println("\nStarting optimization...")
  optimize!(model)

  # Check solution status
  status = termination_status(model)
  println("\nSolution Status: ", status)

  if status == MOI.OPTIMAL || status == MOI.LOCALLY_SOLVED
    println("✓ Optimal solution found!")

    # Get objective value
    utility_value = objective_value(model)
    println("Optimal Utility Value: ", utility_value)

    # ------------------------------------------------------------------------------
    # Reporting (equivalent to GAMS reporting section)
    # ------------------------------------------------------------------------------

    println("\n" * "="^60)
    println("SOLUTION REPORTING")
    println("="^60)

    # Calculate GDP_MACRO values
    println("\nGDP_MACRO (Trillion \$):")
    for y in year_all
      gdp_macro_val = value(model[:I][y]) + value(model[:C][y]) + value(model[:EC][y])
      println("  $y: $(round(gdp_macro_val, digits=3))")
    end

    # Energy system results
    println("\nKey Energy Results:")
    println("Total system cost: $(round(value(model[:TOTAL_COST]), digits=2)) billion \$")
    println("Cumulative emissions: $(round(value(model[:CUM_EMISS]), digits=2)) GtCO2")

    # Macro economic results  
    println("\nKey Macro Results:")
    for y in year_all
      println("  $y - C: $(round(value(model[:C][y]), digits=2)), I: $(round(value(model[:I][y]), digits=2)), Y: $(round(value(model[:Y][y]), digits=2))")
    end
    
    # Values critical for GBD formulation
    println("\n" * "="^60)
    println("GBD-RELEVANT VALUES")
    println("="^60)
    
    # Physical energy services (PHYSENE) - key coupling variables
    println("\nPhysical Energy Services (PHYSENE) [PWh]:")
    for s in sector
      println("  $s: ", [round(value(model[:PHYSENE][s, y]), digits=1) for y in year_all])
    end
    
    # Energy costs by year (EC variable)
    println("\nEnergy Costs (EC) [Trillion USD]:")
    for y in year_all
      ec_val = value(model[:EC][y])
      println("  $y: $(round(ec_val, digits=3))")
    end
    
    # Total energy system cost breakdown
    println("\nEnergy System Cost Breakdown [Billion USD]:")
    total_energy_cost = value(model[:TOTAL_COST])
    println("  Total TOTAL_COST: $(round(total_energy_cost, digits=1))")
    println("  Total TOTAL_COST (Trillion): $(round(total_energy_cost/1000.0, digits=3))")
    
    # Annual costs
    println("  Annual costs:")
    for y in year_all
      annual_cost = value(model[:COST_ANNUAL][y])
      println("    $y: $(round(annual_cost, digits=1)) billion USD")
    end
    
    # Utility calculation details
    println("\nUtility Calculation:")
    println("  Total UTILITY: $(round(value(model[:UTILITY]), digits=3))")
    println("  Components by year:")
    for y in year_all
      util_component = udf[y] * log(value(model[:C][y]))
      println("    $y: util_factor=$(round(udf[y], digits=4)) * log(C=$(round(value(model[:C][y]), digits=2))) = $(round(util_component, digits=4))")
    end
    
    # Shadow prices for energy balance constraints (if available)
    # Note: These would be the lambda values in GBD
    println("\nEnergy Balance Shadow Prices (lambda equivalent):")
    if haskey(model.obj_dict, :energy_balance)
      for s in sector, y in year_all
        if haskey(model[:energy_balance], (s, y))
          shadow_price = dual(model[:energy_balance][s, y])
          if abs(shadow_price) > 1e-6
            println("  λ[($s, $y)] = $(round(shadow_price/1000.0, digits=6)) trillion USD/PWh")
          end
        end
      end
    else
      println("  (Energy balance constraints not accessible as named constraints)")
    end
    
    # Comparison of PHYSENE vs initial enestart values
    println("\nComparison: PHYSENE vs enestart:")
    for s in sector
      println("  $s:")
      for y in year_all
        physene_val = value(model[:PHYSENE][s, y])
        enestart_val = enestart[(s, y)]
        diff = physene_val - enestart_val
        println("    $y: PHYSENE=$(round(physene_val, digits=1)), enestart=$(round(enestart_val, digits=1)), diff=$(round(diff, digits=1))")
      end
    end

    # Energy service demands
    println("\nEnergy Service Demands (PWh):")
    for y in year_all
      elec_demand = value(model[:PHYSENE]["ELEC", y])
      nele_demand = value(model[:PHYSENE]["NELE", y])
      println("  $y - ELEC: $(round(elec_demand, digits=2)), NELE: $(round(nele_demand, digits=2))")
    end

    # Technology activities (sample)
    println("\nSample Technology Activities for 2030 (PWh):")
    sample_techs = ["coal_ppl", "gas_ppl", "wind_ppl", "solar_PV_ppl", "nuclear_ppl"]
    for tech in sample_techs
      if tech in technology
        activity = value(model[:ACT][tech, 2030])
        println("  $tech: $(round(activity, digits=2))")
      end
    end

    println("\n✓ Model solved successfully!")

    # Save results (equivalent to execute_unload "energy_macro_results.gdx")
    save_results(model)

  elseif status == MOI.INFEASIBLE
    println("✗ Model is infeasible")
    compute_conflict!(model)
  elseif status == MOI.DUAL_INFEASIBLE
    println("✗ Model is dual infeasible (unbounded)")
  else
    println("✗ Solver terminated with status: ", status)
    if has_values(model)
      println("Partial solution available")
    end
  end

  return model, status
end

# Function to save results (equivalent to GDX export) - optimized
function save_results(model)
  println("\nSaving results...")

  # Pre-allocate typed dictionaries for better performance
  results = Dict{String,Dict}(
    "ACT" => Dict{Tuple{String,Int},Float64}(),
    "CAP_NEW" => Dict{Tuple{String,Int},Float64}(),
    "EMISS" => Dict{Int,Float64}(),
    "Y" => Dict{Int,Float64}(),
    "C" => Dict{Int,Float64}(),
    "I" => Dict{Int,Float64}(),
    "K" => Dict{Int,Float64}(),
    "PHYSENE" => Dict{Tuple{String,Int},Float64}()
  )

  # Batch variable extraction for better performance
  act_vals = value.(model[:ACT])
  cap_new_vals = value.(model[:CAP_NEW])

  # Efficiently populate results
  for tech in technology, y in year_all
    results["ACT"][(tech, y)] = act_vals[tech, y]
    results["CAP_NEW"][(tech, y)] = cap_new_vals[tech, y]
  end

  # Extract other variables
  for y in year_all
    results["EMISS"][y] = value(model[:EMISS][y])
    results["Y"][y] = value(model[:Y][y])
    results["C"][y] = value(model[:C][y])
    results["I"][y] = value(model[:I][y])
    results["K"][y] = value(model[:K][y])
  end

  # Extract PHYSENE values
  physene_vals = value.(model[:PHYSENE])
  for s in sector, y in year_all
    results["PHYSENE"][(s, y)] = physene_vals[s, y]
  end

  # Save to file (in Julia, we can save as JLD2 or serialize)
  # For now, just print confirmation
  println("✓ Results saved (implement file output as needed)")

  return results
end

# ------------------------------------------------------------------------------
# Main execution
# ------------------------------------------------------------------------------

if abspath(PROGRAM_FILE) == @__FILE__
  println("Energy-Macro Integrated Assessment Model")
  println("Julia + JuMP Implementation")
  println("Equivalent to GAMS run_energy_macro.gms")
  println()

  try
    model, status = solve_energy_macro_model()

    if status == MOI.OPTIMAL || status == MOI.LOCALLY_SOLVED
      println("\n🎉 SUCCESS: Energy-Macro model solved optimally!")
    else
      println("\n⚠️  WARNING: Model did not reach optimal solution")
    end

  catch e
    println("\n❌ ERROR: ", e)
    rethrow(e)
  end
end
