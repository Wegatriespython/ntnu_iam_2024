# run_energy_macro.jl - Equivalent to run_energy_macro.gms
# Integrating energy systems model and macro-economic model

using JuMP
using Ipopt
using LinearAlgebra

println("Starting Energy-Macro Integrated Model...")

# Include all model components in order (equivalent to $INCLUDE statements)
include("shared.jl")
println("✓ Loaded shared definitions")

include("energy_model_world.jl") 
println("✓ Loaded energy model")

include("macro_data_load.jl")
println("✓ Loaded macro data")

include("macro_core.jl")
println("✓ Loaded macro core")

include("macro_presolve.jl")
println("✓ Loaded macro preprocessing")

# ------------------------------------------------------------------------------
# Model definition and creation
# ------------------------------------------------------------------------------

function create_integrated_model()
    # Create the main optimization model
    model = Model(Ipopt.Optimizer)
    
    # Set solver options
    set_optimizer_attribute(model, "print_level", 5)
    set_optimizer_attribute(model, "max_iter", 3000)
    set_optimizer_attribute(model, "tol", 1e-6)
    
    println("Creating integrated energy-macro model...")
    
    # Define shared variables first (from shared.jl)
    @variable(model, PHYSENE[sector, year_all] >= 0)  # Physical end-use service use
    @variable(model, COST_ANNUAL[year_all])           # Annual costs
    
    println("✓ Created shared variables")
    
    # Create energy model components
    create_energy_model!(model)
    println("✓ Created energy model equations")
    
    # Create macro model components  
    create_macro_model!(model)
    println("✓ Created macro model equations")
    
    # Add energy model cost equations that depend on shared variables
    # Pre-compute indices and lifetime ranges for performance
    year_indices = Dict(y => i for (i, y) in enumerate(year_all))
    tech_lifetimes = Dict(tech => get(lifetime, tech, 0) for tech in technology if haskey(lifetime, tech))
    
    # Pre-compute lifetime ranges for each technology and year combination
    lifetime_ranges = Dict{Tuple{String, Int}, Vector{Int}}()
    for tech in keys(tech_lifetimes), y in year_all
        y_idx = year_indices[y]
        lt = tech_lifetimes[tech]
        valid_years = [y2 for y2 in year_all 
                       if year_indices[y2] <= y_idx && 
                          (y_idx - year_indices[y2] + 1) * period_length <= lt]
        lifetime_ranges[(tech, y)] = valid_years
    end
    
    # Pre-compute discount factors (matching GAMS formula)
    discount_factors = Dict(y => (1 - discount_rate)^(period_length * (year_indices[y] - 1)) for y in year_all)
    
    # EQ_COST_ANNUAL - costs per year (optimized)
    for y in year_all
        @constraint(model,
            sum(model[:ACT][tech, y] * get(vom, (tech, y), 0) for tech in technology) +
            sum(sum(model[:CAP_NEW][tech, y2] * get(cost_capacity, (tech, y2), 0)
                    for y2 in lifetime_ranges[(tech, y)])
                for tech in keys(tech_lifetimes)) == COST_ANNUAL[y]
        )
    end
    
    # EQ_COST - total discounted system costs (optimized)
    @constraint(model,
        sum(COST_ANNUAL[y] * period_length * discount_factors[y] for y in year_all) == model[:TOTAL_COST]
    )
    
    println("✓ Added energy cost equations")
    
    # Set energy model calibration constraints
    for (tech_year, value) in energy_calibration
        tech, year = tech_year
        fix(model[:ACT][tech, year], value; force=true)
    end
    
    for (tech_year, value) in energy_calibration_lo
        tech, year = tech_year
        set_lower_bound(model[:ACT][tech, year], value)
    end
    
    println("✓ Applied energy model calibration")
    
    # Set macro bounds and initial values
    set_macro_bounds_and_initial_values!(model)
    println("✓ Set macro bounds and initial values")
    
    # Additional preprocessing
    additional_macro_preprocessing!(model)
    
    # Set objective function to maximize utility
    @objective(model, Max, model[:UTILITY])
    
    println("✓ Model creation completed")
    println("Model has $(num_variables(model)) variables and $(num_constraints(model; count_variable_in_set_constraints=false)) constraints")
    
    return model
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
    results = Dict{String, Dict}(
        "ACT" => Dict{Tuple{String,Int}, Float64}(),
        "CAP_NEW" => Dict{Tuple{String,Int}, Float64}(),
        "EMISS" => Dict{Int, Float64}(),
        "Y" => Dict{Int, Float64}(),
        "C" => Dict{Int, Float64}(),
        "I" => Dict{Int, Float64}(),
        "K" => Dict{Int, Float64}(),
        "PHYSENE" => Dict{Tuple{String,Int}, Float64}()
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