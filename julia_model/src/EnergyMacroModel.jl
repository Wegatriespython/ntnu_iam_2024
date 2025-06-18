module EnergyMacroModel

using JuMP
using Ipopt
using LinearAlgebra

# Include all model components in the module scope (only once)
include("../shared.jl")
include("../energy_model_world.jl") 
include("../macro_data_load.jl")
include("../macro_core.jl")
include("../macro_presolve.jl")

# Define the solve function directly here to avoid including run_energy_macro.jl
# which would cause duplicate includes
function solve_energy_macro_model()
    println("\n" * "="^60)
    println("SOLVING ENERGY-MACRO INTEGRATED MODEL")
    println("="^60)
    
    # Create the model
    model = create_integrated_model()
    
    # Solve the optimization problem
    println("\nStarting optimization...")
    solve_start = time()
    optimize!(model)
    solve_time = time() - solve_start
    println("IPOPT solve time: $(round(solve_time, digits=2)) seconds")
    
    # Check solution status
    status = termination_status(model)
    println("\nSolution Status: ", status)
    
    if status == MOI.OPTIMAL || status == MOI.LOCALLY_SOLVED
        println("✓ Optimal solution found!")
        
        # Get objective value
        utility_value = objective_value(model)
        println("Optimal Utility Value: ", utility_value)
        
        # Solution reporting
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
        
        # Save results
        save_results(model)
        
        println("\n✓ Model solved successfully!")
        
    else
        println("✗ Solver terminated with status: ", status)
    end
    
    return model, status
end

# Model creation function (copy from run_energy_macro.jl to avoid duplicate includes)
function create_integrated_model()
    model = Model(Ipopt.Optimizer)
    set_optimizer_attribute(model, "print_level", 5)
    set_optimizer_attribute(model, "max_iter", 3000)
    set_optimizer_attribute(model, "tol", 1e-6)
    
    println("Creating integrated energy-macro model...")
    
    @variable(model, PHYSENE[sector, year_all] >= 0)
    @variable(model, COST_ANNUAL[year_all])
    
    create_energy_model!(model)
    create_macro_model!(model)
    
    # Add cost constraints (optimized version)
    year_indices = Dict(y => i for (i, y) in enumerate(year_all))
    tech_lifetimes = Dict(tech => get(lifetime, tech, 0) for tech in technology if haskey(lifetime, tech))
    
    lifetime_ranges = Dict{Tuple{String, Int}, Vector{Int}}()
    for tech in keys(tech_lifetimes), y in year_all
        y_idx = year_indices[y]
        lt = tech_lifetimes[tech]
        valid_years = [y2 for y2 in year_all 
                       if year_indices[y2] <= y_idx && 
                          (y_idx - year_indices[y2] + 1) * period_length <= lt]
        lifetime_ranges[(tech, y)] = valid_years
    end
    
    discount_factors = Dict(y => (1 - discount_rate)^(period_length * (year_indices[y] - 1)) for y in year_all)
    
    for y in year_all
        @constraint(model,
            sum(model[:ACT][tech, y] * get(vom, (tech, y), 0) for tech in technology) +
            sum(sum(model[:CAP_NEW][tech, y2] * get(cost_capacity, (tech, y2), 0)
                    for y2 in lifetime_ranges[(tech, y)])
                for tech in keys(tech_lifetimes)) == COST_ANNUAL[y]
        )
    end
    
    @constraint(model,
        sum(COST_ANNUAL[y] * period_length * discount_factors[y] for y in year_all) == model[:TOTAL_COST]
    )
    
    # Set calibration and bounds
    for (tech_year, value) in energy_calibration
        tech, year = tech_year
        fix(model[:ACT][tech, year], value; force=true)
    end
    
    for (tech_year, value) in energy_calibration_lo
        tech, year = tech_year
        set_lower_bound(model[:ACT][tech, year], value)
    end
    
    set_macro_bounds_and_initial_values!(model)
    additional_macro_preprocessing!(model)
    
    @objective(model, Max, model[:UTILITY])
    
    return model
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
    
    println("✓ Results saved (implement file output as needed)")
    
    return results
end

# Export the main function
export solve_energy_macro_model

end # module