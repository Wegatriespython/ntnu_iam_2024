# generate_cost_dataset.jl
# Generate a dataset of cost sequences from the energy model with varying demand levels

using JuMP, GAMS
using CSV, DataFrames
using Random
using Printf
using Statistics
using ArgParse

# Include necessary files
include("shared.jl")
include("energy/energy_model_world.jl")

"""
    vary_demand_levels(base_demand::Dict, variation_range::Float64, num_samples::Int) -> Vector{Dict}

Generate multiple demand scenarios by varying the base demand levels.
- base_demand: Original demand dictionary 
- variation_range: ±percentage variation (e.g., 0.3 for ±30%)
- num_samples: Number of demand scenarios to generate
"""
function vary_demand_levels(base_demand::Dict, variation_range::Float64, num_samples::Int)
    demand_scenarios = []
    
    for i in 1:num_samples
        # Create a copy of base demand
        scenario_demand = deepcopy(base_demand)
        
        # Apply random variations to each demand entry
        for key in keys(scenario_demand)
            # Generate a random multiplier within the specified range
            multiplier = 1.0 + (2 * rand() - 1) * variation_range
            # Ensure demand remains positive
            multiplier = max(0.1, multiplier)  # Minimum 10% of base demand
            scenario_demand[key] = scenario_demand[key] * multiplier
        end
        
        push!(demand_scenarios, scenario_demand)
    end
    
    return demand_scenarios
end

"""
    run_energy_with_demand(demand_scenario::Dict, scenario_id::Int) -> DataFrame

Run the energy model with a specific demand scenario and return results.
"""
function run_energy_with_demand(demand_scenario::Dict, scenario_id::Int)
    # Create model
    model = Model(GAMS.Optimizer)
    set_optimizer_attribute(model, "LP", "CPLEX")
    
    # Temporarily override the global demand dictionary
    global demand
    original_demand = deepcopy(demand)
    demand = demand_scenario
    
    try
        # Recalculate costs with potentially perturbed parameters
        cost_capacity, cost_activity, cost = calculate_costs()
        
        # Create the energy model with the modified demand
        create_energy_model!(model)
        
        # Solve the model
        optimize!(model)
        
        # Extract results if optimal
        if termination_status(model) ∈ (MOI.OPTIMAL, MOI.LOCALLY_SOLVED)
            results = DataFrame()
            
            # Get total cost
            total_cost = value(model[:TOTAL_COST])
            
            # Get annual costs for each year
            annual_costs = []
            for year in year_all
                annual_cost = value(model[:COST_ANNUAL][year])
                push!(annual_costs, annual_cost)
            end
            
            # Get demand values used
            elec_demand = demand_scenario[("electricity", "useful")]
            nele_demand = demand_scenario[("nonelectric", "useful")]
            
            # Get key technology activities for 2020 and 2050
            tech_activities_2020 = Dict()
            tech_activities_2050 = Dict()
            
            key_techs = ["coal_ppl", "gas_ppl", "wind_ppl", "solar_PV_ppl", "nuclear_ppl", 
                        "coal_nele", "gas_nele", "oil_nele"]
            
            for tech in key_techs
                if tech in technology
                    tech_activities_2020[tech] = value(model[:ACT][tech, 2020])
                    tech_activities_2050[tech] = value(model[:ACT][tech, 2050])
                end
            end
            
            # Get emissions
            emissions_2020 = value(model[:EMISS][2020])
            emissions_2050 = value(model[:EMISS][2050])
            cum_emissions = value(model[:CUM_EMISS])
            
            # Calculate time-varying demand projections using GDP scaling
            # GDP scaling factors from energy_parameters.yaml
            gdp_scaling = Dict(
                2020 => 1.0, 2030 => 1.2, 2040 => 1.4, 2050 => 1.7,
                2060 => 2.0, 2070 => 2.4, 2080 => 2.9
            )
            beta = 0.7  # elasticity parameter
            
            # Calculate projected demands for each year
            demand_projections = Dict()
            for year in [2020, 2030, 2040, 2050, 2060, 2070, 2080]
                scaling_factor = (gdp_scaling[year] / gdp_scaling[2020])^beta
                demand_projections["elec_demand_$year"] = elec_demand * scaling_factor
                demand_projections["nele_demand_$year"] = nele_demand * scaling_factor
                demand_projections["total_demand_$year"] = (elec_demand + nele_demand) * scaling_factor
            end
            
            # Create result row with time-varying demands
            result_row = DataFrame(
                scenario_id = scenario_id,
                total_cost = total_cost,
                elec_demand = elec_demand,
                nele_demand = nele_demand,
                cost_2020 = annual_costs[1],
                cost_2030 = annual_costs[2], 
                cost_2040 = annual_costs[3],
                cost_2050 = annual_costs[4],
                cost_2060 = annual_costs[5],
                cost_2070 = annual_costs[6],
                cost_2080 = annual_costs[7],
                emissions_2020 = emissions_2020,
                emissions_2050 = emissions_2050,
                cum_emissions = cum_emissions,
                coal_ppl_2020 = get(tech_activities_2020, "coal_ppl", 0.0),
                gas_ppl_2020 = get(tech_activities_2020, "gas_ppl", 0.0),
                wind_ppl_2020 = get(tech_activities_2020, "wind_ppl", 0.0),
                solar_PV_ppl_2020 = get(tech_activities_2020, "solar_PV_ppl", 0.0),
                nuclear_ppl_2020 = get(tech_activities_2020, "nuclear_ppl", 0.0),
                coal_ppl_2050 = get(tech_activities_2050, "coal_ppl", 0.0),
                gas_ppl_2050 = get(tech_activities_2050, "gas_ppl", 0.0),
                wind_ppl_2050 = get(tech_activities_2050, "wind_ppl", 0.0),
                solar_PV_ppl_2050 = get(tech_activities_2050, "solar_PV_ppl", 0.0),
                nuclear_ppl_2050 = get(tech_activities_2050, "nuclear_ppl", 0.0),
                coal_nele_2020 = get(tech_activities_2020, "coal_nele", 0.0),
                gas_nele_2020 = get(tech_activities_2020, "gas_nele", 0.0),
                oil_nele_2020 = get(tech_activities_2020, "oil_nele", 0.0),
                coal_nele_2050 = get(tech_activities_2050, "coal_nele", 0.0),
                gas_nele_2050 = get(tech_activities_2050, "gas_nele", 0.0),
                oil_nele_2050 = get(tech_activities_2050, "oil_nele", 0.0),
                solve_status = string(termination_status(model))
            )
            
            # Add demand projection columns
            for (col_name, value) in demand_projections
                result_row[!, col_name] = [value]
            end
            
            return result_row
            
        else
            # Return row with error status - still include demand projections
            elec_demand = demand_scenario[("electricity", "useful")]
            nele_demand = demand_scenario[("nonelectric", "useful")]
            
            # Calculate time-varying demand projections for error case too
            gdp_scaling = Dict(
                2020 => 1.0, 2030 => 1.2, 2040 => 1.4, 2050 => 1.7,
                2060 => 2.0, 2070 => 2.4, 2080 => 2.9
            )
            beta = 0.7
            
            error_row = DataFrame(
                scenario_id = scenario_id,
                total_cost = NaN,
                elec_demand = elec_demand,
                nele_demand = nele_demand,
                solve_status = string(termination_status(model))
            )
            
            # Add demand projection columns for error case
            for year in [2020, 2030, 2040, 2050, 2060, 2070, 2080]
                scaling_factor = (gdp_scaling[year] / gdp_scaling[2020])^beta
                error_row[!, "elec_demand_$year"] = [elec_demand * scaling_factor]
                error_row[!, "nele_demand_$year"] = [nele_demand * scaling_factor]
                error_row[!, "total_demand_$year"] = [(elec_demand + nele_demand) * scaling_factor]
            end
            
            return error_row
        end
        
    finally
        # Restore original demand
        demand = original_demand
    end
end

"""
    generate_cost_dataset(num_samples::Int=100, demand_variation::Float64=0.3, output_file::String="cost_dataset.csv")

Main function to generate the cost dataset.
"""
function generate_cost_dataset(num_samples::Int=100, demand_variation::Float64=0.3, output_file::String="cost_dataset.csv")
    println("=== Generating Energy Model Cost Dataset ===")
    println("Number of samples: $num_samples")
    println("Demand variation range: ±$(demand_variation*100)%")
    println("Output file: $output_file")
    
    # Set random seed for reproducibility
    Random.seed!(42)
    
    # Get base demand from the loaded parameters
    base_demand = deepcopy(demand)
    println("\nBase demand levels:")
    for (key, value) in base_demand
        println("  $key: $value PWh")
    end
    
    # Generate demand scenarios
    println("\nGenerating demand scenarios...")
    demand_scenarios = vary_demand_levels(base_demand, demand_variation, num_samples)
    
    # Initialize results DataFrame
    all_results = DataFrame()
    
    # Run model for each scenario
    println("\nRunning energy model for each scenario:")
    successful_runs = 0
    failed_runs = 0
    
    for (i, demand_scenario) in enumerate(demand_scenarios)
        print("Scenario $i/$num_samples... ")
        
        try
            result = run_energy_with_demand(demand_scenario, i)
            
            if nrow(result) > 0 && !isnan(result.total_cost[1])
                println("✅ Success (Total cost: $(round(result.total_cost[1], digits=1)) billion USD)")
                successful_runs += 1
            else
                println("❌ Failed to solve")
                failed_runs += 1
            end
            
            append!(all_results, result)
            
        catch e
            println("❌ Error: $e")
            failed_runs += 1
            
            # Add error row
            error_row = DataFrame(
                scenario_id = i,
                total_cost = NaN,
                elec_demand = demand_scenario[("electricity", "useful")],
                nele_demand = demand_scenario[("nonelectric", "useful")],
                solve_status = "ERROR: $e"
            )
            append!(all_results, error_row)
        end
    end
    
    # Save results
    println("\n=== Results Summary ===")
    println("Successful runs: $successful_runs")
    println("Failed runs: $failed_runs") 
    println("Success rate: $(round(successful_runs/num_samples*100, digits=1))%")
    
    if successful_runs > 0
        successful_results = filter(row -> !isnan(row.total_cost), all_results)
        println("\nCost statistics (successful runs only):")
        println("  Mean total cost: $(round(mean(successful_results.total_cost), digits=1)) billion USD")
        println("  Min total cost: $(round(minimum(successful_results.total_cost), digits=1)) billion USD")
        println("  Max total cost: $(round(maximum(successful_results.total_cost), digits=1)) billion USD")
        println("  Std dev: $(round(std(successful_results.total_cost), digits=1)) billion USD")
    end
    
    # Save to CSV
    CSV.write(output_file, all_results)
    println("\nDataset saved to: $output_file")
    println("Total records: $(nrow(all_results))")
    
    return all_results
end

# Command line execution
if abspath(PROGRAM_FILE) == @__FILE__
    function parse_commandline()
        s = ArgParseSettings()
        
        @add_arg_table s begin
            "--num-samples", "-n"
                help = "Number of demand scenarios to generate"
                arg_type = Int
                default = 100
            "--demand-variation", "-v"
                help = "Demand variation as percentage (e.g., 0.3 for ±30%)"
                arg_type = Float64
                default = 0.3
            "--output", "-o"
                help = "Output CSV filename"
                arg_type = String
                default = "cost_dataset.csv"
        end
        
        return parse_args(s)
    end
    
    # Parse arguments
    parsed_args = parse_commandline()
    num_samples = parsed_args["num-samples"]
    demand_variation = parsed_args["demand-variation"]
    output_file = parsed_args["output"]
    
    # Generate the dataset
    results = generate_cost_dataset(num_samples, demand_variation, output_file)
end