using YAML
using Random
using Statistics
using LinearAlgebra

# Set random seed for reproducibility
Random.seed!(42)

# Define perturbation ranges and parameter groups
const PERTURBATION_RANGES = Dict(
    "inv" => 0.30,          # ±30% for investment costs
    "vom" => 0.30,          # ±30% for variable O&M
    "hours" => 0.20,        # ±20% for capacity factors
    "CO2_emission" => 0.10, # ±10% for emissions
)

# Demand scenarios
const DEMAND_SCENARIOS = Dict(
    "low" => 0.5,
    "medium" => 0.7,
    "high" => 0.9
)

# Parameters to perturb (grouped for correlation)
const PARAMETER_GROUPS = [
    # Group 1: Renewable investment costs (correlated)
    [
        ("inv", "wind_ppl"),
        ("inv", "solar_PV_ppl"),
        ("inv", "hydro_ppl")
    ],
    # Group 2: Fossil investment costs (correlated)
    [
        ("inv", "coal_ppl"),
        ("inv", "gas_ppl"),
        ("inv", "oil_ppl")
    ],
    # Group 3: Fossil fuel extraction costs (correlated)
    [
        ("vom", "coal_extr"),
        ("vom", "gas_extr"),
        ("vom", "oil_extr")
    ],
    # Group 4: Renewable capacity factors (independent)
    [
        ("hours", "wind_ppl"),
        ("hours", "solar_PV_ppl")
    ],
    # Group 5: Fossil emissions (correlated)
    [
        ("CO2_emission", "coal_ppl"),
        ("CO2_emission", "gas_ppl"),
        ("CO2_emission", "oil_ppl")
    ]
]

# Latin Hypercube Sampling
function latin_hypercube_sample(n_samples, n_dims)
    samples = zeros(n_samples, n_dims)
    for j in 1:n_dims
        perm = randperm(n_samples)
        for i in 1:n_samples
            lower = (perm[i] - 1) / n_samples
            upper = perm[i] / n_samples
            samples[i, j] = lower + rand() * (upper - lower)
        end
    end
    return samples
end

# Transform uniform [0,1] samples to perturbation factors
function uniform_to_perturbation(u, range)
    # Convert uniform [0,1] to multiplicative factor
    # u=0 -> 1-range, u=0.5 -> 1.0, u=1 -> 1+range
    return 1.0 + range * (2.0 * u - 1.0)
end

# Apply economic sensibility constraints
function apply_constraints!(params)
    # Ensure gas costs > coal costs for extraction
    if haskey(params, "vom") && haskey(params["vom"], "gas_extr") && haskey(params["vom"], "coal_extr")
        if params["vom"]["gas_extr"] < params["vom"]["coal_extr"]
            # Swap them
            params["vom"]["gas_extr"], params["vom"]["coal_extr"] = 
                params["vom"]["coal_extr"], params["vom"]["gas_extr"]
        end
    end
    
    # Ensure oil costs > gas costs for extraction
    if haskey(params, "vom") && haskey(params["vom"], "oil_extr") && haskey(params["vom"], "gas_extr")
        if params["vom"]["oil_extr"] < params["vom"]["gas_extr"]
            params["vom"]["oil_extr"] = params["vom"]["gas_extr"] * 1.2
        end
    end
    
    # Ensure renewable capacity factors are reasonable
    if haskey(params, "hours")
        # Wind: 1500-3000 hours
        if haskey(params["hours"], "wind_ppl")
            params["hours"]["wind_ppl"] = clamp(params["hours"]["wind_ppl"], 1500, 3000)
        end
        # Solar: 800-1800 hours
        if haskey(params["hours"], "solar_PV_ppl")
            params["hours"]["solar_PV_ppl"] = clamp(params["hours"]["solar_PV_ppl"], 800, 1800)
        end
    end
    
    # Ensure emissions maintain relative ordering: coal > oil > gas
    if haskey(params, "CO2_emission")
        emissions = []
        if haskey(params["CO2_emission"], "coal_ppl")
            push!(emissions, ("coal_ppl", params["CO2_emission"]["coal_ppl"]))
        end
        if haskey(params["CO2_emission"], "oil_ppl")
            push!(emissions, ("oil_ppl", params["CO2_emission"]["oil_ppl"]))
        end
        if haskey(params["CO2_emission"], "gas_ppl")
            push!(emissions, ("gas_ppl", params["CO2_emission"]["gas_ppl"]))
        end
        
        if length(emissions) == 3
            # Sort by value
            sort!(emissions, by=x->x[2])
            # Reassign to maintain coal > oil > gas
            params["CO2_emission"]["gas_ppl"] = emissions[1][2]
            params["CO2_emission"]["oil_ppl"] = emissions[2][2]
            params["CO2_emission"]["coal_ppl"] = emissions[3][2]
        end
    end
end

# Generate perturbed parameters
function generate_perturbed_configs(base_params, n_samples=50)
    configs = []
    
    # Count total dimensions needed
    n_groups = length(PARAMETER_GROUPS)
    
    # Generate LHS samples
    lhs_samples = latin_hypercube_sample(n_samples, n_groups)
    
    for i in 1:n_samples
        for demand_scenario in keys(DEMAND_SCENARIOS)
            # Deep copy base parameters
            new_params = deepcopy(base_params)
            
            # Apply demand scenario
            new_params["model"]["beta"] = DEMAND_SCENARIOS[demand_scenario]
            
            # Apply perturbations for each parameter group
            for (group_idx, param_group) in enumerate(PARAMETER_GROUPS)
                # Get the perturbation factor for this group
                u = lhs_samples[i, group_idx]
                
                # Determine the range based on the first parameter type in the group
                param_type = param_group[1][1]
                range = PERTURBATION_RANGES[param_type]
                perturbation = uniform_to_perturbation(u, range)
                
                # Apply correlated perturbation to all parameters in the group
                for (ptype, pname) in param_group
                    if haskey(new_params, ptype) && haskey(new_params[ptype], pname)
                        value = new_params[ptype][pname]
                        if isa(value, Number)
                            new_params[ptype][pname] = value * perturbation
                        elseif isa(value, Dict)
                            # For time-varying parameters
                            for year in keys(value)
                                new_params[ptype][pname][year] = value[year] * perturbation
                            end
                        end
                    end
                end
            end
            
            # Apply economic sensibility constraints
            apply_constraints!(new_params)
            
            # Store configuration
            push!(configs, (
                params = new_params,
                sample_id = i,
                demand = demand_scenario,
                filename = "config_sample_$(i)_demand_$(demand_scenario).yaml"
            ))
        end
    end
    
    return configs
end

# Main execution
function main()
    # Load base parameters
    base_params = YAML.load_file("energy_parameters.yaml")
    
    # Generate perturbed configurations
    configs = generate_perturbed_configs(base_params, 50)
    
    # Save configurations
    for config in configs
        filepath = joinpath("config_samples", config.filename)
        YAML.write_file(filepath, config.params)
    end
    
    # Generate summary file
    summary = Dict(
        "n_samples" => 50,
        "demand_scenarios" => collect(keys(DEMAND_SCENARIOS)),
        "total_configs" => length(configs),
        "parameter_groups" => PARAMETER_GROUPS,
        "perturbation_ranges" => PERTURBATION_RANGES
    )
    
    YAML.write_file(joinpath("config_samples", "sampling_summary.yaml"), summary)
    
    println("Generated $(length(configs)) configuration files in config_samples/")
    println("Demand scenarios: ", join(keys(DEMAND_SCENARIOS), ", "))
    println("Sample IDs: 1-50")
end

# Run if called directly
if abspath(PROGRAM_FILE) == @__FILE__
    main()
end