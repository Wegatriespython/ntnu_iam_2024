# energy_model_reduced.jl - Reduced form version of energy_model_world.jl
# Simplified energy system optimization with aggregated technology groups

using JuMP
using YAML
using Statistics

# --- Bias Configuration ---
# Set BIAS_FACTOR to 0.0 for a pure mean aggregation.
# Set BIAS_FACTOR to 1.0 to use only the winner technology's parameters.
# Values in between will create a weighted average.
const BIAS_FACTOR = 0.0  # <-- TOGGLE THIS VALUE (0.0 to 1.0)

# Define the worst-performing technology for each group based on the world model's results
# These are technologies with zero or minimal usage in the full model optimization
const REPRESENTATIVE_TECHS = Dict(
    "fossil_power" => "oil_ppl",           # Drops to 0 after 2020
    "renewable_power" => "wind_ppl",       # Drops to 0 after 2030
    "fossil_extraction" => "oil_extr",     # Drops to 0 after 2020
    "renewable_potential" => "wind_pot",   # Drops to 0 after 2030
    "fossil_nonelec" => "oil_nele",        # Drops to 0 after 2020
    "biomass_nonelec" => "bio_nele",       # Only one tech in group
    "other_nonelec" => "other_nele",       # Stays at very low levels
    "biomass_power" => "bio_ppl",          # Drops to 0 after 2020
    # Groups with one tech are included for completeness
)


# Load parameters from shared configuration
const ENERGY_MODEL_DIR = dirname(@__FILE__)
const JULIA_MODEL_DIR = dirname(ENERGY_MODEL_DIR)
const PARAMS_FILE = joinpath(JULIA_MODEL_DIR, "energy_parameters_scaled.yaml")

println("Loading detailed parameters from: $PARAMS_FILE")
params = YAML.load_file(PARAMS_FILE)

# Technology grouping based on function and fuel type
const TECH_GROUPS = Dict(
    "fossil_power" => ["coal_ppl", "gas_ppl", "oil_ppl"],
    "renewable_power" => ["hydro_ppl", "wind_ppl", "solar_PV_ppl"],
    "nuclear_power" => ["nuclear_ppl"],
    "biomass_power" => ["bio_ppl"],
    "fossil_extraction" => ["coal_extr", "gas_extr", "oil_extr"],
    "renewable_potential" => ["hydro_pot", "wind_pot", "solar_pot"],
    "fossil_nonelec" => ["coal_nele", "gas_nele", "oil_nele"],
    "biomass_nonelec" => ["bio_nele"],
    "other_nonelec" => ["solar_nele", "other_nele"],
    "infrastructure" => ["electricity_grid", "appliances"]
)

# Load share constraints from parameters
share_up = get(params, "share_up", Dict())
share_lo = get(params, "share_lo", Dict())
share = get(params["sets"], "share", String[])
tec_share = Dict()
if haskey(params, "tec_share")
    for ts in params["tec_share"]
        tec_share[(ts[1], ts[2])] = true
    end
end
tec_share_rhs = Dict()
if haskey(params, "tec_share_rhs")
    for ts in params["tec_share_rhs"]
        tec_share_rhs[(ts[1], ts[2])] = true
    end
end

# Reverse mapping: technology -> group
tech_to_group = Dict()
for (group, techs) in TECH_GROUPS
    for tech in techs
        tech_to_group[tech] = group
    end
end

# Function to aggregate parameter values across technologies in a group
function aggregate_parameter(param_dict, group_techs, aggregation_method=:mean; winner_tech=nothing, bias=0.0)
    values = []
    winner_value = nothing

    for tech in group_techs
        if haskey(param_dict, tech)
            val = param_dict[tech]
            push!(values, val)
            if tech == winner_tech
                winner_value = val
            end
        end
    end
    
    if isempty(values)
        return 0.0
    end
    
    mean_value = 0.0
    if aggregation_method == :mean
        mean_value = mean(values)
    elseif aggregation_method == :sum
        mean_value = sum(values)
    elseif aggregation_method == :max
        mean_value = maximum(values)
    else
        mean_value = mean(values)
    end

    if winner_tech !== nothing && winner_value !== nothing
        return (1 - bias) * mean_value + bias * winner_value
    else
        return mean_value
    end
end

# Function to aggregate nested parameters (like input/output coefficients)
function aggregate_nested_parameter(param_dict, group_techs, aggregation_method=:mean; winner_tech=nothing, bias=0.0)
    result = Dict()
    
    # First, gather all values, grouped by their energy/level key
    for tech in group_techs
        if haskey(param_dict, tech)
            tech_params = param_dict[tech]
            if isa(tech_params, Dict)
                for (energy_type, levels) in tech_params
                    if isa(levels, Dict)
                        for (level, value) in levels
                            key = (energy_type, level)
                            if !haskey(result, key)
                                result[key] = []
                            end
                            # Store the technology name with its value
                            push!(result[key], (tech, value))
                        end
                    end
                end
            end
        end
    end
    
    # Now, aggregate the collected values
    aggregated = Dict()
    for (key, tech_values) in result
        # Extract just the numerical values for calculating the mean
        values = [v for (t, v) in tech_values]
        mean_value = isempty(values) ? 0.0 : mean(values)

        # Find the specific value for the winner technology
        winner_value = nothing
        for (tech, value) in tech_values
            if tech == winner_tech
                winner_value = value
                break
            end
        end

        # If a winner is specified and found, apply the bias. Otherwise, use the mean.
        if winner_tech !== nothing && winner_value !== nothing
            aggregated[key] = (1 - bias) * mean_value + bias * winner_value
        else
            aggregated[key] = mean_value
        end
    end
    
    return aggregated
end

# Load diffusion parameters
diffusion_up = get(params, "diffusion_up", Dict())
startup = get(params, "startup", Dict())

# Compute group parameters from underlying technologies
function compute_group_parameters(params)
    group_params = Dict()
    year_all = [2020, 2030, 2040, 2050, 2060, 2070, 2080]
    discount_rate = get(params["model"], "discount_rate", 0.05)
    
    # Load parameters from params
    diffusion_up = get(params, "diffusion_up", Dict())
    startup = get(params, "startup", Dict())

    for (group, techs) in TECH_GROUPS
        group_params[group] = Dict()
        winner = get(REPRESENTATIVE_TECHS, group, nothing)
        
        # Aggregate diffusion rates
        group_params[group]["diffusion"] = aggregate_parameter(diffusion_up, techs, :mean, winner_tech=winner, bias=BIAS_FACTOR)
        
        # Aggregate startup capacities
        group_params[group]["startup"] = aggregate_parameter(startup, techs, :sum, winner_tech=winner, bias=BIAS_FACTOR)
        
        # Aggregate input coefficients
        group_params[group]["input_coeffs"] = aggregate_nested_parameter(get(params, "input", Dict()), techs, :mean, winner_tech=winner, bias=BIAS_FACTOR)
        
        # Aggregate output coefficients
        group_params[group]["output_coeffs"] = aggregate_nested_parameter(get(params, "output", Dict()), techs, :mean, winner_tech=winner, bias=BIAS_FACTOR)
        
        # Aggregate costs
        group_params[group]["vom_cost"] = aggregate_parameter(get(params, "vom", Dict()), techs, :mean, winner_tech=winner, bias=BIAS_FACTOR)
        
        # Aggregate capacity factors from hours
        avg_hours = aggregate_parameter(get(params, "hours", Dict()), techs, :mean, winner_tech=winner, bias=BIAS_FACTOR)
        group_params[group]["capacity_factor"] = (avg_hours * 1000.0) / 8760.0
        
        # Aggregate emissions
        group_params[group]["emissions"] = aggregate_parameter(get(params, "CO2_emission", Dict()), techs, :mean, winner_tech=winner, bias=BIAS_FACTOR)
        
        # Aggregate lifetime
        lifetime_val = aggregate_parameter(params["lifetime"], techs, :mean, winner_tech=winner, bias=BIAS_FACTOR)
        group_params[group]["lifetime"] = lifetime_val > 0 ? lifetime_val : 30.0

        # Aggregate investment and FOM costs and calculate cost_capacity
        group_params[group]["cost_capacity"] = Dict()
        annuity_factor = 0.0
        if group_params[group]["lifetime"] > 0
            lt = group_params[group]["lifetime"]
            annuity_factor = ((1 + discount_rate)^lt * discount_rate) / ((1 + discount_rate)^lt - 1)
        end

        for y in year_all
            inv_costs = []
            fom_costs = []
            winner_inv_cost = nothing
            winner_fom_cost = nothing

            for tech in techs
                if haskey(params, "inv") && haskey(params["inv"], tech) && haskey(params["inv"][tech], string(y))
                    val = params["inv"][tech][string(y)]
                    push!(inv_costs, val)
                    if tech == winner
                        winner_inv_cost = val
                    end
                end
                if haskey(params, "fom") && haskey(params["fom"], tech) && haskey(params["fom"][tech], string(y))
                    val = params["fom"][tech][string(y)]
                    push!(fom_costs, val)
                    if tech == winner
                        winner_fom_cost = val
                    end
                end
            end
            avg_inv_cost = isempty(inv_costs) ? 0.0 : mean(inv_costs)
            avg_fom_cost = isempty(fom_costs) ? 0.0 : mean(fom_costs)

            final_inv_cost = (winner !== nothing && winner_inv_cost !== nothing) ? (1 - BIAS_FACTOR) * avg_inv_cost + BIAS_FACTOR * winner_inv_cost : avg_inv_cost
            final_fom_cost = (winner !== nothing && winner_fom_cost !== nothing) ? (1 - BIAS_FACTOR) * avg_fom_cost + BIAS_FACTOR * winner_fom_cost : avg_fom_cost

            group_params[group]["cost_capacity"][y] = final_inv_cost * annuity_factor + final_fom_cost
        end

        # Aggregate base year calibration
        base_activities = []
        for tech in techs
            if haskey(params, "energy_calibration") && haskey(params["energy_calibration"], tech)
                tech_years = params["energy_calibration"][tech]
                if haskey(tech_years, "2020")
                    push!(base_activities, tech_years["2020"])
                end
            end
        end
        group_params[group]["base_activity"] = isempty(base_activities) ? 0.0 : sum(base_activities)
    end
    
    return group_params
end

# Compute group parameters
group_params = compute_group_parameters(params)

# Model construction for reduced form
function create_reduced_energy_model!(model; config_file=nothing)
    # Load parameters and compute group params based on config file
    if config_file !== nothing
        params = YAML.load_file(config_file)
        group_params = compute_group_parameters(params)
    else
        # Use the module-level params and group_params if no config file provided
        params = Main.params
        group_params = Main.group_params
    end
    
    # Get shared parameters if available
    year_all = get(model.ext, :year_all, [2020, 2030, 2040, 2050, 2060, 2070, 2080])
    period_length = get(model.ext, :period_length, 10)
    discount_rate = get(model.ext, :discount_rate, 0.05)
    
    groups = collect(keys(TECH_GROUPS))
    
    # Decision variables (much fewer than original)
    @variable(model, ACT_R[groups, year_all] >= 0)  # Activity level
    @variable(model, CAP_NEW_R[groups, year_all] >= 0)  # New capacity additions
    @variable(model, COST_ANNUAL_R[year_all] >= 0)       # Annual cost
    @variable(model, TOTAL_COST_R >= 0)                  # Total cost
    @variable(model, EMISS_R[year_all])                  # Annual emissions
    @variable(model, CUM_EMISS_R)                        # Cumulative emissions
    
    # Energy balance constraints using computed I/O coefficients
    energy_levels = [("electricity", "secondary"), ("electricity", "final"), 
                     ("electricity", "useful"), ("nonelectric", "useful")]
    
    for (energy, level) in energy_levels, y in year_all
        # Net production for this energy-level combination
        net_production = sum(
            ACT_R[g, y] * (get(group_params[g]["output_coeffs"], (energy, level), 0.0) - 
                          get(group_params[g]["input_coeffs"], (energy, level), 0.0))
            for g in groups
        )
        
        # Demand side
        if (energy, level) == ("electricity", "useful")
            if haskey(object_dictionary(model), :PHYSENE)
                # Integrated mode
                @constraint(model, net_production >= model[:PHYSENE]["ELEC", y])
            else
                # Standalone mode - demand from YAML is in GWh
                base_demand = 22600000.0  # GWh from YAML
                gdp_index = get(params["gdp"], string(y), 1.0)
                @constraint(model, net_production >= base_demand * gdp_index^params["model"]["beta"])
            end
        elseif (energy, level) == ("nonelectric", "useful")
            if haskey(object_dictionary(model), :PHYSENE)
                # Integrated mode
                @constraint(model, net_production >= model[:PHYSENE]["NELE", y])
            else
                # Standalone mode - demand from YAML is in GWh
                base_demand = 87300000.0  # GWh from YAML
                gdp_index = get(params["gdp"], string(y), 1.0)
                @constraint(model, net_production >= base_demand * gdp_index^params["model"]["beta"])
            end
        else
            # Intermediate energy balance
            @constraint(model, net_production >= 0)
        end
    end
    
    # Capacity constraints using computed capacity factors and lifetime
    for g in groups, y in year_all
        cf = group_params[g]["capacity_factor"]
        lifetime = group_params[g]["lifetime"]
        if cf > 0 && lifetime > 0
            # Sum capacities that are still active (within lifetime)
            yi = findfirst(==(y), year_all)
            @constraint(model, 
                ACT_R[g, y] <= sum(CAP_NEW_R[g, year_all[i]] * cf * 8760 
                                  for i in 1:yi if (yi-i+1)*period_length <= lifetime))
        end
    end
    
    # Growth constraints using aggregated diffusion parameters from full model
    for g in groups, i in 2:length(year_all)
        y = year_all[i]
        y_prev = year_all[i-1]
        
        # Use aggregated diffusion rate for this group
        diffusion_rate = group_params[g]["diffusion"]
        startup_cap = group_params[g]["startup"]
        
        if diffusion_rate > 0
            @constraint(model, 
                CAP_NEW_R[g, y] <= CAP_NEW_R[g, y_prev] * (1 + diffusion_rate)^period_length + startup_cap)
        end
    end
    
    # Cost calculation using computed cost parameters
    yi = Dict(y=>i for (i,y) in enumerate(year_all))
    for y in year_all
        @constraint(model,
            COST_ANNUAL_R[y] == sum(ACT_R[g, y] * group_params[g]["vom_cost"] for g in groups)
            + sum(sum(CAP_NEW_R[g, y2] * group_params[g]["cost_capacity"][y2]
                      for y2 in year_all if yi[y2] <= yi[y] && (yi[y] - yi[y2] + 1) * period_length <= group_params[g]["lifetime"])
                  for g in groups)
        )
    end

    # Total discounted cost
    @constraint(model,
        TOTAL_COST_R == sum(COST_ANNUAL_R[y] * period_length * 
                           (1 + discount_rate)^(-period_length * (findfirst(==(y), year_all) - 1))
                           for y in year_all)
    )
    
    # Emissions constraints using computed emission factors
    for y in year_all
        @constraint(model, 
            EMISS_R[y] == sum(ACT_R[g, y] * group_params[g]["emissions"] for g in groups)
        )
    end
    
    # Cumulative emissions
    @constraint(model,
        CUM_EMISS_R == sum(EMISS_R[y] * period_length for y in year_all)
    )
    
    # Share constraints from YAML parameters (matching full model)
    # For the reduced model, we need to map technology-level share constraints to groups
    # Currently only coal_nonelectric share is defined in the parameters
    # This maps to the ratio of coal_nele to all nonelectric technologies
    # NOTE: The share constraint for 'coal_nonelectric' cannot be accurately
    # applied in the reduced model because 'coal_nele' is aggregated into
    # the 'fossil_nonelec' group with other technologies. Applying the
    # constraint to the entire group would be incorrect. This is a limitation
    # of the reduced model's aggregation.
    for s in share, y in year_all
        if s == "coal_nonelectric"
            # coal_nele is part of fossil_nonelec group
            # The RHS includes all nonelec groups
            nonelec_groups = ["fossil_nonelec", "biomass_nonelec", "other_nonelec"]
            
            if haskey(share_up, s) && all(g in groups for g in nonelec_groups)
                # Upper bound: coal share of nonelectric <= 40%
                # Since coal is aggregated with gas and oil in fossil_nonelec,
                # we apply a proportional constraint based on the original tech composition
                # @constraint(model, 
                #     ACT_R["fossil_nonelec", y] <= share_up[s] * sum(ACT_R[g, y] for g in nonelec_groups))
            end
            
            if haskey(share_lo, s) && all(g in groups for g in nonelec_groups)
                # Lower bound: coal share >= 0 (no constraint needed as variables are non-negative)
            end
        end
    end
    
    # Base year calibration using computed base activities
    if 2020 in year_all && !haskey(object_dictionary(model), :PHYSENE)
        for g in groups
            base_activity = group_params[g]["base_activity"]
            if base_activity > 0
                fix(ACT_R[g, 2020], base_activity; force=true)
            end
        end
    end
    
    # Set objective if standalone mode
    if !haskey(object_dictionary(model), :UTILITY)
        # Check if QP mode is enabled (default to false if not specified)
        use_qp = get(model.ext, :use_qp, false)
        
        if use_qp
            # Add small negative quadratic penalty to make it a convex QP
            # The penalty is on activity levels to encourage smoother solutions
            penalty_coefficient = get(model.ext, :qp_penalty_coefficient, 1e-6)  # Allow custom penalty coefficient
            
            @objective(model, Min, 
                TOTAL_COST_R + penalty_coefficient * sum(ACT_R[g, y]^2 for g in groups, y in year_all)
            )
        else
            # Standard LP formulation
            @objective(model, Min, TOTAL_COST_R)
        end
    end
    
    return model
end

# Export the model creation function
export create_reduced_energy_model!