# energy_model_reduced.jl - Reduced form version of energy_model_world.jl
# Simplified energy system optimization with aggregated technology groups

using JuMP
using YAML
using Statistics

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
function aggregate_parameter(param_dict, group_techs, aggregation_method=:mean)
    values = []
    for tech in group_techs
        if haskey(param_dict, tech)
            push!(values, param_dict[tech])
        end
    end
    
    if isempty(values)
        return 0.0
    end
    
    if aggregation_method == :mean
        return mean(values)
    elseif aggregation_method == :sum
        return sum(values)
    elseif aggregation_method == :max
        return maximum(values)
    else
        return mean(values)
    end
end

# Function to aggregate nested parameters (like input/output coefficients)
function aggregate_nested_parameter(param_dict, group_techs, aggregation_method=:mean)
    result = Dict()
    
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
                            push!(result[key], value)
                        end
                    end
                end
            end
        end
    end
    
    # Aggregate collected values
    aggregated = Dict()
    for (key, values) in result
        if aggregation_method == :mean
            aggregated[key] = mean(values)
        elseif aggregation_method == :sum
            aggregated[key] = sum(values)
        else
            aggregated[key] = mean(values)
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

    for (group, techs) in TECH_GROUPS
        group_params[group] = Dict()
        
        # Aggregate diffusion rates
        diffusion_rates = []
        for tech in techs
            if haskey(diffusion_up, tech)
                push!(diffusion_rates, diffusion_up[tech])
            end
        end
        group_params[group]["diffusion"] = isempty(diffusion_rates) ? 0.1 : mean(diffusion_rates)
        
        # Aggregate startup capacities
        startup_caps = []
        for tech in techs
            if haskey(startup, tech)
                push!(startup_caps, startup[tech])
            end
        end
        group_params[group]["startup"] = isempty(startup_caps) ? 0.0 : sum(startup_caps)
        
        # Aggregate input coefficients
        if haskey(params, "input")
            group_params[group]["input_coeffs"] = aggregate_nested_parameter(params["input"], techs, :mean)
        else
            group_params[group]["input_coeffs"] = Dict()
        end
        
        # Aggregate output coefficients
        if haskey(params, "output")
            group_params[group]["output_coeffs"] = aggregate_nested_parameter(params["output"], techs, :mean)
        else
            group_params[group]["output_coeffs"] = Dict()
        end
        
        # Aggregate costs
        if haskey(params, "vom")
            group_params[group]["vom_cost"] = aggregate_parameter(params["vom"], techs, :mean)
        else
            group_params[group]["vom_cost"] = 50.0  # Default cost
        end
        
        # Aggregate capacity factors from hours
        if haskey(params, "hours")
            avg_hours = aggregate_parameter(params["hours"], techs, :mean)
            # Hours are in thousands, so multiply by 1000 before dividing by 8760
            group_params[group]["capacity_factor"] = (avg_hours * 1000.0) / 8760.0
        else
            group_params[group]["capacity_factor"] = 0.8  # Default
        end
        
        # Aggregate emissions
        if haskey(params, "CO2_emission")
            group_params[group]["emissions"] = aggregate_parameter(params["CO2_emission"], techs, :mean)
        else
            group_params[group]["emissions"] = 0.0
        end
        
        # Aggregate lifetime
        lifetime_val = aggregate_parameter(params["lifetime"], techs, :mean)
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
            for tech in techs
                if haskey(params, "inv") && haskey(params["inv"], tech) && haskey(params["inv"][tech], string(y))
                    push!(inv_costs, params["inv"][tech][string(y)])
                end
                if haskey(params, "fom") && haskey(params["fom"], tech) && haskey(params["fom"][tech], string(y))
                    push!(fom_costs, params["fom"][tech][string(y)])
                end
            end
            avg_inv_cost = isempty(inv_costs) ? 0.0 : mean(inv_costs)
            avg_fom_cost = isempty(fom_costs) ? 0.0 : mean(fom_costs)
            group_params[group]["cost_capacity"][y] = avg_inv_cost * annuity_factor + avg_fom_cost
        end

        # Aggregate base year calibration
        if haskey(params, "energy_calibration")
            base_activities = []
            for tech in techs
                if haskey(params["energy_calibration"], tech)
                    tech_years = params["energy_calibration"][tech]
                    if haskey(tech_years, "2020")
                        push!(base_activities, tech_years["2020"])
                    end
                end
            end
            group_params[group]["base_activity"] = isempty(base_activities) ? 0.0 : sum(base_activities)
        else
            group_params[group]["base_activity"] = 0.0
        end
    end
    
    return group_params
end

# Compute group parameters
group_params = compute_group_parameters(params)

# Model construction for reduced form
function create_reduced_energy_model!(model)
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
        @objective(model, Min, TOTAL_COST_R)
    end
    
    return model
end

# Export the model creation function
export create_reduced_energy_model!
