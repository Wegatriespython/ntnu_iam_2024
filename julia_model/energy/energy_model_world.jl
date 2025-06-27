# energy_model_world.jl - Equivalent to energy/energy_model_world.gms
# Energy system optimization model

using JuMP
using YAML

# Load parameters from YAML file
# Get the directory of the current file
const ENERGY_MODEL_DIR = dirname(@__FILE__)
const JULIA_MODEL_DIR = dirname(ENERGY_MODEL_DIR)
const PARAMS_FILE = joinpath(JULIA_MODEL_DIR, "energy_parameters.yaml")

# Load parameters from YAML
println("Loading parameters from: $PARAMS_FILE")
params = YAML.load_file(PARAMS_FILE)

# Set definitions
# Definition of indices (sets in GAMS) for technologies, energy carriers, energy levels and periods
technology = params["sets"]["technology"]
energy = params["sets"]["energy"]
level = params["sets"]["level"]

# Energy carrier and level combinations
energy_level = Dict()
for el in params["energy_level"]
    energy_level[(el[1], el[2])] = true
end

# Technology share mappings
share = params["sets"]["share"]
tec_share = Dict()
for ts in params["tec_share"]
    tec_share[(ts[1], ts[2])] = true
end
tec_share_rhs = Dict()
for ts in params["tec_share_rhs"]
    tec_share_rhs[(ts[1], ts[2])] = true
end

# Energy-sector mapping
map_energy_sector = Dict()
for mes in params["map_energy_sector"]
    map_energy_sector[(mes[1], mes[2], mes[3])] = true
end

# Parameter definitions
period_length = params["model"]["period_length"]

# Input coefficients
input = Dict()
for (tech, inputs) in params["input"]
    for (energy_type, levels) in inputs
        for (lvl, val) in levels
            input[(tech, energy_type, lvl)] = val
        end
    end
end

# Output coefficients
output = Dict()
for (tech, outputs) in params["output"]
    for (energy_type, levels) in outputs
        for (lvl, val) in levels
            output[(tech, energy_type, lvl)] = val
        end
    end
end

# CO2 emission coefficients [tCO2/MWh]
CO2_emission = params["CO2_emission"]
# Maximum annual technology capacity growth rate
diffusion_up = params["diffusion_up"]
# Maximum technology capacity constant addition per period
startup = params["startup"]

# Demand in base year [PWh]
demand = Dict()
for (energy_type, levels) in params["demand"]
    for (lvl, val) in levels
        demand[(energy_type, lvl)] = val
    end
end

# GDP [index]
gdp = Dict()
for (yr_str, val) in params["gdp"]
    gdp[parse(Int, yr_str)] = val
end

discount_rate = params["model"]["discount_rate"]
beta = params["model"]["beta"]

# Share constraints
share_up = params["share_up"]
share_lo = params["share_lo"]

# Investment cost [$/kW]
inv_cost = Dict()
if haskey(params, "inv")
    for (tech, yrs) in params["inv"]
        for (yr_str, val) in yrs
            inv_cost[(tech, parse(Int, yr_str))] = val
        end
    end
end

# Fixed O&M cost [$/(kW/yr)]
fom = Dict()
if haskey(params, "fom")
    for (tech, yrs) in params["fom"]
        for (yr_str, val) in yrs
            fom[(tech, parse(Int, yr_str))] = val
        end
    end
end

# Variable cost [$/MWh]
vom = Dict()
# init zeros
for tech in technology, yr in year_all
    vom[(tech, yr)] = 0.0
end
# set from params
if haskey(params, "vom")
    for (tech, val) in params["vom"]
        for yr in year_all
            vom[(tech, yr)] = val
        end
    end
end

# Technical lifetime
lifetime = params["lifetime"]
# Full load hours
hours = params["hours"]

# Cost calculation
function calculate_costs()
    cost_capacity = Dict{Tuple{String,Int},Float64}()
    cost_activity = Dict{Tuple{String,Int},Float64}()
    cost = Dict{Tuple{String,Int},Float64}()
    annuity_factors = Dict{String,Float64}()
    for tech in technology
        if haskey(lifetime, tech) && lifetime[tech] > 0
            lt = lifetime[tech]
            annuity_factors[tech] = ((1+discount_rate)^lt * discount_rate)/((1+discount_rate)^lt -1)
        end
    end
    valid_techs = [t for t in technology if haskey(lifetime,t) && lifetime[t]>0 && haskey(hours,t) && hours[t]>0]
    for tech in valid_techs, yr in year_all
        if haskey(inv_cost,(tech,yr)) && haskey(fom,(tech,yr))
            cost_capacity[(tech,yr)] = (inv_cost[(tech,yr)]*annuity_factors[tech] + fom[(tech,yr)])*1000
        end
        if haskey(vom,(tech,yr))
            cost_activity[(tech,yr)] = vom[(tech,yr)]
        end
        if haskey(cost_capacity,(tech,yr)) && haskey(cost_activity,(tech,yr))
            cost[(tech,yr)] = cost_capacity[(tech,yr)]/hours[tech] + cost_activity[(tech,yr)]
        end
    end
    return cost_capacity, cost_activity, cost
end
cost_capacity, cost_activity, cost = calculate_costs()

# Model construction
function create_energy_model!(model)
    @variable(model, ACT[technology, year_all] >= 0)
    @variable(model, CAP_NEW[technology, year_all] >= 0)
    @variable(model, EMISS[year_all])
    @variable(model, CUM_EMISS)
    @variable(model, TOTAL_COST)
    if !haskey(object_dictionary(model), :COST_ANNUAL)
        @variable(model, COST_ANNUAL[year_all])
    end
    COST_ANNUAL = model[:COST_ANNUAL]

    # Energy balance
    # Check if we're in integrated mode (PHYSENE exists) or standalone mode
    if haskey(object_dictionary(model), :PHYSENE)
        # Integrated mode: use PHYSENE from macro model
        for e in energy, l in level, y in year_all
            if haskey(energy_level, (e,l))
                # Find sectors that map to this energy-level combination
                sectors_for_el = [s for ((e2,l2,s),_) in map_energy_sector if e2==e && l2==l]
                if !isempty(sectors_for_el)
                    @constraint(model,
                        sum(ACT[t,y]*(get(output,(t,e,l),0)-get(input,(t,e,l),0)) for t in technology)
                        - sum(model[:PHYSENE][s,y] for s in sectors_for_el) >= 0)
                else
                    # No sector mapping, use standalone version
                    @constraint(model,
                        sum(ACT[t,y]*(get(output,(t,e,l),0)-get(input,(t,e,l),0)) for t in technology)
                        - get(demand,(e,l),0)*(gdp[y]/gdp[2020])^beta >= 0)
                end
            end
        end
    else
        # Standalone mode: use fixed demand with GDP scaling
        for e in energy, l in level, y in year_all
            if haskey(energy_level, (e,l))
                @constraint(model,
                    sum(ACT[t,y]*(get(output,(t,e,l),0)-get(input,(t,e,l),0)) for t in technology)
                    - get(demand,(e,l),0)*(gdp[y]/gdp[2020])^beta >= 0)
            end
        end
    end

    # Capacity balance
    for t in technology, y in year_all
        if haskey(hours,t) && haskey(lifetime,t)
            yi = findfirst(==(y),year_all)
            @constraint(model,
                ACT[t,y] <= sum(CAP_NEW[t,year_all[i]]*hours[t] for i in 1:yi if (yi-i+1)*period_length<=lifetime[t]))
        end
    end

    # Emissions
    @constraint(model, [y in year_all], sum(ACT[t,y]*get(CO2_emission,t,0) for t in technology)==EMISS[y])
    @constraint(model, sum(EMISS[y]*period_length for y in year_all)==CUM_EMISS)

    # Diffusion
    for t in technology, y in year_all[2:end]
        if haskey(diffusion_up,t)
            y0 = year_all[findfirst(==(y),year_all)-1]
            @constraint(model,
                CAP_NEW[t,y] <= CAP_NEW[t,y0]*(1+diffusion_up[t])^period_length + get(startup,t,0))
        end
    end

    # Share constraints
    for s in share, y in year_all
        lhs = [t for t in technology if haskey(tec_share,(s,t))]
        rhs = [t for t in technology if haskey(tec_share_rhs,(s,t))]
        if haskey(share_up,s)
            @constraint(model, sum(ACT[t,y] for t in lhs) <= share_up[s]*sum(ACT[t,y] for t in rhs))
        end
        if haskey(share_lo,s)
            @constraint(model, sum(ACT[t,y] for t in lhs) >= share_lo[s]*sum(ACT[t,y] for t in rhs))
        end
    end

    # Cost equations
    yi = Dict(y=>i for (i,y) in enumerate(year_all))
    for y in year_all
        @constraint(model,
            sum(ACT[t,y]*get(vom,(t,y),0) for t in technology)
            + sum(sum(CAP_NEW[t,y2]*get(cost_capacity,(t,y2),0) for y2 in year_all
                      if yi[y2]<=yi[y] && (yi[y]-yi[y2]+1)*period_length<=get(lifetime,t,0) && haskey(lifetime,t))
                  for t in technology)
            == COST_ANNUAL[y])
    end
    @constraint(model,
        sum(COST_ANNUAL[y]*period_length*(1+discount_rate)^(-period_length*(yi[y]-1)) for y in year_all)
        == TOTAL_COST)

    # Historical calibration - only apply if not in integrated mode
    # (integrated mode applies these constraints separately)
    if !haskey(object_dictionary(model), :PHYSENE)
        fix(ACT["coal_ppl",2020],9.462;force=true)
        fix(ACT["oil_ppl",2020],0.7;force=true)
        fix(ACT["solar_PV_ppl",2020],0.839;force=true)
        fix(ACT["gas_ppl",2020],6.36;force=true)
        fix(ACT["nuclear_ppl",2020],2.68;force=true)
        fix(ACT["hydro_ppl",2020],4.36;force=true)
        fix(ACT["wind_ppl",2020],1.6;force=true)
        fix(ACT["bio_ppl",2020],0.69;force=true)
        set_lower_bound(ACT["other_ppl",2020],0.127)
        fix(ACT["coal_nele",2020],10.7;force=true)
        fix(ACT["oil_nele",2020],43.0;force=true)
        fix(ACT["gas_nele",2020],18.7;force=true)
        fix(ACT["bio_nele",2020],10.6;force=true)
        set_lower_bound(ACT["other_nele",2020],0.28)
    end

    # Objective
    # Only set objective if not in integrated mode
    if !haskey(object_dictionary(model), :UTILITY)
        @objective(model, Min, TOTAL_COST)
    end

    return model
end

# Model calibration and resource constraints
energy_calibration = Dict()
if haskey(params, "energy_calibration")
    for (tech, yrs) in params["energy_calibration"]
        for (yr_str, val) in yrs
            energy_calibration[(tech, parse(Int, yr_str))] = val
        end
    end
end
energy_calibration_lo = Dict()
if haskey(params, "energy_calibration_lo")
    for (tech, yrs) in params["energy_calibration_lo"]
        for (yr_str, val) in yrs
            energy_calibration_lo[(tech, parse(Int, yr_str))] = val
        end
    end
end

# end of script

