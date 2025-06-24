# energy_model_world.jl - Equivalent to energy/energy_model_world.gms
# Energy system optimization model

using JuMP

# ------------------------------------------------------------------------------
# Set definitions
# ------------------------------------------------------------------------------

# Definition of indices (sets in GAMS) for technologies, energy carriers, energy levels and periods

technology = [
    "coal_extr",         # coal extraction
    "gas_extr",          # natural gas extraction  
    "oil_extr",          # crude oil extraction
    "nuclear_fuel",      # nuclear fuel for light water reactors
    "bio_pot",           # bioenergy potential
    "hydro_pot",         # hydropower potential
    "wind_pot",          # wind potential
    "solar_pot",         # solar energy potential
    "coal_ppl",          # coal power plant
    "gas_ppl",           # natural gas combined cycle power plant
    "oil_ppl",           # fuel oil power plant
    "bio_ppl",           # biomass power plant
    "hydro_ppl",         # hydroelectric power plant
    "wind_ppl",          # wind turbine
    "solar_PV_ppl",      # solar photovoltaics power plant
    "nuclear_ppl",       # nuclear power plant
    "other_ppl",         # other power plants
    "coal_nele",         # non-electric coal
    "oil_nele",          # non-electric oil
    "gas_nele",          # non-electric gas
    "bio_nele",          # non-electric biomass
    "solar_nele",        # non-electric solar
    "other_nele",        # non-electric other
    "electricity_grid",  # electricity grid
    "appliances"         # electric appliances (other electricity consumption)
]

energy = [
    "coal", "gas", "oil", "biomass", "nuclear", "hydro", 
    "wind", "solar", "electricity", "nonelectric"
]

level = ["primary", "secondary", "final", "useful"]

# Energy carrier and level combinations
energy_level = Dict(
    ("hydro", "primary") => true,
    ("wind", "primary") => true,
    ("nuclear", "primary") => true,
    ("coal", "final") => true,
    ("oil", "final") => true,
    ("gas", "final") => true,
    ("biomass", "final") => true,
    ("solar", "final") => true,
    ("electricity", "secondary") => true,
    ("electricity", "final") => true,
    ("electricity", "useful") => true,
    ("nonelectric", "useful") => true
)

# Technology share mappings
share = ["coal_nonelectric"]

tec_share = Dict(
    ("coal_nonelectric", "coal_nele") => true
)

tec_share_rhs = Dict(
    ("coal_nonelectric", "coal_nele") => true,
    ("coal_nonelectric", "gas_nele") => true,
    ("coal_nonelectric", "oil_nele") => true,
    ("coal_nonelectric", "bio_nele") => true,
    ("coal_nonelectric", "solar_nele") => true,
    ("coal_nonelectric", "other_nele") => true
)

# Energy-sector mapping
map_energy_sector = Dict(
    ("electricity", "useful", "ELEC") => true,
    ("nonelectric", "useful", "NELE") => true
)

# ------------------------------------------------------------------------------
# Parameter definitions
# ------------------------------------------------------------------------------

period_length = 10

# Input coefficients
input = Dict(
    ("electricity_grid", "electricity", "secondary") => 1,
    ("appliances", "electricity", "final") => 1,
    ("coal_ppl", "coal", "final") => 2.9,
    ("gas_ppl", "gas", "final") => 2.3,
    ("oil_ppl", "oil", "final") => 2.88,
    ("bio_ppl", "biomass", "final") => 3.59,
    ("nuclear_ppl", "nuclear", "primary") => 1,
    ("hydro_ppl", "hydro", "primary") => 1,
    ("wind_ppl", "wind", "primary") => 1,
    ("solar_PV_ppl", "solar", "final") => 1,
    ("coal_nele", "coal", "final") => 1,
    ("oil_nele", "oil", "final") => 1,
    ("gas_nele", "gas", "final") => 1,
    ("bio_nele", "biomass", "final") => 1,
    ("solar_nele", "solar", "final") => 1
)

# Output coefficients
output = Dict(
    ("electricity_grid", "electricity", "final") => 0.842,
    ("appliances", "electricity", "useful") => 1,
    ("coal_ppl", "electricity", "secondary") => 1,
    ("gas_ppl", "electricity", "secondary") => 1,
    ("oil_ppl", "electricity", "secondary") => 1,
    ("bio_ppl", "electricity", "secondary") => 1,
    ("hydro_ppl", "electricity", "secondary") => 1,
    ("wind_ppl", "electricity", "secondary") => 1,
    ("solar_PV_ppl", "electricity", "final") => 1,
    ("nuclear_ppl", "electricity", "secondary") => 1,
    ("other_ppl", "electricity", "secondary") => 1,
    ("coal_nele", "nonelectric", "useful") => 1,
    ("oil_nele", "nonelectric", "useful") => 1,
    ("gas_nele", "nonelectric", "useful") => 1,
    ("bio_nele", "nonelectric", "useful") => 1,
    ("solar_nele", "nonelectric", "useful") => 1,
    ("other_nele", "nonelectric", "useful") => 1,
    ("coal_extr", "coal", "final") => 1,
    ("gas_extr", "gas", "final") => 1,
    ("oil_extr", "oil", "final") => 1,
    ("nuclear_fuel", "nuclear", "primary") => 1,
    ("bio_pot", "biomass", "final") => 1,
    ("hydro_pot", "hydro", "primary") => 1,
    ("wind_pot", "wind", "primary") => 1,
    ("solar_pot", "solar", "final") => 1
)

# CO2 emission coefficients [tCO2/MWh]
CO2_emission = Dict(
    "gas_ppl" => 0.367,
    "coal_ppl" => 0.854,
    "oil_ppl" => 0.57,
    "coal_nele" => 0.342,
    "oil_nele" => 0.202,
    "gas_nele" => 0.26
)

# Maximum annual technology capacity growth rate
diffusion_up = Dict(
    "coal_ppl" => 0.075,
    "gas_ppl" => 0.10,
    "oil_ppl" => 0.075,
    "bio_ppl" => 0.075,
    "hydro_ppl" => 0.05,
    "wind_ppl" => 0.10,
    "solar_PV_ppl" => 0.15,
    "nuclear_ppl" => 0.05,
    "other_ppl" => 0.075,
    "coal_nele" => 0.075,
    "oil_nele" => 0.075,
    "gas_nele" => 0.10,
    "bio_nele" => 0.075,
    "solar_nele" => 0.15,
    "other_nele" => 0.075,
    "coal_extr" => 0.05,
    "gas_extr" => 0.05,
    "oil_extr" => 0.05,
    "nuclear_fuel" => 0.05,
    "bio_pot" => 0.05
)

# Maximum technology capacity constant addition per period (exact from GAMS lines 179-200)
startup = Dict(
    "coal_ppl" => 1, "gas_ppl" => 1, "oil_ppl" => 1, "bio_ppl" => 1,
    "hydro_ppl" => 1, "wind_ppl" => 1, "solar_PV_ppl" => 1, "nuclear_ppl" => 1, "other_ppl" => 1,
    "coal_nele" => 1, "oil_nele" => 1, "gas_nele" => 1, "bio_nele" => 1,
    "solar_nele" => 1, "other_nele" => 1,
    "coal_extr" => 1, "gas_extr" => 1, "oil_extr" => 1, "nuclear_fuel" => 1, "bio_pot" => 1,
    "hydro_pot" => 1, "wind_pot" => 1, "solar_pot" => 1, "electricity_grid" => 1, "appliances" => 1
)

# Demand in base year [PWh]
demand = Dict(
    ("electricity", "useful") => 22.60,
    ("nonelectric", "useful") => 87.3
)

# GDP [index]
gdp = Dict(
    2020 => 1,
    2030 => 1.2,
    2040 => 1.4,
    2050 => 1.7,
    2060 => 2.0,
    2070 => 2.4,
    2080 => 2.9
)

discount_rate = 0.05
beta = 0.7

# Share constraints
share_up = Dict("coal_nonelectric" => 0.4)
share_lo = Dict("coal_nonelectric" => 0)

# Investment cost [$/kW] (exact from GAMS table lines 229-240 - ONLY power plants have costs)
inv = Dict(
    # Power plants ONLY (exactly as in GAMS)
    ("coal_ppl", 2020) => 1500, ("coal_ppl", 2030) => 1500, ("coal_ppl", 2040) => 1500,
    ("coal_ppl", 2050) => 1500, ("coal_ppl", 2060) => 1500, ("coal_ppl", 2070) => 1500, ("coal_ppl", 2080) => 1500,
    ("gas_ppl", 2020) => 800, ("gas_ppl", 2030) => 800, ("gas_ppl", 2040) => 800,
    ("gas_ppl", 2050) => 800, ("gas_ppl", 2060) => 800, ("gas_ppl", 2070) => 800, ("gas_ppl", 2080) => 800,
    ("oil_ppl", 2020) => 950, ("oil_ppl", 2030) => 950, ("oil_ppl", 2040) => 950,
    ("oil_ppl", 2050) => 950, ("oil_ppl", 2060) => 950, ("oil_ppl", 2070) => 950, ("oil_ppl", 2080) => 950,
    ("hydro_ppl", 2020) => 3000, ("hydro_ppl", 2030) => 3000, ("hydro_ppl", 2040) => 3000,
    ("hydro_ppl", 2050) => 3000, ("hydro_ppl", 2060) => 3000, ("hydro_ppl", 2070) => 3000, ("hydro_ppl", 2080) => 3000,
    ("bio_ppl", 2020) => 1600, ("bio_ppl", 2030) => 1600, ("bio_ppl", 2040) => 1600,
    ("bio_ppl", 2050) => 1600, ("bio_ppl", 2060) => 1600, ("bio_ppl", 2070) => 1600, ("bio_ppl", 2080) => 1600,
    ("wind_ppl", 2020) => 1000, ("wind_ppl", 2030) => 1000, ("wind_ppl", 2040) => 1000,
    ("wind_ppl", 2050) => 1000, ("wind_ppl", 2060) => 1000, ("wind_ppl", 2070) => 1000, ("wind_ppl", 2080) => 1000,
    ("solar_PV_ppl", 2020) => 4000, ("solar_PV_ppl", 2030) => 4000, ("solar_PV_ppl", 2040) => 4000,
    ("solar_PV_ppl", 2050) => 4000, ("solar_PV_ppl", 2060) => 4000, ("solar_PV_ppl", 2070) => 4000, ("solar_PV_ppl", 2080) => 4000,
    ("nuclear_ppl", 2020) => 6000, ("nuclear_ppl", 2030) => 6000, ("nuclear_ppl", 2040) => 6000,
    ("nuclear_ppl", 2050) => 6000, ("nuclear_ppl", 2060) => 6000, ("nuclear_ppl", 2070) => 6000, ("nuclear_ppl", 2080) => 6000,
    ("other_ppl", 2020) => 4000, ("other_ppl", 2030) => 4000, ("other_ppl", 2040) => 4000,
    ("other_ppl", 2050) => 4000, ("other_ppl", 2060) => 4000, ("other_ppl", 2070) => 4000, ("other_ppl", 2080) => 4000
    # ALL OTHER TECHNOLOGIES HAVE ZERO INVESTMENT COSTS (as in GAMS)
)

# Fixed operation and maintenance cost [$/(kW/yr)] (exact from GAMS table lines 243-254 - ONLY power plants have costs)
fom = Dict(
    # Power plants ONLY (exactly as in GAMS)
    ("coal_ppl", 2020) => 40, ("coal_ppl", 2030) => 40, ("coal_ppl", 2040) => 40,
    ("coal_ppl", 2050) => 40, ("coal_ppl", 2060) => 40, ("coal_ppl", 2070) => 40, ("coal_ppl", 2080) => 40,
    ("gas_ppl", 2020) => 25, ("gas_ppl", 2030) => 25, ("gas_ppl", 2040) => 25,
    ("gas_ppl", 2050) => 25, ("gas_ppl", 2060) => 25, ("gas_ppl", 2070) => 25, ("gas_ppl", 2080) => 25,
    ("oil_ppl", 2020) => 25, ("oil_ppl", 2030) => 25, ("oil_ppl", 2040) => 25,
    ("oil_ppl", 2050) => 25, ("oil_ppl", 2060) => 25, ("oil_ppl", 2070) => 25, ("oil_ppl", 2080) => 25,
    ("bio_ppl", 2020) => 60, ("bio_ppl", 2030) => 60, ("bio_ppl", 2040) => 60,
    ("bio_ppl", 2050) => 60, ("bio_ppl", 2060) => 60, ("bio_ppl", 2070) => 60, ("bio_ppl", 2080) => 60,
    ("hydro_ppl", 2020) => 30, ("hydro_ppl", 2030) => 30, ("hydro_ppl", 2040) => 30,
    ("hydro_ppl", 2050) => 30, ("hydro_ppl", 2060) => 30, ("hydro_ppl", 2070) => 30, ("hydro_ppl", 2080) => 30,
    ("wind_ppl", 2020) => 40, ("wind_ppl", 2030) => 40, ("wind_ppl", 2040) => 40,
    ("wind_ppl", 2050) => 40, ("wind_ppl", 2060) => 40, ("wind_ppl", 2070) => 20, ("wind_ppl", 2080) => 20,
    ("solar_PV_ppl", 2020) => 25, ("solar_PV_ppl", 2030) => 25, ("solar_PV_ppl", 2040) => 25,
    ("solar_PV_ppl", 2050) => 25, ("solar_PV_ppl", 2060) => 25, ("solar_PV_ppl", 2070) => 25, ("solar_PV_ppl", 2080) => 25,
    ("nuclear_ppl", 2020) => 100, ("nuclear_ppl", 2030) => 100, ("nuclear_ppl", 2040) => 100,
    ("nuclear_ppl", 2050) => 100, ("nuclear_ppl", 2060) => 100, ("nuclear_ppl", 2070) => 100, ("nuclear_ppl", 2080) => 100,
    ("other_ppl", 2020) => 25, ("other_ppl", 2030) => 25, ("other_ppl", 2040) => 25,
    ("other_ppl", 2050) => 25, ("other_ppl", 2060) => 25, ("other_ppl", 2070) => 25, ("other_ppl", 2080) => 25
    # ALL OTHER TECHNOLOGIES HAVE ZERO FIXED O&M COSTS (as in GAMS)
)

# Variable cost [$/MWh] (exact from GAMS table lines 257-279)
vom = Dict()
# Initialize all to zero first
for tech in technology, year in year_all
    vom[(tech, year)] = 0.0
end
# Set non-zero values exactly as in GAMS
for year in year_all
    vom[("electricity_grid", year)] = 47.8
    vom[("coal_extr", year)] = 7.2
    vom[("gas_extr", year)] = 14.4
    vom[("oil_extr", year)] = 40.0
    vom[("nuclear_fuel", year)] = 10.0
    vom[("bio_pot", year)] = 18.0
    # All other technologies have 0.0 VOM costs (as in GAMS)
end

# Technical lifetime
lifetime = Dict(
    "coal_extr" => 30, "gas_extr" => 30, "oil_extr" => 30, "nuclear_fuel" => 40,
    "bio_pot" => 20, "hydro_pot" => 80, "wind_pot" => 20, "solar_pot" => 20,
    "coal_ppl" => 30, "gas_ppl" => 20, "oil_ppl" => 20, "bio_ppl" => 20,
    "hydro_ppl" => 80, "wind_ppl" => 20, "solar_PV_ppl" => 20, "nuclear_ppl" => 40,
    "other_ppl" => 20, "coal_nele" => 20, "oil_nele" => 20, "gas_nele" => 20,
    "bio_nele" => 20, "solar_nele" => 20, "other_nele" => 20,
    "electricity_grid" => 40, "appliances" => 15
)

# Full load hours
hours = Dict(
    "coal_extr" => 7000, "gas_extr" => 7000, "oil_extr" => 7000, "nuclear_fuel" => 7000,
    "bio_pot" => 7000, "hydro_pot" => 4500, "wind_pot" => 2000, "solar_pot" => 1200,
    "coal_ppl" => 7000, "gas_ppl" => 6000, "oil_ppl" => 6000, "bio_ppl" => 6000,
    "hydro_ppl" => 4500, "wind_ppl" => 2000, "solar_PV_ppl" => 1200, "nuclear_ppl" => 7500,
    "other_ppl" => 4000, "coal_nele" => 7000, "oil_nele" => 7000, "gas_nele" => 7000,
    "bio_nele" => 7000, "solar_nele" => 7000, "other_nele" => 7000,
    "electricity_grid" => 8760, "appliances" => 8760
)

# Calculate cost parameters with type stability optimizations
function calculate_costs()
    # Type-stable dictionaries
    cost_capacity = Dict{Tuple{String,Int}, Float64}()
    cost_activity = Dict{Tuple{String,Int}, Float64}()
    cost = Dict{Tuple{String,Int}, Float64}()
    
    # Pre-calculate annuity factors to avoid repeated computation
    annuity_factors = Dict{String, Float64}()
    for tech in technology
        if haskey(lifetime, tech) && lifetime[tech] > 0
            lt = lifetime[tech]
            annuity_factors[tech] = ((1 + discount_rate)^lt * discount_rate) / ((1 + discount_rate)^lt - 1)
        end
    end
    
    # Pre-filter valid technologies to avoid repeated checks
    valid_techs = [tech for tech in technology 
                   if haskey(lifetime, tech) && lifetime[tech] > 0 && 
                      haskey(hours, tech) && hours[tech] > 0]
    
    for tech in valid_techs, year in year_all
        # Capacity-related costs (annuity of investment + FOM) in $/kW
        if haskey(inv, (tech, year)) && haskey(fom, (tech, year))
            # GAMS multiplies by 1000 (line 348-349) for unit conversion
            cost_capacity[(tech, year)] = (inv[(tech, year)] * annuity_factors[tech] + fom[(tech, year)]) * 1000
        end
        
        # Activity-related costs in $/MWh
        if haskey(vom, (tech, year))
            cost_activity[(tech, year)] = vom[(tech, year)]
        end
        
        # Total costs in $/MWh
        if haskey(cost_capacity, (tech, year)) && haskey(cost_activity, (tech, year))
            cost[(tech, year)] = cost_capacity[(tech, year)] / hours[tech] + cost_activity[(tech, year)]
        end
    end
    
    return cost_capacity, cost_activity, cost
end

cost_capacity, cost_activity, cost = calculate_costs()

# Function to create energy model equations
function create_energy_model!(model)
    # Variables (exact from GAMS lines 367-373)
    @variable(model, ACT[technology, year_all] >= 0)  # technology activity
    @variable(model, CAP_NEW[technology, year_all] >= 0)  # new technology capacity
    @variable(model, EMISS[year_all])  # CO2 emissions
    @variable(model, CUM_EMISS)  # cumulative CO2 emissions
    @variable(model, TOTAL_COST)  # total discounted systems costs
    
    # COST_ANNUAL might already be defined as a shared variable
    if !haskey(object_dictionary(model), :COST_ANNUAL)
        @variable(model, COST_ANNUAL[year_all])  # costs per year
    end
    
    # Get references to variables for use in constraints
    COST_ANNUAL = model[:COST_ANNUAL]
    
    # Energy balance equations (exact from GAMS lines 396-398)
    for e in energy, l in level, y in year_all
        if haskey(energy_level, (e, l))
            @constraint(model,
                sum(ACT[t, y] * (get(output, (t, e, l), 0) - get(input, (t, e, l), 0)) for t in technology) -
                get(demand, (e, l), 0) * (gdp[y] / gdp[2020])^beta >= 0
            )
        end
    end
    
    # Capacity balance (exact from GAMS lines 400-401)
    for t in technology, y in year_all
        if haskey(hours, t) && haskey(lifetime, t)
            year_idx = findfirst(==(y), year_all)
            @constraint(model,
                ACT[t, y] <= sum(CAP_NEW[t, year_all[i]] * hours[t]
                                 for i in 1:year_idx if (year_idx - i + 1) * period_length <= lifetime[t])
            )
        end
    end
    
    # Emissions (exact from GAMS lines 408-412)
    @constraint(model, [y in year_all],
        sum(ACT[t, y] * get(CO2_emission, t, 0) for t in technology) == EMISS[y])
    @constraint(model,
        sum(EMISS[y] * period_length for y in year_all) == CUM_EMISS)
    
    # Diffusion constraints (exact from GAMS lines 414-415)
    for t in technology, y in year_all[2:end]
        if haskey(diffusion_up, t)
            y_prev = year_all[findfirst(==(y), year_all)-1]
            @constraint(model,
                CAP_NEW[t, y] <= CAP_NEW[t, y_prev] * (1 + diffusion_up[t])^period_length + get(startup, t, 0)
            )
        end
    end
    
    # Share constraints (exact from GAMS lines 417-421)
    for s in share, y in year_all
        if haskey(share_up, s)
            lhs = [t for t in technology if haskey(tec_share, (s, t))]
            rhs = [t for t in technology if haskey(tec_share_rhs, (s, t))]
            @constraint(model, sum(ACT[t, y] for t in lhs) <= share_up[s] * sum(ACT[t, y] for t in rhs))
        end
        if haskey(share_lo, s)
            lhs = [t for t in technology if haskey(tec_share, (s, t))]
            rhs = [t for t in technology if haskey(tec_share_rhs, (s, t))]
            @constraint(model, sum(ACT[t, y] for t in lhs) >= share_lo[s] * sum(ACT[t, y] for t in rhs))
        end
    end
    
    # Cost equations (exact from GAMS lines 423-436)
    year_idx = Dict(y => i for (i, y) in enumerate(year_all))
    for y in year_all
        @constraint(model,
            sum(ACT[t, y] * get(vom, (t, y), 0) for t in technology) +
            sum(sum(CAP_NEW[t, y2] * get(cost_capacity, (t, y2), 0)
                    for y2 in year_all
                    if year_idx[y2] <= year_idx[y] && (year_idx[y] - year_idx[y2] + 1) * period_length <= get(lifetime, t, 0) && haskey(lifetime, t))
                for t in technology)
            ==
            COST_ANNUAL[y]
        )
    end
    
    @constraint(model,
        sum(COST_ANNUAL[y] * period_length * (1 + discount_rate)^(-period_length * (year_idx[y] - 1)) for y in year_all)
        ==
        TOTAL_COST
    )
    
    # Historical calibration (exact from GAMS lines 461-474)
    fix(ACT["coal_ppl", 2020], 9.462; force=true)
    fix(ACT["oil_ppl", 2020], 0.7; force=true)
    fix(ACT["solar_PV_ppl", 2020], 0.839; force=true)
    fix(ACT["gas_ppl", 2020], 6.36; force=true)
    fix(ACT["nuclear_ppl", 2020], 2.68; force=true)
    fix(ACT["hydro_ppl", 2020], 4.36; force=true)
    fix(ACT["wind_ppl", 2020], 1.6; force=true)
    fix(ACT["bio_ppl", 2020], 0.69; force=true)
    set_lower_bound(ACT["other_ppl", 2020], 0.127)
    fix(ACT["coal_nele", 2020], 10.7; force=true)
    fix(ACT["oil_nele", 2020], 43.0; force=true)
    fix(ACT["gas_nele", 2020], 18.7; force=true)
    fix(ACT["bio_nele", 2020], 10.6; force=true)
    set_lower_bound(ACT["other_nele", 2020], 0.28)
    
    # Objective (minimize TOTAL_COST)
    @objective(model, Min, TOTAL_COST)
    
    return model
end

# Model calibration and resource constraints
energy_calibration = Dict(
    ("coal_ppl", 2020) => 9.462,
    ("oil_ppl", 2020) => 0.7,
    ("solar_PV_ppl", 2020) => 0.839,
    ("gas_ppl", 2020) => 6.36,
    ("nuclear_ppl", 2020) => 2.68,
    ("hydro_ppl", 2020) => 4.36,
    ("wind_ppl", 2020) => 1.6,
    ("bio_ppl", 2020) => 0.69,
    ("coal_nele", 2020) => 10.7,
    ("oil_nele", 2020) => 43.0,
    ("gas_nele", 2020) => 18.7,
    ("bio_nele", 2020) => 10.6
)

energy_calibration_lo = Dict(
    ("other_ppl", 2020) => 0.127,
    ("other_nele", 2020) => 0.28
)