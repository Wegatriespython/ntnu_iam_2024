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

# Maximum technology capacity constant addition per period
# Maximum technology capacity constant addition per period
# Set to small values or zero for most technologies
startup = Dict(
    "coal_ppl" => 0.001,      # 1 GW
    "gas_ppl" => 0.001,       # 1 GW
    "oil_ppl" => 0.0001,      # 0.1 GW
    "bio_ppl" => 0.0001,      # 0.1 GW
    "hydro_ppl" => 0.001,     # 1 GW
    "wind_ppl" => 0.001,      # 1 GW
    "solar_PV_ppl" => 0.001,  # 1 GW
    "nuclear_ppl" => 0.001,   # 1 GW
    "other_ppl" => 0.0,       # No startup
    # Extraction technologies shouldn't have startup capacity
    "coal_extr" => 0.0,
    "gas_extr" => 0.0,
    "oil_extr" => 0.0,
    "nuclear_fuel" => 0.0,
    "bio_pot" => 0.0,
    "hydro_pot" => 0.0,
    "wind_pot" => 0.0,
    "solar_pot" => 0.0,
    # Non-electric and grid
    "coal_nele" => 0.0,
    "oil_nele" => 0.0,
    "gas_nele" => 0.0,
    "bio_nele" => 0.0,
    "solar_nele" => 0.0,
    "other_nele" => 0.0,
    "electricity_grid" => 0.0,
    "appliances" => 0.0
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

# Investment cost [$/kW]
inv = Dict(
    # Power plants
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
    ("other_ppl", 2050) => 4000, ("other_ppl", 2060) => 4000, ("other_ppl", 2070) => 4000, ("other_ppl", 2080) => 4000,
    
    # Fossil fuel extraction technologies (mining/extraction infrastructure)
    ("coal_extr", 2020) => 1000, ("coal_extr", 2030) => 1000, ("coal_extr", 2040) => 1000,
    ("coal_extr", 2050) => 1000, ("coal_extr", 2060) => 1000, ("coal_extr", 2070) => 1000, ("coal_extr", 2080) => 1000,
    ("gas_extr", 2020) => 500, ("gas_extr", 2030) => 500, ("gas_extr", 2040) => 500,
    ("gas_extr", 2050) => 500, ("gas_extr", 2060) => 500, ("gas_extr", 2070) => 500, ("gas_extr", 2080) => 500,
    ("oil_extr", 2020) => 800, ("oil_extr", 2030) => 800, ("oil_extr", 2040) => 800,
    ("oil_extr", 2050) => 800, ("oil_extr", 2060) => 800, ("oil_extr", 2070) => 800, ("oil_extr", 2080) => 800,
    
    # Nuclear fuel processing
    ("nuclear_fuel", 2020) => 2000, ("nuclear_fuel", 2030) => 2000, ("nuclear_fuel", 2040) => 2000,
    ("nuclear_fuel", 2050) => 2000, ("nuclear_fuel", 2060) => 2000, ("nuclear_fuel", 2070) => 2000, ("nuclear_fuel", 2080) => 2000,
    
    # Renewable resource potentials (minimal infrastructure costs)
    ("bio_pot", 2020) => 100, ("bio_pot", 2030) => 100, ("bio_pot", 2040) => 100,
    ("bio_pot", 2050) => 100, ("bio_pot", 2060) => 100, ("bio_pot", 2070) => 100, ("bio_pot", 2080) => 100,
    ("hydro_pot", 2020) => 50, ("hydro_pot", 2030) => 50, ("hydro_pot", 2040) => 50,
    ("hydro_pot", 2050) => 50, ("hydro_pot", 2060) => 50, ("hydro_pot", 2070) => 50, ("hydro_pot", 2080) => 50,
    ("wind_pot", 2020) => 50, ("wind_pot", 2030) => 50, ("wind_pot", 2040) => 50,
    ("wind_pot", 2050) => 50, ("wind_pot", 2060) => 50, ("wind_pot", 2070) => 50, ("wind_pot", 2080) => 50,
    ("solar_pot", 2020) => 50, ("solar_pot", 2030) => 50, ("solar_pot", 2040) => 50,
    ("solar_pot", 2050) => 50, ("solar_pot", 2060) => 50, ("solar_pot", 2070) => 50, ("solar_pot", 2080) => 50,
    
    # Infrastructure
    ("electricity_grid", 2020) => 300, ("electricity_grid", 2030) => 300, ("electricity_grid", 2040) => 300,
    ("electricity_grid", 2050) => 300, ("electricity_grid", 2060) => 300, ("electricity_grid", 2070) => 300, ("electricity_grid", 2080) => 300,
    ("appliances", 2020) => 200, ("appliances", 2030) => 200, ("appliances", 2040) => 200,
    ("appliances", 2050) => 200, ("appliances", 2060) => 200, ("appliances", 2070) => 200, ("appliances", 2080) => 200,
    
    # Non-electric technologies
    ("coal_nele", 2020) => 500, ("coal_nele", 2030) => 500, ("coal_nele", 2040) => 500,
    ("coal_nele", 2050) => 500, ("coal_nele", 2060) => 500, ("coal_nele", 2070) => 500, ("coal_nele", 2080) => 500,
    ("oil_nele", 2020) => 400, ("oil_nele", 2030) => 400, ("oil_nele", 2040) => 400,
    ("oil_nele", 2050) => 400, ("oil_nele", 2060) => 400, ("oil_nele", 2070) => 400, ("oil_nele", 2080) => 400,
    ("gas_nele", 2020) => 350, ("gas_nele", 2030) => 350, ("gas_nele", 2040) => 350,
    ("gas_nele", 2050) => 350, ("gas_nele", 2060) => 350, ("gas_nele", 2070) => 350, ("gas_nele", 2080) => 350,
    ("bio_nele", 2020) => 600, ("bio_nele", 2030) => 600, ("bio_nele", 2040) => 600,
    ("bio_nele", 2050) => 600, ("bio_nele", 2060) => 600, ("bio_nele", 2070) => 600, ("bio_nele", 2080) => 600,
    ("solar_nele", 2020) => 1500, ("solar_nele", 2030) => 1500, ("solar_nele", 2040) => 1500,
    ("solar_nele", 2050) => 1500, ("solar_nele", 2060) => 1500, ("solar_nele", 2070) => 1500, ("solar_nele", 2080) => 1500,
    ("other_nele", 2020) => 500, ("other_nele", 2030) => 500, ("other_nele", 2040) => 500,
    ("other_nele", 2050) => 500, ("other_nele", 2060) => 500, ("other_nele", 2070) => 500, ("other_nele", 2080) => 500
)

# Fixed operation and maintenance cost [$/(kW/yr)]
fom = Dict(
    # Power plants
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
    ("other_ppl", 2050) => 25, ("other_ppl", 2060) => 25, ("other_ppl", 2070) => 25, ("other_ppl", 2080) => 25,
    
    # Fossil fuel extraction technologies (3-4% of investment cost)
    ("coal_extr", 2020) => 30, ("coal_extr", 2030) => 30, ("coal_extr", 2040) => 30,
    ("coal_extr", 2050) => 30, ("coal_extr", 2060) => 30, ("coal_extr", 2070) => 30, ("coal_extr", 2080) => 30,
    ("gas_extr", 2020) => 20, ("gas_extr", 2030) => 20, ("gas_extr", 2040) => 20,
    ("gas_extr", 2050) => 20, ("gas_extr", 2060) => 20, ("gas_extr", 2070) => 20, ("gas_extr", 2080) => 20,
    ("oil_extr", 2020) => 25, ("oil_extr", 2030) => 25, ("oil_extr", 2040) => 25,
    ("oil_extr", 2050) => 25, ("oil_extr", 2060) => 25, ("oil_extr", 2070) => 25, ("oil_extr", 2080) => 25,
    
    # Nuclear fuel processing (2% of investment cost)
    ("nuclear_fuel", 2020) => 40, ("nuclear_fuel", 2030) => 40, ("nuclear_fuel", 2040) => 40,
    ("nuclear_fuel", 2050) => 40, ("nuclear_fuel", 2060) => 40, ("nuclear_fuel", 2070) => 40, ("nuclear_fuel", 2080) => 40,
    
    # Renewable resource potentials (2% of investment cost)
    ("bio_pot", 2020) => 2, ("bio_pot", 2030) => 2, ("bio_pot", 2040) => 2,
    ("bio_pot", 2050) => 2, ("bio_pot", 2060) => 2, ("bio_pot", 2070) => 2, ("bio_pot", 2080) => 2,
    ("hydro_pot", 2020) => 1, ("hydro_pot", 2030) => 1, ("hydro_pot", 2040) => 1,
    ("hydro_pot", 2050) => 1, ("hydro_pot", 2060) => 1, ("hydro_pot", 2070) => 1, ("hydro_pot", 2080) => 1,
    ("wind_pot", 2020) => 1, ("wind_pot", 2030) => 1, ("wind_pot", 2040) => 1,
    ("wind_pot", 2050) => 1, ("wind_pot", 2060) => 1, ("wind_pot", 2070) => 1, ("wind_pot", 2080) => 1,
    ("solar_pot", 2020) => 1, ("solar_pot", 2030) => 1, ("solar_pot", 2040) => 1,
    ("solar_pot", 2050) => 1, ("solar_pot", 2060) => 1, ("solar_pot", 2070) => 1, ("solar_pot", 2080) => 1,
    
    # Infrastructure (2.5% of investment cost)
    ("electricity_grid", 2020) => 8, ("electricity_grid", 2030) => 8, ("electricity_grid", 2040) => 8,
    ("electricity_grid", 2050) => 8, ("electricity_grid", 2060) => 8, ("electricity_grid", 2070) => 8, ("electricity_grid", 2080) => 8,
    ("appliances", 2020) => 5, ("appliances", 2030) => 5, ("appliances", 2040) => 5,
    ("appliances", 2050) => 5, ("appliances", 2060) => 5, ("appliances", 2070) => 5, ("appliances", 2080) => 5,
    
    # Non-electric technologies (3% of investment cost)
    ("coal_nele", 2020) => 15, ("coal_nele", 2030) => 15, ("coal_nele", 2040) => 15,
    ("coal_nele", 2050) => 15, ("coal_nele", 2060) => 15, ("coal_nele", 2070) => 15, ("coal_nele", 2080) => 15,
    ("oil_nele", 2020) => 12, ("oil_nele", 2030) => 12, ("oil_nele", 2040) => 12,
    ("oil_nele", 2050) => 12, ("oil_nele", 2060) => 12, ("oil_nele", 2070) => 12, ("oil_nele", 2080) => 12,
    ("gas_nele", 2020) => 10, ("gas_nele", 2030) => 10, ("gas_nele", 2040) => 10,
    ("gas_nele", 2050) => 10, ("gas_nele", 2060) => 10, ("gas_nele", 2070) => 10, ("gas_nele", 2080) => 10,
    ("bio_nele", 2020) => 18, ("bio_nele", 2030) => 18, ("bio_nele", 2040) => 18,
    ("bio_nele", 2050) => 18, ("bio_nele", 2060) => 18, ("bio_nele", 2070) => 18, ("bio_nele", 2080) => 18,
    ("solar_nele", 2020) => 45, ("solar_nele", 2030) => 45, ("solar_nele", 2040) => 45,
    ("solar_nele", 2050) => 45, ("solar_nele", 2060) => 45, ("solar_nele", 2070) => 45, ("solar_nele", 2080) => 45,
    ("other_nele", 2020) => 15, ("other_nele", 2030) => 15, ("other_nele", 2040) => 15,
    ("other_nele", 2050) => 15, ("other_nele", 2060) => 15, ("other_nele", 2070) => 15, ("other_nele", 2080) => 15
)

# Variable cost [$/MWh]
vom = Dict(
    ("coal_ppl", 2020) => 0.0, ("coal_ppl", 2030) => 0.0, ("coal_ppl", 2040) => 0.0,
    ("coal_ppl", 2050) => 0.0, ("coal_ppl", 2060) => 0.0, ("coal_ppl", 2070) => 0.0, ("coal_ppl", 2080) => 0.0,
    ("gas_ppl", 2020) => 0.0, ("gas_ppl", 2030) => 0.0, ("gas_ppl", 2040) => 0.0,
    ("gas_ppl", 2050) => 0.0, ("gas_ppl", 2060) => 0.0, ("gas_ppl", 2070) => 0.0, ("gas_ppl", 2080) => 0.0,
    ("oil_ppl", 2020) => 0.0, ("oil_ppl", 2030) => 0.0, ("oil_ppl", 2040) => 0.0,
    ("oil_ppl", 2050) => 0.0, ("oil_ppl", 2060) => 0.0, ("oil_ppl", 2070) => 0.0, ("oil_ppl", 2080) => 0.0,
    ("bio_ppl", 2020) => 0.0, ("bio_ppl", 2030) => 0.0, ("bio_ppl", 2040) => 0.0,
    ("bio_ppl", 2050) => 0.0, ("bio_ppl", 2060) => 0.0, ("bio_ppl", 2070) => 0.0, ("bio_ppl", 2080) => 0.0,
    ("nuclear_ppl", 2020) => 0.0, ("nuclear_ppl", 2030) => 0.0, ("nuclear_ppl", 2040) => 0.0,
    ("nuclear_ppl", 2050) => 0.0, ("nuclear_ppl", 2060) => 0.0, ("nuclear_ppl", 2070) => 0.0, ("nuclear_ppl", 2080) => 0.0,
    ("electricity_grid", 2020) => 47.8, ("electricity_grid", 2030) => 47.8, ("electricity_grid", 2040) => 47.8,
    ("electricity_grid", 2050) => 47.8, ("electricity_grid", 2060) => 47.8, ("electricity_grid", 2070) => 47.8, ("electricity_grid", 2080) => 47.8,
    ("coal_nele", 2020) => 25.0, ("coal_nele", 2030) => 27.0, ("coal_nele", 2040) => 29.0,
    ("coal_nele", 2050) => 31.0, ("coal_nele", 2060) => 33.0, ("coal_nele", 2070) => 35.0, ("coal_nele", 2080) => 37.0,
    ("gas_nele", 2020) => 35.0, ("gas_nele", 2030) => 37.0, ("gas_nele", 2040) => 39.0,
    ("gas_nele", 2050) => 41.0, ("gas_nele", 2060) => 43.0, ("gas_nele", 2070) => 45.0, ("gas_nele", 2080) => 47.0,
    ("oil_nele", 2020) => 55.0, ("oil_nele", 2030) => 58.0, ("oil_nele", 2040) => 61.0,
    ("oil_nele", 2050) => 64.0, ("oil_nele", 2060) => 67.0, ("oil_nele", 2070) => 70.0, ("oil_nele", 2080) => 73.0,
    ("bio_nele", 2020) => 20.0, ("bio_nele", 2030) => 21.0, ("bio_nele", 2040) => 22.0,
    ("bio_nele", 2050) => 23.0, ("bio_nele", 2060) => 24.0, ("bio_nele", 2070) => 25.0, ("bio_nele", 2080) => 26.0,
    ("solar_nele", 2020) => 5.0, ("solar_nele", 2030) => 4.0, ("solar_nele", 2040) => 3.0,
    ("solar_nele", 2050) => 2.0, ("solar_nele", 2060) => 1.5, ("solar_nele", 2070) => 1.0, ("solar_nele", 2080) => 0.5,
    ("other_nele", 2020) => 30.0, ("other_nele", 2030) => 32.0, ("other_nele", 2040) => 34.0,
    ("other_nele", 2050) => 36.0, ("other_nele", 2060) => 38.0, ("other_nele", 2070) => 40.0, ("other_nele", 2080) => 42.0,
    ("coal_extr", 2020) => 7.2, ("coal_extr", 2030) => 7.2, ("coal_extr", 2040) => 7.2,
    ("coal_extr", 2050) => 7.2, ("coal_extr", 2060) => 7.2, ("coal_extr", 2070) => 7.2, ("coal_extr", 2080) => 7.2,
    ("gas_extr", 2020) => 14.4, ("gas_extr", 2030) => 14.4, ("gas_extr", 2040) => 14.4,
    ("gas_extr", 2050) => 14.4, ("gas_extr", 2060) => 14.4, ("gas_extr", 2070) => 14.4, ("gas_extr", 2080) => 14.4,
    ("oil_extr", 2020) => 40.0, ("oil_extr", 2030) => 40.0, ("oil_extr", 2040) => 40.0,
    ("oil_extr", 2050) => 40.0, ("oil_extr", 2060) => 40.0, ("oil_extr", 2070) => 40.0, ("oil_extr", 2080) => 40.0,
    ("nuclear_fuel", 2020) => 10.0, ("nuclear_fuel", 2030) => 10.0, ("nuclear_fuel", 2040) => 10.0,
    ("nuclear_fuel", 2050) => 10.0, ("nuclear_fuel", 2060) => 10.0, ("nuclear_fuel", 2070) => 10.0, ("nuclear_fuel", 2080) => 10.0,
    ("bio_pot", 2020) => 18.0, ("bio_pot", 2030) => 18.0, ("bio_pot", 2040) => 18.0,
    ("bio_pot", 2050) => 18.0, ("bio_pot", 2060) => 18.0, ("bio_pot", 2070) => 18.0, ("bio_pot", 2080) => 18.0,
    ("hydro_pot", 2020) => 0.0, ("hydro_pot", 2030) => 0.0, ("hydro_pot", 2040) => 0.0,
    ("hydro_pot", 2050) => 0.0, ("hydro_pot", 2060) => 0.0, ("hydro_pot", 2070) => 0.0, ("hydro_pot", 2080) => 0.0,
    ("wind_pot", 2020) => 0.0, ("wind_pot", 2030) => 0.0, ("wind_pot", 2040) => 0.0,
    ("wind_pot", 2050) => 0.0, ("wind_pot", 2060) => 0.0, ("wind_pot", 2070) => 0.0, ("wind_pot", 2080) => 0.0,
    ("solar_pot", 2020) => 0.0, ("solar_pot", 2030) => 0.0, ("solar_pot", 2040) => 0.0,
    ("solar_pot", 2050) => 0.0, ("solar_pot", 2060) => 0.0, ("solar_pot", 2070) => 0.0, ("solar_pot", 2080) => 0.0
)

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
            # Fix: cost_capacity should be in $/kW-year, not $/W-year
            cost_capacity[(tech, year)] = inv[(tech, year)] * annuity_factors[tech] + fom[(tech, year)]
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
    # Variables
    @variable(model, ACT[technology, year_all] >= 0)  # technology activity
    @variable(model, CAP_NEW[technology, year_all] >= 0)  # new technology capacity
    @variable(model, EMISS[year_all])  # CO2 emissions
    @variable(model, CUM_EMISS)  # cumulative CO2 emissions
    @variable(model, TOTAL_COST)  # total discounted systems costs
    
    # Energy balance equation
    for e in energy, l in level, y in year_all
        if haskey(energy_level, (e, l)) && haskey(demand, (e, l))
            # Standard energy balance constraint
            @constraint(model,
                sum(ACT[tech, y] * (get(output, (tech, e, l), 0) - get(input, (tech, e, l), 0)) 
                    for tech in technology) - 
                demand[(e, l)] * (gdp[y]/gdp[2020])^beta >= 0
            )
        end
    end
    
    # Capacity balance equation
    for tech in technology, y in year_all
        if haskey(hours, tech) && hours[tech] > 0 && haskey(lifetime, tech) && lifetime[tech] > 0
            year_idx = findfirst(==(y), year_all)
            @constraint(model,
                ACT[tech, y] <= 
                sum(CAP_NEW[tech, year_all[i]] * hours[tech] 
                    for i in 1:year_idx 
                    if (year_idx - i + 1) * period_length <= lifetime[tech])
            )
        end
    end
    
    # Emission equation
    for y in year_all
        @constraint(model,
            sum(ACT[tech, y] * get(CO2_emission, tech, 0) for tech in technology) == EMISS[y]
        )
    end
    
    # Cumulative emissions
    @constraint(model,
        sum(EMISS[y] * period_length for y in year_all) == CUM_EMISS
    )
    
    # Diffusion constraints
    for tech in technology, y in year_all[2:end]
        if haskey(diffusion_up, tech)
            y_prev = year_all[findfirst(==(y), year_all) - 1]
            @constraint(model,
                CAP_NEW[tech, y] <= 
                CAP_NEW[tech, y_prev] * (1 + diffusion_up[tech])^period_length + 
                get(startup, tech, 0)
            )
        end
    end
    
    # Share constraints
    for s in share, y in year_all
        if haskey(share_up, s)
            lhs_techs = [tech for tech in technology if haskey(tec_share, (s, tech))]
            rhs_techs = [tech for tech in technology if haskey(tec_share_rhs, (s, tech))]
            
            @constraint(model,
                sum(ACT[tech, y] for tech in lhs_techs) <= 
                sum(ACT[tech, y] for tech in rhs_techs) * share_up[s]
            )
        end
        
        if haskey(share_lo, s)
            lhs_techs = [tech for tech in technology if haskey(tec_share, (s, tech))]
            rhs_techs = [tech for tech in technology if haskey(tec_share_rhs, (s, tech))]
            
            @constraint(model,
                sum(ACT[tech, y] for tech in lhs_techs) >= 
                sum(ACT[tech, y] for tech in rhs_techs) * share_lo[s]
            )
        end
    end
    
    # Cost equations
    # Will be defined later when COST_ANNUAL is available from shared variables
    
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