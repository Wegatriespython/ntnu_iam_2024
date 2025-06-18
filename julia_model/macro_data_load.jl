# macro_data_load.jl - Equivalent to macro/macro_data_load.gms
# Set and parameter definitions for MACRO model

using JuMP

# ------------------------------------------------------------------------------
# Set and parameter definitions
# ------------------------------------------------------------------------------

# Sets specific to MACRO
year = year_all  # years included in model instance
macro_horizon = year_all  # set of periods included in MACRO model horizon
macro_base_period = [2020]  # base year period
first_period = [2020]  # first period in model horizon
last_period = [2080]   # last period in model horizon

# Sequential period mapping
seq_period = Dict()
for i in 1:(length(year_all)-1)
    seq_period[year_all[i]] = year_all[i+1]
end

epsilon = 0.01  # small number to avoid divergences

# Parameters
i0 = 0.0          # Initial investment in base year
c0 = 0.0          # Initial consumption in base year  
k0 = 0.0          # Initial capital in base year
y0 = 0.0          # Initial output in base year

gdp_base = 71.0   # Initial GDP (Trillion $) in base year
kgdp = 2.8        # Initial capital to GDP ratio in base year
depr = 0.05       # Annual percent depreciation
drate = 0.05      # Social discount rate
esub = 0.3        # Elasticity between capital-labor (K-L) and energy (Sum E)
rho = 0.0         # Production function exponent (calculated from esub)
kpvs = 0.28       # Capital value share parameter
elvs = 0.42       # Electricity value share parameter
ecst0 = 0.0       # Energy costs in base year

# Demand base for energy services
demand_base = Dict()

# Energy service parameters from MESSAGE
enestart = Dict()  # Consumption level of energy services from MESSAGE
eneprice = Dict()  # Shadow prices of energy services from MESSAGE  
total_cost_energy = Dict()  # Total energy system costs from MESSAGE

# Growth rates of potential GDP
grow = Dict(
    2020 => 0.03, 2030 => 0.03, 2040 => 0.03, 2050 => 0.03,
    2060 => 0.03, 2070 => 0.02, 2080 => 0.01
)

# Autonomous energy efficiency improvement (AEEI)
aeei = Dict(
    ("ELEC", 2020) => 0.02, ("ELEC", 2030) => 0.02, ("ELEC", 2040) => 0.02,
    ("ELEC", 2050) => 0.02, ("ELEC", 2060) => 0.02, ("ELEC", 2070) => 0.02, ("ELEC", 2080) => 0.02,
    ("NELE", 2020) => 0.02, ("NELE", 2030) => 0.02, ("NELE", 2040) => 0.02,
    ("NELE", 2050) => 0.02, ("NELE", 2060) => 0.02, ("NELE", 2070) => 0.02, ("NELE", 2080) => 0.02
)

# Time parameters
duration_period = 10  # duration of one multi-year period (in years)

# Utility and labor parameters
udf = Dict()          # Utility discount factor
labor = Dict()        # Labor force (efficiency units)
newlab = Dict()       # New vintage of labor force
aeei_factor = Dict()  # Cumulative effect of AEEI
finite_time_corr = Dict()  # finite time horizon correction factor
lotol = 0.05          # Tolerance factor for lower bounds

# Start values for variables
SVKN = Dict()         # start values for new capital variable KN
SVNEWE = Dict()       # start values for new energy variable

# MESSAGE data tables
demand_MESSAGE = Dict(
    ("ELEC", 2020) => 22.6, ("ELEC", 2030) => 25.7, ("ELEC", 2040) => 28.60,
    ("ELEC", 2050) => 32.80, ("ELEC", 2060) => 36.7, ("ELEC", 2070) => 41.7, ("ELEC", 2080) => 47.6,
    ("NELE", 2020) => 87.3, ("NELE", 2030) => 99.0, ("NELE", 2040) => 110.00,
    ("NELE", 2050) => 117.00, ("NELE", 2060) => 142.0, ("NELE", 2070) => 161.0, ("NELE", 2080) => 184.0
)

price_MESSAGE = Dict(
    ("ELEC", 2020) => 0.0567, ("ELEC", 2030) => 0.0567, ("ELEC", 2040) => 0.0567,
    ("ELEC", 2050) => 0.0567, ("ELEC", 2060) => 0.0567, ("ELEC", 2070) => 0.0567, ("ELEC", 2080) => 0.0567,
    ("NELE", 2020) => 0.020, ("NELE", 2030) => 0.020, ("NELE", 2040) => 0.020,
    ("NELE", 2050) => 0.020, ("NELE", 2060) => 0.020, ("NELE", 2070) => 0.020, ("NELE", 2080) => 0.020
)

cost_MESSAGE = Dict(
    2020 => 5.053, 2030 => 3.045, 2040 => 3.080, 2050 => 3.516,
    2060 => 3.852, 2070 => 4.376, 2080 => 4.996
)

# Map MESSAGE data to MACRO structure
for s in sector, y in year_all
    enestart[(s, y)] = demand_MESSAGE[(s, y)]
    eneprice[(s, y)] = price_MESSAGE[(s, y)]
end

for y in year_all
    total_cost_energy[y] = cost_MESSAGE[y]
end

# Base year demand levels
for s in sector
    demand_base[s] = enestart[(s, 2020)]
end

# Calculate growth factors and potential GDP
function calculate_growth_factors()
    growth_factor = Dict()
    growth_factor[2020] = 1.0
    
    for i in 2:length(year_all)
        y_curr = year_all[i]
        y_prev = year_all[i-1]
        growth_factor[y_curr] = growth_factor[y_prev] * (1 + grow[y_curr])^duration_period
    end
    
    return growth_factor
end

growth_factor = calculate_growth_factors()

# Potential GDP
potential_gdp = Dict()
for y in year_all
    potential_gdp[y] = gdp_base * growth_factor[y]
end

# Calculate production function exponent
rho = (esub - 1) / esub

# Calculate cumulative AEEI factors
function calculate_aeei_factors()
    # Initialize base year
    for s in sector
        aeei_factor[(s, 2020)] = 1.0
    end
    
    # Calculate for subsequent years
    for i in 2:length(year_all)
        y_curr = year_all[i]
        y_prev = year_all[i-1]
        for s in sector
            aeei_factor[(s, y_curr)] = ((1 - aeei[(s, y_curr)])^duration_period) * aeei_factor[(s, y_prev)]
        end
    end
end

calculate_aeei_factors()

# Calculate labor supply and utility discount factors
function calculate_labor_and_discount()
    # Initialize base year
    udf[2020] = 1.0
    labor[2020] = 1.0
    
    # Calculate for subsequent years
    for i in 2:length(year_all)
        y_curr = year_all[i]
        y_prev = year_all[i-1]
        
        # Labor supply growth
        labor[y_curr] = labor[y_prev] * (1 + grow[y_curr])^duration_period
        
        # New labor supply
        newlab[y_curr] = max(labor[y_curr] - labor[y_prev] * (1 - depr)^duration_period, 0) + epsilon
        
        # Utility discount factor
        udf[y_curr] = udf[y_prev] * (1 - (drate - grow[y_curr]))^duration_period
    end
end

calculate_labor_and_discount()

# Calculate base year economic components
ecst0 = cost_MESSAGE[2020]
k0 = kgdp * gdp_base
i0 = max(k0 * (grow[2020] + depr), 0) + epsilon
c0 = gdp_base - i0 - ecst0
y0 = gdp_base

# Calculate production function coefficients
function calculate_production_coefficients()
    global a, b
    
    # Calculate coefficient b for energy
    b = price_MESSAGE[("NELE", 2020)] / ((1 - elvs) * y0^(1 - rho) * 
        demand_MESSAGE[("ELEC", 2020)]^(rho * elvs) * 
        demand_MESSAGE[("NELE", 2020)]^(rho * (1 - elvs) - 1))
    
    # Calculate coefficient a for capital-labor
    a = (y0^rho - b * demand_MESSAGE[("ELEC", 2020)]^(rho * elvs) * 
         demand_MESSAGE[("NELE", 2020)]^(rho * (1 - elvs))) / (k0^(kpvs * rho))
end

calculate_production_coefficients()

# Finite time horizon correction
for y in year_all
    finite_time_corr[y] = abs(drate - grow[y])
end

# Display key parameters (equivalent to GAMS DISPLAY statements)
println("Growth factors: ", growth_factor)
println("Potential GDP: ", potential_gdp)  
println("AEEI factors: ", aeei_factor)
println("Labor: ", labor)
println("New labor: ", newlab)
println("Utility discount factors: ", udf)
println("Base year values - ecst0: $ecst0, k0: $k0, i0: $i0, c0: $c0, y0: $y0")
println("Production coefficients - a: $a, b: $b")
println("Rho: $rho")