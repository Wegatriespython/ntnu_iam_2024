# model_config.jl - User-definable configuration for the Energy-Macro IAM

using JuMP

"""
Configuration structure for the Energy-Macro Integrated Assessment Model.
All model parameters are centralized here for easy modification.
"""
Base.@kwdef struct ModelConfig
    # Time horizon parameters
    base_year::Int = 2020
    final_year::Int = 2080
    period_step::Int = 10
    duration_period::Int = 10  # MACRO model time discretization (years)
    period_length::Int = 10    # MESSAGE model time discretization (years)
    
    # Economic parameters (base year values)
    gdp_base::Float64 = 71.0      # Initial GDP (Trillion USD)
    kgdp::Float64 = 2.8           # Capital to GDP ratio
    depr::Float64 = 0.05          # Annual depreciation rate
    drate::Float64 = 0.05         # Social discount rate
    esub::Float64 = 0.3           # Elasticity of substitution (K-L vs Energy)
    
    # Production function value shares
    kpvs::Float64 = 0.28          # Capital value share parameter
    elvs::Float64 = 0.42          # Electricity value share parameter
    
    # Technical parameters
    epsilon::Float64 = 0.01       # Small number to avoid divergences
    lotol::Float64 = 0.05         # Tolerance factor for lower bounds
    discount_rate::Float64 = 0.05 # Discount rate for cost calculations
    
    # Growth rates by year (Dict{Int, Float64})
    growth_rates::Dict{Int, Float64} = Dict(
        2020 => 0.03, 2030 => 0.03, 2040 => 0.03, 2050 => 0.03,
        2060 => 0.03, 2070 => 0.02, 2080 => 0.01
    )
    
    # Autonomous Energy Efficiency Improvement (AEEI) rates
    aeei_rates::Dict{Tuple{String, Int}, Float64} = Dict(
        ("ELEC", 2020) => 0.02, ("ELEC", 2030) => 0.02, ("ELEC", 2040) => 0.02,
        ("ELEC", 2050) => 0.02, ("ELEC", 2060) => 0.02, ("ELEC", 2070) => 0.02, ("ELEC", 2080) => 0.02,
        ("NELE", 2020) => 0.02, ("NELE", 2030) => 0.02, ("NELE", 2040) => 0.02,
        ("NELE", 2050) => 0.02, ("NELE", 2060) => 0.02, ("NELE", 2070) => 0.02, ("NELE", 2080) => 0.02
    )
    
    # MESSAGE model data (can be overridden with external data)
    demand_message::Dict{Tuple{String, Int}, Float64} = Dict(
        ("ELEC", 2020) => 22.6, ("ELEC", 2030) => 25.7, ("ELEC", 2040) => 28.60,
        ("ELEC", 2050) => 32.80, ("ELEC", 2060) => 36.7, ("ELEC", 2070) => 41.7, ("ELEC", 2080) => 47.6,
        ("NELE", 2020) => 87.3, ("NELE", 2030) => 99.0, ("NELE", 2040) => 110.00,
        ("NELE", 2050) => 117.00, ("NELE", 2060) => 142.0, ("NELE", 2070) => 161.0, ("NELE", 2080) => 184.0
    )
    
    price_message::Dict{Tuple{String, Int}, Float64} = Dict(
        ("ELEC", 2020) => 0.0567, ("ELEC", 2030) => 0.0567, ("ELEC", 2040) => 0.0567,
        ("ELEC", 2050) => 0.0567, ("ELEC", 2060) => 0.0567, ("ELEC", 2070) => 0.0567, ("ELEC", 2080) => 0.0567,
        ("NELE", 2020) => 0.020, ("NELE", 2030) => 0.020, ("NELE", 2040) => 0.020,
        ("NELE", 2050) => 0.020, ("NELE", 2060) => 0.020, ("NELE", 2070) => 0.020, ("NELE", 2080) => 0.020
    )
    
    cost_message::Dict{Int, Float64} = Dict(
        2020 => 5.053, 2030 => 3.045, 2040 => 3.080, 2050 => 3.516,
        2060 => 3.852, 2070 => 4.376, 2080 => 4.996
    )
    
    # Solver parameters
    solver_print_level::Int = 5
    solver_max_iter::Int = 3000
    solver_tolerance::Float64 = 1e-6
    
    # GBD algorithm parameters
    gbd_max_iterations::Int = 50
    gbd_tolerance::Float64 = 1e-4
    gbd_relative_gap_tolerance::Float64 = 1e-3
end

"""
Generate year sequence based on configuration.
"""
function generate_year_sequence(config::ModelConfig)
    return collect(config.base_year:config.period_step:config.final_year)
end

"""
Calculate derived parameters from configuration.
"""
function calculate_derived_parameters(config::ModelConfig)
    year_all = generate_year_sequence(config)
    
    # Production function exponent
    rho = (config.esub - 1) / config.esub
    
    # Growth factors
    growth_factor = Dict{Int, Float64}()
    growth_factor[config.base_year] = 1.0
    
    for i in 2:length(year_all)
        y_curr = year_all[i]
        y_prev = year_all[i-1]
        growth_factor[y_curr] = growth_factor[y_prev] * (1 + config.growth_rates[y_curr])^config.duration_period
    end
    
    # Potential GDP
    potential_gdp = Dict(y => config.gdp_base * growth_factor[y] for y in year_all)
    
    # AEEI factors
    aeei_factor = Dict{Tuple{String, Int}, Float64}()
    for s in ["ELEC", "NELE"]
        aeei_factor[(s, config.base_year)] = 1.0
    end
    
    for i in 2:length(year_all)
        y_curr = year_all[i]
        y_prev = year_all[i-1]
        for s in ["ELEC", "NELE"]
            aeei_factor[(s, y_curr)] = ((1 - config.aeei_rates[(s, y_curr)])^config.duration_period) * aeei_factor[(s, y_prev)]
        end
    end
    
    # Labor and utility discount factors
    udf = Dict{Int, Float64}()
    labor = Dict{Int, Float64}()
    newlab = Dict{Int, Float64}()
    
    udf[config.base_year] = 1.0
    labor[config.base_year] = 1.0
    
    for i in 2:length(year_all)
        y_curr = year_all[i]
        y_prev = year_all[i-1]
        
        labor[y_curr] = labor[y_prev] * (1 + config.growth_rates[y_curr])^config.duration_period
        newlab[y_curr] = max(labor[y_curr] - labor[y_prev] * (1 - config.depr)^config.duration_period, 0) + config.epsilon
        udf[y_curr] = udf[y_prev] * (1 - (config.drate - config.growth_rates[y_curr]))^config.duration_period
    end
    
    # Base year economic components
    ecst0 = config.cost_message[config.base_year]
    k0 = config.kgdp * config.gdp_base
    i0 = max(k0 * (config.growth_rates[config.base_year] + config.depr), 0) + config.epsilon
    c0 = config.gdp_base - i0 - ecst0
    y0 = config.gdp_base
    
    # Production function coefficients
    b = config.price_message[("NELE", config.base_year)] / ((1 - config.elvs) * y0^(1 - rho) * 
        config.demand_message[("ELEC", config.base_year)]^(rho * config.elvs) * 
        config.demand_message[("NELE", config.base_year)]^(rho * (1 - config.elvs) - 1))
    
    a = (y0^rho - b * config.demand_message[("ELEC", config.base_year)]^(rho * config.elvs) * 
         config.demand_message[("NELE", config.base_year)]^(rho * (1 - config.elvs))) / (k0^(config.kpvs * rho))
    
    # Finite time horizon correction
    finite_time_corr = Dict(y => abs(config.drate - config.growth_rates[y]) for y in year_all)
    
    # Demand base
    demand_base = Dict(s => config.demand_message[(s, config.base_year)] for s in ["ELEC", "NELE"])
    
    return (
        year_all = year_all,
        rho = rho,
        growth_factor = growth_factor,
        potential_gdp = potential_gdp,
        aeei_factor = aeei_factor,
        udf = udf,
        labor = labor,
        newlab = newlab,
        ecst0 = ecst0,
        k0 = k0,
        i0 = i0,
        c0 = c0,
        y0 = y0,
        a = a,
        b = b,
        finite_time_corr = finite_time_corr,
        demand_base = demand_base
    )
end

"""
Create a default configuration instance.
"""
function default_config()
    return ModelConfig()
end

"""
Create a custom configuration with user overrides.
Example usage:
    config = custom_config(
        gdp_base = 80.0,
        drate = 0.04,
        growth_rates = Dict(2020 => 0.025, 2030 => 0.025, ...)
    )
"""
function custom_config(; kwargs...)
    return ModelConfig(; kwargs...)
end

# Export main symbols
export ModelConfig, default_config, custom_config
export generate_year_sequence, calculate_derived_parameters