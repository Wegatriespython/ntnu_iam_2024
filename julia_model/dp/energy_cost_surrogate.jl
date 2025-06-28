# energy_cost_surrogate_julia.jl - Pure Julia implementation of energy cost surrogate
# Uses trained coefficients from AR(2) and time-aware models

# Energy cost surrogate models
struct AR2Surrogate
    # Model: cost = β₀ + β₁*cost_lag1 + β₂*cost_lag2 + β₃*elec_demand + β₄*nele_demand
    β::Vector{Float64}  # Coefficients [intercept, cost_lag1, cost_lag2, elec, nele]
end

struct TimeAwareSurrogate
    # Model: cost = β₀ + β₁*cost_lag1 + β₂*elec + β₃*nele + β₄*time_idx + β₅*demand_time
    β::Vector{Float64}  # Coefficients [intercept, cost_lag1, elec, nele, time, demand*time]
end

# Create surrogates with trained coefficients
function create_ar2_surrogate()
    # From Python model output
    coeffs = [-0.133469, 0.0348350, 113.871257, 19.1798985]
    intercept = 59.41852
    
    return AR2Surrogate([intercept; coeffs])
end

function create_time_aware_surrogate()
    # From Python model output
    coeffs = [-0.0282839, 106.980967, 17.4960798, -58.0962599, 0.157242743]
    intercept = 268.27843
    
    return TimeAwareSurrogate([intercept; coeffs])
end

# Get time-varying demand with GDP growth
function get_demand_projection(base_demand::Float64, year::Int; β::Float64 = 0.7)
    # GDP growth factors (from macro model)
    gdp_growth = Dict(
        2020 => 1.000,
        2030 => 1.344,
        2040 => 1.806,
        2050 => 2.427,
        2060 => 3.262,
        2070 => 4.022,
        2080 => 4.435
    )
    
    return base_demand * gdp_growth[year]^β
end

# Predict using AR(2) model
function predict_ar2(surrogate::AR2Surrogate, cost_lag1::Float64, cost_lag2::Float64,
                    elec_demand::Float64, nele_demand::Float64, year::Int)
    
    # Get time-varying demands
    elec_proj = get_demand_projection(elec_demand, year)
    nele_proj = get_demand_projection(nele_demand, year)
    
    # Linear prediction (in billion USD)
    cost_billion = surrogate.β[1] + 
                   surrogate.β[2] * cost_lag1 +
                   surrogate.β[3] * cost_lag2 +
                   surrogate.β[4] * elec_proj +
                   surrogate.β[5] * nele_proj
    
    # Convert to trillion USD
    return cost_billion / 1000.0
end

# Predict using time-aware model
function predict_time_aware(surrogate::TimeAwareSurrogate, cost_lag1::Float64,
                           elec_demand::Float64, nele_demand::Float64, 
                           year::Int, time_index::Int)
    
    # Get time-varying demands
    elec_proj = get_demand_projection(elec_demand, year)
    nele_proj = get_demand_projection(nele_demand, year)
    total_demand = elec_proj + nele_proj
    
    # Linear prediction (in billion USD)
    cost_billion = surrogate.β[1] + 
                   surrogate.β[2] * cost_lag1 +
                   surrogate.β[3] * elec_proj +
                   surrogate.β[4] * nele_proj +
                   surrogate.β[5] * time_index +
                   surrogate.β[6] * (total_demand * time_index)
    
    # Convert to trillion USD
    return cost_billion / 1000.0
end

# Main energy cost function for DP solver
function energy_cost_surrogate(elec_demand::Float64, nele_demand::Float64, year::Int;
                              model_type::Symbol = :ar2,
                              cost_lag1::Union{Float64,Nothing} = nothing,
                              cost_lag2::Union{Float64,Nothing} = nothing)
    
    # Get time index
    years = [2020, 2030, 2040, 2050, 2060, 2070, 2080]
    time_index = findfirst(==(year), years) - 1  # 0-indexed
    
    # Handle missing lags
    if isnothing(cost_lag1)
        # Use typical cost trajectory (billion USD)
        cost_estimates = Dict(
            2020 => 9264.0,
            2030 => 5268.0,
            2040 => 4961.0,
            2050 => 5671.0,
            2060 => 6266.0,
            2070 => 7119.0,
            2080 => 8128.0
        )
        
        year_idx = findfirst(==(year), years)
        
        if year_idx === nothing || year_idx == 1
            # Base year - use MESSAGE cost
            return 5.053  # Trillion USD
        elseif year_idx == 2
            # 2030 - only one lag
            cost_lag1 = cost_estimates[2020]
            if model_type == :ar2
                cost_lag2 = cost_estimates[2020]  # Use same value
            end
        else
            # Full lags available
            cost_lag1 = cost_estimates[years[year_idx - 1]]
            if model_type == :ar2
                cost_lag2 = cost_estimates[years[year_idx - 2]]
            end
        end
    end
    
    # Predict based on model type
    if model_type == :ar2 && !isnothing(cost_lag2)
        surrogate = create_ar2_surrogate()
        return predict_ar2(surrogate, cost_lag1, cost_lag2, elec_demand, nele_demand, year)
    else
        # Use time-aware model
        surrogate = create_time_aware_surrogate()
        return predict_time_aware(surrogate, cost_lag1, elec_demand, nele_demand, year, time_index)
    end
end

# Predict full trajectory
function predict_cost_trajectory(elec_path::Vector{Float64}, nele_path::Vector{Float64},
                                years::Vector{Int}; model_type::Symbol = :ar2)
    
    T = length(years)
    costs = zeros(T)
    
    # Track lagged costs (in billion USD)
    cost_lag1 = 9264.0  # 2020 cost
    cost_lag2 = 9264.0  # No 2010 data
    
    for t in 1:T
        if t == 1
            # Base year - fixed cost
            costs[t] = 5.053  # Trillion USD
        else
            # Predict using surrogate
            costs[t] = energy_cost_surrogate(
                elec_path[t], nele_path[t], years[t],
                model_type = model_type,
                cost_lag1 = cost_lag1,
                cost_lag2 = (model_type == :ar2 ? cost_lag2 : nothing)
            )
            
            # Update lags (convert back to billion)
            if model_type == :ar2 && t > 1
                cost_lag2 = cost_lag1
            end
            cost_lag1 = costs[t] * 1000.0
        end
    end
    
    return costs
end

# Simple linear approximation (fallback)
function simple_energy_cost(elec_demand::Float64, nele_demand::Float64, year::Int)
    # Linear approximation based on MESSAGE data
    base_cost = Dict(
        2020 => 5.053,
        2030 => 3.045,
        2040 => 3.080,
        2050 => 3.516,
        2060 => 3.852,
        2070 => 4.376,
        2080 => 4.996
    )
    
    # Reference demands
    ref_elec = Dict(
        2020 => 22.6, 2030 => 25.7, 2040 => 28.60,
        2050 => 32.80, 2060 => 36.7, 2070 => 41.7, 2080 => 47.6
    )
    
    ref_nele = Dict(
        2020 => 87.3, 2030 => 99.0, 2040 => 110.00,
        2050 => 117.00, 2060 => 142.0, 2070 => 161.0, 2080 => 184.0
    )
    
    # Get time-varying demands
    elec_proj = get_demand_projection(elec_demand, year)
    nele_proj = get_demand_projection(nele_demand, year)
    
    # Linear approximation with sensitivities
    base = base_cost[year]
    elec_sensitivity = 0.116  # ~116 billion per PWh in trillion = 0.116
    nele_sensitivity = 0.019  # ~19 billion per PWh in trillion = 0.019
    
    cost = base + elec_sensitivity * (elec_proj - ref_elec[year]) + 
                  nele_sensitivity * (nele_proj - ref_nele[year])
    
    return max(cost, 0.5)  # Ensure positive cost
end