# macro_core.jl - Equivalent to macro/macro_core.gms
# MACRO core formulation with utility maximization and CES production

using JuMP

# ------------------------------------------------------------------------------
# MACRO model variables and equations
# ------------------------------------------------------------------------------

# Function to create MACRO model components
function create_macro_model!(model)
    # ------------------------------------------------------------------------------
    # Variable definitions
    # ------------------------------------------------------------------------------
    
    # POSITIVE VARIABLES
    @variable(model, K[year_all] >= 0)                    # Capital stock in period year
    @variable(model, KN[year_all] >= 0)                   # New Capital vintage in period year  
    @variable(model, Y[year_all] >= 0)                    # Production in period year
    @variable(model, YN[year_all] >= 0)                   # New production vintage in period year
    
    @variable(model, PRODENE[sector, year_all] >= 0)      # Value of end-use services in production function
    @variable(model, NEWENE[sector, year_all] >= 0)       # New end-use service (production function value)
    
    @variable(model, C[year_all] >= 0)                    # Consumption (Trillion $)
    @variable(model, I[year_all] >= 0)                    # Investment (Trillion $)
    
    # VARIABLES (unrestricted)
    @variable(model, UTILITY)                             # Utility function (discounted log of consumption)
    @variable(model, EC[year_all])                        # System costs (Trillion $) based on MESSAGE model run
    @variable(model, GDP_MACRO[year_all])                 # GDP for reporting when MESSAGE is run with MACRO
    
    # ------------------------------------------------------------------------------
    # Equation definitions
    # ------------------------------------------------------------------------------
    
    # UTILITY_FUNCTION
    # Utility function which is maximized sums up the discounted logarithm of consumption
    # over the entire time horizon of the model
    utility_terms = []
    
    # Regular periods (not base period, not last period)
    for y in year_all
        if y != 2020 && y != 2080  # not base period and not last period
            push!(utility_terms, udf[y] * log(C[y]) * duration_period)
        end
    end
    
    # Last period with finite time correction
    if 2080 in year_all
        push!(utility_terms, udf[2080] * log(C[2080]) * (duration_period + 1/finite_time_corr[2080]))
    end
    
    @constraint(model, UTILITY == sum(utility_terms))
    
    # CAPITAL_CONSTRAINT 
    # Allocation of total production among consumption, investment and energy system costs
    for y in year_all
        @constraint(model, Y[y] == C[y] + I[y] + EC[y])
    end
    
    # NEW_CAPITAL
    # Accumulation of capital - net capital formation from gross investments
    for y in year_all
        if y != 2020  # not base period
            @constraint(model, KN[y] == duration_period * I[y])
        end
    end
    
    # NEW_PRODUCTION  
    # Nested CES production function with capital, labor and energy services as inputs
    for y in year_all
        if y != 2020  # not base period
            @constraint(model,
                YN[y] == (a * KN[y]^(rho * kpvs) * newlab[y]^(rho * (1 - kpvs)) +
                         b * NEWENE["ELEC", y]^(rho * elvs) * NEWENE["NELE", y]^(rho * (1 - elvs)))^(1/rho)
            )
        end
    end
    
    # TOTAL_CAPITAL
    # Total capital stock is sum of depreciated previous capital and new capital
    for y in year_all
        if y != 2020  # not base period
            y_prev_idx = findfirst(==(y), year_all) - 1
            y_prev = year_all[y_prev_idx]
            @constraint(model, K[y] == K[y_prev] * (1 - depr)^duration_period + KN[y])
        end
    end
    
    # TOTAL_PRODUCTION
    # Total production is sum of depreciated previous production and new production
    for y in year_all
        if y != 2020  # not base period
            y_prev_idx = findfirst(==(y), year_all) - 1
            y_prev = year_all[y_prev_idx]
            @constraint(model, Y[y] == Y[y_prev] * (1 - depr)^duration_period + YN[y])
        end
    end
    
    # NEW_ENERGY
    # Total energy production is sum of depreciated previous energy and new energy
    for s in sector, y in year_all
        if y != 2020  # not base period
            y_prev_idx = findfirst(==(y), year_all) - 1
            y_prev = year_all[y_prev_idx]
            @constraint(model, 
                PRODENE[s, y] == PRODENE[s, y_prev] * (1 - depr)^duration_period + NEWENE[s, y]
            )
        end
    end
    
    # ENERGY_SUPPLY
    # Link between physical energy (PHYSENE) and monetary value (PRODENE) in production function
    for s in sector, y in year_all
        if y != 2020  # not base period
            @constraint(model, model[:PHYSENE][s, y] >= PRODENE[s, y] * aeei_factor[(s, y)])
        end
    end
    
    # COST_ENERGY
    # Energy system costs approximation based on MESSAGE model run (Taylor expansion)
    for y in year_all
        if y != 2020  # not base period
            @constraint(model,
                EC[y] == cost_MESSAGE[y] +
                sum(eneprice[(s, y)] * (model[:PHYSENE][s, y] - enestart[(s, y)]) for s in sector) +
                sum(eneprice[(s, y)] / enestart[(s, y)] * 
                    (model[:PHYSENE][s, y] - enestart[(s, y)]) * (model[:PHYSENE][s, y] - enestart[(s, y)]) 
                    for s in sector)
            )
        end
    end
    
    # COST_ENERGY_LINKED
    # Hard-linked version using actual energy system costs
    for y in year_all
        if y != 2020  # not base period
            @constraint(model, EC[y] == model[:COST_ANNUAL][y] / 1000)
        end
    end
    
    # TERMINAL_CONDITION
    # Terminal constraint to ensure appropriate investment levels at end of horizon
    if 2080 in year_all
        @constraint(model, I[2080] >= K[2080] * (grow[2080] + depr))
    end
    
    return model
end

# Function to set variable bounds and initial values (equivalent to macro_presolve.gms)
function set_macro_bounds_and_initial_values!(model)
    # Calculate start values for variables
    function calculate_start_values()
        SVKN_calc = Dict()
        SVNEWE_calc = Dict()
        
        for y in year_all
            if y != 2020  # not base period
                y_prev_idx = findfirst(==(y), year_all) - 1
                y_prev = year_all[y_prev_idx]
                
                # Start values for new capital
                SVKN_calc[y] = max((potential_gdp[y] - potential_gdp[y_prev] * (1 - depr)^duration_period) * kgdp, 0)
                
                # Start values for new energy
                for s in sector
                    duration_sum = (findfirst(==(y), year_all) - 1) * duration_period
                    SVNEWE_calc[(s, y)] = max(
                        demand_base[s] * growth_factor[y] - 
                        demand_base[s] * (1 - depr)^duration_sum, 0
                    )
                end
            end
        end
        
        return SVKN_calc, SVNEWE_calc
    end
    
    SVKN_calc, SVNEWE_calc = calculate_start_values()
    
    # Set initial values for variables
    for s in sector, y in year_all
        if y != 2020
            # Set starting values
            set_start_value(model[:NEWENE][s, y], max(SVNEWE_calc[(s, y)], epsilon))
        end
        # Physical energy starting values
        set_start_value(model[:PHYSENE][s, y], enestart[(s, y)])
    end
    
    for y in year_all
        if y != 2020
            set_start_value(model[:KN][y], max(SVKN_calc[y], epsilon))
        end
    end
    
    # Set lower bounds on variables to avoid singularities
    for y in year_all
        set_lower_bound(model[:K][y], lotol * k0)
        set_lower_bound(model[:Y][y], lotol * y0)
        set_lower_bound(model[:C][y], lotol * c0)
        set_lower_bound(model[:I][y], lotol * i0)
        
        if y != 2020
            set_lower_bound(model[:KN][y], lotol * i0 * duration_period)
            set_lower_bound(model[:YN][y], lotol * y0 * newlab[y])
        end
        
        for s in sector
            set_lower_bound(model[:PRODENE][s, y], lotol * enestart[(s, y)] / aeei_factor[(s, y)])
            if y != 2020
                set_lower_bound(model[:NEWENE][s, y], lotol * enestart[(s, y)] / aeei_factor[(s, y)])
            end
        end
    end
    
    # Fix base year values to historical values
    fix(model[:Y][2020], y0; force=true)
    fix(model[:K][2020], k0; force=true)
    fix(model[:C][2020], c0; force=true)
    fix(model[:I][2020], i0; force=true)
    fix(model[:EC][2020], y0 - i0 - c0; force=true)
    
    for s in sector
        fix(model[:PRODENE][s, 2020], demand_base[s] / aeei_factor[(s, 2020)]; force=true)
    end
    
    return model
end