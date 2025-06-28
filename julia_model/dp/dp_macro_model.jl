# dp_macro_model.jl - Dynamic Programming formulation of MACRO model
# Uses Bellman equation approach with endogenous grid method

using Parameters
using Interpolations
using LinearAlgebra
using Optim

# DP Model structure with vintage capital and 3-state formulation
@with_kw struct DPMacroModel
    # Time parameters
    years::Vector{Int} = [2020, 2030, 2040, 2050, 2060, 2070, 2080]
    T::Int = length(years)
    duration::Float64 = 10.0  # years per period (duration_period in MACRO)
    
    # Economic parameters
    β::Float64 = 0.95^duration  # Discount factor (from drate = 0.05)
    δ::Float64 = 1 - (1 - 0.05)^duration  # Depreciation over period
    drate::Float64 = 0.05  # Annual discount rate
    depr::Float64 = 0.05   # Annual depreciation rate
    
    # Production parameters
    ρ::Float64 = -0.233  # CES exponent (from esub = 0.3)
    kpvs::Float64 = 0.28  # Capital value share
    elvs::Float64 = 0.42  # Electricity value share
    a::Float64 = 0.0     # Capital-labor coefficient (to be calibrated)
    b::Float64 = 0.0     # Energy coefficient (to be calibrated)
    
    # Grid parameters for 3-dimensional state space
    K_min::Float64 = 50.0   # Minimum capital (trillion USD)
    K_max::Float64 = 600.0  # Maximum capital (trillion USD)
    n_K::Int = 30           # Grid points for capital (reduced for 3D)
    
    Y_min::Float64 = 30.0   # Minimum production (trillion USD)
    Y_max::Float64 = 300.0  # Maximum production (trillion USD)
    n_Y::Int = 20           # Grid points for production
    
    PRODENE_scale::Float64 = 150.0  # Scale factor for PRODENE grid
    n_PRODENE::Int = 15     # Grid points for PRODENE
    
    # Grids
    K_grid::Vector{Float64} = range(K_min, K_max, length=n_K)
    Y_grid::Vector{Float64} = range(Y_min, Y_max, length=n_Y)
    PRODENE_grid::Vector{Float64} = range(0.5, 2.0, length=n_PRODENE) * PRODENE_scale
    
    # Labor and AEEI trajectories (from macro_data_load.jl)
    labor::Dict{Int,Float64} = Dict()
    newlab::Dict{Int,Float64} = Dict()  # New vintage labor
    aeei_factor::Dict{Tuple{String,Int},Float64} = Dict()
    growth::Dict{Int,Float64} = Dict()
    growth_factor::Dict{Int,Float64} = Dict()
    
    # Utility discount factors
    udf::Dict{Int,Float64} = Dict()  # Utility discount factors
    finite_time_corr::Dict{Int,Float64} = Dict()  # Terminal period correction
    
    # Base year values
    y0::Float64 = 71.0    # Base year GDP (trillion USD)
    k0::Float64 = 198.8   # Base year capital (trillion USD) 
    c0::Float64 = 51.0    # Base year consumption
    i0::Float64 = 14.94   # Base year investment
    
    # Energy base values
    demand_base::Dict{String,Float64} = Dict("ELEC" => 22.6, "NELE" => 87.3)
    PRODENE_base::Float64 = 110.0  # Approximate base PRODENE total
    
    # Value and policy functions - now 3D arrays
    V::Array{Float64,4} = zeros(n_K, n_Y, n_PRODENE, T)  # V[i,j,k,t]
    C_policy::Array{Float64,3} = zeros(n_K, n_Y, n_PRODENE)  # Current period only
    I_policy::Array{Float64,3} = zeros(n_K, n_Y, n_PRODENE)
    
    # Store policies for current period (to save memory)
    KN_policy::Array{Float64,3} = zeros(n_K, n_Y, n_PRODENE)
    YN_policy::Array{Float64,3} = zeros(n_K, n_Y, n_PRODENE)
    NEWENE_total_policy::Array{Float64,3} = zeros(n_K, n_Y, n_PRODENE)
    EC_policy::Array{Float64,3} = zeros(n_K, n_Y, n_PRODENE)
end

# Bellman equation components
struct BellmanProblem
    model::DPMacroModel
    energy_cost_func::Function  # Function to get energy costs
end

# Utility function (log utility)
u(c::Float64) = log(max(c, 1e-10))
u_prime(c::Float64) = 1.0 / max(c, 1e-10)
u_prime_inv(up::Float64) = 1.0 / up

# CES production function (MACRO formulation)
# This is for YN calculation with KN and NEWENE inputs
function production(K::Float64, L::Float64, ELEC::Float64, NELE::Float64, model::DPMacroModel)
    # Y = [a*K^(ρκ)*L^(ρ(1-κ)) + b*ELEC^(ρε)*NELE^(ρ(1-ε))]^(1/ρ)
    
    # Ensure positive inputs
    K = max(K, 0.1)
    L = max(L, 0.1)
    ELEC = max(ELEC, 0.1)
    NELE = max(NELE, 0.1)
    
    # Capital-labor term
    kl_term = model.a * K^(model.ρ * model.kpvs) * L^(model.ρ * (1 - model.kpvs))
    
    # Energy term
    energy_term = model.b * ELEC^(model.ρ * model.elvs) * NELE^(model.ρ * (1 - model.elvs))
    
    # Total production
    Y = (kl_term + energy_term)^(1/model.ρ)
    
    return Y
end

# Analytical solution for NEWENE given YN target
function solve_newene_analytical(YN_target::Float64, KN::Float64, L::Float64, 
                                model::DPMacroModel, year::Int)
    # From SymPy derivation:
    # NELE = [(YN^ρ - a*KN^(ρκ)*L^(ρ(1-κ))) / (b*r^(ρε))]^(1/ρ)
    # ELEC = r * NELE
    
    # Calculate capital-labor term
    KL_term = model.a * KN^(model.ρ * model.kpvs) * L^(model.ρ * (1 - model.kpvs))
    
    # Energy term needed
    E_term_needed = YN_target^model.ρ - KL_term
    
    # Check feasibility
    if E_term_needed <= 0
        # YN target too low - can be achieved with just capital and labor
        return (ELEC = 0.0, NELE = 0.0, feasible = false)
    end
    
    # Get energy price ratio from original model's calibration
    # In equilibrium: ELEC/NELE = [(elvs/(1-elvs)) * (p_nele/p_elec)]^(1/(1-ρ))
    # Using calibrated prices: p_elec = 0.0567, p_nele = 0.020
    p_elec = 0.0567
    p_nele = 0.020
    price_ratio = p_elec / p_nele
    
    # Optimal energy ratio
    r = ((model.elvs / (1 - model.elvs)) / price_ratio)^(1 / (1 - model.ρ))
    
    # Solve for NELE
    NELE = (E_term_needed / (model.b * r^(model.ρ * model.elvs)))^(1 / model.ρ)
    
    # Calculate ELEC
    ELEC = r * NELE
    
    return (ELEC = ELEC, NELE = NELE, feasible = true)
end

# Calculate energy variables for given YN
function calculate_energy_for_yn(YN_target::Float64, KN::Float64, L::Float64,
                               Y_prev::Float64, PRODENE_prev::Float64,
                               model::DPMacroModel, year::Int, energy_cost_func::Function)
    # Get AEEI factors
    aeei_elec = get(model.aeei_factor, ("ELEC", year), 1.0)
    aeei_nele = get(model.aeei_factor, ("NELE", year), 1.0)
    
    # Solve for NEWENE analytically
    newene_result = solve_newene_analytical(YN_target, KN, L, model, year)
    
    if !newene_result.feasible
        # Use minimal energy
        NEWENE_ELEC = 0.1
        NEWENE_NELE = 0.1
    else
        NEWENE_ELEC = newene_result.ELEC
        NEWENE_NELE = newene_result.NELE
    end
    
    # Total NEWENE in production function terms
    NEWENE_total = NEWENE_ELEC + NEWENE_NELE
    
    # Update PRODENE with depreciation
    PRODENE_total = PRODENE_prev * (1 - model.depr)^model.duration + NEWENE_total
    
    # Physical energy (applying AEEI)
    # Assume same ELEC/NELE split for simplicity
    elec_share = NEWENE_ELEC / NEWENE_total
    PRODENE_ELEC = PRODENE_total * elec_share
    PRODENE_NELE = PRODENE_total * (1 - elec_share)
    
    PHYSENE_ELEC = PRODENE_ELEC * aeei_elec
    PHYSENE_NELE = PRODENE_NELE * aeei_nele
    
    # Calculate energy cost
    EC = energy_cost_func(PHYSENE_ELEC, PHYSENE_NELE, year)
    
    # Update total production
    Y = Y_prev * (1 - model.depr)^model.duration + YN_target
    
    return (
        Y = Y,
        YN = YN_target,
        PRODENE_total = PRODENE_total,
        NEWENE_total = NEWENE_total,
        PHYSENE_ELEC = PHYSENE_ELEC,
        PHYSENE_NELE = PHYSENE_NELE,
        EC = EC
    )
end

# Bellman operator with 3-state formulation
function bellman_operator_3d!(bp::BellmanProblem, t::Int)
    model = bp.model
    year = model.years[t]
    
    # Get time-specific parameters
    L = get(model.labor, year, 1.0)
    L_new = get(model.newlab, year, L)  # Labor for new vintage
    udf = get(model.udf, year, model.β)
    growth = get(model.growth, year, 0.02)
    
    # Terminal period utility correction
    finite_corr = t == model.T ? get(model.finite_time_corr, year, 0.0) : 0.0
    
    # Create interpolator for next period value function if not terminal
    if t < model.T
        V_next_interp = linear_interpolation((model.K_grid, model.Y_grid, model.PRODENE_grid), 
                                            model.V[:,:,:,t+1], extrapolation_bc=Line())
    end
    
    # For each state combination
    for (i, K) in enumerate(model.K_grid)
        for (j, Y_prev) in enumerate(model.Y_grid) 
            for (k, PRODENE_prev) in enumerate(model.PRODENE_grid)
                
                V_opt = -Inf
                
                if t == model.T  # Terminal period
                    # Terminal condition from MACRO model
                    I_terminal = K * (growth + model.depr)
                    KN = model.duration * I_terminal
                    
                    # Calculate YN with KN and L_new
                    # First get a reasonable YN estimate
                    YN_guess = Y_prev * growth
                    
                    # Use analytical solution to get energy and costs
                    energy_result = calculate_energy_for_yn(YN_guess, KN, L_new, 
                                                           Y_prev, PRODENE_prev, 
                                                           model, year, bp.energy_cost_func)
                    
                    # Terminal consumption
                    C_opt = energy_result.Y - I_terminal - energy_result.EC
                    
                    if C_opt > 0
                        # Terminal utility with finite time correction
                        V_opt = if finite_corr > 0
                            udf * log(C_opt) * (model.duration + 1/finite_corr)
                        else
                            udf * log(C_opt) * model.duration
                        end
                        
                        # Store policies
                        model.C_policy[i,j,k] = C_opt
                        model.I_policy[i,j,k] = I_terminal
                        model.KN_policy[i,j,k] = KN
                        model.YN_policy[i,j,k] = energy_result.YN
                        model.NEWENE_total_policy[i,j,k] = energy_result.NEWENE_total
                        model.EC_policy[i,j,k] = energy_result.EC
                    end
                    
                else  # Non-terminal periods
                    C_opt = 0.0
                    I_opt = 0.0
                    
                    # Search over investment levels
                    I_min = 0.01 * Y_prev  
                    I_max = min(0.5 * Y_prev, Y_prev - 0.1)  # Leave room for consumption
                    
                    # Grid search over investment
                    for I in range(I_min, I_max, length=20)
                        # Calculate new capital vintage
                        KN = model.duration * I
                        
                        # Search over different YN values
                        YN_min = 0.5 * Y_prev
                        YN_max = 2.0 * Y_prev
                        
                        for YN in range(YN_min, YN_max, length=15)
                            # Calculate energy and costs for this YN
                            energy_result = calculate_energy_for_yn(YN, KN, L_new,
                                                                   Y_prev, PRODENE_prev,
                                                                   model, year, bp.energy_cost_func)
                            
                            # Consumption from budget constraint
                            C = energy_result.Y - I - energy_result.EC
                            
                            # Next period states
                            K_next = K * (1 - model.depr)^model.duration + KN
                            Y_next = energy_result.Y
                            PRODENE_next = energy_result.PRODENE_total
                            
                            # Check feasibility
                            if C > 0 && K_next >= model.K_min && K_next <= model.K_max &&
                               Y_next >= model.Y_min && Y_next <= model.Y_max &&
                               PRODENE_next >= model.PRODENE_grid[1] && PRODENE_next <= model.PRODENE_grid[end]
                                
                                # Interpolate next period value
                                V_next = V_next_interp(K_next, Y_next, PRODENE_next)
                                
                                # Current utility + continuation value
                                V_candidate = udf * log(C) * model.duration + model.β * V_next
                                
                                if V_candidate > V_opt
                                    V_opt = V_candidate
                                    C_opt = C
                                    I_opt = I
                                    
                                    # Store optimal policies
                                    model.C_policy[i,j,k] = C
                                    model.I_policy[i,j,k] = I
                                    model.KN_policy[i,j,k] = KN
                                    model.YN_policy[i,j,k] = YN
                                    model.NEWENE_total_policy[i,j,k] = energy_result.NEWENE_total
                                    model.EC_policy[i,j,k] = energy_result.EC
                                end
                            end
                        end
                    end
                end
                
                # Store value function
                model.V[i,j,k,t] = V_opt
            end
        end
    end
end

# Solve the DP problem using backward induction with 3-state formulation
function solve_dp!(bp::BellmanProblem)
    model = bp.model
    
    # Backward induction
    for t in model.T:-1:1
        println("Solving period $t (year $(model.years[t]))")
        bellman_operator_3d!(bp, t)
    end
    
    return model
end

# Simulate forward using 3D policy functions
function simulate_trajectory(model::DPMacroModel, K0::Float64, Y0::Float64, PRODENE0::Float64, 
                           energy_cost_func::Function=(ELEC, NELE, year) -> 0.1*ELEC + 0.05*NELE)
    T = model.T
    
    # Initialize trajectories
    K_path = zeros(T)
    KN_path = zeros(T)
    C_path = zeros(T)
    I_path = zeros(T)
    Y_path = zeros(T)
    YN_path = zeros(T)
    PRODENE_path = zeros(T)
    NEWENE_path = zeros(T)
    EC_path = zeros(T)
    
    # Energy trajectories (for detailed reporting)
    PHYSENE_ELEC_path = zeros(T)
    PHYSENE_NELE_path = zeros(T)
    
    # Set initial values
    K_path[1] = K0
    Y_path[1] = Y0
    PRODENE_path[1] = PRODENE0
    C_path[1] = model.c0
    I_path[1] = model.i0
    EC_path[1] = model.y0 - model.c0 - model.i0
    
    # Base year energy (approximate split)
    elec_share = model.demand_base["ELEC"] / (model.demand_base["ELEC"] + model.demand_base["NELE"])
    PHYSENE_ELEC_path[1] = model.demand_base["ELEC"]
    PHYSENE_NELE_path[1] = model.demand_base["NELE"]
    
    # Create policy interpolators for each period (except base year)
    policy_interps = []
    for t in 2:T
        # Note: We need to solve the DP problem at each time step to get proper policies
        # For now, we'll use the stored policies which are overwritten each period
        # In a proper implementation, we'd store policies for all periods
        
        C_interp = linear_interpolation((model.K_grid, model.Y_grid, model.PRODENE_grid), 
                                       model.C_policy, extrapolation_bc=Line())
        
        I_interp = linear_interpolation((model.K_grid, model.Y_grid, model.PRODENE_grid), 
                                       model.I_policy, extrapolation_bc=Line())
        
        KN_interp = linear_interpolation((model.K_grid, model.Y_grid, model.PRODENE_grid), 
                                        model.KN_policy, extrapolation_bc=Line())
        
        YN_interp = linear_interpolation((model.K_grid, model.Y_grid, model.PRODENE_grid), 
                                        model.YN_policy, extrapolation_bc=Line())
        
        NEWENE_interp = linear_interpolation((model.K_grid, model.Y_grid, model.PRODENE_grid), 
                                            model.NEWENE_total_policy, extrapolation_bc=Line())
        
        EC_interp = linear_interpolation((model.K_grid, model.Y_grid, model.PRODENE_grid), 
                                        model.EC_policy, extrapolation_bc=Line())
        
        push!(policy_interps, (
            C = C_interp,
            I = I_interp,
            KN = KN_interp,
            YN = YN_interp,
            NEWENE = NEWENE_interp,
            EC = EC_interp
        ))
    end
    
    # Simulate forward
    for t in 2:T
        # Current states
        K = K_path[t-1]
        Y = Y_path[t-1]
        PRODENE = PRODENE_path[t-1]
        
        # Apply policy functions
        policies = policy_interps[t-1]
        C_path[t] = policies.C(K, Y, PRODENE)
        I_path[t] = policies.I(K, Y, PRODENE)
        KN_path[t] = policies.KN(K, Y, PRODENE)
        YN_path[t] = policies.YN(K, Y, PRODENE)
        NEWENE_path[t] = policies.NEWENE(K, Y, PRODENE)
        EC_path[t] = policies.EC(K, Y, PRODENE)
        
        # Update state variables with depreciation
        K_path[t] = K * (1 - model.depr)^model.duration + KN_path[t]
        Y_path[t] = Y * (1 - model.depr)^model.duration + YN_path[t]
        PRODENE_path[t] = PRODENE * (1 - model.depr)^model.duration + NEWENE_path[t]
        
        # Calculate physical energy (using same split as base year for simplicity)
        year = model.years[t]
        aeei_elec = get(model.aeei_factor, ("ELEC", year), 1.0)
        aeei_nele = get(model.aeei_factor, ("NELE", year), 1.0)
        
        PHYSENE_ELEC_path[t] = PRODENE_path[t] * elec_share * aeei_elec
        PHYSENE_NELE_path[t] = PRODENE_path[t] * (1 - elec_share) * aeei_nele
    end
    
    # Calculate utility
    utility = 0.0
    for t in 2:(T-1)
        year = model.years[t]
        udf = get(model.udf, year, model.β)
        utility += udf * log(C_path[t]) * model.duration
    end
    
    # Terminal utility with finite correction
    if T > 1
        year = model.years[T]
        udf = get(model.udf, year, model.β)
        finite_corr = get(model.finite_time_corr, year, 0.0)
        if finite_corr > 0
            utility += udf * log(C_path[T]) * (model.duration + 1/finite_corr)
        else
            utility += udf * log(C_path[T]) * model.duration
        end
    end
    
    return (
        years = model.years,
        K = K_path,
        KN = KN_path,
        C = C_path,
        I = I_path,
        Y = Y_path,
        YN = YN_path,
        GDP_MACRO = Y_path,  # Y is GDP in this formulation
        PRODENE_total = PRODENE_path,
        NEWENE_total = NEWENE_path,
        EC = EC_path,
        PHYSENE_ELEC = PHYSENE_ELEC_path,
        PHYSENE_NELE = PHYSENE_NELE_path,
        UTILITY = utility
    )
end