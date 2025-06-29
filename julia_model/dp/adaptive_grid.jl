"""
Adaptive grid refinement for dynamic programming MACRO model.
Concentrates grid points where value function curvature is highest.
"""

using LinearAlgebra
using Statistics

"""
    compute_value_curvature(V, grid_K, grid_Y, grid_P)

Compute local curvature of value function using second derivatives.
Returns curvature measure for each grid point.
"""
function compute_value_curvature(V::Array{Float64, 4}, grid_K, grid_Y, grid_P, t::Int)
    n_K, n_Y, n_P = length(grid_K), length(grid_Y), length(grid_P)
    curvature = zeros(n_K, n_Y, n_P)
    
    # Compute second derivatives using finite differences
    for i in 2:n_K-1, j in 2:n_Y-1, k in 2:n_P-1
        # Second derivatives in each dimension
        d2V_dK2 = (V[i+1,j,k,t] - 2*V[i,j,k,t] + V[i-1,j,k,t]) / ((grid_K[i+1] - grid_K[i-1])/2)^2
        d2V_dY2 = (V[i,j+1,k,t] - 2*V[i,j,k,t] + V[i,j-1,k,t]) / ((grid_Y[j+1] - grid_Y[j-1])/2)^2
        d2V_dP2 = (V[i,j,k+1,t] - 2*V[i,j,k,t] + V[i,j,k-1,t]) / ((grid_P[k+1] - grid_P[k-1])/2)^2
        
        # Mixed derivatives
        d2V_dKdY = ((V[i+1,j+1,k,t] - V[i+1,j-1,k,t]) - (V[i-1,j+1,k,t] - V[i-1,j-1,k,t])) / 
                   (4 * (grid_K[i+1] - grid_K[i-1])/2 * (grid_Y[j+1] - grid_Y[j-1])/2)
        
        # Total curvature (Frobenius norm of Hessian)
        curvature[i,j,k] = sqrt(d2V_dK2^2 + d2V_dY2^2 + d2V_dP2^2 + 2*d2V_dKdY^2)
    end
    
    return curvature
end

"""
    adaptive_grid_refinement(model, energy_cost_func, initial_solve_periods=3)

Perform adaptive grid refinement based on value function curvature.
"""
function adaptive_grid_refinement(model::DPMacroModel, energy_cost_func::Function, initial_solve_periods::Int=3)
    println("\n=== ADAPTIVE GRID REFINEMENT ===")
    
    # Step 1: Solve with coarse grid for initial periods
    println("Step 1: Initial solve with coarse grid...")
    bp = BellmanProblem(model, energy_cost_func)
    solve_dp!(bp)
    
    # Step 2: Analyze value function curvature
    println("Step 2: Analyzing value function curvature...")
    total_curvature = zeros(model.n_K, model.n_Y, model.n_PRODENE)
    
    # Average curvature across solved periods
    for t in 1:min(initial_solve_periods, length(model.years))
        curvature_t = compute_value_curvature(model.V, model.grid_K_vec, model.grid_Y_vec, model.grid_PRODENE_vec, t)
        total_curvature .+= curvature_t
    end
    total_curvature ./= min(initial_solve_periods, length(model.years))
    
    # Step 3: Identify high-curvature regions
    positive_curvature = vec(total_curvature[total_curvature .> 0])
    
    if isempty(positive_curvature)
        println("Warning: No positive curvature values found. Using uniform refinement.")
        # Use uniform refinement pattern
        high_curvature_mask = trues(size(total_curvature))
        curvature_threshold = 0.0
    else
        curvature_threshold = quantile(positive_curvature, 0.75)
        high_curvature_mask = total_curvature .> curvature_threshold
        
        println("Curvature statistics:")
        println("  Mean curvature: $(mean(positive_curvature))")
        println("  75th percentile: $(curvature_threshold)")
        println("  High curvature points: $(sum(high_curvature_mask)) / $(length(high_curvature_mask))")
    end
    
    # Step 4: Create refined grid
    refined_model = create_refined_grid(model, Array{Bool,3}(high_curvature_mask), 1.5)
    
    # Step 5: Solve with refined grid
    println("Step 3: Solving with refined grid...")
    refined_bp = BellmanProblem(refined_model, energy_cost_func)
    solve_dp!(refined_bp)
    
    return refined_model, total_curvature
end

"""
    create_refined_grid(model, high_curvature_mask, refinement_factor=1.5)

Create a new model with refined grid in high-curvature regions.
"""
function create_refined_grid(model::DPMacroModel, high_curvature_mask::Array{Bool, 3}, refinement_factor::Float64=1.5)
    # Identify regions needing refinement for each dimension
    needs_K_refinement = any(high_curvature_mask, dims=(2,3))[:, 1, 1]
    needs_Y_refinement = any(high_curvature_mask, dims=(1,3))[1, :, 1]  
    needs_P_refinement = any(high_curvature_mask, dims=(1,2))[1, 1, :]
    
    # Create refined grids
    new_grid_K = create_adaptive_1d_grid(model.grid_K_vec, needs_K_refinement, refinement_factor)
    new_grid_Y = create_adaptive_1d_grid(model.grid_Y_vec, needs_Y_refinement, refinement_factor)
    new_grid_P = create_adaptive_1d_grid(model.grid_PRODENE_vec, needs_P_refinement, refinement_factor)
    
    println("Grid refinement:")
    println("  K: $(length(model.grid_K_vec)) → $(length(new_grid_K)) points")
    println("  Y: $(length(model.grid_Y_vec)) → $(length(new_grid_Y)) points")
    println("  P: $(length(model.grid_PRODENE_vec)) → $(length(new_grid_P)) points")
    println("  Total: $(length(model.grid_K_vec) * length(model.grid_Y_vec) * length(model.grid_PRODENE_vec)) → $(length(new_grid_K) * length(new_grid_Y) * length(new_grid_P)) states")
    
    # Create new model with refined grids
    refined_model = DPMacroModel(
        # Copy all parameters from original model
        years = model.years,
        ρ = model.ρ,
        kpvs = model.kpvs,
        elvs = model.elvs,
        a = model.a,
        b = model.b,
        
        # Use refined grids
        K_min = new_grid_K[1],
        K_max = new_grid_K[end],
        n_K = length(new_grid_K),
        
        Y_min = new_grid_Y[1],
        Y_max = new_grid_Y[end],
        n_Y = length(new_grid_Y),
        
        PRODENE_scale = model.PRODENE_scale,
        n_PRODENE = length(new_grid_P),
        
        use_linear_interpolation = model.use_linear_interpolation,
        
        # Copy time-dependent parameters
        labor = model.labor,
        newlab = model.newlab,
        aeei_factor = model.aeei_factor,
        growth_factor = model.growth_factor,
        udf = model.udf,
        finite_time_corr = model.finite_time_corr,
        y0 = model.y0,
        k0 = model.k0,
        c0 = model.c0,
        i0 = model.i0
    )
    
    # Override the grids with our adaptive ones - convert to ranges for interpolation compatibility
    refined_model.K_grid = range(new_grid_K[1], new_grid_K[end], length=length(new_grid_K))
    refined_model.Y_grid = range(new_grid_Y[1], new_grid_Y[end], length=length(new_grid_Y))
    refined_model.PRODENE_grid = range(new_grid_P[1], new_grid_P[end], length=length(new_grid_P))
    
    # Store as vectors for our calculations
    refined_model.grid_K_vec = new_grid_K
    refined_model.grid_Y_vec = new_grid_Y  
    refined_model.grid_PRODENE_vec = new_grid_P
    
    # Reinitialize arrays with new dimensions
    refined_model.V = zeros(refined_model.n_K, refined_model.n_Y, refined_model.n_PRODENE, refined_model.T)
    refined_model.V_K = zeros(refined_model.n_K, refined_model.n_Y, refined_model.n_PRODENE, refined_model.T)
    refined_model.V_Y = zeros(refined_model.n_K, refined_model.n_Y, refined_model.n_PRODENE, refined_model.T)
    refined_model.V_P = zeros(refined_model.n_K, refined_model.n_Y, refined_model.n_PRODENE, refined_model.T)
    
    refined_model.C_policy = zeros(refined_model.n_K, refined_model.n_Y, refined_model.n_PRODENE)
    refined_model.I_policy = zeros(refined_model.n_K, refined_model.n_Y, refined_model.n_PRODENE)
    refined_model.KN_policy = zeros(refined_model.n_K, refined_model.n_Y, refined_model.n_PRODENE)
    refined_model.YN_policy = zeros(refined_model.n_K, refined_model.n_Y, refined_model.n_PRODENE)
    refined_model.NEWENE_total_policy = zeros(refined_model.n_K, refined_model.n_Y, refined_model.n_PRODENE)
    refined_model.EC_policy = zeros(refined_model.n_K, refined_model.n_Y, refined_model.n_PRODENE)
    
    return refined_model
end

"""
    create_adaptive_1d_grid(original_grid, needs_refinement, refinement_factor)

Create 1D adaptive grid with higher density in flagged regions.
"""
function create_adaptive_1d_grid(original_grid::Vector{Float64}, needs_refinement::Vector{Bool}, refinement_factor::Float64)
    refined_points = Float64[]
    
    for i in 1:length(original_grid)-1
        # Always include the current point
        push!(refined_points, original_grid[i])
        
        # If this region needs refinement, add intermediate points
        if needs_refinement[i]
            n_extra = max(1, floor(Int, refinement_factor))
            for j in 1:n_extra
                extra_point = original_grid[i] + j * (original_grid[i+1] - original_grid[i]) / (n_extra + 1)
                push!(refined_points, extra_point)
            end
        end
    end
    
    # Always include the last point
    push!(refined_points, original_grid[end])
    
    return sort(unique(refined_points))
end

"""
    compare_grid_performance(coarse_model, refined_model, energy_cost_func)

Compare performance metrics between coarse and refined grids.
"""
function compare_grid_performance(coarse_model::DPMacroModel, refined_model::DPMacroModel, energy_cost_func::Function)
    println("\n=== GRID PERFORMANCE COMPARISON ===")
    
    # Simulate both models
    coarse_results = simulate_trajectory(coarse_model, coarse_model.k0, coarse_model.y0, coarse_model.PRODENE_base, energy_cost_func)
    refined_results = simulate_trajectory(refined_model, refined_model.k0, refined_model.y0, refined_model.PRODENE_base, energy_cost_func)
    
    println("Utility comparison:")
    coarse_utility = coarse_results.UTILITY
    refined_utility = refined_results.UTILITY
    println("  Coarse grid utility: $(round(coarse_utility, digits=2))")
    println("  Refined grid utility: $(round(refined_utility, digits=2))")
    println("  Improvement: $(round(100*(refined_utility - coarse_utility)/coarse_utility, digits=2))%")
    
    println("\nGrid size comparison:")
    coarse_states = coarse_model.n_K * coarse_model.n_Y * coarse_model.n_PRODENE
    refined_states = refined_model.n_K * refined_model.n_Y * refined_model.n_PRODENE
    println("  Coarse grid: $(coarse_states) states")
    println("  Refined grid: $(refined_states) states")
    println("  Size increase: $(round(100*(refined_states - coarse_states)/coarse_states, digits=1))%")
    
    return coarse_results, refined_results
end