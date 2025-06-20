"""
Numerical Stability Tests
========================

Tests for numerical stability, scaling, and robustness of the GBD implementation.
These tests identify potential numerical issues that could cause infeasibility.
"""

using Test
using JuMP, Ipopt
using LinearAlgebra

# Import from module
using EnergyMacroModel

@testset "Variable Scaling Analysis" begin
    # Analyze variable scales across the model
    
    # Energy subproblem scaling
    S_base = Dict((s, y) => enestart[(s, y)] for s in sector for y in year_all)
    sp, _, slack = create_energy_subproblem(S_base)
    optimize!(sp)
    
    if termination_status(sp) in [MOI.OPTIMAL, MOI.LOCALLY_SOLVED]
        # Analyze variable magnitudes
        cost_val = value(sp[:TOTAL_COST])
        slack_vals = [value(slack[s, y]) for s in sector, y in year_all]
        
        println("Energy subproblem variable scaling:")
        println("  TOTAL_COST: $(round(cost_val, digits=3)) billion USD")
        println("  Max slack: $(round(maximum(slack_vals), digits=6)) PWh")
        println("  Total slack: $(round(sum(slack_vals), digits=6)) PWh")
        
        # Test for extreme values
        @test 0.1 < cost_val < 10000
        @test maximum(slack_vals) < 100
    end
    
    # Master problem scaling  
    master = create_master_problem()
    set_gbd_master_bounds_and_initial_values!(master)
    optimize!(master)
    
    if termination_status(master) in [MOI.OPTIMAL, MOI.LOCALLY_SOLVED]
        theta_val = value(master[:theta])
        utility_val = value(master[:UTILITY])
        
        K_vals = [value(master[:K][y]) for y in year_all]
        Y_vals = [value(master[:Y][y]) for y in year_all]
        S_vals = [value(master[:S][s, y]) for s in sector, y in year_all]
        
        println("Master problem variable scaling:")
        println("  theta: $(round(theta_val, digits=3)) trillion USD")
        println("  UTILITY: $(round(utility_val, digits=3))")
        println("  K range: $(round(minimum(K_vals), digits=1)) - $(round(maximum(K_vals), digits=1)) trillion USD")
        println("  Y range: $(round(minimum(Y_vals), digits=1)) - $(round(maximum(Y_vals), digits=1)) trillion USD") 
        println("  S range: $(round(minimum(S_vals), digits=1)) - $(round(maximum(S_vals), digits=1)) PWh")
        
        # Test for reasonable magnitudes
        @test 0.1 < theta_val < 100
        @test 0 < utility_val < 10000
        @test all(0.1 .< K_vals .< 10000)
        @test all(0.1 .< Y_vals .< 10000)
        @test all(0.01 .< S_vals .< 10000)
    end
end

@testset "Condition Number Analysis" begin
    # Analyze condition numbers of key constraints (simplified)
    
    S_test = Dict((s, y) => enestart[(s, y)] for s in sector for y in year_all)
    sp, bal, _ = create_energy_subproblem(S_test)
    
    # Check if we can access constraint matrix (this is solver-dependent)
    # For now, just test that the problem structure is reasonable
    
    @test num_variables(sp) > 0
    @test num_constraints(sp) > 0
    
    var_count = num_variables(sp)
    con_count = num_constraints(sp)
    
    println("Problem size analysis:")
    println("  Energy subproblem: $var_count variables, $con_count constraints")
    
    @test var_count < 100000
    @test con_count < 100000
    @test con_count / var_count < 10
end

@testset "Solver Tolerance Sensitivity" begin
    # Test sensitivity to solver tolerances
    
    S_base = Dict((s, y) => enestart[(s, y)] for s in sector for y in year_all)
    
    tolerances = [1e-4, 1e-5, 1e-6, 1e-7]
    results = []
    
    for tol in tolerances
        sp, _, _ = create_energy_subproblem(S_base)
        set_optimizer_attribute(sp, "tol", tol)
        set_optimizer_attribute(sp, "print_level", 0)
        
        optimize!(sp)
        
        if termination_status(sp) in [MOI.OPTIMAL, MOI.LOCALLY_SOLVED]
            obj_val = objective_value(sp)
            push!(results, (tol, obj_val, true))
        else
            push!(results, (tol, NaN, false))
        end
    end
    
    # Analyze results
    successful_results = [r for r in results if r[3]]
    @test length(successful_results) >= 2
    
    if length(successful_results) >= 2
        obj_values = [r[2] for r in successful_results]
        obj_range = maximum(obj_values) - minimum(obj_values)
        obj_mean = mean(obj_values)
        
        @test obj_range / obj_mean < 0.01
        
        println("Tolerance sensitivity test:")
        for (tol, obj, success) in results
            if success
                println("  tol=$tol: obj=$(round(obj, digits=6))")
            else
                println("  tol=$tol: FAILED")
            end
        end
    end
end

@testset "Parameter Perturbation" begin
    # Test sensitivity to small parameter changes
    
    # Store original values
    orig_drate = drate
    orig_depr = depr
    orig_esub = esub
    
    perturbations = [-0.01, 0.0, 0.01]  # ±1% perturbations
    
    for δ in perturbations
        # Temporarily modify parameters
        global drate = orig_drate * (1 + δ)
        global depr = orig_depr * (1 + δ)
        global esub = orig_esub * (1 + δ)
        
        # Recalculate derived parameters
        global rho = (esub - 1) / esub
        
        try
            # Test energy subproblem
            S_test = Dict((s, y) => enestart[(s, y)] for s in sector for y in year_all)
            sp, _, _ = create_energy_subproblem(S_test)
            set_optimizer_attribute(sp, "print_level", 0)
            optimize!(sp)
            
            @test termination_status(sp) in [MOI.OPTIMAL, MOI.LOCALLY_SOLVED]
            
            # Test master problem
            master = create_master_problem()
            set_gbd_master_bounds_and_initial_values!(master)
            set_optimizer_attribute(master, "print_level", 0)
            optimize!(master)
            
            @test termination_status(master) in [MOI.OPTIMAL, MOI.LOCALLY_SOLVED]
            
        catch e
            @test false
        end
    end
    
    # Restore original values
    global drate = orig_drate
    global depr = orig_depr
    global esub = orig_esub
    global rho = (esub - 1) / esub
end

@testset "Initialization Sensitivity" begin
    # Test sensitivity to different initial values
    
    master = create_master_problem()
    set_gbd_master_bounds_and_initial_values!(master)
    
    # Test different theta starting values
    theta_starts = [1.0, 10.0, 40.0, 49.0]
    
    for theta_start in theta_starts
        set_start_value(master[:theta], theta_start)
        
        # Reset other start values to defaults
        for y in year_all
            set_start_value(master[:C][y], 0.8 * c0)
            set_start_value(master[:I][y], 0.8 * i0)
            set_start_value(master[:Y][y], 0.8 * y0)
        end
        
        set_optimizer_attribute(master, "print_level", 0)
        optimize!(master)
        
        @test termination_status(master) in [MOI.OPTIMAL, MOI.LOCALLY_SOLVED]
        
        if termination_status(master) in [MOI.OPTIMAL, MOI.LOCALLY_SOLVED]
            final_theta = value(master[:theta])
            println("  theta_start=$theta_start -> final_theta=$(round(final_theta, digits=3))")
        end
    end
end

@testset "Bounds Sensitivity" begin
    # Test sensitivity to variable bounds
    
    # Test theta upper bound sensitivity
    master = create_master_problem()
    set_gbd_master_bounds_and_initial_values!(master)
    
    theta_bounds = [10.0, 25.0, 50.0, 100.0]
    
    for ub in theta_bounds
        set_upper_bound(master[:theta], ub)
        set_optimizer_attribute(master, "print_level", 0)
        optimize!(master)
        
        status = termination_status(master)
        
        if status in [MOI.OPTIMAL, MOI.LOCALLY_SOLVED]
            theta_val = value(master[:theta])
            @test theta_val <= ub + 1e-6
            
            # Check if bound is active
            is_active = abs(theta_val - ub) < 1e-4
            println("  theta_bound=$ub: theta=$(round(theta_val, digits=3)), active=$is_active")
            
            if ub < 25.0 && is_active
                @warn "Theta upper bound $ub may be too restrictive"
            end
        else
            println("  theta_bound=$ub: INFEASIBLE")
            if ub < 25.0
                @test status != MOI.INFEASIBLE
            end
        end
    end
end

@testset "Scaling Factor Impact" begin
    # Test impact of the billion->trillion scaling factor
    
    S_base = Dict((s, y) => enestart[(s, y)] for s in sector for y in year_all)
    sp, bal, _ = create_energy_subproblem(S_base)
    optimize!(sp)
    
    if termination_status(sp) in [MOI.OPTIMAL, MOI.LOCALLY_SOLVED]
        cost_billion = objective_value(sp)
        
        # Test different scaling factors
        scaling_factors = [100.0, 1000.0, 10000.0]  # Different billion->trillion conversions
        
        for scale in scaling_factors
            cost_scaled = cost_billion / scale
            
            # Extract and scale dual values
            lambda_scaled = Dict{Tuple{String,Int},Float64}()
            for ((e, l, y), con) in bal
                s_idx = findfirst(t -> haskey(map_energy_sector, (e, l, t)), sector)
                if s_idx !== nothing
                    s = sector[s_idx]
                    dual_val = dual(con)
                    
                    if !haskey(lambda_scaled, (s, y))
                        lambda_scaled[(s, y)] = 0.0
                    end
                    lambda_scaled[(s, y)] += -dual_val / scale
                end
            end
            
            # Test that scaled values are reasonable
            @test 0.0001 < cost_scaled < 100
            
            for (key, lambda_val) in lambda_scaled
                @test abs(lambda_val) < 10
            end
            
            println("Scaling factor $scale: cost=$(round(cost_scaled, digits=6)), max_dual=$(round(maximum(abs.(values(lambda_scaled))), digits=6))")
        end
    end
end

@testset "Numerical Precision" begin
    # Test for numerical precision issues
    
    S_base = Dict((s, y) => enestart[(s, y)] for s in sector for y in year_all)
    
    # Solve twice with identical inputs
    sp1, _, _ = create_energy_subproblem(S_base)
    sp2, _, _ = create_energy_subproblem(S_base)
    
    set_optimizer_attribute(sp1, "print_level", 0)
    set_optimizer_attribute(sp2, "print_level", 0)
    
    optimize!(sp1)
    optimize!(sp2)
    
    if termination_status(sp1) in [MOI.OPTIMAL, MOI.LOCALLY_SOLVED] && 
       termination_status(sp2) in [MOI.OPTIMAL, MOI.LOCALLY_SOLVED]
        
        obj1 = objective_value(sp1)
        obj2 = objective_value(sp2)
        
        @test abs(obj1 - obj2) < 1e-8
        
        println("Numerical precision test:")
        println("  Run 1: $(round(obj1, digits=10))")
        println("  Run 2: $(round(obj2, digits=10))")
        println("  Difference: $(abs(obj1 - obj2))")
    end
end

@testset "Memory and Performance" begin
    # Basic performance and memory tests
    
    println("Performance test:")
    
    # Time energy subproblem construction and solve
    S_test = Dict((s, y) => enestart[(s, y)] for s in sector for y in year_all)
    
    construction_time = @elapsed begin
        sp, _, _ = create_energy_subproblem(S_test)
    end
    
    solve_time = @elapsed begin
        set_optimizer_attribute(sp, "print_level", 0)
        optimize!(sp)
    end
    
    @test construction_time < 10.0
    @test solve_time < 30.0
    
    println("  Construction time: $(round(construction_time, digits=3))s")
    println("  Solve time: $(round(solve_time, digits=3))s")
    
    # Time master problem
    master_construction_time = @elapsed begin
        master = create_master_problem()
        set_gbd_master_bounds_and_initial_values!(master)
    end
    
    master_solve_time = @elapsed begin
        set_optimizer_attribute(master, "print_level", 0)
        optimize!(master)
    end
    
    @test master_construction_time < 5.0
    @test master_solve_time < 30.0
    
    println("  Master construction time: $(round(master_construction_time, digits=3))s")
    println("  Master solve time: $(round(master_solve_time, digits=3))s")
end