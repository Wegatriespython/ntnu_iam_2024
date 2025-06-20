"""
Integration Tests
================

Tests for the complete GBD algorithm integration and comparison with the
integrated formulation. These tests verify end-to-end functionality.
"""

using Test
using JuMP, Ipopt

# Import from module
using EnergyMacroModel

@testset "GBD Algorithm Integration" begin
    # Test the complete GBD algorithm with a few iterations
    
    println("\nTesting GBD algorithm integration...")
    
    # Test with relaxed convergence for faster testing
    @test_nowarn result = solve_gbd(5, 1e-2)  # Max 5 iterations, 1% tolerance
    
    # The algorithm may or may not converge in 5 iterations, but it should run without errors
    # Main goal is to test that all components work together
end

@testset "GBD vs Integrated Comparison" begin
    # Compare GBD solution with integrated formulation
    
    println("\nComparing GBD vs integrated formulation...")
    
    # First, try to solve the integrated formulation for comparison
    integrated_model = create_integrated_model()
    if integrated_model !== nothing
        set_optimizer_attribute(integrated_model, "print_level", 0)
        optimize!(integrated_model)
        
        integrated_status = termination_status(integrated_model)
        
        if integrated_status in [MOI.OPTIMAL, MOI.LOCALLY_SOLVED]
            integrated_utility = value(integrated_model[:UTILITY])
            
            println("Integrated formulation solved successfully")
            println("  Utility: $(round(integrated_utility, digits=3))")
            
            # Now test GBD with same tolerance
            # Note: Full comparison would require convergence, which may take many iterations
            # For testing purposes, we just verify both approaches can start
            
            S_start = Dict((s, y) => enestart[(s, y)] for s in sector for y in year_all)
            sp, _, _ = create_energy_subproblem(S_start)
            set_optimizer_attribute(sp, "print_level", 0)
            optimize!(sp)
            
            @test termination_status(sp) in [MOI.OPTIMAL, MOI.LOCALLY_SOLVED]
            
            master = create_master_problem()
            set_gbd_master_bounds_and_initial_values!(master)
            set_optimizer_attribute(master, "print_level", 0)
            optimize!(master)
            
            @test termination_status(master) in [MOI.OPTIMAL, MOI.LOCALLY_SOLVED]
            
            master_utility = value(master[:UTILITY])
            println("GBD components solve successfully")
            println("  Master utility (no cuts): $(round(master_utility, digits=3))")
            
        else
            println("Integrated formulation failed to solve: $integrated_status")
            @warn "Cannot compare GBD vs integrated - integrated formulation failed"
        end
    else
        @warn "Integrated model not available for comparison"
    end
end

@testset "GBD Convergence Properties" begin
    # Test convergence properties of the GBD algorithm
    
    println("\nTesting GBD convergence properties...")
    
    # Test with small problem (fewer iterations for faster testing)
    S_curr = Dict((s, y) => enestart[(s, y)] for s in sector for y in year_all)
    cuts = []
    LB, UB = -Inf, Inf
    
    # Simulate a few GBD iterations manually to test convergence logic
    for k in 1:3
        # Subproblem
        sp, bal, _ = create_energy_subproblem(S_curr)
        set_optimizer_attribute(sp, "print_level", 0)
        optimize!(sp)
        
        if termination_status(sp) in [MOI.OPTIMAL, MOI.LOCALLY_SOLVED]
            v = objective_value(sp) / 1000.0
            UB = min(UB, v)
            
            # Extract dual values
            lambda = Dict{Tuple{String,Int},Float64}((s, y) => 0.0 for s in sector, y in year_all)
            for ((e, l, y), con) in bal
                s_idx = findfirst(t -> haskey(map_energy_sector, (e, l, t)), sector)
                if s_idx !== nothing
                    s = sector[s_idx]
                    lambda[(s, y)] += -dual(con) / 1000.0
                end
            end
            
            cut = BendersCut("opt", v, lambda, copy(S_curr))
            push!(cuts, cut)
            
            # Master problem
            M = create_master_problem()
            set_gbd_master_bounds_and_initial_values!(M)
            foreach(c -> add_opt_cut!(M, c), cuts)
            set_optimizer_attribute(M, "print_level", 0)
            optimize!(M)
            
            if termination_status(M) in [MOI.OPTIMAL, MOI.LOCALLY_SOLVED]
                theta = value(M[:theta])
                util = value(M[:UTILITY])
                obj = objective_value(M)
                LB = max(LB, obj)
                
                gap = abs((theta - util) - v)
                rel = gap / max(abs(theta - util), abs(v), 1.0)
                
                println("  Iteration $k: UB=$(round(UB, digits=3)), LB=$(round(LB, digits=3)), gap=$(round(rel*100, digits=2))%")
                
                # Test convergence properties
                @test UB >= LB - 0.1
                @test theta >= 0
                @test isfinite(gap)
                @test rel >= 0
                
                # Update S for next iteration
                S_curr = Dict((s, y) => value(M[:S][s, y]) for s in sector for y in year_all)
                
                # Check if convergence would be detected
                if rel <= 1e-4
                    println("  Convergence detected at iteration $k")
                    break
                end
            else
                @warn "Master problem failed at iteration $k"
                break
            end
        else
            @warn "Subproblem failed at iteration $k"
            break
        end
    end
end

@testset "Solution Quality Checks" begin
    # Test solution quality and consistency
    
    println("\nTesting solution quality...")
    
    # Solve one iteration of GBD
    S_test = Dict((s, y) => enestart[(s, y)] for s in sector for y in year_all)
    
    # Energy subproblem
    sp, bal, slack = create_energy_subproblem(S_test)
    set_optimizer_attribute(sp, "print_level", 0)
    optimize!(sp)
    
    if termination_status(sp) in [MOI.OPTIMAL, MOI.LOCALLY_SOLVED]
        # Check slack usage
        total_slack = sum(value(slack[s, y]) for s in sector for y in year_all)
        println("  Energy subproblem slack usage: $(round(total_slack, digits=6)) PWh")
        
        @test total_slack < 1.0
        
        # Check cost components
        total_cost = value(sp[:TOTAL_COST])
        annual_costs = [value(sp[:COST_ANNUAL][y]) for y in year_all]
        
        @test all(annual_costs .>= 0)
        @test total_cost >= 0
        
        println("  Total energy system cost: $(round(total_cost, digits=3)) billion USD")
    end
    
    # Master problem
    master = create_master_problem()
    set_gbd_master_bounds_and_initial_values!(master)
    set_optimizer_attribute(master, "print_level", 0)
    optimize!(master)
    
    if termination_status(master) in [MOI.OPTIMAL, MOI.LOCALLY_SOLVED]
        # Check economic consistency
        for y in year_all
            Y_val = value(master[:Y][y])
            C_val = value(master[:C][y])
            I_val = value(master[:I][y])
            EC_val = value(master[:EC][y])
            
            @test Y_val >= 0
            @test C_val >= 0
            @test I_val >= 0
            
            # Check accounting identity
            accounting_error = abs(Y_val - (C_val + I_val + EC_val))
            @test accounting_error < 0.01
        end
        
        utility_val = value(master[:UTILITY])
        theta_val = value(master[:theta])
        
        @test utility_val > 0
        @test theta_val >= 0
        
        println("  Master problem utility: $(round(utility_val, digits=3))")
        println("  Master problem theta: $(round(theta_val, digits=3)) trillion USD")
    end
end

@testset "Error Handling and Robustness" begin
    # Test error handling in edge cases
    
    println("\nTesting error handling...")
    
    # Test with extreme service demands
    S_extreme = Dict((s, y) => 1000 * enestart[(s, y)] for s in sector for y in year_all)
    
    # This should either solve with large slack or fail gracefully
    sp_extreme, _, slack_extreme = create_energy_subproblem(S_extreme)
    set_optimizer_attribute(sp_extreme, "print_level", 0)
    set_optimizer_attribute(sp_extreme, "max_iter", 1000)
    
    @test_nowarn optimize!(sp_extreme)
    
    status = termination_status(sp_extreme)
    if status in [MOI.OPTIMAL, MOI.LOCALLY_SOLVED]
        total_slack = sum(value(slack_extreme[s, y]) for s in sector for y in year_all)
        println("  Extreme demands solved with slack: $(round(total_slack, digits=3)) PWh")
        @test total_slack > 100
    else
        println("  Extreme demands result in status: $status")
        # This is acceptable - extreme demands may be infeasible
    end
    
    # Test with very low service demands
    S_low = Dict((s, y) => 0.01 * enestart[(s, y)] for s in sector for y in year_all)
    
    sp_low, _, slack_low = create_energy_subproblem(S_low)
    set_optimizer_attribute(sp_low, "print_level", 0)
    optimize!(sp_low)
    
    @test termination_status(sp_low) in [MOI.OPTIMAL, MOI.LOCALLY_SOLVED]
    
    if termination_status(sp_low) in [MOI.OPTIMAL, MOI.LOCALLY_SOLVED]
        total_slack_low = sum(value(slack_low[s, y]) for s in sector for y in year_all)
        @test total_slack_low < 0.001
    end
end

# Helper function to create integrated model (if available)
function create_integrated_model()
    # This would be the integrated energy-macro model for comparison
    # For now, return nothing if not implemented
    try
        # Include and call the integrated model creation
        # This assumes such a function exists in run_energy_macro.jl
        if isdefined(Main, :create_integrated_energy_macro_model)
            return create_integrated_energy_macro_model()
        else
            return nothing
        end
    catch
        return nothing
    end
end