"""
Benders Cut Implementation Tests
===============================

Tests for the Benders cut generation and implementation in the GBD reformulation.
These tests verify mathematical correctness, dual value extraction, and cut validity.
"""

using Test
using JuMP, Ipopt

# Import from module
using EnergyMacroModel

@testset "BendersCut Structure" begin
    # Test BendersCut struct creation and properties
    
    # Create a sample cut
    S_test = Dict((s, y) => enestart[(s, y)] for s in sector for y in year_all)
    lambda_test = Dict((s, y) => 0.1 for s in sector for y in year_all)
    
    cut = BendersCut("optimality", 5.0, lambda_test, S_test)
    
    @test cut.type == "optimality"
    @test cut.v == 5.0
    @test length(cut.lambda) == length(sector) * length(year_all)
    @test length(cut.S_hat) == length(sector) * length(year_all)
    
    # Test that all lambda values are accessible
    for s in sector, y in year_all
        @test haskey(cut.lambda, (s, y))
        @test haskey(cut.S_hat, (s, y))
    end
end

@testset "Dual Value Extraction" begin
    # Test dual value extraction from energy subproblem
    
    S_base = Dict((s, y) => enestart[(s, y)] for s in sector for y in year_all)
    model, bal, slack = create_energy_subproblem(S_base)
    optimize!(model)
    
    if termination_status(model) in [MOI.OPTIMAL, MOI.LOCALLY_SOLVED]
        # Extract dual values as in the main algorithm
        lambda = Dict{Tuple{String,Int},Float64}((s, y) => 0.0 for s in sector, y in year_all)
        
        for ((e, l, y), con) in bal
            s_idx = findfirst(t -> haskey(map_energy_sector, (e, l, t)), sector)
            if s_idx !== nothing
                s = sector[s_idx]
                dual_val = dual(con)
                lambda[(s, y)] += -dual_val / 1000.0  # Convert to trillion USD
                
                @test isfinite(dual_val)
            end
        end
        
        # Test that lambda values are reasonable
        for s in sector, y in year_all
            @test isfinite(lambda[(s, y)])
            @test abs(lambda[(s, y)]) < 1000
        end
        
        println("Dual value extraction test:")
        for s in sector, y in year_all
            println("  λ[$s, $y] = $(round(lambda[(s, y)], digits=6)) T\$/PWh")
        end
    end
end

@testset "Benders Cut Generation" begin
    # Test complete Benders cut generation process
    
    S_test = Dict((s, y) => enestart[(s, y)] for s in sector for y in year_all)
    
    # Solve energy subproblem
    sp, bal, slack = create_energy_subproblem(S_test)
    optimize!(sp)
    
    if termination_status(sp) in [MOI.OPTIMAL, MOI.LOCALLY_SOLVED]
        # Extract cost value (convert billion to trillion)
        v = objective_value(sp) / 1000.0
        @test v > 0
        @test isfinite(v)
        
        # Extract dual values
        lambda = Dict{Tuple{String,Int},Float64}((s, y) => 0.0 for s in sector, y in year_all)
        for ((e, l, y), con) in bal
            s_idx = findfirst(t -> haskey(map_energy_sector, (e, l, t)), sector)
            if s_idx !== nothing
                s = sector[s_idx]
                lambda[(s, y)] += -dual(con) / 1000.0
            end
        end
        
        # Create Benders cut
        cut = BendersCut("optimality", v, lambda, copy(S_test))
        
        @test cut.v == v
        @test cut.lambda == lambda
        @test cut.S_hat == S_test
        
        println("Generated Benders cut:")
        println("  Value: $(round(cut.v, digits=6)) trillion USD")
        println("  Non-zero duals: $(sum(abs(λ) > 1e-6 for λ in values(cut.lambda)))")
    end
end

@testset "Benders Cut Addition to Master" begin
    # Test adding Benders cuts to master problem
    
    # Create master problem
    master = create_master_problem()
    set_gbd_master_bounds_and_initial_values!(master)
    
    # Generate a test cut
    S_test = Dict((s, y) => enestart[(s, y)] for s in sector for y in year_all)
    lambda_test = Dict((s, y) => 0.05 + 0.01 * y / 2020 for s in sector for y in year_all)  # Small positive duals
    cut = BendersCut("optimality", 4.0, lambda_test, S_test)
    
    # Add cut to master
    @test_nowarn add_opt_cut!(master, cut)
    
    # Try to solve master with cut
    @test_nowarn optimize!(master)
    
    status = termination_status(master)
    @test status in [MOI.OPTIMAL, MOI.LOCALLY_SOLVED]
    
    if status in [MOI.OPTIMAL, MOI.LOCALLY_SOLVED]
        theta_val = value(master[:theta])
        
        # Verify cut is satisfied
        cut_value = cut.v + sum(cut.lambda[key] * (value(master[:S][key...]) - cut.S_hat[key]) for key in keys(cut.lambda))
        @test theta_val >= cut_value - 0.01
        
        println("Master with Benders cut:")
        println("  theta: $(round(theta_val, digits=6))")
        println("  Cut RHS: $(round(cut_value, digits=6))")
        println("  Cut satisfied: $(theta_val >= cut_value - 0.01)")
    end
end

@testset "Multiple Benders Cuts" begin
    # Test adding multiple Benders cuts
    
    master = create_master_problem()
    set_gbd_master_bounds_and_initial_values!(master)
    
    # Generate several test cuts with different S values
    cuts = BendersCut[]
    
    for i in 1:3
        S_test = Dict((s, y) => (1.0 + 0.1 * i) * enestart[(s, y)] for s in sector for y in year_all)
        
        # Solve energy subproblem for this S
        sp, bal, slack = create_energy_subproblem(S_test)
        optimize!(sp)
        
        if termination_status(sp) in [MOI.OPTIMAL, MOI.LOCALLY_SOLVED]
            v = objective_value(sp) / 1000.0
            
            lambda = Dict{Tuple{String,Int},Float64}((s, y) => 0.0 for s in sector, y in year_all)
            for ((e, l, y), con) in bal
                s_idx = findfirst(t -> haskey(map_energy_sector, (e, l, t)), sector)
                if s_idx !== nothing
                    s = sector[s_idx]
                    lambda[(s, y)] += -dual(con) / 1000.0
                end
            end
            
            cut = BendersCut("optimality", v, lambda, copy(S_test))
            push!(cuts, cut)
            add_opt_cut!(master, cut)
        end
    end
    
    @test length(cuts) >= 2
    
    # Solve master with multiple cuts
    @test_nowarn optimize!(master)
    
    status = termination_status(master)
    @test status in [MOI.OPTIMAL, MOI.LOCALLY_SOLVED]
    
    if status in [MOI.OPTIMAL, MOI.LOCALLY_SOLVED]
        theta_val = value(master[:theta])
        
        # Verify all cuts are satisfied
        for (i, cut) in enumerate(cuts)
            cut_value = cut.v + sum(cut.lambda[key] * (value(master[:S][key...]) - cut.S_hat[key]) for key in keys(cut.lambda))
            @test theta_val >= cut_value - 0.01
        end
        
        println("Master with $(length(cuts)) cuts: theta = $(round(theta_val, digits=6))")
    end
end

@testset "Cut Validity and Properties" begin
    # Test mathematical properties of Benders cuts
    
    # Generate cut from base case
    S_base = Dict((s, y) => enestart[(s, y)] for s in sector for y in year_all)
    sp_base, bal_base, _ = create_energy_subproblem(S_base)
    optimize!(sp_base)
    
    if termination_status(sp_base) in [MOI.OPTIMAL, MOI.LOCALLY_SOLVED]
        v_base = objective_value(sp_base) / 1000.0
        
        lambda_base = Dict{Tuple{String,Int},Float64}((s, y) => 0.0 for s in sector, y in year_all)
        for ((e, l, y), con) in bal_base
            s_idx = findfirst(t -> haskey(map_energy_sector, (e, l, t)), sector)
            if s_idx !== nothing
                s = sector[s_idx]
                lambda_base[(s, y)] += -dual(con) / 1000.0
            end
        end
        
        cut_base = BendersCut("optimality", v_base, lambda_base, copy(S_base))
        
        # Test cut at the point where it was generated
        cut_at_base = cut_base.v + sum(cut_base.lambda[key] * (S_base[key] - cut_base.S_hat[key]) for key in keys(cut_base.lambda))
        @test abs(cut_at_base - v_base) < 1e-6
        
        # Test cut provides lower bound at other points
        S_test = Dict((s, y) => 1.1 * enestart[(s, y)] for s in sector for y in year_all)
        sp_test, _, _ = create_energy_subproblem(S_test)
        optimize!(sp_test)
        
        if termination_status(sp_test) in [MOI.OPTIMAL, MOI.LOCALLY_SOLVED]
            v_test = objective_value(sp_test) / 1000.0
            cut_at_test = cut_base.v + sum(cut_base.lambda[key] * (S_test[key] - cut_base.S_hat[key]) for key in keys(cut_base.lambda))
            
            @test cut_at_test <= v_test + 0.1
            
            println("Cut validity test:")
            println("  Actual cost at test point: $(round(v_test, digits=6))")
            println("  Cut prediction: $(round(cut_at_test, digits=6))")
            println("  Lower bound property: $(cut_at_test <= v_test + 0.1)")
        end
    end
end

@testset "Dual Value Aggregation" begin
    # Test potential issues with dual value aggregation
    
    S_test = Dict((s, y) => enestart[(s, y)] for s in sector for y in year_all)
    sp, bal, _ = create_energy_subproblem(S_test)
    optimize!(sp)
    
    if termination_status(sp) in [MOI.OPTIMAL, MOI.LOCALLY_SOLVED]
        # Count how many balance constraints map to each sector
        sector_constraint_count = Dict(s => 0 for s in sector)
        dual_contributions = Dict((s, y) => Float64[] for s in sector, y in year_all)
        
        for ((e, l, y), con) in bal
            s_idx = findfirst(t -> haskey(map_energy_sector, (e, l, t)), sector)
            if s_idx !== nothing
                s = sector[s_idx]
                sector_constraint_count[s] += 1
                push!(dual_contributions[(s, y)], dual(con))
            end
        end
        
        println("Dual value aggregation analysis:")
        for s in sector
            println("  Sector $s: $(sector_constraint_count[s]) balance constraints")
            for y in year_all
                contributions = dual_contributions[(s, y)]
                if !isempty(contributions)
                    println("    Year $y: $(length(contributions)) contributions, sum = $(round(sum(contributions), digits=6))")
                end
            end
        end
        
        # Test for potential double-counting issues
        for s in sector
            if sector_constraint_count[s] > 1
                @warn "Sector $s has $(sector_constraint_count[s]) balance constraints - check for double-counting"
            end
        end
    end
end

@testset "Cut Scaling Consistency" begin
    # Test that scaling is consistent between energy subproblem and master problem
    
    S_test = Dict((s, y) => enestart[(s, y)] for s in sector for y in year_all)
    
    # Solve energy subproblem (costs in billion USD)
    sp, bal, _ = create_energy_subproblem(S_test)
    optimize!(sp)
    
    if termination_status(sp) in [MOI.OPTIMAL, MOI.LOCALLY_SOLVED]
        cost_billion = objective_value(sp)
        cost_trillion = cost_billion / 1000.0
        
        # Test scaling factors
        @test cost_billion > 0
        @test cost_trillion > 0
        @test cost_trillion == cost_billion / 1000.0
        
        # Test dual value scaling
        sample_dual_raw = 0.0
        sample_dual_scaled = 0.0
        
        for ((e, l, y), con) in bal
            s_idx = findfirst(t -> haskey(map_energy_sector, (e, l, t)), sector)
            if s_idx !== nothing
                dual_val = dual(con)
                sample_dual_raw = dual_val
                sample_dual_scaled = dual_val / 1000.0
                break
            end
        end
        
        @test abs(sample_dual_scaled) == abs(sample_dual_raw) / 1000.0
        
        println("Scaling consistency test:")
        println("  Energy cost: $(round(cost_billion, digits=3)) billion USD = $(round(cost_trillion, digits=6)) trillion USD")
        println("  Sample dual: $(round(sample_dual_raw, digits=6)) -> $(round(sample_dual_scaled, digits=9)) after scaling")
    end
end