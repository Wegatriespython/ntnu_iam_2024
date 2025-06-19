"""
Energy Subproblem Tests
======================

Tests for the energy subproblem component of the GBD reformulation.
These tests verify feasibility, optimality, and robustness of the energy system optimization.
"""

using Test
using JuMP, Ipopt

# Import from module
using EnergyMacroModel

@testset "Energy Subproblem Construction" begin
    # Test with base case service demands
    S_base = Dict((s, y) => enestart[(s, y)] for s in sector for y in year_all)
    
    @test_nowarn model, bal, slack = create_energy_subproblem(S_base)
    
    model, bal, slack = create_energy_subproblem(S_base)
    
    @test isa(model, Model)
    @test isa(bal, Dict)
    @test isa(slack, Dict)
    
    # Check that model has expected variables
    @test haskey(model.obj_dict, :ACT)
    @test haskey(model.obj_dict, :CAP_NEW)
    @test haskey(model.obj_dict, :TOTAL_COST)
    @test haskey(model.obj_dict, :demand_slack)
    
    # Check variable dimensions
    ACT = model[:ACT]
    @test size(ACT) == (length(technology), length(year_all))
    
    CAP_NEW = model[:CAP_NEW]
    @test size(CAP_NEW) == (length(technology), length(year_all))
    
    demand_slack = model[:demand_slack]
    @test size(demand_slack) == (length(sector), length(year_all))
end

@testset "Energy Subproblem Feasibility - Base Case" begin
    # Test feasibility with base case demands
    S_base = Dict((s, y) => enestart[(s, y)] for s in sector for y in year_all)
    
    model, bal, slack = create_energy_subproblem(S_base)
    
    # Try to solve
    @test_nowarn optimize!(model)
    
    status = termination_status(model)
    @test status in [MOI.OPTIMAL, MOI.LOCALLY_SOLVED]
    
    if status in [MOI.OPTIMAL, MOI.LOCALLY_SOLVED]
        @test objective_value(model) > 0
        @test isfinite(objective_value(model))
        
        # Check slack variables - should be minimal for base case
        total_slack = sum(value(slack[(s, y)]) for s in sector for y in year_all)
        @test total_slack < 1.0
        
        # Log the solution for debugging
        println("Base case energy subproblem solved successfully")
        println("  Objective value: $(round(objective_value(model), digits=3)) billion USD")
        println("  Total slack: $(round(total_slack, digits=6)) PWh")
        
        # Check individual slack values
        for s in sector, y in year_all
            slack_val = value(slack[(s, y)])
            if slack_val > 0.01  # Report significant slack
                println("  Slack for $s in $y: $(round(slack_val, digits=3)) PWh")
            end
        end
    end
end

@testset "Energy Subproblem Feasibility - Stress Tests" begin
    # Test 1: Slightly increased demands (10% above base)
    S_high = Dict((s, y) => 1.1 * enestart[(s, y)] for s in sector for y in year_all)
    
    model_high, _, slack_high = create_energy_subproblem(S_high)
    @test_nowarn optimize!(model_high)
    
    status_high = termination_status(model_high)
    @test status_high in [MOI.OPTIMAL, MOI.LOCALLY_SOLVED]
    
    if status_high in [MOI.OPTIMAL, MOI.LOCALLY_SOLVED]
        total_slack_high = sum(value(slack_high[(s, y)]) for s in sector for y in year_all)
        println("High demand case: objective = $(round(objective_value(model_high), digits=3)), slack = $(round(total_slack_high, digits=6))")
    end
    
    # Test 2: Decreased demands (50% of base - should be very feasible)
    S_low = Dict((s, y) => 0.5 * enestart[(s, y)] for s in sector for y in year_all)
    
    model_low, _, slack_low = create_energy_subproblem(S_low)
    @test_nowarn optimize!(model_low)
    
    status_low = termination_status(model_low)
    @test status_low in [MOI.OPTIMAL, MOI.LOCALLY_SOLVED]
    
    if status_low in [MOI.OPTIMAL, MOI.LOCALLY_SOLVED]
        total_slack_low = sum(value(slack_low[(s, y)]) for s in sector for y in year_all)
        @test total_slack_low < 0.001
        println("Low demand case: objective = $(round(objective_value(model_low), digits=3)), slack = $(round(total_slack_low, digits=6))")
    end
end

@testset "Energy Subproblem Infeasibility Detection" begin
    # Test 3: Extremely high demands (should require significant slack or be infeasible)
    S_extreme = Dict((s, y) => 100 * enestart[(s, y)] for s in sector for y in year_all)
    
    model_extreme, _, slack_extreme = create_energy_subproblem(S_extreme)
    
    # Increase solver limits for this test
    set_optimizer_attribute(model_extreme, "max_iter", 3000)
    set_optimizer_attribute(model_extreme, "tol", 1e-4)
    
    @test_nowarn optimize!(model_extreme)
    
    status_extreme = termination_status(model_extreme)
    
    if status_extreme in [MOI.OPTIMAL, MOI.LOCALLY_SOLVED]
        total_slack_extreme = sum(value(slack_extreme[(s, y)]) for s in sector for y in year_all)
        println("Extreme demand case solved with slack = $(round(total_slack_extreme, digits=3)) PWh")
        
        # With 100x demands, we expect significant slack usage
        @test total_slack_extreme > 10
        
        # Cost should be dominated by slack penalty
        obj_val = objective_value(model_extreme)
        slack_penalty = 1e6 * total_slack_extreme
        @test slack_penalty / obj_val > 0.8
        
    else
        println("Extreme demand case status: $status_extreme")
        @test status_extreme in [MOI.INFEASIBLE, MOI.DUAL_INFEASIBLE, MOI.INFEASIBLE_OR_UNBOUNDED]
    end
end

@testset "Energy Subproblem Solution Properties" begin
    S_base = Dict((s, y) => enestart[(s, y)] for s in sector for y in year_all)
    model, bal, slack = create_energy_subproblem(S_base)
    optimize!(model)
    
    if termination_status(model) in [MOI.OPTIMAL, MOI.LOCALLY_SOLVED]
        # Test monotonicity: higher demands should lead to higher costs
        S_higher = Dict((s, y) => 1.2 * enestart[(s, y)] for s in sector for y in year_all)
        model_higher, _, _ = create_energy_subproblem(S_higher)
        optimize!(model_higher)
        
        if termination_status(model_higher) in [MOI.OPTIMAL, MOI.LOCALLY_SOLVED]
            @test objective_value(model_higher) >= objective_value(model)
        end
        
        # Test that TOTAL_COST matches sum of discounted annual costs
        total_cost = value(model[:TOTAL_COST])
        annual_costs = model[:COST_ANNUAL]
        
        # Calculate expected total from annual costs
        year_idx = Dict(y => i for (i, y) in enumerate(year_all))
        disc = Dict(y => (1 - drate)^(duration_period * (year_idx[y] - 1)) for y in year_all)
        expected_total = sum(value(annual_costs[y]) * duration_period * disc[y] for y in year_all)
        
        @test abs(total_cost - expected_total) / total_cost < 0.01
    end
end

@testset "Energy Balance Constraint Analysis" begin
    S_test = Dict((s, y) => enestart[(s, y)] for s in sector for y in year_all)
    model, bal, slack = create_energy_subproblem(S_test)
    optimize!(model)
    
    if termination_status(model) in [MOI.OPTIMAL, MOI.LOCALLY_SOLVED]
        # Analyze balance constraints and their dual values
        println("\nBalance Constraint Analysis:")
        
        for ((e, l, y), constraint) in bal
            # Find corresponding sector
            s_idx = findfirst(s -> haskey(map_energy_sector, (e, l, s)), sector)
            if s_idx !== nothing
                s = sector[s_idx]
                dual_val = dual(constraint)
                slack_val = value(slack[(s, y)])
                
                println("  ($e, $l, $y) -> $s: dual = $(round(dual_val, digits=6)), slack = $(round(slack_val, digits=6))")
                
                # Test dual value properties
                @test isfinite(dual_val)
                
                # If slack is positive, dual should be large (due to penalty)
                if slack_val > 0.01
                    @test abs(dual_val) > 1000
                end
            end
        end
    end
end

@testset "Energy Subproblem Variable Bounds" begin
    S_base = Dict((s, y) => enestart[(s, y)] for s in sector for y in year_all)
    model, _, _ = create_energy_subproblem(S_base)
    
    # Check that all variables have appropriate bounds
    ACT = model[:ACT]
    CAP_NEW = model[:CAP_NEW]
    demand_slack = model[:demand_slack]
    
    # All activity and capacity variables should be non-negative
    for t in technology, y in year_all
        @test lower_bound(ACT[t, y]) == 0
        @test lower_bound(CAP_NEW[t, y]) == 0
    end
    
    # All slack variables should be non-negative
    for s in sector, y in year_all
        @test lower_bound(demand_slack[s, y]) == 0
    end
    
    # TOTAL_COST should be unrestricted (no explicit bounds typically)
    @test !has_lower_bound(model[:TOTAL_COST]) || lower_bound(model[:TOTAL_COST]) <= 0
end

@testset "Energy Subproblem Scaling" begin
    S_base = Dict((s, y) => enestart[(s, y)] for s in sector for y in year_all)
    model, _, _ = create_energy_subproblem(S_base)
    optimize!(model)
    
    if termination_status(model) in [MOI.OPTIMAL, MOI.LOCALLY_SOLVED]
        obj_val = objective_value(model)
        
        # Test that objective value is in reasonable range (billion USD)
        @test 0.1 < obj_val < 1000
        
        # Test scaling consistency - convert to trillion
        obj_trillion = obj_val / 1000
        @test 0.0001 < obj_trillion < 1
        
        println("Energy subproblem objective: $(round(obj_val, digits=3)) billion USD = $(round(obj_trillion, digits=6)) trillion USD")
    end
end