"""
Master Problem (MACRO) Tests
============================

Tests for the master problem component of the GBD reformulation.
These tests verify the MACRO model formulation, feasibility, and economic consistency.
"""

using Test
using JuMP, Ipopt

# Import from module
using EnergyMacroModel

@testset "Master Problem Construction" begin
    @test_nowarn model = create_master_problem()
    
    model = create_master_problem()
    
    @test isa(model, Model)
    
    # Check that model has expected macro variables
    @test haskey(model.obj_dict, :S)
    @test haskey(model.obj_dict, :theta)
    @test haskey(model.obj_dict, :K)
    @test haskey(model.obj_dict, :Y)
    @test haskey(model.obj_dict, :C)
    @test haskey(model.obj_dict, :I)
    @test haskey(model.obj_dict, :UTILITY)
    
    # Check variable dimensions
    S = model[:S]
    @test size(S) == (length(sector), length(year_all))
    
    K = model[:K]
    @test length(K) == length(year_all)
    
    # Check theta bounds
    theta = model[:theta]
    @test lower_bound(theta) == 0
    @test upper_bound(theta) == 50
end

@testset "Master Problem Without Benders Cuts" begin
    # Test master problem feasibility without any Benders cuts
    model = create_master_problem()
    set_gbd_master_bounds_and_initial_values!(model)
    
    @test_nowarn optimize!(model)
    
    status = termination_status(model)
    @test status in [MOI.OPTIMAL, MOI.LOCALLY_SOLVED]
    
    if status in [MOI.OPTIMAL, MOI.LOCALLY_SOLVED]
        @test objective_value(model) < 0
        @test isfinite(objective_value(model))
        
        # Check that theta is positive (energy costs)
        theta_val = value(model[:theta])
        @test theta_val > 0
        @test theta_val < 50
        
        # Check utility value
        utility_val = value(model[:UTILITY])
        @test utility_val > 0
        
        println("Master problem without cuts solved successfully")
        println("  Objective value: $(round(objective_value(model), digits=3))")
        println("  theta: $(round(theta_val, digits=3)) trillion USD")
        println("  Utility: $(round(utility_val, digits=3))")
    end
end

@testset "Master Problem Variable Values" begin
    model = create_master_problem()
    set_gbd_master_bounds_and_initial_values!(model)
    optimize!(model)
    
    if termination_status(model) in [MOI.OPTIMAL, MOI.LOCALLY_SOLVED]
        # Test economic consistency of solution
        
        # Capital stock should be reasonable
        for y in year_all
            K_val = value(model[:K][y])
            @test K_val > 0
            @test K_val > 0.8 * k0
            @test K_val < 10 * k0
        end
        
        # Production should be reasonable
        for y in year_all
            Y_val = value(model[:Y][y])
            @test Y_val > 0
            @test Y_val > 0.8 * y0
        end
        
        # Consumption should be reasonable
        for y in year_all
            C_val = value(model[:C][y])
            @test C_val > 0
            @test C_val > 0.5 * c0
        end
        
        # Investment should be reasonable
        for y in year_all
            I_val = value(model[:I][y])
            @test I_val > 0
        end
        
        # Service demands should be reasonable
        for s in sector, y in year_all
            S_val = value(model[:S][s, y])
            @test S_val > 0
            @test S_val > 0.01
            @test S_val < 1000
        end
    end
end

@testset "Master Problem Economic Accounting" begin
    model = create_master_problem()
    set_gbd_master_bounds_and_initial_values!(model)
    optimize!(model)
    
    if termination_status(model) in [MOI.OPTIMAL, MOI.LOCALLY_SOLVED]
        # Test accounting identities
        
        for y in year_all
            Y_val = value(model[:Y][y])
            C_val = value(model[:C][y])
            I_val = value(model[:I][y])
            EC_val = value(model[:EC][y])
            
            # Capital constraint: Y = C + I + EC
            @test abs(Y_val - (C_val + I_val + EC_val)) < 0.01
        end
        
        # Test capital accumulation
        for i in 2:length(year_all)
            y = year_all[i]
            y_prev = year_all[i-1]
            
            K_val = value(model[:K][y])
            K_prev = value(model[:K][y_prev])
            KN_val = value(model[:KN][y])
            
            expected_K = K_prev * (1 - depr)^duration_period + KN_val
            @test abs(K_val - expected_K) < 0.01
        end
        
        # Test new capital relationship
        for i in 2:length(year_all)
            y = year_all[i]
            I_val = value(model[:I][y])
            KN_val = value(model[:KN][y])
            
            expected_KN = duration_period * I_val
            @test abs(KN_val - expected_KN) < 0.01
        end
    end
end

@testset "Master Problem CES Production Function" begin
    model = create_master_problem()
    set_gbd_master_bounds_and_initial_values!(model)
    optimize!(model)
    
    if termination_status(model) in [MOI.OPTIMAL, MOI.LOCALLY_SOLVED]
        # Test CES production function for non-base years
        
        for i in 2:length(year_all)
            y = year_all[i]
            
            YN_val = value(model[:YN][y])
            KN_val = value(model[:KN][y])
            NEWENE_ELEC = value(model[:NEWENE]["ELEC", y])
            NEWENE_NELE = value(model[:NEWENE]["NELE", y])
            
            # Calculate expected YN from CES function
            ces_term = (a * KN_val^(rho * kpvs) * newlab[y]^(rho * (1 - kpvs)) +
                       b * NEWENE_ELEC^(rho * elvs) * NEWENE_NELE^(rho * (1 - elvs)))^(1/rho)
            
            @test abs(YN_val - ces_term) / max(YN_val, ces_term) < 0.01
        end
    end
end

@testset "Master Problem Energy-Economy Linkage" begin
    model = create_master_problem()
    set_gbd_master_bounds_and_initial_values!(model)
    optimize!(model)
    
    if termination_status(model) in [MOI.OPTIMAL, MOI.LOCALLY_SOLVED]
        # Test energy-economy linkage
        
        for s in sector
            for i in 2:length(year_all)
                y = year_all[i]
                
                S_val = value(model[:S][s, y])
                PRODENE_val = value(model[:PRODENE][s, y])
                
                # S = PRODENE * aeei_factor
                expected_S = PRODENE_val * aeei_factor[(s, y)]
                @test abs(S_val - expected_S) < 0.01
            end
        end
        
        # Test energy accumulation
        for s in sector
            for i in 2:length(year_all)
                y = year_all[i]
                y_prev = year_all[i-1]
                
                PRODENE_val = value(model[:PRODENE][s, y])
                PRODENE_prev = value(model[:PRODENE][s, y_prev])
                NEWENE_val = value(model[:NEWENE][s, y])
                
                expected_PRODENE = PRODENE_prev * (1 - depr)^duration_period + NEWENE_val
                @test abs(PRODENE_val - expected_PRODENE) < 0.01
            end
        end
    end
end

@testset "Master Problem Energy Cost Expansion" begin
    model = create_master_problem()
    set_gbd_master_bounds_and_initial_values!(model)
    optimize!(model)
    
    if termination_status(model) in [MOI.OPTIMAL, MOI.LOCALLY_SOLVED]
        theta_val = value(model[:theta])
        
        # Test energy cost expansion
        totW = sum(cost_MESSAGE[y] for y in year_all if y != 2020)
        
        for i in 2:length(year_all)
            y = year_all[i]
            EC_val = value(model[:EC][y])
            
            disc_y = (1 - drate)^(duration_period * (i - 1))
            expected_EC = theta_val * (cost_MESSAGE[y] / totW) / (disc_y * duration_period)
            
            @test abs(EC_val - expected_EC) < 0.01
        end
    end
end

@testset "Master Problem Bounds and Initialization" begin
    model = create_master_problem()
    
    # Test before setting bounds
    @test lower_bound(model[:theta]) == 0
    @test upper_bound(model[:theta]) == 50
    
    # Apply bounds and initial values
    set_gbd_master_bounds_and_initial_values!(model)
    
    # Test K bounds
    @test has_upper_bound(model[:K][2020])
    @test has_lower_bound(model[:K][2020])
    @test is_fixed(model[:K][2020])
    @test abs(fix_value(model[:K][2020]) - k0) < 1e-10
    
    # Test S bounds
    for s in sector, y in year_all
        @test lower_bound(model[:S][s, y]) == 0.01
        @test upper_bound(model[:S][s, y]) == 1000.0
    end
    
    # Test that bounds are feasible
    @test_nowarn optimize!(model)
    status = termination_status(model)
    @test status in [MOI.OPTIMAL, MOI.LOCALLY_SOLVED]
end

@testset "Master Problem Terminal Condition" begin
    model = create_master_problem()
    set_gbd_master_bounds_and_initial_values!(model)
    optimize!(model)
    
    if termination_status(model) in [MOI.OPTIMAL, MOI.LOCALLY_SOLVED]
        # Test terminal investment condition for 2080
        if 2080 in year_all
            I_2080 = value(model[:I][2080])
            K_2080 = value(model[:K][2080])
            
            required_investment = K_2080 * (grow[2080] + depr)
            @test I_2080 >= required_investment - 0.01
        end
    end
end

@testset "Master Problem Scaling and Magnitude" begin
    model = create_master_problem()
    set_gbd_master_bounds_and_initial_values!(model)
    optimize!(model)
    
    if termination_status(model) in [MOI.OPTIMAL, MOI.LOCALLY_SOLVED]
        # Test that all variables are in reasonable ranges (trillion USD scale)
        
        theta_val = value(model[:theta])
        @test 0.1 < theta_val < 50
        
        utility_val = value(model[:UTILITY])
        @test 0 < utility_val < 1000
        
        obj_val = objective_value(model)
        @test -1000 < obj_val < 100
        
        for y in year_all
            Y_val = value(model[:Y][y])
            @test 10 < Y_val < 1000
            
            C_val = value(model[:C][y])
            @test 5 < C_val < 500
        end
        
        println("Master problem scaling check:")
        println("  theta: $(round(theta_val, digits=3)) trillion USD")
        println("  Utility: $(round(utility_val, digits=3))")
        println("  Objective: $(round(obj_val, digits=3))")
    end
end