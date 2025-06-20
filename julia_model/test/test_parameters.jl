"""
Parameter Validation Tests
=========================

Tests for validating model parameters, dimensions, and data consistency.
These tests ensure all parameters are properly loaded and dimensionally consistent.
"""

using Test
using JuMP

# Import from module - all components already loaded
using EnergyMacroModel

@testset "Basic Parameter Loading" begin
    @test length(year_all) == 7
    @test year_all[1] == 2020
    @test year_all[end] == 2080
    @test length(sector) == 2
    @test "ELEC" ∈ sector
    @test "NELE" ∈ sector
    @test duration_period == 10
end

@testset "Configuration Structure" begin
    # Test default configuration
    config = default_config()
    @test config.base_year == 2020
    @test config.final_year == 2080
    @test config.gdp_base == 71.0
    @test config.drate == 0.05
    
    # Test custom configuration
    custom = custom_config(gdp_base = 80.0, drate = 0.04)
    @test custom.gdp_base == 80.0
    @test custom.drate == 0.04
    @test custom.base_year == 2020  # default values preserved
    
    # Test year sequence generation
    years = generate_year_sequence(config)
    @test years == [2020, 2030, 2040, 2050, 2060, 2070, 2080]
    
    # Test derived parameters calculation
    derived = calculate_derived_parameters(config)
    @test haskey(derived, :year_all)
    @test haskey(derived, :rho)
    @test haskey(derived, :k0)
    @test derived.rho ≈ (config.esub - 1) / config.esub
end

@testset "Economic Parameters" begin
    @test gdp_base > 0
    @test kgdp > 0  
    @test 0 < depr < 1
    @test 0 < drate < 1
    @test 0 < esub < 10
    @test 0 < kpvs < 1
    @test 0 < elvs < 1
    
    # Test derived parameters
    @test abs(rho - (esub - 1) / esub) < 1e-10
    @test k0 > 0
    @test y0 > 0
    @test c0 > 0
    @test i0 > 0
end

@testset "Energy Parameters" begin
    @test length(enestart) == length(sector) * length(year_all)
    @test length(eneprice) == length(sector) * length(year_all)
    @test length(cost_MESSAGE) == length(year_all)
    
    # Check all energy start values are positive
    for s in sector, y in year_all
        @test enestart[(s, y)] > 0
        @test eneprice[(s, y)] > 0
    end
    
    # Check MESSAGE costs are positive
    for y in year_all
        @test cost_MESSAGE[y] > 0
    end
end

@testset "AEEI Factors" begin
    @test length(aeei_factor) == length(sector) * length(year_all)
    
    # Base year should be 1.0
    for s in sector
        @test abs(aeei_factor[(s, 2020)] - 1.0) < 1e-10
    end
    
    # Subsequent years should be decreasing (efficiency improvement)
    for s in sector
        for i in 2:length(year_all)
            y_curr = year_all[i]
            y_prev = year_all[i-1]
            @test aeei_factor[(s, y_curr)] <= aeei_factor[(s, y_prev)]
            @test aeei_factor[(s, y_curr)] > 0
        end
    end
end

@testset "Labor and Growth Parameters" begin
    @test length(newlab) >= length(year_all) - 1
    @test length(udf) == length(year_all)
    @test length(growth_factor) == length(year_all)
    @test length(potential_gdp) == length(year_all)
    
    # Base year checks
    @test abs(growth_factor[2020] - 1.0) < 1e-10
    @test abs(potential_gdp[2020] - gdp_base) < 1e-10
    
    # Growth over time
    for i in 2:length(year_all)
        y_curr = year_all[i]
        y_prev = year_all[i-1]
        @test growth_factor[y_curr] >= growth_factor[y_prev]
        @test potential_gdp[y_curr] >= potential_gdp[y_prev]
        @test newlab[y_curr] > 0
        @test udf[y_curr] > 0
    end
end

@testset "Production Function Coefficients" begin
    @test a > 0
    @test b > 0
    @test isfinite(a)
    @test isfinite(b)
    
    # Test that coefficients satisfy the calibration equation
    # This is a crucial test - the production function should match base year data
    y_base = a * k0^(rho * kpvs) * 1.0^(rho * (1 - kpvs)) + 
             b * enestart[("ELEC", 2020)]^(rho * elvs) * enestart[("NELE", 2020)]^(rho * (1 - elvs))
    y_base = y_base^(1/rho)
    
    @test abs(y_base - y0) / y0 < 0.01
end

@testset "Dimensional Consistency" begin
    # Energy dimensions: PWh (petawatt-hours)
    # Money dimensions: Trillion USD for macro, Billion USD for MESSAGE
    
    @test typeof(gdp_base) <: Real
    @test typeof(k0) <: Real
    @test typeof(y0) <: Real
    
    # Check energy-economy scale compatibility
    total_energy_value = sum(enestart[(s, 2020)] * eneprice[(s, 2020)] for s in sector)
    @test total_energy_value < y0
    @test cost_MESSAGE[2020] / 1000 < y0
end

@testset "Parameter Ranges and Bounds" begin
    # Test that all parameters are within reasonable economic ranges
    
    @test 50 < gdp_base < 200
    @test 1 < kgdp < 5
    @test 0.01 < depr < 0.2
    @test 0.01 < drate < 0.2
    @test 0.1 < esub < 2
    
    # Energy parameter ranges
    for s in sector, y in year_all
        @test 0.1 < enestart[(s, y)] < 1000
        @test 0.001 < eneprice[(s, y)] < 1
    end
    
    for y in year_all
        @test 0.1 < cost_MESSAGE[y] < 50
    end
end

@testset "Finite Time Correction" begin
    @test length(finite_time_corr) == length(year_all)
    
    for y in year_all
        @test finite_time_corr[y] > 0
        @test isfinite(finite_time_corr[y])
        @test finite_time_corr[y] == abs(drate - grow[y])
    end
end

@testset "Data Consistency Checks" begin
    # Cross-validate different parameter calculations
    
    # Check that base year accounting identity holds
    @test abs(y0 - (c0 + i0 + ecst0)) < 0.01
    
    # Check that capital accumulation is consistent
    expected_i0 = k0 * (grow[2020] + depr)
    @test abs(i0 - expected_i0) < 0.1 || i0 > expected_i0
    
    # Check energy-economy linkage consistency
    for s in sector
        @test demand_base[s] == enestart[(s, 2020)]
    end
    
    # Check time consistency
    @test year_all == sort(year_all)
    for i in 2:length(year_all)
        @test year_all[i] - year_all[i-1] == duration_period
    end
end