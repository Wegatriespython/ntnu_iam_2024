#!/usr/bin/env julia
"""
GBD Energy-Macro Model Test Suite
=================================

Comprehensive test suite for the Generalized Benders Decomposition reformulation
of the Energy-Macro (MESSAGE+MACRO) model.

This test suite systematically verifies:
1. Mathematical consistency of the formulation
2. Component-wise functionality (energy subproblem, master problem)
3. Benders decomposition implementation
4. Numerical stability and parameter validation
5. Convergence properties and solution quality

Usage: julia runtests.jl
Or: julia -e "using Pkg; Pkg.test()"
"""

using Test
using EnergyMacroModel
using JuMP, Ipopt
using LinearAlgebra
using JSON

@testset "GBD Energy-Macro Model Test Suite" begin
    
    @testset "Parameter Validation" begin
        include("test_parameters.jl")
    end
    
    @testset "Energy Subproblem" begin
        include("test_energy_subproblem.jl")
    end
    
    @testset "Master Problem (MACRO)" begin
        include("test_master_problem.jl")
    end
    
    @testset "Benders Cut Implementation" begin
        include("test_benders_cuts.jl")
    end
    
    @testset "Numerical Stability" begin
        include("test_numerical_stability.jl")
    end
    
    @testset "Integration Tests" begin
        include("test_integration.jl")
    end
    
end

println("\n" * "="^60)
println("GBD Test Suite Complete")
println("="^60)