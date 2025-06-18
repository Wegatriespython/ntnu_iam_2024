#!/usr/bin/env julia --project=.

"""
Simple AOT System Image Creator

Creates a system image for the current Julia project.
Either works or fails - no fallbacks.

Usage: julia --project=. create_sysimage.jl [name]
"""

using Pkg
using PackageCompiler
using TOML

function main()
    # Get sysimage name from command line or use default
    sysimage_name = length(ARGS) > 0 ? ARGS[1] : "OptimizationModel.so"
    
    println("Creating system image: $sysimage_name")
    
    # Read packages from Project.toml
    project_data = TOML.parsefile("Project.toml")
    packages = collect(keys(project_data["deps"]))
    
    println("Packages: $(join(packages, ", "))")
    
    # Create system image
    create_sysimage(
        packages;
        sysimage_path=sysimage_name,
        precompile_execution_file="precompile_workload.jl",
        project="."
    )
    
    println("✓ System image created: $sysimage_name")
    println("Usage: julia --sysimage=$sysimage_name --project=.")
end

if abspath(PROGRAM_FILE) == @__FILE__
    main()
end