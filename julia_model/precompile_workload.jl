#!/usr/bin/env julia --project=.

# Minimal precompile workload - just load everything available

using Pkg
project_deps = keys(Pkg.project().dependencies)

# Load all packages from Project.toml
for pkg in project_deps
    @eval using $(Symbol(pkg))
end

# Load all model files
for file in readdir(".")
    if endswith(file, ".jl") && file ∉ ["precompile_workload.jl", "create_sysimage.jl"]
        try
            include(file)
        catch
            # Skip files that can't be included
        end
    end
end

# Load modules from src/
if isdir("src")
    for file in readdir("src") 
        if endswith(file, ".jl")
            module_name = replace(file, ".jl" => "")
            try
                @eval using $(Symbol(module_name))
            catch
                # Skip modules that can't be loaded
            end
        end
    end
end

println("Precompile workload completed")