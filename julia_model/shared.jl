# shared.jl - Equivalent to shared.gms
# Set definitions, parameters and variables shared across models

using JuMP
include("model_config.jl")

# SETS
# Time periods (10 years per period) in FaIR
t = 1:41

# Periods in energy and macro models (generated from config)
year_all = generate_year_sequence(default_config())

# Energy Sectors for macro-economic analysis in MACRO
sector = ["ELEC", "NELE"]

# Function to create shared variables in a JuMP model
function create_shared_variables!(model::Model; gbd_mode::Bool=false)
    # POSITIVE VARIABLES
    @variable(model, PHYSENE[sector, year_all] >= 0)  # Physical end-use service or commodity use
    
    # VARIABLES  
    # In GBD mode, skip COST_ANNUAL so macro_core.jl uses Taylor expansion instead of hard-linking
    if !gbd_mode
        @variable(model, COST_ANNUAL[year_all])  # costs per year_all
    end
    
    return model
end