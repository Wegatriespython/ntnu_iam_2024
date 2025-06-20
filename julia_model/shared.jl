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

# POSITIVE VARIABLES
# Will be defined in JuMP model as:
# @variable(model, PHYSENE[sector, year_all] >= 0)  # Physical end-use service or commodity use

# VARIABLES  
# Will be defined in JuMP model as:
# @variable(model, COST_ANNUAL[year_all])  # costs per year_all