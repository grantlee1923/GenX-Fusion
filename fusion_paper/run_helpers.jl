using JuMP
using OrderedCollections
using DataFrames
using CSV

function get_settings_path(case::AbstractString)
    return joinpath(case, "Settings")
end

function get_settings_path(case::AbstractString, filename::AbstractString)
    return joinpath(get_settings_path(case), filename)
end

function get_default_output_folder(case::AbstractString)
    return joinpath(case, "Results")
end

function time_domain_reduced_files_exist(tdrpath)
    tdr_load = isfile(joinpath(tdrpath,"Load_data.csv"))
    tdr_genvar = isfile(joinpath(tdrpath,"Generators_variability.csv"))
    tdr_fuels = isfile(joinpath(tdrpath,"Fuels_data.csv"))
    return (tdr_load && tdr_genvar && tdr_fuels)
end