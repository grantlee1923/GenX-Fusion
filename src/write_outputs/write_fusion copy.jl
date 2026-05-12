"""
GenX: An Configurable Capacity Expansion Model
Copyright (C) 2021,  Massachusetts Institute of Technology
This program is free software; you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation; either version 2 of the License, or
(at your option) any later version.
This program is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.
A complete copy of the GNU General Public License v2 (GPLv2) is available
in LICENSE.txt.  Users uncompressing this from an archive may not have
received this license file.  If not, see <http://www.gnu.org/licenses/>.
"""

function write_fusion(path::AbstractString, inputs::Dict, setup::Dict, EP::Model)
    # dfFusion = inputs["dfFusion"]
    FUSION = inputs["FUSION"]
    T = inputs["T"]     # Number of time steps (hours)

    # Make ~/fusion directory if it doesn't exist
    fusion_results_path = joinpath(path, "fusion")
    if !isdir(fusion_results_path)
        mkpath(fusion_results_path)
    end

    for y in FUSION
        result_filename = string("fusion_time_", y, ".csv")
        fusion_timeseries = OrderedDict{String, Any}(
            "Reactor Commitment State" => :vFusionReactorCommit,
            "Turbine Commitment State" => :eFusionTurbCommit,
            "Reactor Thermal Output MWht" => :vThermOutput,
            "Turbine Thermal Input MWht" => :eTurbThermal,
            "Turbine Gross Electric Output MWhe" => :eTurbGrossElec,
            "Turbine Net Electric Output MWhe" => :eFusionNetElec,
            "FPP Imports MWhe" => :vFusionElecImports,
            "Total Recirc Power MWhe" => :eRecircPwr,
            "Recirc Salt Heating MWhe" => :vSaltElecHeating,
            "Fixed Recirc Power MWhe" => :ePlantFix,
            "Var Recirc Power MWhe" => :ePlantVar,
            "Tritium Inventory kg" => :vTritInventory,
            "Tritium Breeding kg" => :eTritBreeding,
            "Tritium Consumption kg" => :eTritConsumption,
            "Tritium Exports kg" => :vTritExports,
            "Tritium Decay kg" => :eTritDecay,
            "Tritium Leakage kg" => :eTritLeakage,
            "Deuterium Inventory kg" => :vDeuInventory,
            "Deuterium Imports kg" => :vDeuImports,
            "Deuterium Consumption kg" => :eDeuConsumption,
            "Deuterium Leakage kg" => :eDeuLeakage,
            "Thermal Storage Inventory MWht" => :vThermStor,
            "Thermal Storage Charging MWht" => :vThermChar,
            "Thermal Storage Discharging MWht" => :vThermDis,
            "Thermal Storage Net Discharge MWht" => :eThermStorNetDischarge
        )
        for (name, model_key) in fusion_timeseries
            if haskey(EP, model_key)
                fusion_timeseries[name] = get_fusion_data_vector(EP[model_key], y, T)
            else
                println("Can't find key: ", name, " in EP")
                fusion_timeseries[name] = zeros(T)
            end
        end
        fusion_df = DataFrame(fusion_timeseries)    
        CSV.write(joinpath(fusion_results_path, result_filename), fusion_df)
    end
end

function get_fusion_data_vector(model_data::JuMP.Containers.DenseAxisArray{AffExpr, 2, Tuple{Vector{Int64}, Base.OneTo{Int64}}, Tuple{JuMP.Containers._AxisLookup{Dict{Int64, Int64}}, JuMP.Containers._AxisLookup{Base.OneTo{Int64}}}}, R_ID::Int64, maxtime::Int64, mintime::Int64=1)
    return value.(model_data[R_ID, mintime:maxtime]).data
end

function get_fusion_data_vector(model_data::JuMP.Containers.DenseAxisArray{VariableRef, 2, Tuple{Vector{Int64}, Base.OneTo{Int64}}, Tuple{JuMP.Containers._AxisLookup{Dict{Int64, Int64}}, JuMP.Containers._AxisLookup{Base.OneTo{Int64}}}}, R_ID::Int64, maxtime::Int64, mintime::Int64=1)
    return value.(model_data[R_ID, mintime:maxtime]).data
end

function get_fusion_fleet_data(EP::JuMP.Model, model_key::Symbol, R_ID::Int64)
    return value(EP[model_key][R_ID])
end

function write_fusion_summary(path::AbstractString, inputs::Dict, setup::Dict, EP::Model)
    # dfFusion = inputs["dfFusion"]
    FUSION = inputs["FUSION"]
    # dfGen = inputs["dfGen"]

    fusion_results_path = joinpath(path, "fusion")
    if !isdir(fusion_results_path)
        mkpath(fusion_results_path)
    end

    fusion_var_names = OrderedDict{String, Any}(
        "Reactor Fusion Capacity MWt" => :eFusionFusionCap,
        "Reactor Thermal Capacity MWt" => :eFusionThermCap,
        "Turbine Gross Capacity MWe" => :eFusionTurbGrossCap,
        "Turbine Net Capacity MWe" => :eFusionTurbNetCap,
        "Num Reactors" => :vFusionNumReactors,
        "Num Turbines" => :eFusionNumTurbines,
        "Single Reactor Fusion Capacity MWt" => :eFusionFusionCapSize,
        "Single Reactor Thermal Capacity MWt" => :eFusionThermCapSize,
        "Single Turbine Gross Capacity MWe" => :eFusionTurbGrossCapSize,
        "Single Turbine Net Capacity MWe" => :eFusionTurbNetCapSize,
        "Total Thermal Output MWht" => :eThermOutputTot,
        "Total Net Electric Output MWhe" => :vP,
        "Tritium Storage Capacity kg" => :vTritCap,
        "Deuterium Storage Capacity kg" => :vDeuCap,
        "Thermal Storage Energy Capacity MWht" => :vThermStorCap,
        "Thermal Storage Discharge Capacity MWt" => :vThermDisCap,
        "Vessel Fixed Cost \$ per period" => :eVesselFixCosts,
        "Vessel Variable Cost \$ per period" => :eVesselVarCosts,
        "Turbine Fixed Cost \$ per period" => :eTurbFinalCost
    )

    fusion_variables = OrderedDict{String, Vector{Float64}}()
    for (name, _) in fusion_var_names
        fusion_variables[name] = []
    end

    file_name = string("fusion_summary.csv")

    for y in FUSION
        for (name, model_key) in fusion_var_names
            if haskey(EP, model_key)
                append!(fusion_variables[name], get_fusion_fleet_data(EP, model_key, y))
            else
                append!(fusion_variables[name], 0.0)
                println("Can't find key: ", name, " in EP")
            end
        end
        fusion_var_df = DataFrame(fusion_variables)    
        CSV.write(joinpath(fusion_results_path, file_name), fusion_var_df)
    end
end
