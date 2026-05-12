using GenX
using JuMP
using OrderedCollections
using DataFrames
using CSV

input_name = "updated_basecase_cost8500"
case_name = "updated_basecase_cost8500"

case_path = @__DIR__
results_path = joinpath(case_path, "Results")

function gethomedir(case_path::String)
    path_split = splitpath(case_path)
    home_dir = ""
    for s in path_split
        if s == "fusion_paper"
            home_dir = joinpath(home_dir, s)
            break
        end
        home_dir = joinpath(home_dir, s)
    end

    return home_dir
end

# Find the home directory, to let us load the run_helpers.jl file
home_dir = gethomedir(case_path)
println(home_dir) 

## Load helper functions
include(joinpath(home_dir,"run_helpers.jl"))

## Define input and output paths
inputs_path = @__DIR__

## Load settings
genx_settings = get_settings_path(inputs_path, "genx_settings.yml") #Settings YAML file path
mysetup = configure_settings(genx_settings) # mysetup dictionary stores settings and GenX-specific parameters
settings_path = get_settings_path(inputs_path)

## Cluster time series inputs if necessary and if specified by the user
TDRpath = joinpath(inputs_path, mysetup["TimeDomainReductionFolder"])
if mysetup["TimeDomainReduction"] == 1
    if !time_domain_reduced_files_exist(TDRpath)
        println("Clustering Time Series Data (Grouped)...")
        cluster_inputs(inputs_path, settings_path, mysetup)
    else
        println("Time Series Data Already Clustered.")
    end
end

## Configure solver
println("Configuring Solver")
OPTIMIZER = configure_solver(mysetup["Solver"], settings_path)

# Turn this setting on if you run into numerical stability issues
# set_optimizer_attribute(OPTIMIZER, "BarHomogeneous", 1)

#### Running a case

## Load inputs
println("Loading Inputs")
myinputs = load_inputs(mysetup, inputs_path)

# Total load = 4,827,887,023 MWh across all 20 scenarios
# The scenarios vary in load by ~1%, so we'll treat the equally
# The emission intensity limits are: [4.0, 12.0, 50.0] gCO2 / kWh
# GenX requires the limit to be in millions tonnes (metric), so we'll convert by:
# g / tonne = 1e6
# kWh / MWh = 1e3
# total = 4,827,887,023[MWh] * limit[g/kWh] * 1e3[kWh/MWh] / 1e6[tonne/g] / 1e6[MMT / tonne]
# total = 4,827,887,023 * limit / 1e9
# emiss_lim_list = [2.5, 4.0, 12.0, 50.0] .* 4827887023.0 ./ 1e6
emiss_lim_list = [4.0, 50.0]


mysetup["CO2Cap"] = 2
scale_factor = mysetup["ParameterScale"] == 1 ? ModelScalingFactor : 1

# fusion_cap_list = vcat([0.0, 500.0, 1000.0], range(start=2500.0, stop=30000.0, step=2500.0))
fusion_cap_list = vcat(range(start=0.0, stop=2000.0, step=500.0), range(start=2500.0, stop=30000.0, step=2500.0))

mkpath(results_path)

task_id = parse(Int,ARGS[1])
num_tasks = parse(Int,ARGS[2])
num_threads = parse(Int,ARGS[3])

set_optimizer_attribute(OPTIMIZER, "Threads", num_threads)

# Get all cases as tuples of (emiss_lim, fusion_cap)
all_cases = vcat(collect(Iterators.product(emiss_lim_list, fusion_cap_list))...)

reduced_cases = []

# Go through the cases and add any where !isfile(joinpath(outputs_path, "costs.csv"))
for idx in task_id+1:num_tasks:length(all_cases)
    emiss_lim = all_cases[idx][1]
    fusion_cap = all_cases[idx][2]
    outputs_path = joinpath(results_path, "Dual_$(fusion_cap)mw_EmissLevel_$(emiss_lim)")
    if !isfile(joinpath(outputs_path, "costs.csv"))
        println("Including Case for emiss limit = $emiss_lim, fusion cap = $fusion_cap")
        push!(reduced_cases, (emiss_lim, fusion_cap))
        rm(outputs_path, force=true, recursive=true)
    end
end

for idx in task_id+1:num_tasks:length(all_cases)
    emiss_lim = all_cases[idx][1]
    fusion_cap = all_cases[idx][2]

    println("Emiss Limit: $emiss_lim, Fusion Cap: $fusion_cap")
    GC.gc()

    myinputs["dfMaxCO2Rate"][2] = emiss_lim / scale_factor ./ 1e3
    outputs_path = joinpath(results_path, "Dual_$(fusion_cap)mw_EmissLevel_$(emiss_lim)")

    # Find all the fusion resources in the model
    # and set their investment and fixed O&M costs to zero
    dfGen = myinputs["dfGen"]
    fusion_rid = findall(x -> startswith(x, "fusion"), dfGen[!,:Resource])
    for y in fusion_rid
        dfGen[y,:Inv_Cost_per_MWyr] = 0.0
        dfGen[y,:Fixed_OM_Cost_per_MWyr] = 0.0
    end

    # This check will cause the case to be skipped if the results already exist
    if isfile(joinpath(outputs_path, "costs.csv"))
        println("Skipping Case for emiss limit = " * string(emiss_lim) * " because it already exists.")
        continue
    end

    mkpath(dirname(outputs_path))
    
    ## Generate model
    println("Generating the Optimization Model")
    EP = generate_model(mysetup, myinputs, OPTIMIZER)

    ########################
    #### Add any additional constraints
    HYDRO_RES = myinputs["HYDRO_RES"]

    # Empty arrays for indexing
    jan1_idxs = Int[]
    may1_idxs = Int[]

    # 20 year indexing
    for year_num in 1:1
        # Calculate the index for the beginning of the years
        start_year = (year_num-1) * 8760 + 1
        push!(jan1_idxs, start_year)

        # Calculate the index for the middle of the years
        mid_year = (year_num-1) * 8760 + 2879
        push!(may1_idxs, mid_year)
    end

    ## Hydro storage == 0.70 * Existing Capacity at the start of the year
    @constraint(EP, cHydroJan[y in HYDRO_RES, jan1_idx in jan1_idxs], EP[:vS_HYDRO][y, jan1_idx]  .== 0.70 .* EP[:eTotalCap][y] .* dfGen[y,:Hydro_Energy_to_Power_Ratio])
    
    ## Hydro storage <= 0.55 * Existing Capacity at start of May 1st 
    @constraint(EP, cHydroSpring[y in HYDRO_RES, may1_idx in may1_idxs], EP[:vS_HYDRO][y, may1_idx] .<= 0.55 .* EP[:eTotalCap][y] .* dfGen[y,:Hydro_Energy_to_Power_Ratio])
        
    ## Maine -> Quebec transmission limited to 2170MWe.
    # The line is defined as Quebec -> Maine in Network.csv, so these flows will be negative
    # Make sure to correc the line index if the order is changed in Network.csv
    @constraint(EP, cMaine2Quebec[t=1:myinputs["T"]], EP[:vFLOW][2, t] >= -170.0)

    ## Fusion <= fusion_cap
    @constraint(EP, cFusionCap, sum(EP[:eTotalCap][y] for y in fusion_rid) <= fusion_cap)
    
    ########################

    ## Solve model
    println("Solving Model")
    EP, solve_time = solve_model(EP, mysetup)
    myinputs["solve_time"] = solve_time # Store the model solve time in myinputs

    ## Run MGA if the MGA flag is set to 1 else only save the least cost solution
    println("Writing Output")
    # outputs_path = get_default_output_folder(outputs_path)

    ## Write outputs
    write_outputs(EP, outputs_path, mysetup, myinputs)

    result_summ = DataFrame(Cost=objective_value(EP), Dual=dual(EP[:cFusionCap]))
    CSV.write(joinpath(outputs_path, "fpp_results.csv"), result_summ)

end

