using GenX
using JuMP
using OrderedCollections
using DataFrames
using CSV

input_name = "updated_basecase_90NGCCS"
case_name = "updated_basecase_90NGCCS"

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
dfGen = myinputs["dfGen"]
dfFusion = myinputs["dfFusion"]
FUSION = myinputs["FUSION"]
MMBTU_PER_MWH = 3.412 ### hard coded MMBTU to MWh conversion

# Total load = 4,827,887,023 MWh across all 20 scenarios
# The scenarios vary in load by ~1%, so we'll treat the equally
# The emission intensity limits are: [4.0, 12.0, 50.0] gCO2 / kWh
# GenX requires the limit to be in millions tonnes (metric), so we'll convert by:
# g / tonne = 1e6
# kWh / MWh = 1e3
# total = 4,827,887,023[MWh] * limit[g/kWh] * 1e3[kWh/MWh] / 1e6[tonne/g] / 1e6[MMT / tonne]
# total = 4,827,887,023 * limit / 1e9
# emiss_lim_list = [2.5, 4.0, 12.0, 50.0]
emiss_lim_list = [13.6] #[4.0, 12.0, 50.0]

fusion_cost_list = [6000.0] #[8500.0, 3000.0, 6000.0, 12000.0]
# fusion_cost_list = [8500.0, 6000.0]

mysetup["CO2Cap"] = 2
scale_factor = mysetup["ParameterScale"] == 1 ? ModelScalingFactor : 1

mkpath(results_path)

# task_id = parse(Int,ARGS[1])
task_id = 0
# num_tasks = parse(Int,ARGS[2])
num_tasks = 1
# num_threads = parse(Int,ARGS[3])
num_threads = 16

set_optimizer_attribute(OPTIMIZER, "Threads", num_threads)

all_cases = vcat(collect(Iterators.product(emiss_lim_list, fusion_cost_list))...)

reduced_cases = []

# Find all the fusion resources in the model
fusion_rid = findall(x -> startswith(x, "fusion"), dfGen[!,:Resource])

TES_option = 0
TES_cha = 0
TES_dis = 0

for y in fusion_rid
    global TES_option = Int(dfFusion[y,:Add_Therm_Stor])
    global TES_cha = Int(dfFusion[y,:Therm_Stor_Cha_position])
    global TES_dis = Int(dfFusion[y,:Therm_Stor_Dis_position])
end

# Go through the cases and add any where !isfile(joinpath(outputs_path, "costs.csv"))
for idx in eachindex(all_cases)
    emiss_lim = all_cases[idx][1]
    fusion_cost = all_cases[idx][2]
    if mysetup["CO2Cap"] > 0
        outputs_path = joinpath(results_path, "TES_$(TES_option)$(TES_cha)$(TES_dis)_Cost_$(fusion_cost)_EmissLevel_$(emiss_lim)_gCO2perkWh")
    else
        outputs_path = joinpath(results_path, "TES_$(TES_option)$(TES_cha)$(TES_dis)_Cost_$(fusion_cost)")
    end
    if !isfile(joinpath(outputs_path, "costs.csv"))
        println("Including Case: $outputs_path")
        push!(reduced_cases, (emiss_lim, fusion_cost))
        rm(outputs_path, force=true, recursive=true)
    else
        println("Skipping Case (already exists): $outputs_path")
    end
end

for idx in task_id+1:num_tasks:length(reduced_cases)
    GC.gc()
    
    emiss_lim = reduced_cases[idx][1]
    fusion_cost = reduced_cases[idx][2]

    if mysetup["CO2Cap"] > 0
        outputs_path = joinpath(results_path, "TES_$(TES_option)$(TES_cha)$(TES_dis)_Cost_$(fusion_cost)_EmissLevel_$(emiss_lim)_gCO2perkWh")
        myinputs["dfMaxCO2Rate"][2] = emiss_lim / scale_factor ./ 1e3
        println("Emiss Limit: $emiss_lim, Fusion Cost: $fusion_cost, TES : $(TES_option)$(TES_cha)$(TES_dis)")
    else
        outputs_path = joinpath(results_path, "TES_$(TES_option)$(TES_cha)$(TES_dis)_Cost_$(fusion_cost)")
        println("Fusion Cost: $fusion_cost, TES : $(TES_option)$(TES_cha)$(TES_dis)")
    end

    discount_factor = 0.06
    FPP_lifetime = 40.0
    FPP_annuity = discount_factor / (1.0 - (1.0 + discount_factor)^(-FPP_lifetime))

    TES_Original_turb_cost = 1270.0   # SEE fusion.jl:145. Turbine cost is annualized with plant lifetime. Hence, this cost is inserted in Fusion_data.csv
    FPP_Original_turb_cost = 1700.0   # SEE fusion.jl:88. Turbine cost is annualized with plant lifetime. Hence, this cost is inserted in Fusion_data.csv
    TES_New_turb_cost = 1070.0   # SEE fusion.jl:146. Turbine cost is annualized with plant lifetime. Hence, this cost is inserted in Fusion_data.csv
    FPP_New_turb_cost = TES_New_turb_cost * FPP_Original_turb_cost/TES_Original_turb_cost   # SEE fusion.jl:81. Turbine cost is annualized with plant lifetime. Hence, this cost is inserted in Fusion_data.csv. We assume the same gross/net ratio for the turbine, but this remains to be verified.

    vessel_cost = 150.0          # SEE fusion.jl:........ Vessel cost is annualized with vessel lifetime and the guess on reactor utilization rate. Hence, this cost is inserted in Fusion_data.csv
    
    TES_energy_cost_kWh_e = 45.0*1 + 11.0*0 + 1716/1000/.4*0
    TES_power_cost_kW_e = 375.0*1
    indirect_TES_cha_position = 2           # indirect TES charges on the steam loop by default
    
    Original_Turbine_Eff_Frac = 0.40
    New_Turbine_Eff_Frac = 0.4931
    
    location_adjustment = 1.12
    fixed_cost_ratio = 0.15
    num_years = 1.0
    
    fusion_annual_cost = (fusion_cost .- FPP_Original_turb_cost .- vessel_cost) .* FPP_annuity .* num_years .* 1000 .* location_adjustment
    fusion_fixed_cost = fusion_cost .* FPP_annuity .* num_years .* 1000 .* fixed_cost_ratio
    
    for y in fusion_rid
        println("Running FPP at R_ID $y with investment costs: $(fusion_annual_cost) and fixed costs: $(fusion_fixed_cost)")

        if dfFusion[y,:Add_Therm_Stor] ∉ [0,1,2]
            error("Invalid Add_Therm_Stor value: $(dfFusion[y,:Add_Therm_Stor]). Must be <0, 1, or 2.")
        end
        
        if dfFusion[y,:Add_Therm_Stor] > 0
            # Checking TES structural assumptions
            if dfFusion[y,:Therm_Stor_Cha_position] ∉ [0,1,2]
                error("Invalid Therm_Stor_Cha_position value: $(dfFusion[y,:Therm_Stor_Cha_position]). Must be 0, 1, or 2.")
            elseif (dfFusion[y,:Therm_Stor_Cha_position] == 0)
                @warn "Currently, the TES charge from plant loop is assumed equivalent to TES charge from salt loop. You won't avoid the FPP salt losses."
            elseif (dfFusion[y,:Therm_Stor_Cha_position] == 1) .& (dfFusion[y,:Add_Therm_Stor]== 2)
                error("Indirect TES cannot charge into existing salt loop (Therm_Stor_Cha_position == 1)")
            elseif (dfFusion[y,:Therm_Stor_Cha_position] == 2)
                @warn "Charging TES on steam loop (Therm_Stor_Cha_position == 2) is suboptimal in most cases."
            end

            if dfFusion[y,:Therm_Stor_Dis_position] ∉ [0,1,2]
                error("Invalid Therm_Stor_Dis_position value: $(dfFusion[y,:Therm_Stor_Dis_position]). Must be 0, 1, or 2.")
            elseif dfFusion[y,:Therm_Stor_Dis_position] ∈ [1,2]
                @warn "The turbine technology has been changed due to TES discharge settings. New turbine efficiency applied: $(round(100 * dfFusion[y,:New_Turb_Eff_Frac], digits=2)) %. New turbine CAPEX applied: $(round(dfFusion[y,:FPP_New_Turb_CAPEX], digits=2)) \$/MW_e/yr."
                @warn "The new turbine cost is using the old gross/net ratio. This remains to be verified."
            end

            if dfFusion[y,:Therm_Stor_Dis_position] == 0
                # TES uses the original turbine technology
                TES_turb_eff = Original_Turbine_Eff_Frac
            elseif ((dfFusion[y,:Add_Therm_Stor] == 1) && (dfFusion[y,:Therm_Stor_Dis_position] == 1)) || ((dfFusion[y,:Add_Therm_Stor] == 2) && (dfFusion[y,:Therm_Stor_Dis_position] == 2))
                # TES uses the new technology
                TES_turb_eff = New_Turbine_Eff_Frac
            else
                error("Impossible combination of Add_Therm_Stor and Therm_Stor_Dis_position.")
            end

            if dfFusion[y,:Add_Therm_Stor] == 2
                # Indirect TES: force to charge on the steam loop. 0 can also be a choice (indirect TES on the plant fluid)
                dfFusion[y,:Therm_Stor_Cha_position] = indirect_TES_cha_position
            end
        else
            TES_turb_eff = 0
        end
        
        TES_energy_annual_cost = TES_energy_cost_kWh_e .* FPP_annuity .* num_years .* 1000 .* TES_turb_eff #.* location_adjustment
        TES_power_annual_cost = TES_power_cost_kW_e .* FPP_annuity .* num_years .* 1000 .* TES_turb_eff #.* location_adjustment

        println(TES_energy_annual_cost)
        println(TES_power_annual_cost)

        dfGen[y,:Inv_Cost_per_MWyr] = fusion_annual_cost 
        dfGen[y,:Fixed_OM_Cost_per_MWyr] = fusion_fixed_cost

        dfGen[y,:Heat_Rate_MMBTU_per_MWh] = MMBTU_PER_MWH / Original_Turbine_Eff_Frac

        dfFusion[y,:Plant_Life] = FPP_lifetime
        dfFusion[y,:Dis_Rate] = discount_factor

        dfFusion[y,:New_Turb_Eff_Frac] = New_Turbine_Eff_Frac
                
        dfFusion[y,:TES_Original_Turb_CAPEX] = TES_Original_turb_cost .* 1000 .* location_adjustment # See line 125
        dfFusion[y,:FPP_Original_Turb_CAPEX] = FPP_New_turb_cost .* 1000 .* location_adjustment # See line 126
        dfFusion[y,:TES_New_Turb_CAPEX] = TES_New_turb_cost .* 1000 .* location_adjustment # See line 127
        dfFusion[y,:FPP_New_Turb_CAPEX] = FPP_Original_turb_cost .* 1000 .* location_adjustment # See line 128

        dfFusion[y,:Inv_Vessel_per_MWe] = vessel_cost .* 1000 .* location_adjustment # See line 126

        dfFusion[y,:TES_Stor_Cost_per_MWht] = TES_energy_annual_cost
        dfFusion[y,:TES_HEX_Plant_Salt_annCAPEX] = TES_power_annual_cost ./ 2
        dfFusion[y,:TES_HEX_Steam_Salt_annCAPEX] = TES_power_annual_cost ./ 2
        dfFusion[y,:TES_HEX_Salt_Steam_annCAPEX] = TES_power_annual_cost ./ 2
        dfFusion[y,:TES_HEX_Salt_TESWF_annCAPEX] = TES_power_annual_cost ./ 2
    end

    # This check will cause the case to be skipped if the results already exist
    if isfile(joinpath(outputs_path, "costs.csv"))
        println("Skipping Case (already exists): $outputs_path")
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
    for year_num in 1:20
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


    ########################

    for y in FUSION
        if dfFusion[y, :Min_Reactor_Therm_MWth] > 0.0
            @constraint(EP, EP[:eFusionThermCap][y] >= dfFusion[y, :Min_Reactor_Therm_MWth])
        end
        if dfFusion[y, :Max_Reactor_Therm_MWth] > 0.0
            @constraint(EP, EP[:eFusionThermCap][y] <= dfFusion[y, :Max_Reactor_Therm_MWth])
        end
    end

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

    result_summ = DataFrame(Cost=objective_value(EP), Dual=0.0)
    CSV.write(joinpath(outputs_path, "fpp_results.csv"), result_summ)

end

