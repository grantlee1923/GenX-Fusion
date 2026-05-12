#### This file runs the initial paper runs for the inter-annual variability analysis

####### Run greenfields and corresponding brownfields ########

using GenX
using JuMP
using CSV, DataFrames

task_id = parse(Int, ARGS[1])
num_tasks = parse(Int, ARGS[2])

year_list = 2001:2021

for idx in task_id+1:num_tasks:length(year_list)

    for i in 2001:2021

        include("/home/gridsan/nbhatt1/GenX/fusion_paper/paper_runs/2zone_greenfields/no_fusion/$i/Run.jl")

        endcap_data = CSV.read("/home/gridsan/nbhatt1/GenX/fusion_paper/paper_runs/results/greenfields/no_fusion/$i/capacity.csv", DataFrame)

        yearly_cap = endcap_data[1:end-1, "EndCap"]
        yearly_mwh = endcap_data[1:end-1, "EndEnergyCap"]

        gen_source = CSV.read("/home/gridsan/nbhatt1/GenX/fusion_paper/paper_runs/2zone_brownfields/data/no_fusion/Generators_data.csv", DataFrame)

        gen_source[!, "Existing_Cap_MW"] = yearly_cap
        gen_source[!, "Existing_Cap_MWh"] = yearly_mwh

        CSV.write("/home/gridsan/nbhatt1/GenX/fusion_paper/paper_runs/2zone_brownfields/data/no_fusion/Generators_data.csv", gen_source)

        ##### Modeling starts here #####
        case = "/home/gridsan/nbhatt1/GenX/fusion_paper/paper_runs/2zone_brownfields/data/no_fusion"

        # function get_settings_path(case::AbstractString)
        #     return joinpath(case, "Settings")
        # end

        # function get_settings_path(case::AbstractString, filename::AbstractString)
        #     return joinpath(get_settings_path(case), filename)
        # end

        # function get_default_output_folder(case::AbstractString)
        #     return joinpath(case, "Results")
        # end

        genx_settings = "/home/gridsan/nbhatt1/GenX/fusion_paper/paper_runs/2zone_brownfields/data/no_fusion/Settings/genx_settings.yml"
        # get_settings_path(case, "genx_settings.yml") #Settings YAML file path
        mysetup = configure_settings(genx_settings) # mysetup dictionary stores settings and GenX-specific parameters

        inputs_path = case
        settings_path = "/home/gridsan/nbhatt1/GenX/fusion_paper/paper_runs/2zone_brownfields/data/no_fusion/Settings"
        #get_settings_path(case)

        ### Cluster time series inputs if necessary and if specified by the user
        TDRpath = joinpath(case, mysetup["TimeDomainReductionFolder"])

        function time_domain_reduced_files_exist(tdrpath)
            tdr_load = isfile(joinpath(tdrpath,"Load_data.csv"))
            tdr_genvar = isfile(joinpath(tdrpath,"Generators_variability.csv"))
            tdr_fuels = isfile(joinpath(tdrpath,"Fuels_data.csv"))
            return (tdr_load && tdr_genvar && tdr_fuels)
        end

        if mysetup["TimeDomainReduction"] == 1
            if !time_domain_reduced_files_exist(TDRpath)
                println("Clustering Time Series Data (Grouped)...")
                cluster_inputs(inputs_path, settings_path, mysetup)
            else
                println("Time Series Data Already Clustered.")
            end
        end

        ### Configure solver
        println("Configuring Solver")
        OPTIMIZER = configure_solver(mysetup["Solver"], settings_path)


        #### Running a case

        ### Load inputs
        println("Loading Inputs")
        myinputs = load_inputs(mysetup, case)

        println("Generating the Optimization Model")
        EP = generate_model(mysetup, myinputs, OPTIMIZER)

        ########################
        #### Add any additional constraints
        HYDRO_RES = myinputs["HYDRO_RES"]
        dfGen = myinputs["dfGen"]

        # Empty arrays for indexing
        yr_start = Int[]
        yr_mid = Int[]

        # 20 year indexing
        for index in 1:20
            # Calculate the value for the beginning of the year
            start_year = (index-1) * 8760 + 1
            push!(yr_start, start_year)

            # Calculate the value for the middle of the year
            mid_year = (index-1) * 8760 + 2879
            push!(yr_mid, mid_year)
        end

        # Now define constraints using populated arrays
        @constraint(EP, cHydroSpring[y in HYDRO_RES, i in yr_mid], EP[:vS_HYDRO][y, i] .<= 0.55 .* EP[:eTotalCap][y] .* dfGen[y,:Hydro_Energy_to_Power_Ratio])
        @constraint(EP, cHydroJan[y in HYDRO_RES, i in yr_start], EP[:vS_HYDRO][y, i]  .== 0.70 .* EP[:eTotalCap][y] .* dfGen[y,:Hydro_Energy_to_Power_Ratio])
    

        ## Maine -> Quebec transmission limited to 2170MWe.
        # The line is defined as Quebec -> Maine in Network.csv, so these flows will be negative
        # Make sure to correc the line index if the order is changed in Network.csv
        @constraint(EP, cMaine2Quebec[t=1:myinputs["T"]], EP[:vFLOW][2, t] >= -170.0)

        # ## Solar <= 22GWe
        # solar_rid = findall(x -> startswith(x, "solar"), dfGen[!,:Resource])
        # @constraint(EP, cSolarCap, sum(EP[:eTotalCap][y] for y in solar_rid) <= 22e3)

        # ## Commercial Solar <= 15GWe
        # com_rid = findall(x -> startswith(x, "commercial"), dfGen[!,:Resource])
        # @constraint(EP, cComCap, sum(EP[:eTotalCap][y] for y in com_rid) <= 15e3)

        # ## Residential Solar <= 10GWe
        # res_rid = findall(x -> startswith(x, "residential"), dfGen[!,:Resource])
        # @constraint(EP, cResCap, sum(EP[:eTotalCap][y] for y in res_rid) <= 10e3)

        # ## Onshore wind <= 10GWe
        # onshore_rid = findall(x -> startswith(x, "onshore"), dfGen[!,:Resource])
        # @constraint(EP, cOnshoreCap, sum(EP[:eTotalCap][y] for y in onshore_rid) <= 10e3)

        # ## Offshore_fixed <= 37.5 GWe 
        # offshore_fixed = findall(x -> startswith(x, "fixed_offshore"), dfGen[!,:Resource])
        # @constraint(EP, cOffshoreFixCap, sum(EP[:eTotalCap][y] for y in offshore_fixed) <= 37.5e3)

        # ## Offshore_floating <= 275 GWe 
        # offshore_float = findall(x -> startswith(x, "float_offshore"), dfGen[!,:Resource])
        # @constraint(EP, cOffshoreFloatCap, sum(EP[:eTotalCap][y] for y in offshore_float) <= 275e3)

        ########################


        println("Solving Model")
        EP, solve_time = solve_model(EP, mysetup)
        myinputs["solve_time"] = solve_time # Store the model solve time in myinputs

        # Run MGA if the MGA flag is set to 1 else only save the least cost solution
        println("Writing Output")
        outputs_path = "/home/gridsan/nbhatt1/GenX/fusion_paper/paper_runs/results/brownfields/no_fusion/$i"
        elapsed_time = @elapsed write_outputs(EP, outputs_path, mysetup, myinputs)
        println("Time elapsed for writing is")
        println(elapsed_time)

    end
end