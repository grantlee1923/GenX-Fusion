local_path = @__DIR__
scenarios_path = joinpath(local_path, "Scenarios")

for scenario in readdir("fusion_paper\\paper_runs\\dual_runs\\final_report_runs\\thermstor_1year\\Scenarios")
    println("Running scenario: $scenario")
    if scenario != "Scenario_00"
        scenario_path = joinpath(scenarios_path, "$scenario\\Run.jl")
        include(scenario_path)
    end
end