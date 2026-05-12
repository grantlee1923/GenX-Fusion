#### This file runs the initial paper runs for the inter-annual variability analysis

task_id = parse(Int, ARGS[1])
num_tasks = parse(Int, ARGS[2])

####### Run 20 scenarios ########
include("/home/gridsan/nbhatt1/GenX/fusion_paper/paper_runs/2zone_20yr/fusion/Run.jl")

include("/home/gridsan/nbhatt1/GenX/fusion_paper/paper_runs/2zone_20yr/no_fusion/Run.jl")