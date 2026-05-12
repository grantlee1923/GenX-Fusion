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

function write_lorenz_curve(path::AbstractString, inputs::Dict, setup::Dict, EP::Model)
    dfGen = inputs["dfGen"]
	T = inputs["T"]     			# Number of time steps (hours)
	Z = inputs["Z"]     			# Number of zones
	G = inputs["G"]     			# Number of generators


    gen = 1:G
    FLEX = inputs["FLEX"]
	NONFLEX = setdiff(collect(1:G), FLEX)

    # Calaculate the energy revenue

    energyrevenue = zeros(G, T)
	energyrevenue[NONFLEX, :] = value.(EP[:vP][NONFLEX, :]) .* transpose(dual.(EP[:cPowerBalance]) ./ inputs["omega"])[dfGen[NONFLEX, :Zone], :]
	if !isempty(FLEX)
		energyrevenue[FLEX, :] = value.(EP[:vCHARGE_FLEX][FLEX, :]).data .* transpose(dual.(EP[:cPowerBalance]) ./ inputs["omega"])[dfGen[FLEX, :Zone], :]
	end
	if setup["ParameterScale"] == 1
		energyrevenue *= ModelScalingFactor^2
	end

    # Sort the energy revenue and create two dataframes (one normalized, other not)
    dfer_sort = cumsum(sort(energyrevenue, dims=2), dims=2)
    df_sorted = DataFrame(transpose(dfer_sort), Symbol.(gen))

    dfer_sort_norm = dfer_sort ./ dfer_sort[:,end]
    df_sortnorm = DataFrame(transpose(dfer_sort_norm), Symbol.(gen))

    # Write the two different dataframes into csv files

    # Not normalized lorenz curve  
    CSV.write(joinpath(path, "LorenzCurveData.csv"), df_sorted)

    # Normalized Lorenz Curve data
    CSV.write(joinpath(path, "LorenzCurveDataNorm.csv"), df_sortnorm)
end
