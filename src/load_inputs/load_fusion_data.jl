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

function load_fusion_data!(setup::Dict, path::AbstractString, inputs_ffuel::Dict)

	# All inputs are converted to Float64 by default.
	# The input columns listed here should be converted to Int64
	int_inputs = Vector{Symbol}([
		:Zone,
		:THERM,
		:FUSION,
		:Add_Therm_Stor,
		:Var_Turb_Cap
	])

    filename = "Fusion_data.csv"
	fusion_in = DataFrame(CSV.File(joinpath(path, filename), header=true, typemap=Dict(Int64 => Float64)), copycols=false)

	# Convert specified columns to Int64
	for col in int_inputs
		fusion_in[!,col] = Int64.(fusion_in[!,col])
	end

	# Add Resource IDs after reading to prevent user errors
	fusion_in[!,:R_ID] = 1:length(collect(skipmissing(fusion_in[!,1])))

	# Store DataFrame of generators/resources input data for use in model
	inputs_ffuel["dfFusion"] = fusion_in

	println(string(filename," Successfully Read!"))
end
