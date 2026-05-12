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

function fusioninvestment_fpp!(EP::Model, inputs::Dict, setup::Dict)
    # Adds variables and constraints describing the 
    # installed capacity and cost of fusion power plant and their components.
    # For now, we are assuming a fixed reactor pulse ratio

    ##############################################################################

    ## Relevant inputs
    # dfGen[y,:Cap_Size] is the net capacity of the standard 1-to-1 turbine
    # dfGen[y,:Blanket_Power_Mult] = (FPP thermal output) / (FPP fusion output)
 
    ## Relevant pre-existing variables and expressions
    # EP[:eTotalCap] , set in investment_discharge.jl = total net capacity of the FPP fleet

    ##############################################################################

    dfGen = inputs["dfGen"]
    dfFusion = inputs["dfFusion"]
    FUSION = inputs["FUSION"]
    
    MMBTU_PER_MWH = 3.412 ### hard coded MMBTU to MWh conversion
    # FIXME -> we may need to multiply this by parameter scale

    ##############################################################################

    ## FPP fleet net electic installed capacity (MWe)
    @expression(EP, eFusionTurbNetCap[y in FUSION], EP[:eTotalCap][y])

    @expression(EP, eFusionTurbNetCapSize[y in FUSION], dfGen[y,:Cap_Size])

    ## Number of standard turbines
    # To prevent the problem from being quadratic, 
    # we will model the additional turbine capacity
    # as adding additional turbines of the 'standard' turbine size
    @expression(EP, eFusionNumTurbines[y in FUSION], EP[:eTotalCap][y] / EP[:eFusionTurbNetCapSize][y])

    ## Number of FPP reactors
    @variable(EP, 0 <= vFusionNumReactors[y in FUSION])

    # Minimum and maximum ratios between the number of turbines and the number of reactors
    for y in FUSION
        if dfFusion[y, :Turb_Min_Ratio] > 0
            @constraint(EP, eFusionNumTurbines[y] >= vFusionNumReactors[y] * dfFusion[y, :Turb_Min_Ratio])
        end
        if dfFusion[y, :Turb_Max_Ratio] > 0
            @constraint(EP, eFusionNumTurbines[y] <= vFusionNumReactors[y] * dfFusion[y, :Turb_Max_Ratio])
        end
    end 
    
    # Var_Turb_Cap controls whether the turbine size can be changed.
    # If it is not allowed for any generators, 
    # then constrain the turbine capacity to the standard capacity
    # and make turbine always work in lock-step with the reactor
    # TODO: We could allow them to stay separate in the future
    variable_turbine_flags = dfFusion[!,:Var_Turb_Cap]
    for y in FUSION
        if variable_turbine_flags[y] == 0
            @constraint(EP, vFusionNumReactors[y] == eFusionNumTurbines[y])
        end
    end

    # Sets the FPP turbine efficiency and cost according to TES settings
    @expression(EP, eFPP_turb_eff[y in FUSION], (       (dfFusion[y,:Add_Therm_Stor] == 1) * (dfFusion[y,:Therm_Stor_Dis_position] == 1)) * dfFusion[y,:New_Turb_Eff_Frac] +                        # Direct TES discharging in new turbine: need to change the efficiency of the FPP main turbine
                                                    (1 -(dfFusion[y,:Add_Therm_Stor] == 1) * (dfFusion[y,:Therm_Stor_Dis_position] == 1)) * MMBTU_PER_MWH ./ dfGen[y,:Heat_Rate_MMBTU_per_MWh])     # Other cases: FPP efficiency is the original one
    @expression(EP, eFPP_Turb_Cost[y in FUSION],(       (dfFusion[y,:Add_Therm_Stor] == 1) * (dfFusion[y,:Therm_Stor_Dis_position] == 1)) * dfFusion[y,:FPP_New_Turb_CAPEX] +                       # idem for the FPP turbine costs
                                                    (1 -(dfFusion[y,:Add_Therm_Stor] == 1) * (dfFusion[y,:Therm_Stor_Dis_position] == 1)) * dfFusion[y,:FPP_Original_Turb_CAPEX]
                * dfFusion[y,:Dis_Rate] / (1 - (1 + dfFusion[y,:Dis_Rate])^(-dfFusion[y,:Plant_Life])))                                                                                             # CAPEX annualization

    # We want to calculate the gross electric capacity from the net capacity and other parameters
    # We are assuming that the system installs an unknown numer of  turbines of one size, 
    # rather than a fixed number of turbines of unknown size. So there may be more turbines than reactors.
    # However, we know that each turbine is sized to match one reactor, including its recirculating power
    #
    # Therefore, for one reactor:
    # Gross_1 = Net_1 + Max_Recirc_1
    # Max_Recirc_1 = mag_cool + plant_fixed + plant_var * Thermal_1
    # Gross_1 = turb_eff * Thermal_1
    # Net_1 + Max_Recirc_1 = turb_eff * Thermal_1
    # Net_1 + mag_cool + plant_fixed + plant_var * Thermal_1 = turb_eff * Thermal_1
    # Thermal_1 = (Net_1 + mag_cool + plant_fixed) / (turb_eff - plant_var)

    ## Reactor thermal capacity (MWt)
    @expression(EP, eFusionThermCapSize[y in FUSION],
        (EP[:eFusionTurbNetCapSize][y] + dfFusion[y,:Mag_Cool] + dfFusion[y,:Plant_fixed_power_MWe])
        / 
        (eFPP_turb_eff[y] - dfFusion[y,:Plant_var_power_MWeperMWt])
    )

    ## FPP fleet thermal installed capacity (MWt)
    @expression(EP, eFusionThermCap[y in FUSION], vFusionNumReactors[y] * eFusionThermCapSize[y])

    ## Maximum recirculating power of a single FPP
    @expression(EP, eMaxRecircPwr[y in FUSION], 
        + dfFusion[y,:Mag_Cool] 
        + dfFusion[y,:Plant_fixed_power_MWe]
        + dfFusion[y,:Plant_var_power_MWeperMWt] * eFusionThermCapSize[y]
    )

    ## Maximum recirculating power of the FPP fleet in each zone
    @expression(EP, eFleetMaxRecircPwr[y in FUSION],
        eMaxRecircPwr[y] * vFusionNumReactors[y]
    )

    pulse_ratio = fusionpulseratio.(dfFusion[!,:Pulse_length_min], dfFusion[!,:Pulse_gap_min]) # Fraction of time spent on

    ## Reactor fusion capacity (MWt)
    # This is where we account for fact that the reactor is pulsing
    @expression(EP, eFusionFusionCapSize[y in FUSION], eFusionThermCapSize[y] / dfFusion[y,:Blanket_Power_Mult] / pulse_ratio[y])

    ## FPP fleet fusion power installed capacity (MWt)
    @expression(EP, eFusionFusionCap[y in FUSION], eFusionFusionCapSize[y] * vFusionNumReactors[y])
    
    ## Gross electric capacity of a single FPP turbine
    @expression(EP, eFusionTurbGrossCapSize[y in FUSION],
        eFusionThermCapSize[y] * eFPP_turb_eff[y] 
    )
    ## Alternative formulation:
    # @expression(EP, eFusionTurbGrossCapSize[y in FUSION],
        # EP[:eFusionTurbNetCapSize][y] + eMaxRecircPwr[y]
    # )

    @expression(EP, eFusionTurbGrossCap[y in FUSION], vFusionNumReactors[y] * eFusionTurbGrossCapSize[y])
    # @expression(EP, eTESTurbNetCap[y in FUSION], (eFusionNumTurbines[y] .- vFusionNumReactors[y]) * eFusionTurbNetCapSize[y])

    # Calculate the turbine cost ($/MWe). If there is a TES, the turbine cost depends on the TES discharge setting.
    @expression(EP, eTES_Turb_Cost[y in FUSION], ((dfFusion[y,:Therm_Stor_Dis_position] == 0) * dfFusion[y,:TES_Original_Turb_CAPEX] +              # TES discharging in the original (oversized) turbine
                                                    (dfFusion[y,:Therm_Stor_Dis_position] > 0) * dfFusion[y,:TES_New_Turb_CAPEX])                     # TES discharging in a new technology turbine
                                                    * dfFusion[y,:Dis_Rate] / (1 - (1 + dfFusion[y,:Dis_Rate])^(-dfFusion[y,:Plant_Life])))           # Annualization of the CAPEX selected

    # The overshoots will be zero-vectors on lines without TES
    @expression(EP, eFusionTESInvFPPCostOvershoot[y in FUSION], (dfFusion[y,:Add_Therm_Stor] > 0) * (vFusionNumReactors[y] - eFusionNumTurbines[y]) * eFusionTurbNetCapSize[y] * dfGen[y,:Inv_Cost_per_MWyr])
    @expression(EP, eFusionTESFixedFPPCostOvershoot[y in FUSION], (dfFusion[y,:Add_Therm_Stor] > 0) * (vFusionNumReactors[y] - eFusionNumTurbines[y]) * eFusionTurbNetCapSize[y] * dfGen[y,:Fixed_OM_Cost_per_MWyr])

    # Calculate the total cost of the turbines
    # The cost is input as a net-MWe cost. If we have a TES, we have more turbines than reactors: remove the difference between turbine net cost for FPP and for TES (FPP cost includes recirculating power)
    @warn "There may be something wrong here because the number of reactor do not match exactly the number of turbines without TES so it creates a little delta in TES cases"
    @expression(EP, eTurbFinalCost[y in FUSION], eFPP_Turb_Cost[y] * eFusionNumTurbines[y] * eFusionTurbNetCapSize[y] +                                                                                                 # Cost without TES
                                                (dfFusion[y,:Add_Therm_Stor] > 0) * ((vFusionNumReactors[y] - eFusionNumTurbines[y]) * eFusionTurbNetCapSize[y] * (eFPP_Turb_Cost[y] - eTES_Turb_Cost[y])))             # Cost-adjustement with TES

    for y in FUSION
        update_fixed_cost!(EP, EP[:eTurbFinalCost][y], y)
        update_fixed_cost!(EP, EP[:eFusionFusionCap][y] * 1e-9, y)
        if dfFusion[y,:Add_Therm_Stor] > 0
            # Remove the overshoot FPP costs if there is a TES
            update_fixed_cost!(EP, eFusionTESInvFPPCostOvershoot[y] .+ eFusionTESFixedFPPCostOvershoot[y], y)
        end
    end
end

function fusionthermalstorage_instantiation!(EP::Model, inputs::Dict, setup::Dict)
    # Adds variables and constraints describing the
    # commitment state of fusion reactors and turbines

    ##############################################################################

    dfGen = inputs["dfGen"]
    dfFusion = inputs["dfFusion"]
    FUSION = inputs["FUSION"]

    T = inputs["T"]     # Number of time steps (hours)
    MMBTU_PER_MWH = 3.412 ### hard coded MMBTU to MWh conversion

    ##############################################################################

    ## Instantiate these variables in all cases
    @variable(EP, vThermStorCap[y in FUSION], lower_bound = 0.)
    @variable(EP, vThermDisCap[y in FUSION], lower_bound = 0.)
    @variable(EP, vThermChaCap[y in FUSION], lower_bound = 0.)

    @variable(EP, vThermStor[y in FUSION, t=1:T], lower_bound = 0.)
    @variable(EP, vThermDis[y in FUSION, t=1:T], lower_bound = 0.)
    @variable(EP, vThermChar[y in FUSION, t=1:T], lower_bound = 0.)

    @expression(EP, eTES_turb_eff[y in FUSION], (dfFusion[y,:Therm_Stor_Dis_position] == 0) * MMBTU_PER_MWH / dfGen[y,:Heat_Rate_MMBTU_per_MWh] +
                                                (dfFusion[y,:Therm_Stor_Dis_position] > 0) * dfFusion[y,:New_Turb_Eff_Frac]
                                                )
    @expression(EP, eTES_cha_eff[y in FUSION],  (dfFusion[y,:Therm_Stor_Cha_position] == 0) * dfFusion[y,:TES_HEX_Plant_Salt_Eff] +                    # TES charging from the plant
                                                          (dfFusion[y,:Therm_Stor_Cha_position] == 1) * (dfFusion[y,:Add_Therm_Stor] == 1) * 1 +                 # Direct TES charging from the existing salt loop without losses
                                                          (dfFusion[y,:Therm_Stor_Cha_position] == 2) * dfFusion[y,:TES_HEX_Steam_Salt_Eff]                      # TES charging from the steam loop
                                                          )
    @expression(EP, eTES_dis_eff[y in FUSION],  (dfFusion[y,:Therm_Stor_Dis_position] == 0) * (dfFusion[y,:Add_Therm_Stor] == 1) * dfFusion[y,:TES_HEX_Salt_Steam_Eff_direct] +                    # Direct TES discharging into steam
                                                          (dfFusion[y,:Therm_Stor_Dis_position] == 0) * (dfFusion[y,:Add_Therm_Stor] == 2) * dfFusion[y,:TES_HEX_Salt_Steam_Eff_indirect] +                  # Indirect TES discharging into steam: account for possible losses because of BPV
                                                          (dfFusion[y,:Therm_Stor_Dis_position] > 0) * dfFusion[y,:TES_Salt_TESWF_Eff]
                                                          )
end

function fusioninvestment_thermstor!(EP::Model, inputs::Dict, setup::Dict)

    ##############################################################################

    dfGen = inputs["dfGen"]
    dfFusion = inputs["dfFusion"]

    FUSION = inputs["FUSION"]
    FUSION_ThermStor = intersect(FUSION, dfFusion[dfFusion[!,:Add_Therm_Stor].>0,:R_ID])
    FUSION_noThermStor = intersect(FUSION, dfFusion[dfFusion[!,:Add_Therm_Stor].==0,:R_ID])

    #Add_Therm_Stor = 0 is no TES
    #Add_Therm_Stor = 1 is a direct TES
    #Add_Therm_stor = 2 is an indirect TES

    #Therm_Stor_Cha_position = 0 TES charges from plant loop (requires a salt loop)
    #Therm_Stor_Cha_position = 1 TES is plugged INTO an existing salt loop (only available for direct TES)
    #Therm_Stor_Cha_position = 2 TES charges from steam loop (requires a salt loop)

    #Therm_Stor_Dis_position = 0 TES+FPP discharge via the original (oversized) steam turbine loop
    #Therm_Stor_Dis_position = 1 TES+FPP discharge via a new turbine fed with working fluid TESWF (can be air, sCO2, etc.) (only available for direct TES)
    #Therm_Stor_Dis_position = 2 TES discharges via a new turbine fed with working fluid TESWF (can be air, sCO2, etc.) and FPP discharges via the original steam turbine loop (only available for indirect TES)
    
    MMBTU_PER_MWH = 3.412 ### hard coded MMBTU to MWh conversion

    ##############################################################################

    # Setting the correct values for the TES
    @expression(EP, eTES_cha_cost[y in FUSION], (dfFusion[y,:Therm_Stor_Cha_position] == 0) * dfFusion[y,:TES_HEX_Plant_Salt_annCAPEX] +               # idem with costs
                                                          (dfFusion[y,:Therm_Stor_Cha_position] == 1) * (dfFusion[y,:Add_Therm_Stor] == 1) * 0 +                 # Direct TES charging from the existing salt loop without added cost
                                                          (dfFusion[y,:Therm_Stor_Cha_position] == 2) * dfFusion[y,:TES_HEX_Steam_Salt_annCAPEX]
                                                          )
    
    @warn "Double check if we should consider the HEX cost or the difference between this HEX cost and the original HEX (already included in the turbine cost)"
    @expression(EP, eTES_dis_cost[y in FUSION], (dfFusion[y,:Therm_Stor_Dis_position] == 0) * dfFusion[y,:TES_HEX_Salt_Steam_annCAPEX] +
                                                          (dfFusion[y,:Therm_Stor_Dis_position] > 0) * dfFusion[y,:TES_HEX_Salt_TESWF_annCAPEX]
                                                          )
    @expression(EP, eTES_stor_cost[y in FUSION], dfFusion[y,:TES_Stor_Cost_per_MWht])

    @expression(EP, eDeltaCapexSH[y in FUSION], EP[:eTES_cha_eff][y] * EP[:eTES_dis_eff][y] * EP[:eTES_turb_eff][y] * (eTES_dis_cost[y] + EP[:eTES_Turb_Cost][y]) - dfFusion[y,:HEX_Plant_Steam_Eff] * (MMBTU_PER_MWH / dfGen[y,:Heat_Rate_MMBTU_per_MWh]) * (dfFusion[y,:HEX_Plant_Steam_annCAPEX] + dfFusion[y,:TES_Original_Turb_CAPEX]))
    @expression(EP, eTESSaltLoopCost[y in FUSION], (dfFusion[y, :Add_Therm_Stor] == 1) * (dfFusion[y, :Therm_Stor_Cha_position] == 0) * (eTES_cha_cost[y] + eDeltaCapexSH[y]) * EP[:eFusionThermCap][y])      # Non-null only for direct TES without salt-loop by design

 
    ## Thermal storage energy capacity [MWht]

    ## Thermal storage hourly charge and discharge capacity [MWt]

    for y in FUSION_ThermStor
        # Constrain the thermal storage energy capacity
        if dfFusion[y, :Min_Therm_Stor_MWh] > 0.0
            set_lower_bound(EP[:vThermStorCap][y], dfFusion[y, :Min_Therm_Stor_MWh])
        end
        if dfFusion[y, :Max_Therm_Stor_MWh] > 0.0
            set_upper_bound(EP[:vThermStorCap][y], dfFusion[y, :Max_Therm_Stor_MWh])
        end

        # Constrain the thermal storage discharge capacity
        if dfFusion[y, :Min_Therm_Stor_Dis_MWth] > 0.0
            set_lower_bound(EP[:vThermDisCap][y], dfFusion[y, :Min_Therm_Stor_Dis_MWth])
        end
        if dfFusion[y, :Max_Therm_Stor_Dis_MWth] > 0
            set_upper_bound(EP[:vThermDisCap][y], dfFusion[y, :Max_Therm_Stor_Dis_MWth])
        end
        
        # Constrain the thermal storage charge capacity
        if dfFusion[y, :Min_Therm_Stor_Cha_MWth] > 0.0
            set_lower_bound(EP[:vThermChaCap][y], dfFusion[y, :Min_Therm_Stor_Cha_MWth])
        end
        if dfFusion[y, :Max_Therm_Stor_Cha_MWth] > 0
            set_upper_bound(EP[:vThermChaCap][y], dfFusion[y, :Max_Therm_Stor_Cha_MWth])
        end

        # Constrain the ratio of thermal storage energy capacity to discharge capacity
        if dfFusion[y, :Min_Therm_Dur] > 0
            @constraint(EP, EP[:vThermStorCap][y] >= EP[:vThermDisCap][y] * dfFusion[y, :Min_Therm_Dur])
        end
        if dfFusion[y, :Max_Therm_Dur] > 0
            @constraint(EP, EP[:vThermStorCap][y] <= EP[:vThermDisCap][y] * dfFusion[y, :Max_Therm_Dur])
        end

        # Constrain the symmetry of the TES
        if dfFusion[y, :Sym_Therm_Stor] > 0
            @constraint(EP, EP[:vThermChaCap][y] == EP[:vThermDisCap][y])
            @warn "The TES is constrained to be symmetrical, this is a strong assumption"
        end
    end

    for y in FUSION_noThermStor
        # Constrain the thermal storage energy capacity to 0
        set_lower_bound(EP[:vThermStorCap][y], 0)
        set_upper_bound(EP[:vThermStorCap][y], 0)

        set_lower_bound(EP[:vThermDisCap][y], 0)
        set_upper_bound(EP[:vThermDisCap][y], 0)
        
        set_lower_bound(EP[:vThermChaCap][y], 0)
        set_upper_bound(EP[:vThermChaCap][y], 0)
    end

    ## Fixed costs associated with the thermal storage. We separate the TES cost from the turbine cost. The turbine cost is already included in fusion_investment_fpp! function above
    @expression(EP, eTESCost[y in FUSION],        eTES_stor_cost[y] * EP[:vThermStorCap][y] +                                                                                                                                  # This part of TES cost is architecture-independent
                                                            (dfFusion[y, :Add_Therm_Stor] == 1) * (eTESSaltLoopCost[y] + eTES_dis_cost[y] * (EP[:vThermDisCap][y] - EP[:eTES_cha_eff][y] * EP[:eTES_dis_eff][y] * EP[:eFusionThermCap][y])) +      # This part is specific to direct TES, SaltLoopCost is non-null only for direct TES without salt loop by-design
                                                            (dfFusion[y, :Add_Therm_Stor] == 2) * (eTES_cha_cost[y] * EP[:vThermChaCap][y] + eTES_dis_cost[y] * EP[:vThermDisCap][y])                                                         # This part is specific to indirect TES
                                                            )
    # The TES turbine cost is already accounted in the fusioninvesment_fpp! function. This is only an expression for future use.
    @expression(EP, eTESTurbineCost[y in FUSION], ((dfFusion[y, :Add_Therm_Stor] == 1) * (EP[:vThermDisCap][y] - EP[:eTES_cha_eff][y] * EP[:eTES_dis_eff][y] * EP[:eFusionThermCap][y]) +                                                # This part is specific to direct TES
                                                             (dfFusion[y, :Add_Therm_Stor] == 2) * (EP[:vThermDisCap][y])                                                                                                                # This part is specific to indirect TES
                                                            ) * EP[:eTES_Turb_Cost][y])
    
    ## Add thermal storage fixed costs to the fusion generator costs
    for y in FUSION
        update_fixed_cost!(EP, eTESCost[y], y)
    end
end

function fusioninvestment_fuel!(EP::Model, inputs::Dict, setup::Dict)

    ##############################################################################

    dfFusion = inputs["dfFusion"]
    FUSION = inputs["FUSION"]

    ##############################################################################

    ## Tritium storage capacity
    @variable(EP, vTritCap[y in FUSION] >= 0.0)

    ## Deuterium storage capacity
    @variable(EP, vDeuCap[y in FUSION] .>= 0.0)

    ### Add Tritium and Deuterium Storage Capacity to the Fixed Costs

    # Pull values of cost from dfFusion
    trit_stor_cost = dfFusion[!,:Trit_Stor_Cost]
    deu_stor_cost = dfFusion[!,:Deu_Stor_Cost]

    # Multiply by the storage capacity and add to the fixed costs
    # fuel_fixed_costs = vTritCap[FUSION] .* trit_stor_cost[FUSION] .+ vDeuCap[FUSION] .* deu_stor_cost[FUSION]
    @expression(EP, eFuel_fixed_costs[y in FUSION], vTritCap[y] .* trit_stor_cost[y] .+ vDeuCap[y] .* deu_stor_cost[y])

    for y in FUSION
        update_fixed_cost!(EP, eFuel_fixed_costs[y], y)
    end
end

function fusioninvestment_vessel(EP::Model, inputs::Dict, setup::Dict)

    ##############################################################################

    dfFusion = inputs["dfFusion"]

    # plant_life = dfFusion[!,:Plant_Life]
    vessel_name = dfFusion[!,:Vessel_Life]
    replace_dur = dfFusion[!,:Rep_Dur]
    discount = dfFusion[!,:Dis_Rate]
    vessel_inv = dfFusion[!,:Inv_Vessel_per_MWe]
    # plant_cost = dfGen[!,:Inv_Cost_per_MWyr]
    reactor_util_guess = dfFusion[!,:Reactor_Utilization_Guess]

    FUSION = inputs["FUSION"]
    FUSION_VESSEL = intersect(FUSION, dfFusion[dfFusion[!,:VV_Util_Deg].>0,:R_ID])

    ##############################################################################

    @expression(EP, eC1[y in FUSION_VESSEL], calc_vacvessel_c1(vessel_inv[y], discount[y], vessel_name[y], reactor_util_guess[y], replace_dur[y]))

    # The vessel capex is input based on the net-electric power of a matched turbine + reactor
    # However, the turbine and reactor capacity are not fixed,
    # so we need to calculate the total cost based on the equivalent reactor capacity
    # Reactor gross-elec cap 
    # Units: $/year/MW-net * (# reactors) * (MW-net / reactor) = $/year
    @expression(EP, eVesselFixCosts[y in FUSION_VESSEL], eC1[y] * EP[:eFusionTurbNetCapSize][y] * EP[:vFusionNumReactors][y])

    for y in FUSION_VESSEL
        update_fixed_cost!(EP, eVesselFixCosts[y], y)
    end
end

function fusioninvestment!(EP::Model, inputs::Dict, setup::Dict)

    ##############################################################################
    
    dfFusion = inputs["dfFusion"]

    ##############################################################################

    fusioninvestment_fpp!(EP, inputs, setup)
    fusionthermalstorage_instantiation!(EP, inputs, setup)
    if any(dfFusion[!,:Add_Therm_Stor].>0)
        fusioninvestment_thermstor!(EP, inputs, setup)
    end
    fusioninvestment_fuel!(EP, inputs, setup)
    if any(dfFusion[!,:VV_Util_Deg].>0)
        fusioninvestment_vessel(EP, inputs, setup)
    end
end

function fusioncommit!(EP::Model, inputs::Dict, setup::Dict)
    # Adds variables and constraints describing the
    # commitment state of fusion reactors and turbines

    ##############################################################################

    dfFusion = inputs["dfFusion"]
    FUSION = inputs["FUSION"]

    T = inputs["T"]     # Number of time steps (hours)

    ##############################################################################

    ## FPP turbine commit state == general commitment state
    @expression(EP, eFusionTurbCommit[y in FUSION,t=1:T], EP[:vCOMMIT][y,t])

    ## FPP reactor commit state
    @variable(EP, 0 <= vFusionReactorCommit[y in FUSION,t=1:T])

    # Reactor commit state <= number of reactors
    @constraint(EP, [y in FUSION,t=1:T], vFusionReactorCommit[y,t] <= EP[:vFusionNumReactors][y])

    # Var_Turb_Cap controls whether the turbine size can be changed.
    # If it is not allowed for any generators, 
    # then constrain the turbine capacity to the standard capacity
    # and make turbine always work in lock-step with the reactor
    # TODO: We could allow them to stay separate in the future
    if !any(dfFusion[!,:Var_Turb_Cap] .> 0)
        @constraint(EP, [y in FUSION, t=1:T], vFusionReactorCommit[y,t] == eFusionTurbCommit[y,t])
    end

    ## Commitment state of the FPP
    if setup["UCommit"] == 1
        set_integer.(vFusionReactorCommit[y,:])
    end
    # This is automatically done for eFusionTurbCommit via vCommit
end

function fusionthermalstorage!(EP::Model, inputs::Dict, setup::Dict)

    ##############################################################################

    T = inputs["T"]     # Number of time steps (hours)

    dfFusion = inputs["dfFusion"]

    FUSION = inputs["FUSION"]
    FUSION_ThermStor = intersect(FUSION, dfFusion[dfFusion[!,:Add_Therm_Stor].>0,:R_ID])

    START_SUBPERIODS = inputs["START_SUBPERIODS"]
    INTERIOR_SUBPERIODS = inputs["INTERIOR_SUBPERIODS"]
    hours_per_subperiod = inputs["hours_per_subperiod"]

    ##############################################################################

    ## Energy stored <= Storage Capacity
    @constraint(EP, [y in FUSION,t=1:T], EP[:vThermStor][y,t] <= EP[:vThermStorCap][y])

    ## Limit charge and discharge rates
    @constraint(EP, [y in FUSION,t=1:T], EP[:vThermDis][y,t]  <= EP[:vThermDisCap][y])
    @constraint(EP, [y in FUSION,t=1:T], EP[:vThermChar][y,t] <= EP[:vThermChaCap][y])
    # @constraint(EP, [y in FUSION_ThermStor,t=1:T], EP[:vThermDis][y,t] + EP[:vThermChar][y,t] <= EP[:vThermDisCap][y])

    if "Therm_Stor_Leakage_Daily_Frac" in names(dfFusion)
        hourly_leakage = dfFusion[!,:Therm_Stor_Leakage_Daily_Frac] ./ 24.
    else
        hourly_leakage = zeros(length(FUSION_ThermStor))
    end

    ## Calculate the net discharge rate
    # This calculation is done from the TES point of view,
    # so the discharge energy is after storage -> discharge losses
    # and the charge energy is before charge -> storage losses
    @expression(EP, eThermStorNetDischarge_TES_PoV_thermal[y in FUSION, t=1:T], EP[:vThermDis][y,t] / EP[:eTES_dis_eff][y]  - EP[:vThermChar][y,t]  * EP[:eTES_cha_eff][y])  # Architecture independent
    @expression(EP, eThermStorNetDischarge_Turb_PoV_thermal[y in FUSION, t=1:T], EP[:vThermDis][y,t] - EP[:vThermChar][y,t])
    @expression(EP, eThermStorNetDischarge_Turb_PoV_elec[y in FUSION, t=1:T], (EP[:vThermDis][y,t] - EP[:vThermChar][y,t]) * EP[:eTES_turb_eff][y])                                                           # Architecture dependent: this is for direct
    @expression(EP, eThermStorLosses[y in FUSION, t=1:T], hourly_leakage[y] * EP[:vThermStor][y,t])


    ## Change in stored energy == net discharge rate
    ## Later, we can add losses or other factors which differentiate these two terms
    @constraint(EP, [y in FUSION_ThermStor,t in START_SUBPERIODS], 
        (1 - hourly_leakage[y]) * EP[:vThermStor][y,cyclicindex(t+hours_per_subperiod-1,T)] 
        - EP[:vThermStor][y,t] 
        == 
        EP[:vThermDis][y,t] / EP[:eTES_dis_eff][y] 
        - EP[:vThermChar][y,t] * EP[:eTES_cha_eff][y])

    @constraint(EP, [y in FUSION_ThermStor,t in INTERIOR_SUBPERIODS], 
        (1 - hourly_leakage[y]) * EP[:vThermStor][y,cyclicindex(t-1,T)] 
        - EP[:vThermStor][y,t] 
        == 
        EP[:vThermDis][y,t] / EP[:eTES_dis_eff][y] 
        - EP[:vThermChar][y,t] * EP[:eTES_cha_eff][y])
end

function fusionthermalbalance!(EP::Model, inputs::Dict, setup::Dict)
    # Set up the thermal energy balance for the fusion power plant.
    # This excludes operation of the thermal storage.
    # Thermal storage operation must be run beforehand

    ##############################################################################

    dfFusion = inputs["dfFusion"]
    FUSION = inputs["FUSION"]

    ##############################################################################

    T = inputs["T"]     # Number of time steps (hours)

    ## FPP fleet thermal output
    @variable(EP, vThermOutput[y in FUSION,t=1:T] >= 0)

    # Thermal power <= thermal capacity of non-standby units
    @constraint(EP, [y in FUSION,t=1:T], 
        vThermOutput[y,t] 
        <= 
        EP[:eFusionThermCapSize][y] * EP[:vFusionReactorCommit][y,t]
    )

    ### Calculation for the thermal balance of the salt loop
    saltEff = dfFusion[!,:Salt_Eff]      ## Salt electric heating efficiency (already converts MW to thermal)
    saltLosses = dfFusion[!,:Salt_Loss]  ## Hourly Thermal losses from salt loop (2 MW/hr) for each reactor

    ## Electric power used to heat the salt (MWe)
    # The thermal energy delivered to the salt = vSaltElecHeating * saltEff
    @variable(EP, 0. <= vSaltElecHeating[y in FUSION, t=1:T])
    @constraint(EP, [y in FUSION,t=1:T], vSaltElecHeating[y,t] <= EP[:vFusionNumReactors][y] * saltLosses[y] / saltEff[y])

    ## FPP thermal energy which reaches the plant steam generator > main turbine
    @expression(EP, eFPPThermalOutput[y in FUSION,t=1:T], EP[:vThermOutput][y,t] - EP[:vFusionNumReactors][y] * saltLosses[y] + vSaltElecHeating[y,t] * saltEff[y] -            # This is common to all architectures
                                                          (dfFusion[y, :Add_Therm_Stor] >= 1) * (dfFusion[y, :Therm_Stor_Cha_position] <= 1) * (EP[:vThermChar][y,t])           # If there is a TES charging from the plant/salt loop, it charges from the FPP output
                                                        )
    ## Steam loop thermal energy reaching the main turbine
    @expression(EP, eTurbThermal[y in FUSION,t=1:T], (dfFusion[y, :Add_Therm_Stor] != 1) * (EP[:eFPPThermalOutput][y,t] * dfFusion[y, :HEX_Plant_Steam_Eff]) -                  # No TES/Indirect TES: the FPP Thermal output is connected to the turbine
                                                     (dfFusion[y, :Add_Therm_Stor] >= 1) * (dfFusion[y, :Therm_Stor_Cha_position] == 2) * (EP[:vThermChar][y,t]) +              # If there is a TES charging from the steam loop, it charges from the turbine input
                                                     (dfFusion[y, :Add_Therm_Stor] == 1) * (dfFusion[y,:Therm_Stor_Dis_position] == 0) * (EP[:vThermDis][y,t]) +                # Direct TES discharging via the main turbine
                                                     (dfFusion[y, :Add_Therm_Stor] == 2) * (dfFusion[y,:Therm_Stor_Dis_position] == 0) * (EP[:vThermDis][y,t])                  # Indirect TES discharging via the main turbine
                                                     )

    ## TESWF loop thermal energy which reaches the TES turbine
    @expression(EP, eTESTurbThermal[y in FUSION,t=1:T], (dfFusion[y, :Add_Therm_Stor] == 1) * (dfFusion[y,:Therm_Stor_Dis_position] == 1) * (EP[:vThermDis][y,t]) +             # Direct TES case with a new turbine
                                                        (dfFusion[y, :Add_Therm_Stor] == 2) * (dfFusion[y,:Therm_Stor_Dis_position] == 2) * (EP[:vThermDis][y,t])               # Indirect TES case with a dedicated TES turbine
                                                        )

    @constraint(EP, [y in FUSION, t=1:T], 0.0 <= eFPPThermalOutput[y,t])
    @constraint(EP, [y in FUSION, t=1:T], 0.0 <= eTurbThermal[y,t])
    @constraint(EP, [y in FUSION, t=1:T], 0.0 <= eTESTurbThermal[y,t])
end

### This function will calculate the power provided by the reactor to the grid
function fusionelecbalance!(EP::Model, inputs::Dict, setup::Dict)

    ##############################################################################

    MMBTU_PER_MWH = 3.412 ### hard coded MMBTU to MWh conversion

    T = inputs["T"]     # Number of time steps (hours)
    Z = inputs["Z"]     # Number of zones

    dfGen = inputs["dfGen"]
    dfFusion = inputs["dfFusion"]

    FUSION = inputs["FUSION"]

    ## Define key variables
    magcool = dfFusion[!,:Mag_Cool]  #Amount of power (MW) being delivered to cool the magnets

    ##############################################################################

    @expression(EP, ePlantFix[y in FUSION,t=1:T], dfFusion[y,:Plant_fixed_power_MWe] * EP[:vFusionReactorCommit][y,t])  # Fixed power being delivered to the power plant (Number multiplied with the binary variable for start)
    @expression(EP, ePlantVar[y in FUSION,t=1:T], dfFusion[y,:Plant_var_power_MWeperMWt] * EP[:vThermOutput][y,t])  # Variable power being delivered to the power plant (Function of the ThermalOutput of the power plant)

    ## Combine to make recirculating power expression
    @expression(EP, eRecircPwr[y in FUSION,t=1:T], 
        + EP[:vFusionNumReactors][y] * magcool[y] 
        + ePlantFix[y,t] 
        + EP[:vSaltElecHeating][y,t] 
        + ePlantVar[y,t]
    )
    @constraint(EP, [y in FUSION, t=1:T], 0.0 <= eRecircPwr[y,t])

    # The amount of electric energy exiting the main turbine
    @expression(EP, eTurbGrossElec[y in FUSION,t=1:T], EP[:eTurbThermal][y,t] * EP[:eFPP_turb_eff][y])
    @expression(EP, eTESTurbGrossElec[y in FUSION,t=1:T], EP[:eTESTurbThermal][y,t] * EP[:eTES_turb_eff][y])

    # Note, we have an existing constraint on the grid power:
    # cMax_vP_therm_commit:
        # vP[y,t] <= availability[y,t] * Cap_Size[y] * vCOMMIT[y,t]
    # This is set in thermal_commit.jl

    ## Expression for the amount of net power that is actually going to the grid
    @expression(EP, eFusionNetElec[y in FUSION,t=1:T], eTurbGrossElec[y,t] + eTESTurbGrossElec[y,t] - eRecircPwr[y,t])
    # @constraint(EP, [y in FUSION,t=1:T], 0.0 <= eFusionNetElec[y,t])

    ## Variable for the amount of power that is being imported by the reactor when it is in standby mode
    @variable(EP, vFusionElecImports[y in FUSION,t=1:T] >= 0)
    @constraint(EP, [y in FUSION,t=1:T], vFusionElecImports[y,t] <= eRecircPwr[y,t])

    ## Expression to sum the fusion imports across the zones 
    @expression(EP, eFusionElecImports_byZone[t=1:T, z=1:Z], 
        sum(vFusionElecImports[y,t] for y in intersect(FUSION, dfGen[dfGen[!,:Zone].==z,:R_ID])))
    for z in 1:Z
        for t in 1:T
            add_to_expression!(EP[:ePowerBalance][t,z], -eFusionElecImports_byZone[t,z])
        end
    end

    ## Constraints for the electric power node (Supply of power on the left, Demand for power on right)
    @constraint(EP, PwrNode[y in FUSION, t=1:T], 
        + eTurbGrossElec[y,t]
        + eTESTurbGrossElec[y,t] 
        + vFusionElecImports[y,t] 
        == 
        + EP[:vP][y,t] 
        + EP[:eRecircPwr][y,t]
    )
end

function fusionfuel!(EP::Model, inputs::Dict, setup::Dict)

    ##############################################################################

    T = inputs["T"]     # Number of time steps (hours)

    FUSION = inputs["FUSION"]
    dfFusion = inputs["dfFusion"]

    START_SUBPERIODS = inputs["START_SUBPERIODS"]
    INTERIOR_SUBPERIODS = inputs["INTERIOR_SUBPERIODS"]
    hours_per_subperiod = inputs["hours_per_subperiod"]

    ##############################################################################

    ##### Tritium Balance #####

    ## Tritium inventory
    @variable(EP, vTritInventory[y in FUSION, t=1:T], lower_bound = 0.)
    ## Tritium inventory <= Tritium capacity
    @constraint(EP, [y in FUSION,t=1:T], vTritInventory[y,t] <= EP[:vTritCap][y])

    # Tritium exports
    @variable(EP, vTritExports[y in FUSION, t=1:T] >= 0.0)

    ## Fuel useage
    # eV produced per fusion = 14.1e6
    # Joules per fusion = 14.1e6 * 1.6e-19
    # Fusions per MWh = 3600 * 1e6 / (14.1e6 * 1.6e-19)
    # Avagadro's number = 6.022e23
    # Mols per MWh = (3600 * 1e6 / (14.1e6 * 1.6e-19)) / 6.022e23
    # Tritium molar mass = 3g / mol
    # kg tritium per MWh = 3 / 1000 * (3600 * 1e6 / (14.1e6 * 1.6e-19)) / 6.022e23
    # Extra thermal energy from secondary reactions in the blanket ~ 15%
    # Adjusted kg tritium per MWh = 3 / 1000 * (3600 * 1e6 / (14.1e6 * 1.6e-19)) / 6.022e23 / 1.15 = 6.91267e-6
    # We require the same number of mols of deuterium, but that equals 2/3 the mass

    TRIT_FUEL_RATIO = 3. / 1000. * (3600. * 1e6 / (14.1e6 * 1.6e-19)) / 6.022e23 ./ dfFusion[!,:Blanket_Power_Mult]
    DEU_FUEL_RATIO = TRIT_FUEL_RATIO .* (2. / 3.)

    ## Hourly tritium decay rate [/ hr]
    TRIT_DECAY = 0.00000642259

    @expression(EP, eTritBreeding[y in FUSION, t=1:T],
        TRIT_FUEL_RATIO[y] * EP[:vThermOutput][y,t] * dfFusion[y,:Trit_Breed]
    )

    @expression(EP, eTritConsumption[y in FUSION, t=1:T],
        TRIT_FUEL_RATIO[y] * EP[:vThermOutput][y,t] * dfFusion[y,:Trit_Fuel]
    )

    @expression(EP, eTritDecay[y in FUSION, t=1:T],
        vTritInventory[y, t] * TRIT_DECAY
    )

    @expression(EP, eTritLeakage[y in FUSION, t=1:T],
        vTritInventory[y, t] * dfFusion[y,:Trit_Loss]
    )

    ## Hourly tritium balance
    @constraint(EP, cTrit_balance_start[y in FUSION, t in START_SUBPERIODS], 
        vTritInventory[y,t] 
        == 
        + vTritInventory[y, cyclicindex(t+hours_per_subperiod-1, T)]
        + eTritBreeding[y,t]
        - eTritConsumption[y,t]
        - eTritDecay[y,cyclicindex(t+hours_per_subperiod-1, T)]
        - eTritLeakage[y,cyclicindex(t+hours_per_subperiod-1, T)]
        - vTritExports[y,t]
    )

    @constraint(EP, cTrit_balance_interior[y in FUSION, t in INTERIOR_SUBPERIODS], 
        vTritInventory[y,t] 
        == 
        + vTritInventory[y, cyclicindex(t-1, T)]
        + eTritBreeding[y,t]
        - eTritConsumption[y,t]
        - eTritDecay[y,cyclicindex(t-1, T)]
        - eTritLeakage[y,cyclicindex(t-1, T)]
        - vTritExports[y,t]
    )

    @constraint(EP, [y in FUSION,t in START_SUBPERIODS], 
        vTritInventory[y,cyclicindex(t, T)] 
        >= 
        eTritConsumption[y,cyclicindex(t+hours_per_subperiod-1, T)]
    )

    @constraint(EP, [y in FUSION,t in INTERIOR_SUBPERIODS], 
        vTritInventory[y,cyclicindex(t, T)] 
        >= 
        eTritConsumption[y,cyclicindex(t-1, T)]
    )

    ##### Deuterium Balance #####

    # Deuterium inventory
    @variable(EP, vDeuInventory[y in FUSION, t=1:T], lower_bound = 0.)
    @constraint(EP, [y in FUSION,t=1:T], vDeuInventory[y,t] <= EP[:vDeuCap][y])

    # Deuterium exports
    @variable(EP, 0. <= vDeuImports[y in FUSION, t=1:T])

    # Deuterium Consumption
    @expression(EP, eDeuConsumption[y in FUSION, t=1:T],
        DEU_FUEL_RATIO[y] * EP[:vThermOutput][y,t] * dfFusion[y,:Deu_Fuel]
    )

    # Deuterium Leakage
    @expression(EP, eDeuLeakage[y in FUSION, t=1:T],
        vDeuInventory[y, t] * dfFusion[y,:Deu_Loss]
    )

    # Deuterium balance
    @constraint(EP, cDeu_balance_start[y in FUSION, t in START_SUBPERIODS], 
        vDeuInventory[y,t] 
        == 
        + vDeuInventory[y, cyclicindex(t+hours_per_subperiod-1, T)]
        + vDeuImports[y,t]
        - eDeuConsumption[y,t]
        - eDeuLeakage[y,cyclicindex(t+hours_per_subperiod-1, T)]
    )

    @constraint(EP, cDeu_balance_interior[y in FUSION, t in INTERIOR_SUBPERIODS], 
        vDeuInventory[y,t] 
        == 
        + vDeuInventory[y, cyclicindex(t-1, T)]
        + vDeuImports[y,t]
        - eDeuConsumption[y,t]
        - eDeuLeakage[y,cyclicindex(t-1, T)]
    )

    @constraint(EP, [y in FUSION,t in START_SUBPERIODS], 
        vDeuInventory[y,cyclicindex(t, T)] 
        >= 
        DEU_FUEL_RATIO[y] * EP[:vThermOutput][y,cyclicindex(t+hours_per_subperiod-1, T)] * dfFusion[y,:Deu_Fuel]
    )

    @constraint(EP, [y in FUSION,t in INTERIOR_SUBPERIODS], 
        vDeuInventory[y,cyclicindex(t, T)] 
        >= 
        DEU_FUEL_RATIO[y] * EP[:vThermOutput][y,cyclicindex(t-1, T)] * dfFusion[y,:Deu_Fuel]
    )

    # Adding costs for deuterium imports
    @expression(EP, deuterium_var_costs[y in FUSION, t=1:T],
        vDeuImports[y,t] .* dfFusion[y,:Deu_Import_Cost]
    )

    for t in 1:T
        for y in FUSION
            update_var_cost_t!(EP, deuterium_var_costs[y, t], y, t)
        end
    end
end

# This function incorporates the vessel costs
function fusionvessel!(EP::Model, inputs::Dict, setup::Dict)

    ##############################################################################

    T = inputs["T"]     # Number of time steps (hours)

    dfGen = inputs["dfGen"]
    dfFusion = inputs["dfFusion"]

    # Grab variables from the fusion_data.csv
    vessel_name = dfFusion[!,:Vessel_Life]
    replace_dur = dfFusion[!,:Rep_Dur]
    discount = dfFusion[!,:Dis_Rate]
    vessel_inv = dfFusion[!,:Inv_Vessel_per_MWe]
    reactor_util_guess = dfFusion[!,:Reactor_Utilization_Guess]

    MMBTU_PER_MWH = 3.412

    FUSION = inputs["FUSION"]
    FUSION_VESSEL = intersect(FUSION, dfFusion[dfFusion[!,:VV_Util_Deg].>0,:R_ID])

    ##############################################################################

    ## Calculate annual plant utilization based on the thermal output
    @expression(EP, eThermOutputTot[y in FUSION_VESSEL], sum(EP[:vThermOutput][y,t] for t=1:T))

    ## Calculate Vessel Investment Costs
    @expression(EP, eC2[y in FUSION_VESSEL], calc_vacvessel_c2(vessel_inv[y], discount[y], vessel_name[y], reactor_util_guess[y], replace_dur[y]))

    ## Constrain annual utilization
    if "Cap_Fac_Constraint_Type" in names(dfFusion)
        # FUSION_100 = intersect(FUSION_VESSEL, dfFusion[dfFusion[!,:Cap_Fac_Constraint_Type].==2,:R_ID])
        FUSION_VESSEL_any = intersect(FUSION_VESSEL, dfFusion[dfFusion[!,:Cap_Fac_Constraint_Type].==1,:R_ID])
    else 
        FUSION_VESSEL_any = FUSION_VESSEL
    end
    @constraint(EP, cAnnualUtilMax[y in FUSION_VESSEL_any], 
        eThermOutputTot[y] 
        <= 
        EP[:eFusionThermCap][y] * T * calc_fpp_maxutil(vessel_name[y], replace_dur[y])
    )

    ## Constrain minimum annual utilization to avoid negative vessel costs
    @constraint(EP, cAnnualUtilMin[y in FUSION_VESSEL], 
        eThermOutputTot[y] 
        >= 
        EP[:eFusionThermCap][y] * T * calc_fpp_minutil(EP[:eC1][y], eC2[y])
    )
    
    ## Add Vessel Investment Costs to Var Costs
    @expression(EP, 
        eVesselVarCosts[y in FUSION_VESSEL], 
        eC2[y] / T * eThermOutputTot[y] * EP[:eFusionTurbNetCapSize][y] / EP[:eFusionThermCapSize][y]
    )

    for y in FUSION_VESSEL
        update_var_cost!(EP, eVesselVarCosts[y], y)
    end
end

function cyclicindex(idx::Int, maxidx::Int)
    return mod(idx - 1, maxidx) + 1
end

function fusionpulseratio(pulse_length_min::Float64, pulse_gap_min::Float64)
    return pulse_length_min / (pulse_length_min + pulse_gap_min)
end

function calc_fpp_maxutil(nom_lifetime::Float64, replace_dur::Float64)
    # return (-nom_lifetime + sqrt(nom_lifetime^2 + 4 * nom_lifetime * replace_dur)) / (2 * replace_dur)
    return(1.0 - replace_dur / nom_lifetime)
end

function vessel_degradation_taylor_first(discount_rate::Float64, nom_lifetime::Float64, util_guess::Float64, rep_dur::Float64)
    vessel_lifetime = nom_lifetime / util_guess + rep_dur
    return (
        nom_lifetime * log(1 + discount_rate) * (1 + discount_rate)^vessel_lifetime
        / 
        ((1 + discount_rate)^vessel_lifetime - 1)^2 
        / util_guess^2
    )
end

function calc_vacvessel_c1(capex::Float64, discount_rate::Float64, nom_lifetime::Float64, util_guess::Float64, rep_dur::Float64)
    return capex * discount_rate *
        (
            (1 / (1 - (1 + discount_rate)^-(nom_lifetime / util_guess + rep_dur)))  
            - 
            util_guess * vessel_degradation_taylor_first(discount_rate, nom_lifetime, util_guess, rep_dur)
        )
end

function calc_vacvessel_c2(capex::Float64, discount_rate::Float64, nom_lifetime::Float64, util_guess::Float64, rep_dur::Float64)
    return capex * discount_rate * vessel_degradation_taylor_first(discount_rate, nom_lifetime, util_guess, rep_dur)
end

function update_fixed_cost!(EP::Model, term::Union{AffExpr, Float64}, R_ID::Int64)
    add_to_expression!(EP[:eCFix][R_ID], term)
    add_to_expression!(EP[:eTotalCFix], term)
    add_to_expression!(EP[:eObj], term)
end

function update_var_cost!(EP::Model, term::Union{AffExpr, Float64}, R_ID::Int64)
    # add_to_expression!.(EP[:eCVar_out][R_ID,:], term)
    # add_to_expression!.(EP[:eTotalCVarOutT], term)
    add_to_expression!(EP[:eTotalCVarOut], term)
    add_to_expression!(EP[:eObj], term)
end

function update_var_cost_t!(EP::Model, term::Union{AffExpr, Float64}, R_ID::Int64, t::Int64)
    add_to_expression!(EP[:eCVar_out][R_ID,t], term)
    # add_to_expression!.(EP[:eTotalCVarOutT], term)
    add_to_expression!(EP[:eTotalCVarOut], term)
    add_to_expression!(EP[:eObj], term)
end

function calc_fpp_minutil(c1::Float64, c2::Float64)
    # The periodic vacuum vessel costs are of the form y = c1 + c2 * x
    # Where x is the utilization factor and y is the cost
    # If y = 0, x = -c1 / c2
    return -c1 / c2
end


## This function is the overall function for fusion power 
function fusion!(EP::Model, inputs::Dict, setup::Dict)
    println("Fusion Power Module")

    print("-- Fusion investments ... ")
    fusioninvestment!(EP, inputs, setup)
    print("done")

    print("\n-- Fusion commitment states ... ")
    fusioncommit!(EP, inputs, setup)
    print("done")

    dfFusion = inputs["dfFusion"]
    if any(dfFusion[!,:Add_Therm_Stor].>0)
        print("\n-- Fusion thermal storage ... ")
        fusionthermalstorage!(EP,inputs,setup)
        print("done")
    end

    print("\n-- Fusion thermal balance ... ")
    fusionthermalbalance!(EP, inputs, setup)
    print("done")

    print("\n-- Fusion electric balance ... ")
    fusionelecbalance!(EP, inputs, setup)
    print("done")

    print("\n-- Fusion fuel ... ")
    fusionfuel!(EP, inputs, setup)
    print("done")

    if any(dfFusion[!,:VV_Util_Deg].>0)
        print("\n-- Fusion vessel costs ... ")
        fusionvessel!(EP,inputs,setup)
        print("done\n")
    end
end     