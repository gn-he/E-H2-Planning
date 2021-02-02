# ENV["CPLEX_STUDIO_BINARIES"] = "/home/gridsan/gnhe/opt/ibm/ILOG/CPLEX_Studio129/cplex/bin/x86-64_linux"

# export PATH=/opt/gurobi801/linux64/bin:$PATH
# export LD_LIBRARY_PATH=/opt/gurobi801/linux64/lib:$LD_LIBRARY_PATH

module Hydrogen_Modular
#export Module functions
export load_H2_inputs, adjust_parameters, print_setting, generate_master_model, generate_H2_model, generate_Power_H2_model,generate_Power_model, fix_planning_var,warm_start,solve_model, generate_Power_results, generate_H2_results,save_H2_results,solve_H2_model, solve_Power_H2_model, H2_Gen, H2_PowerGen,H2_Demand, H2_Storage, H2Pipeline, H2Truck,H2Truck_old,Battery_Truck

cd(dirname(@__FILE__))

working_path = pwd()

# GenX path
genx_path = pwd()

# Load GenX modules
push!(LOAD_PATH, genx_path)
println("Loading hydrogen packages")

using PowerCapExp
using JuMP
using DataFrames #This package allows put together data into a matrix
# using Gurobi #Gurobi solver
# using CPLEX #CPLEX solver
using MathProgBase #for fix_integers
using CSV
using StatsBase
using LinearAlgebra
using Random
using JLD2, FileIO
using Dates

using Statistics
using JSON
using DelimitedFiles

cd(dirname(@__FILE__))



function Battery_Truck(HY::Model, dModuleArgs::Dict)

	println("Battery Truck Module")

	inputs_H2 = dModuleArgs["inputs_H2"]
	setup = dModuleArgs["setup"]
	dExpressions_H2 = dModuleArgs["dExpressions_H2"]
	dObjective_H2 = dModuleArgs["dObjective_H2"]
	dBalance_H2 = dModuleArgs["dBalance_H2"]

	Z_total = inputs_H2["Z_total"]
    PT_total = inputs_H2["PT_total"]
    RT_total = inputs_H2["RT_total"]
    K_prod_total = inputs_H2["K_prod_total"]
    K_stor_total = inputs_H2["K_stor_total"]
    G_H2_total = inputs_H2["G_H2_total"]
    T = inputs_H2["T"]
    Z = inputs_H2["Z"]

    BT = 1

    W = inputs_H2["W"]
    Tw = inputs_H2["Tw"]

	weights = inputs_H2["weights"]
	RouteLength = inputs_H2["RouteLength"]
	H2TruckUnitCapex = inputs_H2["H2TruckUnitCapex"]
    TruckCap = inputs_H2["TruckCap"]
	Full_weight = inputs_H2["Full_weight_tonne_per_unit"]
	Empty_weight = inputs_H2["Empty_weight_tonne_per_unit"]
    AvgTruckSpeed = inputs_H2["AvgTruckSpeed"]
    H2TruckUnitOpex_full = inputs_H2["H2TruckUnitOpex_full"]
	H2TruckUnitOpex_empty = inputs_H2["H2TruckUnitOpex_empty"]
    TD = inputs_H2["TD"]
    H2TLoss = inputs_H2["H2TLoss"]
    truck_emission_rate = inputs_H2["truck_emission_rate"]
	H2TruckCompressionUnitCapex = inputs_H2["H2TruckCompressionUnitCapex"]
	H2TruckCompressionEnergy = inputs_H2["H2TruckCompressionEnergy"]
	H2TruckCompressionUnitOpex = inputs_H2["H2TruckCompressionUnitOpex"]

	Truck_integer_model = inputs_H2["Truck_integer_model"]
	truck_model_simp = inputs_H2["truck_model_simp"]

	Transmission_cost_factor = inputs_H2["Transmission_cost_factor"]


	H2TruckUnitCapex = [735000]
	H2TruckUnitOpex = [20/AvgTruckSpeed]
	BatTruckLoss = [0.05]
	TruckCap = [2.7]
	### Variables ###

	@variable(HY, vBatTruckFlow[z=1:Z, r=1:BT, t = 1:T] )
	# @variable(HY, vH2TruckFlow_pos[zz=1:Z,z=1:Z, r=1:BT, t = 1:T] >=0 )
	# @variable(HY, vH2TruckFlow_neg[zz=1:Z,z=1:Z, r=1:BT, t = 1:T] >=0 )
	# @variable(HY, vH2TruckLevel[zz=1:Z,z=1:Z, r=1:BT, t = 1:T] >=0 )


	if inputs_H2["Truck_option_avai"] == 0
		@variable(HY, 0 >= vN_Bat_Truck[r=1:BT] >= 0 )
        @variable(HY, 0 >= vN_Bat_TruckRoute[r=1:BT] >= 0 )
        @variable(HY, 0 >= vN_Bat_avail_full[z=1:Z,r=1:BT,t = 1:T] >= 0  )
        @variable(HY, 0 >= vN_Bat_travel_full[zz=1:Z,z=1:Z,r=1:BT,t = 1:T] >= 0 )
        @variable(HY, 0 >= vN_Bat_arrive_full[zz=1:Z,z=1:Z,r=1:BT,t = 1:T] >= 0 )
        @variable(HY, 0 >= vN_Bat_depart_empty[zz=1:Z,z=1:Z,r=1:BT,t = 1:T] >= 0 )
        @variable(HY, 0 >= vN_Bat_avail_empty[z=1:Z,r=1:BT,t = 1:T] >= 0  )
        @variable(HY, 0 >= vN_Bat_travel_empty[zz=1:Z,z=1:Z,r=1:BT,t = 1:T] >= 0 )
        @variable(HY, 0 >= vN_Bat_arrive_empty[zz=1:Z,z=1:Z,r=1:BT,t = 1:T] >= 0 )
        @variable(HY, 0 >= vN_Bat_depart_full[zz=1:Z,z=1:Z,r=1:BT,t = 1:T] >= 0 )
        @variable(HY, 0 >= vN_Bat_charged[z=1:Z,r=1:BT,t = 1:T] >= 0 )
        @variable(HY, 0 >= vN_Bat_discharged[z=1:Z,r=1:BT,t = 1:T] >= 0 )
        @variable(HY, 0 >= vN_Bat__full[r=1:BT,t = 1:T] >= 0 )
        @variable(HY, 0 >= vN_Bat__empty[r=1:BT,t = 1:T] >= 0 )
	else
		@variable(HY, vN_Bat_Truck[r=1:BT] >= 0 )
        @variable(HY, vN_Bat_TruckRoute[r=1:BT] >= 0 )
        @variable(HY, vN_Bat_avail_full[z=1:Z,r=1:BT,t = 1:T] >= 0  )
        @variable(HY, vN_Bat_travel_full[zz=1:Z,z=1:Z,r=1:BT,t = 1:T] >= 0 )
        @variable(HY, vN_Bat_arrive_full[zz=1:Z,z=1:Z,r=1:BT,t = 1:T] >= 0 )
        @variable(HY, vN_Bat_depart_full[zz=1:Z,z=1:Z,r=1:BT,t = 1:T] >= 0 )
        @variable(HY, vN_Bat_avail_empty[z=1:Z,r=1:BT,t = 1:T] >= 0  )
        @variable(HY, vN_Bat_travel_empty[zz=1:Z,z=1:Z,r=1:BT,t = 1:T] >= 0 )
        @variable(HY, vN_Bat_arrive_empty[zz=1:Z,z=1:Z,r=1:BT,t = 1:T] >= 0 )
        @variable(HY, vN_Bat_depart_empty[zz=1:Z,z=1:Z,r=1:BT,t = 1:T] >= 0 )
        @variable(HY, vN_Bat_charged[z=1:Z,r=1:BT,t = 1:T] >= 0 )
        @variable(HY, vN_Bat_discharged[z=1:Z,r=1:BT,t = 1:T] >= 0 )
        @variable(HY, vN_Bat__full[r=1:BT,t = 1:T] >= 0 )
        @variable(HY, vN_Bat__empty[r=1:BT,t = 1:T] >= 0 )
	end


	### Expressions ###
	## Objective Function Expressions ##
	discount_factor = dExpressions_H2["discount_factor"]

	@expression(HY, CAPEX_Bat_Truck, sum(inputs_H2["discount_factor_H2Truck"][r] * vN_Bat_Truck[r] * H2TruckUnitCapex[r] for r = 1:BT))
	dExpressions_H2["CAPEX_Bat_Truck"] = CAPEX_Bat_Truck
	dObjective_H2["CAPEX_Bat_Truck"] = CAPEX_Bat_Truck

    @expression(HY, OPEX_Bat_Truck, sum(weights[t] *(vN_Bat_arrive_full[zz,z,r,t] + vN_Bat_arrive_empty[zz,z,r,t])*H2TruckUnitOpex[r] * RouteLength[zz,z]  for zz = 1:Z,z=1:Z,r = 1:BT,t = 1:T if zz != z))

	dExpressions_H2["OPEX_Bat_Truck"] = OPEX_Bat_Truck
	dObjective_H2["OPEX_Bat_Truck"] = OPEX_Bat_Truck

	# @expression(HY, CAPEX_Compression_Truck, Transmission_cost_factor * inputs_H2["discount_factor_H2Compression"] *sum(vH2TruckCompressionCap[z,r] * H2TruckCompressionUnitCapex[r] for z=1:Z,r = 1:BT))
	# dExpressions_H2["CAPEX_Compression_Truck"] = CAPEX_Compression_Truck
	# dObjective_H2["CAPEX_Compression_Truck"] = CAPEX_Compression_Truck

	## End Objective Function Expressions ##

	## Balance Expressions ##
	# H2 Power Consumption balance
	# @expression(HY, eH2TruckCompressionPowerConsumption[t=1:T, z=1:Z],
	# sum(vN_Bat_charged[z,r,t] * TruckCap[r]* H2TruckCompressionEnergy[r] for r = 1:BT)
	# )
	# dExpressions_H2["eH2TruckCompressionPowerConsumption"] = eH2TruckCompressionPowerConsumption
	# dBalance_H2["H2PowerBalance"] = dBalance_H2["H2PowerBalance"] + eH2TruckCompressionPowerConsumption

	# H2 balance
	@expression(HY, BatTruckFlow[t=1:T,z=1:Z], sum(vBatTruckFlow[z,r,t] for r = 1:BT ))
	dExpressions_H2["BatTruckFlow"] = BatTruckFlow
	dBalance_H2["H2PowerBalance"] = dBalance_H2["H2PowerBalance"] - BatTruckFlow
	# dBalance_H2["H2Balance"] = dBalance_H2["H2Balance"] + TruckFlow

	# Carbon emission balance
	# @expression(HY, Truck_carbon_emission, sum(weights[t] * (vN_Bat_arrive_full[z,zz,r,t]*Full_weight[r]+vN_Bat_arrive_empty[z,zz,r,t]*Empty_weight[r]) *RouteLength[z,zz] *truck_emission_rate[r,z] for zz = 1:Z, z = 1:Z, r = 1:BT, t = 1:T if zz != z))
	# dExpressions_H2["Truck_carbon_emission"] = Truck_carbon_emission
	# dBalance_H2["H2_Carbon_Emission"] = dBalance_H2["H2_Carbon_Emission"] + Truck_carbon_emission
	## End Balance Expressions ##
	### End Expressions ###

	### Constraints ###
	for w in 1:W
    	for h in 1:Tw
	        t = Tw*(w-1)+h
	        tw_min = Tw*(w-1)+1
	        tw_max = Tw*(w-1)+Tw
			for z in 1:Z
			    # for t in 1:T
		        for zz in 1:Z
		            for r in 1:BT
						# @constraint(HY, vN_Bat_travel[zz,z,r,t] + vN_Bat_travel[z,zz,r,t] + vN_Bat_avail[zz,z,r,t] + vN_Bat_avail[z,zz,r,t] == vN_Bat_TruckRoute[z,zz,r])


					    ## Total number
					    @constraint(HY,
					    vN_Bat__full[r,t] +
					    vN_Bat__empty[r,t] +
					    0
					    ==
					    vN_Bat_Truck[r]
					    )
					    # The number of total full trucks
					    @constraint(HY, vN_Bat__full[r,t] ==
					    sum(vN_Bat_travel_full[zz,z,r,t] for zz = 1:Z, z = 1:Z if zz != z) +
					    sum(vN_Bat_avail_full[z,r,t] for z = 1:Z) +
					    0
					    )
					    ## The number of total empty trucks
						@constraint(HY, vN_Bat__empty[r,t] ==
					   sum(vN_Bat_travel_empty[zz,z,r,t] for zz = 1:Z, z = 1:Z if zz != z) +
					   sum(vN_Bat_avail_empty[z,r,t] for z = 1:Z) +
					   0
					   )


					    # @constraint(HY, vN_Bat_discharged[z,r,t] <= vN_Bat_avail_full[z,r,t] )
					    # @constraint(HY, vN_Bat_charged[zz,z,r,t] <= vN_Bat_TruckRoute[zz,z,r] )

					    # @constraint(HY, vN_Bat_TruckRoute[zz,z,r] == vN_Bat_TruckRoute[z,zz,r])
						@constraint(HY, vN_Bat_TruckRoute[r] == vN_Bat_Truck[r] )
					    if z==zz
						    # @constraint(HY, vN_Bat_TruckRoute[zz,z,r] == 0 )
						    @constraint(HY, vN_Bat_travel_full[zz,z,r,t] == 0 )
						    @constraint(HY, vN_Bat_travel_empty[zz,z,r,t] == 0 )
							@constraint(HY, vN_Bat_arrive_full[zz,z,r,t] == 0 )
						    @constraint(HY, vN_Bat_depart_full[zz,z,r,t] == 0 )
						    @constraint(HY, vN_Bat_arrive_empty[zz,z,r,t] == 0 )
						    @constraint(HY, vN_Bat_depart_empty[zz,z,r,t] == 0 )
							# @constraint(HY, vH2TruckFlow[zz,z,r,t] == 0 )
					    else
					    	# @constraint(HY, vN_Bat_TruckRoute[zz,z,r] >=0 )
					    # @constraint(HY, 5000>= vN_Bat_TruckRoute[zz,z,r] )
					    end

					    t_arrive = 1
					    t_depart = 1

					    # Change of the number of full available trucks
					    if h > 1
						    @constraint(HY, vN_Bat_avail_full[z,r,t] - vN_Bat_avail_full[z,r,t-1] ==
						    vN_Bat_charged[z,r,t] - vN_Bat_discharged[z,r,t] +
							sum(vN_Bat_arrive_full[zz,z,r,t-t_arrive] for zz = 1:Z if zz != z) -
							sum(vN_Bat_depart_full[z,zz,r,t-t_depart] for zz = 1:Z if zz != z) +
						    0
						    )
					    else
					        @constraint(HY, vN_Bat_avail_full[z,r,t] - vN_Bat_avail_full[z,r,tw_max] ==
					        vN_Bat_charged[z,r,t] - vN_Bat_discharged[z,r,t] +
					        sum(vN_Bat_arrive_full[zz,z,r,tw_max] for zz = 1:Z if zz != z) -
					        sum(vN_Bat_depart_full[z,zz,r,tw_max] for zz = 1:Z if zz != z) +
						    0
							)
					    end
					    # Change of the number of empty available trucks
					    if h > 1
						    @constraint(HY, vN_Bat_avail_empty[z,r,t] - vN_Bat_avail_empty[z,r,t-1] ==
						    -vN_Bat_charged[z,r,t] + vN_Bat_discharged[z,r,t] +
						    sum(vN_Bat_arrive_empty[zz,z,r,t-t_arrive] for zz = 1:Z if zz != z) -
						    sum(vN_Bat_depart_empty[z,zz,r,t-t_depart] for zz = 1:Z if zz != z)
							)
					    else
					        @constraint(HY, vN_Bat_avail_empty[z,r,t] - vN_Bat_avail_empty[z,r,tw_max] ==
					        -vN_Bat_charged[z,r,t] + vN_Bat_discharged[z,r,t] +
					        sum(vN_Bat_arrive_empty[zz,z,r,tw_max] for zz = 1:Z if zz != z) -
					        sum(vN_Bat_depart_empty[z,zz,r,tw_max] for zz = 1:Z if zz != z)
							)
					    end

					    # Change of the number of full traveling trucks
					    if h > 1
					    @constraint(HY, vN_Bat_travel_full[z,zz,r,t] - vN_Bat_travel_full[z,zz,r,t-1] ==
					    vN_Bat_depart_full[z,zz,r,t-t_depart] -
					    vN_Bat_arrive_full[z,zz,r,t-t_arrive]
					    )
					    else
					        @constraint(HY, vN_Bat_travel_full[z,zz,r,t] - vN_Bat_travel_full[z,zz,r,tw_max] ==
					        vN_Bat_depart_full[z,zz,r,tw_max]
					        - vN_Bat_arrive_full[z,zz,r,tw_max])
					    end

					    # Change of the number of empty traveling trucks
					    if h > 1
						    @constraint(HY, vN_Bat_travel_empty[z,zz,r,t] - vN_Bat_travel_empty[z,zz,r,t-1] ==
						    vN_Bat_depart_empty[z,zz,r,t-t_depart] -
						    vN_Bat_arrive_empty[z,zz,r,t-t_arrive])
					    else
					        @constraint(HY, vN_Bat_travel_empty[z,zz,r,t] - vN_Bat_travel_empty[z,zz,r,tw_max] ==
					        vN_Bat_depart_empty[z,zz,r,tw_max] -
					        vN_Bat_arrive_empty[z,zz,r,tw_max])
					    end


					    # Travel delay

					    if t-TD[zz,z]+1 >= tw_min && t-TD[zz,z]+1 <= tw_max && t-TD[zz,z]+1 <= t
					    @constraint(HY, vN_Bat_travel_full[zz,z,r,t] >= sum(vN_Bat_depart_full[zz,z,r,tt] for tt = (t-TD[zz,z]+1):t) ) # deaprt from zz to z
					    @constraint(HY, vN_Bat_travel_empty[zz,z,r,t] >= sum(vN_Bat_depart_empty[zz,z,r,tt] for tt = (t-TD[zz,z]+1):t) )
					    end
					    if t+TD[zz,z] >= tw_min && t+TD[zz,z] <= tw_max && t+1 <= t+TD[zz,z]
					    @constraint(HY, vN_Bat_travel_full[zz,z,r,t] >= sum(vN_Bat_arrive_full[zz,z,r,tt] for tt = (t+1):(t+TD[zz,z]) ) ) # arrive to z to zz
					    @constraint(HY, vN_Bat_travel_empty[zz,z,r,t] >= sum(vN_Bat_arrive_empty[zz,z,r,tt] for tt = (t+1):(t+TD[zz,z]) ) )
					    end

					    @constraint(HY, vBatTruckFlow[z,r,t]  == vN_Bat_discharged[z,r,t] * TruckCap[r] * (1 - BatTruckLoss[r]) -
					    vN_Bat_charged[z,r,t] * TruckCap[r]/(1 - BatTruckLoss[r]) )

						# @constraint(HY, vH2TruckFlow[z,zz,r,t]  == vH2TruckFlow_pos[z,zz,r,t] -
                        # vH2TruckFlow_neg[z,zz,r,t]  +
                        # 0
                        # )

					    # @constraint(HY, vN_Bat_charged[z,r,t] * TruckCap[r] <= vH2TruckCompressionCap[z,r]
					    # )



					end
				end
			end
		end
	end
	### End Constraints ###

	dModuleArgs["dExpressions_H2"] = dExpressions_H2
	dModuleArgs["dObjective_H2"] = dObjective_H2
	dModuleArgs["dBalance_H2"] = dBalance_H2
	return HY, dModuleArgs
end # end H2Truck module


function load_H2_inputs(setup::Dict,path::AbstractString)
	## Use appropriate directory separator depending on Mac or Windows config
    if setup["MacOrWindows"]=="Mac"
		sep = "/"
	else
		sep = "\U005c"
	end

	## Read input files
	println("About to Read In CSV Files")

    Par_setting = CSV.read(string(path,sep,"Case_setting.csv"), header=true)

    inputs_H2 = Dict()

	inputs_H2["Power2H2"] = 0

    inputs_H2["Z_total"] = 7
    inputs_H2["PT_total"] = 1
    inputs_H2["RT_total"] = 2
    inputs_H2["K_prod_total"] = 8
    inputs_H2["K_stor_total"] = 2
    inputs_H2["G_H2_total"] = 1

    # Z_set = [1,2,3,4,5,6]


    Z_set = collect(skipmissing(Par_setting[!,:Zone_set]))
    PT_set = collect(skipmissing(Par_setting[!,:PT_set]))
    RT_set = collect(skipmissing(Par_setting[!,:RT_set]))
    K_prod_set = collect(skipmissing(Par_setting[!,:K_prod_set]))
    K_prod_set_int = collect(skipmissing(Par_setting[!,:K_prod_set_int]))
    K_stor_set = collect(skipmissing(Par_setting[!,:K_stor_set]))
    G_H2_set = collect(skipmissing(Par_setting[!,:G_H2_set]))

    inputs_H2["Z_set"] = Z_set
    inputs_H2["PT_set"] = PT_set
    inputs_H2["RT_set"] = RT_set
    inputs_H2["K_prod_set"] = K_prod_set
    inputs_H2["K_prod_set_int"] = K_prod_set_int
    inputs_H2["K_stor_set"] = K_stor_set
    inputs_H2["G_H2_set"] = G_H2_set

    inputs_H2["Z"] =  length(Z_set)    # Number of zones
    inputs_H2["PT"] = length(PT_set)    # Number of pipeline type
    inputs_H2["RT"] = length(RT_set)    # Number of truck type
    inputs_H2["K_prod"] = length(K_prod_set)    # Number of H2 production type
    inputs_H2["K_prod_int"] = length(K_prod_set_int)    # Number of H2 production type
    inputs_H2["K_stor"] = length(K_stor_set)     # Number of H2 storage type
    inputs_H2["G_H2"] = length(G_H2_set)   # Number of H2 power generation type
    inputs_H2["W"] = collect(skipmissing(Par_setting[!,:Num_week]))[1]     # Number of weeks
    inputs_H2["Tw"] = 7*24 # Number of time steps (hours)

    inputs_H2["W"] = setup["W"]
	inputs_H2["Tw"] = setup["Tw"]

	par_list = names(Par_setting)
	for par in par_list
	    if size(collect(skipmissing(Par_setting[!,par])))[1] == 1
	        inputs_H2[par] = collect(skipmissing(Par_setting[!,par]))[1]
	    elseif size(collect(skipmissing(Par_setting[!,par])))[1] > 1
	        inputs_H2[par] = collect(skipmissing(Par_setting[!,par]))
	    end

	end

    # Truck_option_avai = 0
    # no_integer = 1
    if inputs_H2["no_integer"] == 1
        inputs_H2["Truck_integer_model"] = 0
    end

    inputs_H2["T"] = inputs_H2["W"] * inputs_H2["Tw"]    # Number of time steps (hours)

    ## Set indices
    Z_total = inputs_H2["Z_total"]
    PT_total = inputs_H2["PT_total"]
    RT_total = inputs_H2["RT_total"]
    K_prod_total = inputs_H2["K_prod_total"]
    K_stor_total = inputs_H2["K_stor_total"]
    G_H2_total = inputs_H2["G_H2_total"]

    T = inputs_H2["T"]
    Z = inputs_H2["Z"]
    PT = inputs_H2["PT"]
    RT = inputs_H2["RT"]

    K_prod = inputs_H2["K_prod"]
    K_prod_int = inputs_H2["K_prod_int"]
    K_stor = inputs_H2["K_stor"]

    G_H2 = inputs_H2["G_H2"]
    W = inputs_H2["W"]
    Tw = inputs_H2["Tw"]

    zone_distance = CSV.read(string(path,sep,"zone-distances-km.csv"), header=true)
    # inputs_H2["zone_distance"] = zone_distance

    RouteLength = zone_distance[Z_set,Z_set.+1]
    inputs_H2["RouteLength"] = RouteLength

    # Gas_Price_data = collect(skipmissing(Par_setting[!,:GasPrice]))
	#
    # GasPrice = Gas_Price_factor*(Gas_Price_data .* ones(Z_total,T))[Z_set,:] # $/MMBtu
	#
    # inputs_H2["GasPrice"] = GasPrice

	Avg_Hrly_Demand = CSV.read(string(path,sep,"hydrogen_transportation_demand.csv"), header=true)
	total_LDV_demand = collect(skipmissing(Avg_Hrly_Demand[!,:LDV
	])) # tonne-H2/hour
	total_HDV_demand = collect(skipmissing(Avg_Hrly_Demand[!,:HDV
	])) # tonne-H2/hour
	# Avg_Hrly_Demand_50 = [12.475,79.520,28.462,61.375,19.376,27.291]
	inputs_H2["total_LDV_demand"] = total_LDV_demand
	inputs_H2["total_HDV_demand"] = total_HDV_demand
	BEV_convert_factor = 33.3*65/120 ## MW/(tonne/hour)
	inputs_H2["BEV_convert_factor"] = BEV_convert_factor

	if inputs_H2["H2T_constant"] == 1
	    FCEV_avg_demand = (total_LDV_demand+total_HDV_demand)*inputs_H2["FCEV_penetration"]
	    Random.seed!(2019)
	    demand_profile = rand(1,T)
	    demand_lo = 0.6
	    demand_profile = ones(1,T)
	    # H2T = round.(Int,Avg_Hrly_Demand) * (demand_profile.*(1-demand_lo).+demand_lo)
		if inputs_H2["hydrogen_demand_option"] == 0
			demand_profile = zeros(1,T)
		end
		inputs_H2["unit_FCEV_profile"] = demand_profile
		H2T = round.(Int,FCEV_avg_demand) * demand_profile
		inputs_H2["FCEV_demand"] = transpose(H2T)
	else
		FCEV_profile_data = CSV.read(string(path,sep,"FCEV_demand_zone.csv"), copycols=true, header=true)
		# start = findall(s -> s == Symbol("1"),names(FCEV_profile_data))[1]
		start = 2
		FCEV_profile = convert(Matrix, FCEV_profile_data[1:inputs_H2["T"],start:start-1+inputs_H2["Z"]])
		unit_FCEV_profile = convert(Matrix, FCEV_profile_data[1:inputs_H2["T"],start:start-1+inputs_H2["Z"]])./(mean(FCEV_profile, dims=1).+0.0001)
		FCEV_avg_demand = transpose(total_LDV_demand+total_HDV_demand)*inputs_H2["FCEV_penetration"]
		FCEV_demand = transpose(unit_FCEV_profile.*FCEV_avg_demand)
		inputs_H2["FCEV_demand"] = FCEV_demand
		inputs_H2["unit_FCEV_profile"] = unit_FCEV_profile
		H2T= FCEV_demand
	end

	BEV_profile_data = CSV.read(string(path,sep,"BEV_demand_zone.csv"), copycols=true, header=true)
	# start = findall(s -> s == "1",names(BEV_profile_data))[1]
	start = 2
	BEV_profile = convert(Matrix, BEV_profile_data[1:inputs_H2["T"],start:start-1+inputs_H2["Z"]])
	unit_BEV_profile = convert(Matrix, BEV_profile_data[1:inputs_H2["T"],start:start-1+inputs_H2["Z"]])./(mean(BEV_profile, dims=1).+0.0001)
	BEV_avg_demand = transpose((total_LDV_demand+total_HDV_demand)*BEV_convert_factor)*inputs_H2["BEV_penetration"]
	BEV_demand = unit_BEV_profile.*BEV_avg_demand/1000
	inputs_H2["BEV_demand"] = BEV_demand
	inputs_H2["unit_BEV_profile"] = unit_BEV_profile

	pD_extra = zeros(T,Z)
	inputs_H2["pD_extra"] = pD_extra

    inputs_H2["H2T"] = H2T
	# 	prices_grnfld_EI0 = CSV.read(string(path,sep,"prices_grnfld_EI0.csv"), header=true)
	# 	prices_grnfld_EI2 = CSV.read(string(path,sep,"prices_grnfld_EI2.csv"), header=true)
	# 	prices_grnfld_EI50 = CSV.read(string(path,sep,"prices_grnfld_EI50.csv"), header=true)
	# 	prices_grnfld_EI100 = CSV.read(string(path,sep,"prices_grnfld_EI100.csv"), header=true)
	# 	prices_grnfld_EI1000 = CSV.read(string(path,sep,"prices_grnfld_EI1000.csv"), header=true)

	    # electricity_price_EI0 = prices_grnfld_EI0[:,Z_set.+1]
	    # electricity_price_EI2 = prices_grnfld_EI2[:,Z_set.+1]
	    # electricity_price_EI50 = prices_grnfld_EI50[:,Z_set.+1]
	    # electricity_price_EI100 = prices_grnfld_EI100[:,Z_set.+1]
	    # electricity_price_EI1000 = prices_grnfld_EI1000[:,Z_set.+1]

	#     electricity_price_EI0 = prices_grnfld_EI0[Z_set,2:T+1]
	#     electricity_price_EI2 = prices_grnfld_EI2[Z_set,2:T+1]
	#     electricity_price_EI50 = prices_grnfld_EI50[Z_set,2:T+1]
	#     electricity_price_EI100 = prices_grnfld_EI100[Z_set,2:T+1]
	#     electricity_price_EI1000 = prices_grnfld_EI1000[Z_set,2:T+1]

		# electricity_price = ones(Z,T) .* [50,5] # $/MWh
	# 	if emission_scenario == 0
	# 		electricity_price = transpose(convert(Matrix, electricity_price_EI0))
	# 	elseif emission_scenario == 2
	# 		electricity_price = transpose(convert(Matrix, electricity_price_EI2))
	# 	elseif emission_scenario == 50
	# 		electricity_price = transpose(convert(Matrix, electricity_price_EI50))
	# 	elseif emission_scenario == 100
	# 		electricity_price = transpose(convert(Matrix, electricity_price_EI100))
	# 	else
	# 		electricity_price = transpose(convert(Matrix, electricity_price_EI1000))
	# 	end
		# electricity_price = electricity_price_EI1000



	if setup["single_hydrogen_sector"] == 1
		e_prices = CSV.read(string(path,sep,"PRICE.csv"), header=true)
		electricity_price_df = e_prices[Z_set,2:T+1]
		electricity_price = transpose(convert(Matrix, electricity_price_df))
    else
        electricity_price = ones(Z,T) .* 0 # $/MWh
	end


	for pi = 1:size(electricity_price)[1]
		for pj = 1:size(electricity_price)[2]
			electricity_price[pi,pj] = round.(electricity_price[pi,pj],digits = 1)
			if abs(electricity_price[pi,pj]) < 1e-3
				electricity_price[pi,pj] = 0
			end
		end
	end

	inputs_H2["electricity_price"] = electricity_price
	# electricity_price

	weights = setup["omega"]
    inputs_H2["weights"] = weights


	datetimenow = Dates.format(now(),"mmddyy_HHMMSS")

    println("Data Successfully Read!")
    println("Current time:",datetimenow)


    ## Set parameters
    discount_rate = 0.07
    life = 20
	inputs_H2["discount_rate"] = discount_rate
	inputs_H2["life"] = life


    flag = 1
    # H2 Generator
    if flag == 1
        H2GenData = CSV.read(string(path,sep,"H2_generation.csv"), header=true)
        LHV = 33.4 # MWh/tonne-H2
        etaP2G = (ones(K_prod_total,Z_total) .* H2GenData[!,:etaP2G_MWh_per_tonne])[K_prod_set,Z_set]# MWh/tonne-H2
        etaGas = (ones(K_prod_total,Z_total) .* H2GenData[!,:etaGas_MMBtu_per_tonne])[K_prod_set,Z_set] # MMBtu-gas/tonne-H2
        # GasPrice = ones(Z,T) .* [30,3] # $/MMBtu

		Gas_Price_data = (ones(K_prod_total,Z_total,T) .* H2GenData[!,:GasPrice])[K_prod_set,Z_set,:] # $/MMBtu

	    GasPrice = inputs_H2["Gas_Price_factor"]*Gas_Price_data # $/MMBtu

		inputs_H2["GasPrice"] = GasPrice


        rhoH2Gen_min = (ones(K_prod_total,Z_total) .* H2GenData[!,:rhoH2Gen_min])[K_prod_set,Z_set] #
        rhoH2Gen_max = (ones(K_prod_total,Z_total) .* H2GenData[!,:rhoH2Gen_max])[K_prod_set,Z_set] #
        H2GenSize = (ones(K_prod_total,Z_total) .* H2GenData[!,:H2GenSize_tonne_per_hour])[K_prod_set,Z_set] # tonne-H2/hour
        # P2GenUnitCapex = H2GenSize .* [900*53/LHV,910] .* 1000 # $/unit
        P2GenUnitCapex = (H2GenData[!,:P2GenUnitCapex_per_unit])[K_prod_set,:] # $/unit
		start_up_cost = (round.(Int,ones(K_prod_total,Z_total) .* H2GenData[!,:start_up_cost]))[K_prod_set,Z_set] # $
        tau_up = (round.(Int,ones(K_prod_total,Z_total) .* H2GenData[!,:tau_up_hour]))[K_prod_set,Z_set] # hour
        tau_down = (round.(Int,ones(K_prod_total,Z_total) .* H2GenData[!,:tau_down_hour]))[K_prod_set,Z_set] # hour
        min_up_ratio = ((ones(K_prod_total,Z_total) .* H2GenData[!,:min_up_ratio]))[K_prod_set,Z_set]
        gen_emission_rate = ((ones(K_prod_total,Z_total) .* H2GenData[!,:generation_emission_rate_tonne_per_tonne]))[K_prod_set,Z_set]

		P2Genlifetime = ((ones(K_prod_total,Z_total) .* H2GenData[!,:lifetime]))[K_prod_set,Z_set]

        inputs_H2["etaP2G"] = etaP2G
        inputs_H2["etaGas"] = etaGas
        inputs_H2["GasPrice"] = GasPrice
        inputs_H2["rhoH2Gen_min"] = rhoH2Gen_min
        inputs_H2["rhoH2Gen_max"] = rhoH2Gen_max
        inputs_H2["H2GenSize"] = H2GenSize
        inputs_H2["P2GenUnitCapex"] = P2GenUnitCapex
		inputs_H2["start_up_cost"] = start_up_cost
        inputs_H2["tau_up"] = tau_up
        inputs_H2["tau_down"] = tau_down
        inputs_H2["min_up_ratio"] = min_up_ratio
        inputs_H2["gen_emission_rate"] = gen_emission_rate
		inputs_H2["P2Genlifetime"] = P2Genlifetime
		life = inputs_H2["P2Genlifetime"]
		inputs_H2["discount_factor_P2Gen"] = discount_rate./(float(1).-(1+discount_rate).^(-life))

    end
    if flag == 1
        LHV = 33.4 # MWh/tonne-H2
        etaP2G_int = (ones(K_prod_total,Z_total) .* H2GenData[!,:etaP2G_MWh_per_tonne])[K_prod_set_int,Z_set]# MWh/tonne-H2
        etaGas_int = (ones(K_prod_total,Z_total) .* H2GenData[!,:etaGas_MMBtu_per_tonne])[K_prod_set_int,Z_set] # MMBtu-gas/tonne-H2
        # GasPrice = ones(Z,T) .* [30,3] # $/MMBtu

		Gas_Price_data_int = (ones(K_prod_total,Z_total,T) .* H2GenData[!,:GasPrice])[K_prod_set_int,Z_set,:] # $/MMBtu

	    GasPrice_int = inputs_H2["Gas_Price_factor"]*Gas_Price_data_int # $/MMBtu

		inputs_H2["GasPrice_int"] = GasPrice_int

        rhoH2Gen_min_int = (ones(K_prod_total,Z_total) .* H2GenData[!,:rhoH2Gen_min])[K_prod_set_int,Z_set] #
        rhoH2Gen_max_int = (ones(K_prod_total,Z_total) .* H2GenData[!,:rhoH2Gen_max])[K_prod_set_int,Z_set] #
        H2GenSize_int = (ones(K_prod_total,Z_total) .* H2GenData[!,:H2GenSize_tonne_per_hour])[K_prod_set_int,Z_set] # tonne-H2/hour
        # P2GenUnitCapex = H2GenSize .* [900*53/LHV,910] .* 1000 # $/unit
        P2GenUnitCapex_int = (H2GenData[!,:P2GenUnitCapex_per_unit])[K_prod_set_int,:] # $/unit
		start_up_cost_int = (round.(Int,ones(K_prod_total,Z_total) .* H2GenData[!,:start_up_cost]))[K_prod_set_int,Z_set] # $
		tau_up_int = (round.(Int,ones(K_prod_total,Z_total) .* H2GenData[!,:tau_up_hour]))[K_prod_set_int,Z_set] # hour
        tau_down_int = (round.(Int,ones(K_prod_total,Z_total) .* H2GenData[!,:tau_down_hour]))[K_prod_set_int,Z_set] # hour
        min_up_ratio_int = ((ones(K_prod_total,Z_total) .* H2GenData[!,:min_up_ratio]))[K_prod_set_int,Z_set]
        gen_emission_rate_int = ((ones(K_prod_total,Z_total) .* H2GenData[!,:generation_emission_rate_tonne_per_tonne]))[K_prod_set_int,Z_set]

		P2Genlifetime_int = ((ones(K_prod_total,Z_total) .* H2GenData[!,:lifetime]))[K_prod_set_int,Z_set]

		inputs_H2["etaP2G_int"] = etaP2G_int
        inputs_H2["etaGas_int"] = etaGas_int
        inputs_H2["rhoH2Gen_min_int"] = rhoH2Gen_min_int
        inputs_H2["rhoH2Gen_max_int"] = rhoH2Gen_max_int
        inputs_H2["H2GenSize_int"] = H2GenSize_int
        inputs_H2["P2GenUnitCapex_int"] = P2GenUnitCapex_int
		inputs_H2["start_up_cost_int"] = start_up_cost_int
        inputs_H2["tau_up_int"] = tau_up_int
        inputs_H2["tau_down_int"] = tau_down_int
        inputs_H2["min_up_ratio_int"] = min_up_ratio_int
        inputs_H2["gen_emission_rate_int"] = gen_emission_rate_int
		inputs_H2["P2Genlifetime_int"] = P2Genlifetime_int
		life = inputs_H2["P2Genlifetime_int"]
		inputs_H2["discount_factor_P2Gen_int"] = discount_rate./(float(1).-(1+discount_rate).^(-life))
    end

    # H2 to Power
    PfGData = CSV.read(string(path,sep,"H2_to_power.csv"), header=true)

    etaPfG = (ones(G_H2_total,Z_total) * LHV .* PfGData[!,:etaPfG])[G_H2_set,Z_set] # MWh/tonne
    rhoPfG_min = (ones(G_H2_total,Z_total) .* PfGData[!,:rhoPfG_min])[G_H2_set,Z_set] #
    rhoPfG_max = (ones(G_H2_total,Z_total) .* PfGData[!,:rhoPfG_max])[G_H2_set,Z_set] #
    PfGSize = (ones(G_H2_total,Z_total) .* PfGData[!,:PfGSize_MWe])[G_H2_set,Z_set] # MW-e
    PfGUnitCapex = (ones(G_H2_total,Z_total).* PfGData[!,:PfGUnitCapex_per_unit])[G_H2_set,Z_set]# $/unit
	PfGVOM = (ones(G_H2_total,Z_total).* PfGData[!,:PfG_VOM_per_MWh])[G_H2_set,Z_set]# $/MWh
	PfGFOM = (ones(G_H2_total,Z_total).* PfGData[!,:PfG_FOM_per_MW_yr])[G_H2_set,Z_set]# $/MW-yr
    tau_up_g = (round.(Int,ones(G_H2_total,Z_total) .* PfGData[!,:tau_up_hour]))[G_H2_set,Z_set] # hour
    tau_down_g = (round.(Int,ones(G_H2_total,Z_total) .* PfGData[!,:tau_down_hour]))[G_H2_set,Z_set] # hour

	PfGlifetime = (round.(Int,ones(G_H2_total,Z_total) .* PfGData[!,:lifetime]))[G_H2_set,Z_set]

    inputs_H2["etaPfG"] = etaPfG
    inputs_H2["rhoPfG_min"] = rhoPfG_min
    inputs_H2["rhoPfG_max"] = rhoPfG_max
    inputs_H2["PfGSize"] = PfGSize
    inputs_H2["PfGUnitCapex"] = PfGUnitCapex
	inputs_H2["PfGVOM"] = PfGVOM
	inputs_H2["PfGFOM"] = PfGFOM
    inputs_H2["tau_up_g"] = tau_up_g
    inputs_H2["tau_down_g"] = tau_down_g
	inputs_H2["PfGlifetime"] = PfGlifetime

	life = inputs_H2["PfGlifetime"]
	inputs_H2["discount_factor_PfG"] = discount_rate./(float(1).-(1+discount_rate).^(-life))


    # H2 Storage
    H2StorData = CSV.read(string(path,sep,"H2_storage.csv"), header=true)

    etaStor = zeros(K_stor_total,Z_total) #
    rhoH2Stor_min = zeros(K_stor_total,Z_total)
    H2StorUnitCapex = zeros(K_stor_total,Z_total)
    H2StorCap_min = zeros(K_stor_total,Z_total)
    H2StorCap_max = zeros(K_stor_total,Z_total)
    H2StorCap_avai = zeros(K_stor_total,Z_total)
	H2Storlifetime = zeros(K_stor_total,Z_total)
    # for s in unique(H2StorData[!,:Storage_Type])
    #     for z in unique(H2StorData[!,:Zone])
    #     etaStor[s,z] = H2StorData[(H2StorData[!,:Storage_Type] == s) && (H2StorData[!,:Zone] == z),:etaStor]
    #     println(etaStor[s,z])
    # end
    # end
    storage_set = unique(H2StorData[!,:Storage_Type])
    zone_set_unique = unique(H2StorData[!,:Zone])
    for s in 1:size(storage_set)[1]
        # for z in 1:size(zone_set_unique)[1]
            etaStor[s,:] = H2StorData[(H2StorData[!,:Storage_Type].== storage_set[s]),:etaStor]
            rhoH2Stor_min[s,:] = H2StorData[(H2StorData[!,:Storage_Type].== storage_set[s]),:rhoH2Stor_min]
            H2StorUnitCapex[s,:] = H2StorData[(H2StorData[!,:Storage_Type].== storage_set[s]),:H2StorUnitCapex_per_tonne]
            H2StorCap_min[s,:] = H2StorData[(H2StorData[!,:Storage_Type].== storage_set[s]),:H2StorCap_min_tonne]
            H2StorCap_max[s,:] = H2StorData[(H2StorData[!,:Storage_Type].== storage_set[s]),:H2StorCap_max_tonne]
            H2StorCap_avai[s,:] = H2StorData[(H2StorData[!,:Storage_Type].== storage_set[s]),:H2StorCap_avai]
			H2Storlifetime[s,:] = H2StorData[(H2StorData[!,:Storage_Type].== storage_set[s]),:lifetime]
    # end
    end
    inputs_H2["etaStor"] = etaStor
    inputs_H2["rhoH2Stor_min"] = rhoH2Stor_min
    inputs_H2["H2StorUnitCapex"] = H2StorUnitCapex
    inputs_H2["H2StorCap_min"] = H2StorCap_min
    inputs_H2["H2StorCap_max"] = H2StorCap_max
    inputs_H2["H2StorCap_avai"] = H2StorCap_avai
	inputs_H2["H2Storlifetime"] = H2Storlifetime

	life = inputs_H2["H2Storlifetime"]
	inputs_H2["discount_factor_H2Stor"] = discount_rate./(float(1).-(1+discount_rate).^(-life))

    # H2StorRate = ones(K_stor,Z) * 100 # tonne-H2
	#     etaStor = (ones(K_stor_total,Z_total) .* H2StorData[!,:etaStor])[K_stor_set,Z_set] #
	#     rhoH2Stor_min = (ones(K_stor_total,Z_total) .* H2StorData[!,:rhoH2Stor_min
	# ])[K_stor_set,Z_set] #
	#     H2StorUnitCapex = (ones(K_stor_total) .* H2StorData[!,:H2StorUnitCapex
	# ])[K_stor_set,:] # $/tonne-H2
    #
    H2StorCompressionEnergy = (ones(K_stor_total) .* [0.5,2])[K_stor_set,:] # MWh/tonne
    H2StorCompressionUnitCapex = (ones(K_stor_total).*[10000,500000])[K_stor_set,:] # $/(tonne-H2/hour)

	inputs_H2["H2StorCompressionEnergy"] = H2StorCompressionEnergy
	inputs_H2["H2StorCompressionUnitCapex"] = H2StorCompressionUnitCapex

    # H2 Transportation
    # Pipeline
    H2PipeData = CSV.read(string(path,sep,"H2_pipeline.csv"), header=true)
    # H2PipeUnitCapex = ones(Z,Z,PT) * 741162 # $/mile
    H2PipeUnitCapex = (ones(PT_total) .* H2PipeData[!,:H2PipeUnitCapex_per_mile])[PT_set,:] # $/mile
    PipeLength = RouteLength

    H2PipeFlowSize = (ones(PT_total) .* H2PipeData[!,:H2PipeFlowSize_tonne_per_hour])[PT_set,:] # tonne-H2/hour
    H2PipeCap = (ones(PT_total) .* H2PipeData[!,:H2PipeCap_tonne_per_mile])[PT_set,:] # tonne-H2/mile
    rhoH2PipeCap_min = (ones(Z_total,Z_total,PT) .* H2PipeData[!,:rhoH2PipeCap_min])[Z_set,Z_set,PT_set] #

	H2Pipelifetime = (ones(PT_total) .* H2PipeData[!,:lifetime])[PT_set,:]

    H2PipeCompressionUnitCapex = (ones(PT_total) * 746600)[PT_set,:] # $/(tonne/hour)
    H2PipeCompressionUnitOpex = (ones(PT_total) * 0)[PT_set,:]
    H2PipeCompressionEnergy = (ones(PT_total) * 1)[PT_set,:] # MWh/tonne

    Number_online_compression = PipeLength./70
    H2PipeCompressionOnlineUnitCapex = H2PipeCompressionUnitCapex /221000*15000 # $/(tonne/hour)
    H2PipeCompressionOnlineUnitOpex = (ones(PT_total) * 0)[PT_set,:]
    H2PipeCompressionOnlineEnergy = (ones(PT_total) * 1)[PT_set,:] # MWh/tonne

    inputs_H2["H2PipeUnitCapex"] = H2PipeUnitCapex
    inputs_H2["PipeLength"] = PipeLength
    inputs_H2["H2PipeFlowSize"] = H2PipeFlowSize
    inputs_H2["H2PipeCap"] = H2PipeCap
    inputs_H2["rhoH2PipeCap_min"] = rhoH2PipeCap_min
	inputs_H2["H2Pipelifetime"] = H2Pipelifetime

	life = inputs_H2["H2Pipelifetime"]
	inputs_H2["discount_factor_H2Pipe"] = discount_rate./(float(1).-(1+discount_rate).^(-life))

	inputs_H2["Number_online_compression"] = Number_online_compression
	inputs_H2["H2PipeCompressionUnitCapex"] = H2PipeCompressionUnitCapex
	inputs_H2["H2PipeCompressionUnitOpex"] = H2PipeCompressionUnitOpex
	inputs_H2["H2PipeCompressionEnergy"] = H2PipeCompressionEnergy
	inputs_H2["H2PipeCompressionOnlineUnitCapex"] = H2PipeCompressionOnlineUnitCapex
	inputs_H2["H2PipeCompressionOnlineUnitOpex"] = H2PipeCompressionOnlineUnitOpex
	inputs_H2["H2PipeCompressionOnlineEnergy"] = H2PipeCompressionOnlineEnergy

    # Truck
    H2TruckData = CSV.read(string(path,sep,"H2_truck.csv"), header=true)
    H2TruckUnitCapex = (ones(RT_total) .* H2TruckData[!,:H2TruckUnitCapex_per_unit])[RT_set,:] # $/truck
    TruckCap = ones(RT_total) .* H2TruckData[!,:TruckCap_tonne_per_unit][RT_set,:] # tonne-H2

	Full_weight_tonne_per_unit = ones(RT_total) .* H2TruckData[!,:Full_weight_tonne_per_unit][RT_set,:] # tonne-H2
	Empty_weight_tonne_per_unit = ones(RT_total) .* H2TruckData[!,:Empty_weight_tonne_per_unit][RT_set,:] # tonne-H2

	H2Trucklifetime = (ones(RT_total) .* H2TruckData[!,:lifetime])[RT_set,:] # $/truck

    AvgTruckSpeed = H2TruckData[!,:AvgTruckSpeed_mile_per_hour][1] # mph
    H2TruckUnitOpex_full = (ones(RT_total) .* H2TruckData[!,:H2TruckUnitOpex_per_mile_full])[RT_set,:] # $/mile, fuel + driver
	H2TruckUnitOpex_empty = (ones(RT_total) .* H2TruckData[!,:H2TruckUnitOpex_per_mile_empty])[RT_set,:] # $/mile, fuel + driver
    TD = round.(Int,RouteLength./AvgTruckSpeed)

    H2TLoss = (ones(RT_total) .* H2TruckData[!,:H2TLoss_per_mile])[RT_set,:] # Loss %

    truck_emission_rate = (ones(RT_total,Z_total) .* H2TruckData[!,:truck_emission_rate_tonne_per_tonne_mile])[RT_set,Z_set] # Loss/mile

    H2TruckCompressionUnitCapex = (ones(RT_total) .* [15000 / 10 * 1000, 40000000/(30/24)])[RT_set,:] # $/(tonne/hour)
    H2TruckCompressionEnergy = (ones(RT_total) .* [1,11])[RT_set,:] # MWh/tonne


    H2TruckCompressionUnitOpex = (ones(RT_total) * 0)[RT_set,:]
    # TD = round.(Int,ones(Z,Z,RT) * 2) #

    inputs_H2["H2TruckUnitCapex"] = H2TruckUnitCapex
    inputs_H2["TruckCap"] = TruckCap
	inputs_H2["Full_weight_tonne_per_unit"] = Full_weight_tonne_per_unit
	inputs_H2["Empty_weight_tonne_per_unit"] = Empty_weight_tonne_per_unit
    inputs_H2["AvgTruckSpeed"] = AvgTruckSpeed
    inputs_H2["H2TruckUnitOpex_full"] = H2TruckUnitOpex_full
	inputs_H2["H2TruckUnitOpex_empty"] = H2TruckUnitOpex_empty
    inputs_H2["TD"] = TD
    inputs_H2["H2TLoss"] = H2TLoss
    inputs_H2["truck_emission_rate"] = truck_emission_rate
	inputs_H2["H2TruckCompressionUnitCapex"] = H2TruckCompressionUnitCapex
	inputs_H2["H2TruckCompressionEnergy"] = H2TruckCompressionEnergy
	inputs_H2["H2TruckCompressionUnitOpex"] = H2TruckCompressionUnitOpex
	inputs_H2["H2Trucklifetime"] = H2Trucklifetime

	life = inputs_H2["H2Trucklifetime"]
	inputs_H2["discount_factor_H2Truck"] = discount_rate./(float(1).-(1+discount_rate).^(-life))

    # Compression

	inputs_H2["H2Compressionlifetime"] = 30

	life = inputs_H2["H2Compressionlifetime"]
	inputs_H2["discount_factor_H2Compression"] = discount_rate./(float(1).-(1+discount_rate).^(-life))

    H2CompressionData = CSV.read(string(path,sep,"H2_compression.csv"), header=true)

    UnitCompressionGen2Stor = ones(K_prod_total,K_stor_total).*[H2CompressionData[!,:UnitCompressCost_per_tonne][2],H2CompressionData[!,:UnitCompressCost_per_tonne][3]]'
    UnitCompressionGen2Pipe = ones(K_prod_total,PT_total).*[H2CompressionData[!,:UnitCompressCost_per_tonne][2]]
    UnitCompressionGen2Truck = ones(K_prod_total,RT_total).*[H2CompressionData[!,:UnitCompressCost_per_tonne][3],H2CompressionData[!,:UnitCompressCost_per_tonne][4]]'
    UnitCompressionStor2Truck =
    ones(K_stor_total,RT_total).*[H2CompressionData[!,:UnitCompressCost_per_tonne][3] H2CompressionData[!,:UnitCompressCost_per_tonne][4]; H2CompressionData[!,:UnitCompressCost_per_tonne][1] H2CompressionData[!,:UnitCompressCost_per_tonne][4];]
    UnitCompressionPipe2Truck = ones(PT_total,RT_total).*[H2CompressionData[!,:UnitCompressCost_per_tonne][3],H2CompressionData[!,:UnitCompressCost_per_tonne][4]]'
    UnitCompressionPipe2Stor = ones(PT_total,K_stor_total).*[H2CompressionData[!,:UnitCompressCost_per_tonne][1],H2CompressionData[!,:UnitCompressCost_per_tonne][3]]'

    # UnitCompressionGen2Stor = ones(K_prod_total,K_stor_total)*30
    # UnitCompressionGen2Pipe = ones(K_prod_total,PT_total)*30
    # UnitCompressionGen2Truck = ones(K_prod_total,RT_total)*270
    # UnitCompressionStor2Truck = ones(K_stor_total,RT_total)*270
    # UnitCompressionPipe2Truck = ones(PT_total,RT_total)*270
    # UnitCompressionPipe2Stor = ones(PT_total,K_stor_total)*0

    # UnitCompressionGen2Stor = ones(K_prod_total,K_stor_total)*0
    # UnitCompressionGen2Pipe = ones(K_prod_total,PT_total)*0
    # UnitCompressionGen2Truck = ones(K_prod_total,RT_total)*0
    # UnitCompressionStor2Truck = ones(K_stor_total,RT_total)*0
    # UnitCompressionPipe2Truck = ones(PT_total,RT_total)*0
    # UnitCompressionPipe2Stor = ones(PT_total,K_stor_total)*0

    UnitEnergyGen2Stor = ones(K_prod_total,K_stor_total).*[H2CompressionData[!,:UnitCompressEnergy_MWh_per_tonne][2],H2CompressionData[!,:UnitCompressEnergy_MWh_per_tonne][3]]'
    UnitEnergyGen2Pipe = ones(K_prod_total,PT_total).*[H2CompressionData[!,:UnitCompressEnergy_MWh_per_tonne][2]]
    UnitEnergyGen2Truck = ones(K_prod_total,RT_total).*[H2CompressionData[!,:UnitCompressEnergy_MWh_per_tonne][3],H2CompressionData[!,:UnitCompressEnergy_MWh_per_tonne][4]]'
    UnitEnergyStor2Truck =
    ones(K_stor_total,RT_total).*[H2CompressionData[!,:UnitCompressEnergy_MWh_per_tonne][3] H2CompressionData[!,:UnitCompressEnergy_MWh_per_tonne][4]; H2CompressionData[!,:UnitCompressEnergy_MWh_per_tonne][1] H2CompressionData[!,:UnitCompressEnergy_MWh_per_tonne][4];]
    UnitEnergyPipe2Truck = ones(PT_total,RT_total).*[H2CompressionData[!,:UnitCompressEnergy_MWh_per_tonne][3],H2CompressionData[!,:UnitCompressEnergy_MWh_per_tonne][4]]'
    UnitEnergyPipe2Stor = ones(PT_total,K_stor_total).*[H2CompressionData[!,:UnitCompressEnergy_MWh_per_tonne][1],H2CompressionData[!,:UnitCompressEnergy_MWh_per_tonne][3]]'

    # UnitEnergyGen2Stor = ones(K_prod_total,K_stor_total)*0.3
    # UnitEnergyGen2Pipe = ones(K_prod_total,PT_total)*0.3
    # UnitEnergyGen2Truck = ones(K_prod_total,RT_total)*1
    # UnitEnergyStor2Truck = ones(K_stor_total,RT_total)*1
    # UnitEnergyPipe2Truck = ones(PT_total,RT_total)*1
    # UnitEnergyPipe2Stor = ones(PT_total,K_stor_total)*1

    # UnitEnergyGen2Stor = ones(K_prod_total,K_stor_total)*0
    # UnitEnergyGen2Pipe = ones(K_prod_total,PT_total)*0
    # UnitEnergyGen2Truck = ones(K_prod_total,RT_total)*0
    # UnitEnergyStor2Truck = ones(K_stor_total,RT_total)*0
    # UnitEnergyPipe2Truck = ones(PT_total,RT_total)*0
    # UnitEnergyPipe2Stor = ones(PT_total,K_stor_total)*0

    # H2T = [10,0] * [0,0,0,0,10,15,20,25,30,40,70,100,100,120,150,160,120,100,100,70,50,20,10,10]'
    # tonne
    H2H_min = (ones(Z_total,T) * 0)[Z_set,:] # tonne
    H2H_max = (ones(Z_total,T) * 0)[Z_set,:] # tonne
    H2I_min = (ones(Z_total,T) * 0)[Z_set,:] # tonne
    H2I_max = (ones(Z_total,T) * 0)[Z_set,:] # tonne
	# end # END data_input()

    inputs_H2["H2H_min"] = H2H_min
    inputs_H2["H2H_max"] = H2H_max
    inputs_H2["H2I_min"] = H2I_min
    inputs_H2["H2I_max"] = H2I_max
	# load_inputs(mysetup, inpath)
    inputs_H2["mysetup"] = setup
	return inputs_H2
end # END data_input()

function adjust_parameters(mysetup::Dict,inputs::Dict,inputs_H2::Dict,Scenario_Param,c)

	if "CO2_max_ton_MWh" in names(Scenario_Param)
		if mysetup["CO2Cap"] >= 1
			inputs["pMaxCO2Rate"] .= Scenario_Param[!,:CO2_max_ton_MWh][c]
		end
		inputs["CO2_max_ton_MWh"] = Scenario_Param[!,:CO2_max_ton_MWh][c]*1000
	end

	if "CO2_price" in names(Scenario_Param)
		inputs["CO2_price"] = Scenario_Param[!,:CO2_price][c]
	end
	inputs_H2["CO2_price"] = inputs["CO2_price"]

	if mysetup["reliability_test"] == 1
		task_id = c
		working_path = pwd()
		load_in = CSV.read(string(working_path,"/Inputs_Wks_365/Load_data.csv"), copycols=true, header=true)
		gen_var = CSV.read(string(working_path,"/Inputs_Wks_365/Generators_variability.csv"), copycols=true, header=true)
		# start = findall(s -> s == Symbol("Load_MW_z1"), names(load_in))[1]
		start = 6
		if mysetup["ParameterScale"] ==1  # Parameter scaling turned on
			# Demand in MW
			inputs["pD"] =convert(Matrix, load_in[(inputs["T"]*(task_id-1)+1) : inputs["T"]*task_id ,start:start-1+inputs["Z"]])*mysetup["power_load_factor"]/1e+3  # convert to GW

		else # No scaling
			# Demand in MW
			inputs["pD"] =convert(Matrix, load_in[(inputs["T"]*(task_id-1)+1) : inputs["T"]*task_id ,start:start-1+inputs["Z"]])*mysetup["power_load_factor"]
		end
		inputs["pP_Max"] = transpose(convert(Matrix{Float64}, gen_var[(inputs["T"]*(task_id-1)+1) : inputs["T"]*task_id,2:(inputs["G"]+1)]))
	end

	if "H2Gen_avai" in names(Scenario_Param)
		inputs_H2["H2Gen_avai"] = Scenario_Param[!,:H2Gen_avai][c]
	end

	if "H2_emission_rate" in names(Scenario_Param)
		inputs_H2["H2_emission_rate"] = Scenario_Param[!,:H2_emission_rate][c]
	end

	if "H2Gen_integer" in names(Scenario_Param)
		inputs_H2["H2Gen_integer"] = Scenario_Param[!,:H2Gen_integer][c]
	end

	if "hydrogen_demand_option" in names(Scenario_Param)
		inputs_H2["hydrogen_demand_option"] = Scenario_Param[!,:hydrogen_demand_option][c]
	end
	if "hydrogen_demand_factor" in names(Scenario_Param)
		inputs_H2["hydrogen_demand_factor"] = Scenario_Param[!,:hydrogen_demand_factor][c]
	end
	if "BEV_penetration" in names(Scenario_Param)
		inputs_H2["BEV_penetration"] = Scenario_Param[!,:BEV_penetration][c]
	end
	if "FCEV_penetration" in names(Scenario_Param)
		inputs_H2["FCEV_penetration"] = Scenario_Param[!,:FCEV_penetration][c]
	end
	FCEV_avg_demand = transpose(inputs_H2["total_LDV_demand"] + inputs_H2["total_HDV_demand"])*inputs_H2["FCEV_penetration"]
	FCEV_demand = transpose(inputs_H2["unit_FCEV_profile"].*FCEV_avg_demand)
	inputs_H2["FCEV_demand"] = FCEV_demand
	inputs_H2["H2T"]= FCEV_demand* inputs_H2["hydrogen_demand_factor"]
	BEV_avg_demand = transpose((inputs_H2["total_LDV_demand"] + inputs_H2["total_HDV_demand"])*inputs_H2["BEV_convert_factor"])*inputs_H2["BEV_penetration"]
	BEV_demand = inputs_H2["unit_BEV_profile"].*BEV_avg_demand/1000
	inputs_H2["BEV_demand"] = BEV_demand

	if "SMR_option_avai" in names(Scenario_Param)
		inputs_H2["SMR_option_avai"] = Scenario_Param[!,:SMR_option_avai][c]
	end
	if "EL_option_avai" in names(Scenario_Param)
		inputs_H2["EL_option_avai"] = Scenario_Param[!,:EL_option_avai][c]
	end
	if "EL_cost_factor" in names(Scenario_Param)
		inputs_H2["EL_cost_factor"] = Scenario_Param[!,:EL_cost_factor][c]
	end

	if "hydrogen_generation_option" in names(Scenario_Param)
		inputs_H2["hydrogen_generation_option"] = Scenario_Param[!,:hydrogen_generation_option][c]
	end
	if "Truck_option_avai" in names(Scenario_Param)
		inputs_H2["Truck_option_avai"] = Scenario_Param[!,:Truck_option_avai][c]
	end
	if "H2_to_power_option" in names(Scenario_Param)
		inputs_H2["H2_to_power_option"] = Scenario_Param[!,:H2_to_power_option][c]
	end


	if "TD_factor" in names(Scenario_Param)
		inputs_H2["TD_factor"] = Scenario_Param[!,:TD_factor][c]
		inputs_H2["TD"] = round.(Int,inputs_H2["RouteLength"]./inputs_H2["AvgTruckSpeed"].*inputs_H2["TD_factor"])
	end

	if "Pipe_option_avai" in names(Scenario_Param)
		inputs_H2["Pipe_option_avai"] = Scenario_Param[!,:Pipe_option_avai][c]
	end

	if "Pipe_integer" in names(Scenario_Param)
		inputs_H2["Pipe_integer"] = Scenario_Param[!,:Pipe_integer][c]
	end

	if "Pipe_cost_factor" in names(Scenario_Param)
		inputs_H2["Pipe_cost_factor"] = Scenario_Param[!,:Pipe_cost_factor][c]
	end

	if "battery_truck_option" in names(Scenario_Param)
		mysetup["battery_truck_option"] = Scenario_Param[!,:battery_truck_option][c]
	end

	if "LDS" in names(Scenario_Param)
		mysetup["LDS"] = Scenario_Param[!,:LDS][c]
	end

	if "H2_to_power_cost_factor" in names(Scenario_Param)
		inputs_H2["H2_to_power_cost_factor"] = Scenario_Param[!,:H2_to_power_cost_factor][c]
	end

	if "extra_pD_avai" in names(Scenario_Param)
		inputs_H2["extra_pD_avai"] = Scenario_Param[!,:extra_pD_avai][c]
	end
	if inputs_H2["extra_pD_avai"] == 1
		inputs_H2["pD_extra"] = transpose(inputs_H2["H2T"]*inputs_H2["etaP2G"][1]/1000)
	end

	if "Transmission_cost_factor" in names(Scenario_Param)
		inputs_H2["Transmission_cost_factor"] = Scenario_Param[!,:Transmission_cost_factor][c]
	end

	if inputs_H2["hydrogen_generation_option"] == 0
		inputs_H2["hydrogen_demand_option"] = 0
	end
	if inputs_H2["hydrogen_demand_option"] == 0
		inputs_H2["H2T"] = inputs_H2["H2T"] * 0
	end

	if "share_truck_options" in names(Scenario_Param)
		inputs_H2["share_truck_options"] = Scenario_Param[!,:share_truck_options][c]
	end

	if "warm_start" in names(Scenario_Param)
		inputs_H2["warm_start"] = Scenario_Param[!,:warm_start][c]
	end
	if "fix_investment_option" in names(Scenario_Param)
		inputs_H2["fix_investment_option"] = Scenario_Param[!,:fix_investment_option][c]
	end

	return inputs, inputs_H2, mysetup
end

function print_setting(setup::Dict,inputs::Dict,inputs_H2::Dict)
	println("----------------------------------------------------------")
	println("----------------------------------------------------------")
	println("Case Setting:")
	println("----------------------------------------------------------")
	println("Z_set = ", inputs_H2["Z_set"])
	println("PT_set = ", inputs_H2["PT_set"])
	println("RT_set = ", inputs_H2["RT_set"])
	println("K_prod_set = ", inputs_H2["K_prod_set"])
	println("K_prod_set_int = ", inputs_H2["K_prod_set_int"])
	println("K_stor_set = ", inputs_H2["K_stor_set"])
	println("G_H2_set = ", inputs_H2["G_H2_set"])
	println("W = ", inputs_H2["W"])
	println("FCEV_penetration = ", inputs_H2["FCEV_penetration"])
	println("emission = ", inputs_H2["emission_scenario"])
	println("----------------------------------------------------------")
	println("extra_pD_avai = ", inputs_H2["extra_pD_avai"])
	println("hydrogen_demand_option = ", inputs_H2["hydrogen_demand_option"])
	println("hydrogen_demand_factor = ", inputs_H2["hydrogen_demand_factor"])
	println("hydrogen_generation_option = ", inputs_H2["hydrogen_generation_option"])
	println("SMR_option_avai = ", inputs_H2["SMR_option_avai"])
	println("EL_option_avai = ", inputs_H2["EL_option_avai"])
	println("H2Gen_zone_option = ", inputs_H2["H2Gen_zone_option"])
	println("H2Gen_int_zone_option = ", inputs_H2["H2Gen_int_zone_option"])
	println("gen_category = ", inputs_H2["gen_category"])
	println("----------------------------------------------------------")
	println("Pipe_option_avai = ", inputs_H2["Pipe_option_avai"])
	println("Truck_option_avai = ", inputs_H2["Truck_option_avai"])
	println("Truck_integer_model = ", inputs_H2["Truck_integer_model"])
	println("H2_to_power_cost_factor = ", inputs_H2["H2_to_power_cost_factor"])
	println("----------------------------------------------------------")
	println("EL_cost_factor = ", inputs_H2["EL_cost_factor"])
	println("Gas_Price_factor = ", inputs_H2["Gas_Price_factor"])
	println("Transmission_cost_factor = ", inputs_H2["Transmission_cost_factor"])
	println("H2_emission_rate = ", inputs_H2["H2_emission_rate"])
	println("Power emission = ", inputs_H2["emission_scenario"])
	println("H2StorCap_avai = ", inputs_H2["H2StorCap_avai"])
	println("----------------------------------------------------------")
	println("CO2_price:",inputs["CO2_price"])
	println("----------------------------------------------------------")
	println("----------------------------------------------------------")
	return
end

function generate_master_model(setup::Dict,inputs::Dict,OPTIMIZER,modeloutput = nothing)

	MasterModel = Model(solver=OPTIMIZER)

	# inputs["T"] = 8760     # Number of time steps (hours)
	# inputs["Z"] = 1     # Number of zones
	# inputs["G"] = 1
	# inputs["K_prod"] = 1

	T = inputs["T"]     # Number of time steps (hours)
	Z = inputs["Z"]     # Number of zones

	@variable(MasterModel, 1e4 >= vPower2H2_Master[t = 1:T,z=1:Z] >=0)
	@variable(MasterModel, alpha >= -1e15)
	@objective(MasterModel, Min, alpha)

	return MasterModel
end

function fix_planning_var(setup::Dict,inputs::Dict,inputs_H2::Dict,EP::Model, dModuleArgs::Dict)
	## Loading investment values
	sep = "/"
	investment_path = "$working_path/Inputs_Investment"

	searchdir(path,key) = filter(x->occursin(key,x), readdir(path))

	file_string = searchdir(investment_path,"H2")[1]
	output_filename = string(investment_path,sep,file_string)
	@load output_filename inputs_H2 results_H2

	# file_string = searchdir(investment_path,"power")[1]
	# output_filename = string(investment_path,sep,file_string)
	# @load output_filename inputs_power results_power

	Z_total = inputs_H2["Z_total"]
	PT_total = inputs_H2["PT_total"]
	RT_total = inputs_H2["RT_total"]
	K_prod_total = inputs_H2["K_prod_total"]
	K_stor_total = inputs_H2["K_stor_total"]
	G_H2_total = inputs_H2["G_H2_total"]
	T = inputs_H2["T"]
	Z = inputs_H2["Z"]
	PT = inputs_H2["PT"]
	RT = inputs_H2["RT"]
	K_prod = inputs_H2["K_prod"]
	K_prod_int = inputs_H2["K_prod_int"]
	K_stor = inputs_H2["K_stor"]
	G_H2 = inputs_H2["G_H2"]
	W = inputs_H2["W"]
	Tw = inputs_H2["Tw"]


	if setup["single_hydrogen_sector"] == 0

		println("Fixing power sector investment")

		CAP = CSV.read(string(investment_path,sep,"CAP.csv"), header=true)


		G = inputs["G"]

		dfGen = inputs["dfGen"]

		## Fixing power sector planning variables
		# CAP = results_power["CAP"]
		vRETCAP = CAP[1:end-1,:RetCap]./dfGen[!,:Cap_size]
		vCAP = CAP[1:end-1,:NewCap]./dfGen[!,:Cap_size]
		for y=1:G
		    JuMP.fix(EP[:vRETCAP][y],vRETCAP[y])
		    JuMP.fix(EP[:vCAP][y],vCAP[y])
		end

		aux_i = (dfGen[(dfGen[!,:STOR].>=2),:][!,:R_ID])
		vCAPSTORAGE = CAP[1:end-1,:NewEnergyCap][aux_i]
		vRETCAPSTORAGE = CAP[1:end-1,:RetEnergyCap][aux_i]
		for k = 1:size(aux_i)[1]
		    y=aux_i[k]
		    JuMP.fix(EP[:vCAPSTORAGE][y],vCAPSTORAGE[k])
		    JuMP.fix(EP[:vRETCAPSTORAGE][y],vRETCAPSTORAGE[k])
		end

		aux_i = (dfGen[(dfGen[!,:STOR].==3),:][!,:R_ID])
		vCAPCHARGE = CAP[1:end-1,:NewChargeCap][aux_i]
		vRETCAPCHARGE = CAP[1:end-1,:RetChargeCap][aux_i]
		for k = 1:size(aux_i)[1]
		    y=aux_i[k]
		    JuMP.fix(EP[:vCAPCHARGE][y],vCAPCHARGE[k])
		    JuMP.fix(EP[:vRETCAPCHARGE][y],vRETCAPCHARGE[k])
		end
	end

	if setup["single_power_sector"] == 0
		println("Fixing hydrogen sector investment")
		## Fixing hydrogen sector planning variables
		for key in keys(results_H2["H2PLAN"])
		    # JuMP.fix(EP[:key][y],results_H2["H2PLAN"][:key][k])
		    var_size = size(results_H2["H2PLAN"][key])
		    var_dim = length(size(results_H2["H2PLAN"][key]))

		    if var_dim == 1
		        for i = 1:var_size[1]
		            JuMP.fix(EP[Symbol(key)][i],results_H2["H2PLAN"][key][i])
		        end
		    end

		    if var_dim == 2
		        for i = 1:var_size[1]
		            for j = 1:var_size[2]
		                JuMP.fix(EP[Symbol(key)][i,j],results_H2["H2PLAN"][key][i,j])
		            end
		        end
		    end

		    if var_dim == 3
		        for i = 1:var_size[1]
		            for j = 1:var_size[2]
		                for k = 1:var_size[3]
		                    JuMP.fix(EP[Symbol(key)][i,j,k],results_H2["H2PLAN"][key][i,j,k])
		                end
		            end
		        end
		    end
		end
	    # println(results_H2["H2PLAN"][key])
	end

	return EP, dModuleArgs
end

function warm_start(setup::Dict,inputs::Dict,inputs_H2::Dict,EP::Model, dModuleArgs::Dict)
	println("Presetting investment variables")
	## Loading investment values
	sep = "/"
	warm_start_path = "$working_path/Inputs_Warm_start"

	searchdir(path,key) = filter(x->occursin(key,x), readdir(path))

	file_string = searchdir(warm_start_path,"H2")[1]
	output_filename = string(warm_start_path,sep,file_string)
	@load output_filename inputs_H2 results_H2

	# file_string = searchdir(investment_path,"power")[1]
	# output_filename = string(investment_path,sep,file_string)
	# @load output_filename inputs_power results_power

	if setup["single_hydrogen_sector"] == 0
		CAP = CSV.read(string(warm_start_path,sep,"CAP.csv"), header=true)

		G = inputs["G"]

		dfGen = inputs["dfGen"]

		## Fixing power sector planning variables
		# CAP = results_power["CAP"]
		vRETCAP = CAP[1:end-1,:RetCap]./dfGen[!,:Cap_size]
		vCAP = CAP[1:end-1,:NewCap]./dfGen[!,:Cap_size]
		for y=1:G
		    JuMP.setvalue(EP[:vRETCAP][y],vRETCAP[y])
		    JuMP.setvalue(EP[:vCAP][y],vCAP[y])
		end

		aux_i = (dfGen[(dfGen[!,:STOR].>=2),:][!,:R_ID])
		vCAPSTORAGE = CAP[1:end-1,:NewEnergyCap][aux_i]
		vRETCAPSTORAGE = CAP[1:end-1,:RetEnergyCap][aux_i]
		for k = 1:size(aux_i)[1]
		    y=aux_i[k]
		    JuMP.setvalue(EP[:vCAPSTORAGE][y],vCAPSTORAGE[k])
		    JuMP.setvalue(EP[:vRETCAPSTORAGE][y],vRETCAPSTORAGE[k])
		end

		aux_i = (dfGen[(dfGen[!,:STOR].==3),:][!,:R_ID])
		vCAPCHARGE = CAP[1:end-1,:NewChargeCap][aux_i]
		vRETCAPCHARGE = CAP[1:end-1,:RetChargeCap][aux_i]
		for k = 1:size(aux_i)[1]
		    y=aux_i[k]
		    JuMP.setvalue(EP[:vCAPCHARGE][y],vCAPCHARGE[k])
		    JuMP.setvalue(EP[:vRETCAPCHARGE][y],vRETCAPCHARGE[k])
		end
	end


	Z_total = inputs_H2["Z_total"]
	PT_total = inputs_H2["PT_total"]
	RT_total = inputs_H2["RT_total"]
	K_prod_total = inputs_H2["K_prod_total"]
	K_stor_total = inputs_H2["K_stor_total"]
	G_H2_total = inputs_H2["G_H2_total"]
	T = inputs_H2["T"]
	Z = inputs_H2["Z"]
	PT = inputs_H2["PT"]
	RT = inputs_H2["RT"]
	K_prod = inputs_H2["K_prod"]
	K_prod_int = inputs_H2["K_prod_int"]
	K_stor = inputs_H2["K_stor"]
	G_H2 = inputs_H2["G_H2"]
	W = inputs_H2["W"]
	Tw = inputs_H2["Tw"]



	## Fixing hydrogen sector planning variables
	for key in keys(results_H2["H2PLAN"])
	    # JuMP.fix(EP[:key][y],results_H2["H2PLAN"][:key][k])
	    var_size = size(results_H2["H2PLAN"][key])
	    var_dim = length(size(results_H2["H2PLAN"][key]))

	    if var_dim == 1
	        for i = 1:var_size[1]
	            JuMP.setvalue(EP[Symbol(key)][i],results_H2["H2PLAN"][key][i])
	        end
	    end

	    if var_dim == 2
	        for i = 1:var_size[1]
	            for j = 1:var_size[2]
	                JuMP.setvalue(EP[Symbol(key)][i,j],results_H2["H2PLAN"][key][i,j])
	            end
	        end
	    end

	    if var_dim == 3
	        for i = 1:var_size[1]
	            for j = 1:var_size[2]
	                for k = 1:var_size[3]
	                    JuMP.setvalue(EP[Symbol(key)][i,j,k],results_H2["H2PLAN"][key][i,j,k])
	                end
	            end
	        end
	    end
	    # println(results_H2["H2PLAN"][key])
	end

	return EP, dModuleArgs
end

function generate_Power_H2_model(setup::Dict,inputs::Dict,inputs_H2::Dict,OPTIMIZER,modeloutput = nothing)
	println("Generating Model")
	presolver_start_time = time()
	# println("Tock!")

	## Model Definition
	# Define the Energy Portfolio (EP) model
	# Set solver to use Gurobi
	EP=Model(solver=OPTIMIZER)

	power_sector_option = 1
	if power_sector_option == 1

		T = inputs["T"]     # Number of time steps (hours)
		Z = inputs["Z"]     # Number of zones

		dModuleArgs = Dict()
		dModuleArgs["inputs"] = inputs
		dModuleArgs["setup"] = setup
		dModuleArgs["dExpressions"] = Dict()
		dModuleArgs["dObjective"] = Dict()
		dModuleArgs["dPowerBalance"] = Dict()

		# If Unit Commitment is on, we actually modify the values in the generators data
		if (setup["UCommit"]==0)
			inputs["dfGen"][!,:Cap_size].=1
		end

		# Infrastructure
		EP, dModuleArgs =  economic_dispatch(EP, dModuleArgs)
		EP, dModuleArgs = investment(EP, dModuleArgs)
		EP, dModuleArgs = ucommit(EP, dModuleArgs)
		EP, dModuleArgs = reserves(EP, dModuleArgs)
		EP, dModuleArgs = transmission(EP, dModuleArgs)
		# Technologies
		EP, dModuleArgs = dispatchable(EP, dModuleArgs)
		EP, dModuleArgs = nondispatchable(EP, dModuleArgs)
		EP, dModuleArgs = storage(EP, dModuleArgs)
		EP, dModuleArgs = hydro(EP, dModuleArgs)
		EP, dModuleArgs = dr(EP, dModuleArgs)
		EP, dModuleArgs = heat(EP, dModuleArgs)
		EP, dModuleArgs = thermal(EP, dModuleArgs)
		EP, dModuleArgs = hybrid(EP, dModuleArgs)
		# Policies
		EP, dModuleArgs = trading(EP, dModuleArgs)
		EP, dModuleArgs = co2(EP, dModuleArgs)
		EP, dModuleArgs = rps(EP, dModuleArgs)
		EP, dModuleArgs = ces(EP, dModuleArgs)
		EP, dModuleArgs = specificshare(EP, dModuleArgs)

		dExpressions = dModuleArgs["dExpressions"]
		dObjective = dModuleArgs["dObjective"]
		dPowerBalance = dModuleArgs["dPowerBalance"]

		if setup["CO2_price_option"] == 1
			dObjective["eTotalCO2EmissionCost_Power"] = dExpressions["eCO2Emissions_power"] * inputs["CO2_price"]/1e+6
		else
			dObjective["eTotalCO2EmissionCost_Power"] = 0
		end

		eObjective_power = dObjective["eTotalCFix"] + dObjective["eTotalCVar"] + dObjective["eTotalCNSE"] + dObjective["eTotalCStart"] + dObjective["eTotalHeat_Cost"] - dObjective["eTotalHeat_Rev"] + dObjective["eTotalCRsvPen"] + dObjective["eTotalCNetworkExp"] - dObjective["eTotalexportRev"] + dObjective["eTotalImportCost"] + dObjective["eTotalCO2EmissionCost_Power"]

		dExpressions["eObjective"] = eObjective_power

		# dPowerBalance["ePowerBalance"] = dPowerBalance["ePowerBalance"] - inputs["pD"]
		dPowerBalance["ePowerBalance"] = dPowerBalance["ePowerBalance"] - inputs["pD"] - inputs_H2["pD_extra"] - inputs_H2["BEV_demand"]

		dModuleArgs["dPowerBalance_power_sector"] = dPowerBalance["ePowerBalance"]

		dModuleArgs["dExpressions"] = dExpressions
		dModuleArgs["dObjective"] = dObjective
		dModuleArgs["dPowerBalance"] = dPowerBalance
	end

	H2_sector_option = 1
	if H2_sector_option == 1
		T = inputs_H2["T"]     # Number of time steps (hours)
		Z = inputs_H2["Z"]     # Number of zones

		dModuleArgs["inputs_H2"] = inputs_H2
		dModuleArgs["setup"] = setup
		dModuleArgs["dExpressions_H2"] = Dict()
		dModuleArgs["dObjective_H2"] = Dict()
		dModuleArgs["dBalance_H2"] = Dict()

		dExpressions_H2 = dModuleArgs["dExpressions_H2"]
		dObjective_H2 = dModuleArgs["dObjective_H2"]
		dBalance_H2 = dModuleArgs["dBalance_H2"]


		@expression(EP, H2Balance[t=1:T, z=1:Z], 0)
		dBalance_H2["H2Balance"] = H2Balance
		@expression(EP, H2PowerBalance[t=1:T, z=1:Z], 0)
		dBalance_H2["H2PowerBalance"] = H2PowerBalance
		@expression(EP, H2_Carbon_Emission, 0)
		dBalance_H2["H2_Carbon_Emission"] = H2_Carbon_Emission

		discount_rate = inputs_H2["discount_rate"]
		life = inputs_H2["life"]
		@expression(EP, discount_factor, discount_rate/(1-(1+discount_rate)^(-life)))
		dExpressions_H2["discount_factor"] = discount_factor

		dModuleArgs["dExpressions_H2"] = dExpressions_H2
		dModuleArgs["dObjective_H2"] = dObjective_H2
		dModuleArgs["dBalance_H2"] = dBalance_H2

		EP, dModuleArgs =  H2_Gen(EP, dModuleArgs)
        EP, dModuleArgs =  H2_PowerGen(EP, dModuleArgs)
		EP, dModuleArgs =  H2_Demand(EP, dModuleArgs)
		EP, dModuleArgs =  H2_Storage(EP, dModuleArgs)
		EP, dModuleArgs =  H2Pipeline(EP, dModuleArgs)
		EP, dModuleArgs =  H2Truck(EP, dModuleArgs)


		if setup["battery_truck_option"] == 1
			EP, dModuleArgs =  Battery_Truck(EP, dModuleArgs)
		end

		dExpressions_H2 = dModuleArgs["dExpressions_H2"]
		dObjective_H2 = dModuleArgs["dObjective_H2"]
		dBalance_H2 = dModuleArgs["dBalance_H2"]
		dPowerBalance = dModuleArgs["dPowerBalance"]
		inputs_H2 = dModuleArgs["inputs_H2"]

		weights = inputs_H2["weights"]
		@expression(EP, H2_Power_Comsumption[t=1:T, z=1:Z], dBalance_H2["H2PowerBalance"][t,z])
		dExpressions_H2["H2_Power_Comsumption"] = H2_Power_Comsumption

		dPowerBalance["ePowerBalance"] = dPowerBalance["ePowerBalance"] - dExpressions_H2["H2_Power_Comsumption"]/1e+3

		# electricity_price = inputs_H2["electricity_price"]
		# @expression(HY, Cost_electricity, sum(weights[t] * electricity_price[t,z] * dBalance_H2["H2PowerBalance"][t,z] for z=1:Z, t = 1:T))
		# dExpressions_H2["Cost_electricity"] = Cost_electricity

		# dPowerBalance["ePowerBalance"][t,z] - dExpressions_H2["H2_Power_Comsumption"][t,z] == 0)

		@constraint(EP, cH2Balance[t=1:T, z=1:Z], dBalance_H2["H2Balance"][t,z] == 0)

		Total_H2_Demand = dExpressions_H2["Total_H2_Demand"]
		H2_emission_rate = inputs_H2["H2_emission_rate"]

		dExpressions_H2["H2_Carbon_Emission"] = dBalance_H2["H2_Carbon_Emission"]
		dExpressions_H2["Total_CO2_Emission_Cost_H2"] = dBalance_H2["H2_Carbon_Emission"] * inputs_H2["CO2_price"]

		if setup["CO2_price_option"] == 1
			dObjective_H2["eTotalCO2EmissionCost_H2"] = dExpressions_H2["Total_CO2_Emission_Cost_H2"]
		else
			dObjective_H2["Total_CO2_Emission_Cost_H2"] = 0
			@constraint(EP, cCarbonEmission, dBalance_H2["H2_Carbon_Emission"] <= Total_H2_Demand * H2_emission_rate)
		end

		global eObjective_H2 = 0
		for key in keys(dObjective_H2)
		    global eObjective_H2
		    eObjective_H2 = eObjective_H2 + dObjective_H2[key]
		end
		dExpressions_H2["eObjective"] = eObjective_H2

		dExpressions_H2["Total_Compression_Cost"] = dExpressions_H2["CAPEX_Compression_Truck"] + dExpressions_H2["CAPEX_Compression_Pipe"] + dExpressions_H2["CAPEX_Compression_Stor"]
		dExpressions_H2["Total_H2_Cost"] = dExpressions_H2["eObjective"]

	end

	@objective(EP,Min, dExpressions["eObjective"] + dExpressions_H2["eObjective"]/1e+6)
	@constraint(EP, cPowerBalance[t=1:T, z=1:Z], dPowerBalance["ePowerBalance"][t,z] == 0)

	dModuleArgs["dExpressions_H2"] = dExpressions_H2
	dModuleArgs["dObjective_H2"] = dObjective_H2
	dModuleArgs["dBalance_H2"] = dBalance_H2
	dModuleArgs["dExpressions"] = dExpressions
	dModuleArgs["dObjective"] = dObjective
	dModuleArgs["dPowerBalance"] = dPowerBalance
	return EP, dModuleArgs
end

function generate_Power_model(setup::Dict,inputs::Dict,inputs_H2::Dict,OPTIMIZER,modeloutput = nothing)
	println("Generating Model")
	presolver_start_time = time()
	# println("Tock!")

	## Model Definition
	# Define the Energy Portfolio (EP) model
	# Set solver to use Gurobi
	EP=Model(solver=OPTIMIZER)

	power_sector_option = 1
	if power_sector_option == 1

		T = inputs["T"]     # Number of time steps (hours)
		Z = inputs["Z"]     # Number of zones

		dModuleArgs = Dict()
		dModuleArgs["inputs"] = inputs
		dModuleArgs["setup"] = setup
		dModuleArgs["dExpressions"] = Dict()
		dModuleArgs["dObjective"] = Dict()
		dModuleArgs["dPowerBalance"] = Dict()

		# If Unit Commitment is on, we actually modify the values in the generators data
		if (setup["UCommit"]==0)
			inputs["dfGen"][!,:Cap_size].=1
		end

		# Infrastructure
		EP, dModuleArgs =  economic_dispatch(EP, dModuleArgs)
		EP, dModuleArgs = investment(EP, dModuleArgs)
		EP, dModuleArgs = ucommit(EP, dModuleArgs)
		EP, dModuleArgs = reserves(EP, dModuleArgs)
		EP, dModuleArgs = transmission(EP, dModuleArgs)
		# Technologies
		EP, dModuleArgs = dispatchable(EP, dModuleArgs)
		EP, dModuleArgs = nondispatchable(EP, dModuleArgs)
		EP, dModuleArgs = storage(EP, dModuleArgs)
		EP, dModuleArgs = hydro(EP, dModuleArgs)
		EP, dModuleArgs = dr(EP, dModuleArgs)
		EP, dModuleArgs = heat(EP, dModuleArgs)
		EP, dModuleArgs = thermal(EP, dModuleArgs)
		EP, dModuleArgs = hybrid(EP, dModuleArgs)
		# Policies
		EP, dModuleArgs = trading(EP, dModuleArgs)
		EP, dModuleArgs = co2(EP, dModuleArgs)
		EP, dModuleArgs = rps(EP, dModuleArgs)
		EP, dModuleArgs = ces(EP, dModuleArgs)
		EP, dModuleArgs = specificshare(EP, dModuleArgs)

		dExpressions = dModuleArgs["dExpressions"]
		dObjective = dModuleArgs["dObjective"]
		dPowerBalance = dModuleArgs["dPowerBalance"]

		if setup["CO2_price_option"] == 1
			dObjective["eTotalCO2EmissionCost_Power"] = dExpressions["eCO2Emissions_power"] * inputs["CO2_price"]/1e+6
		else
			dObjective["eTotalCO2EmissionCost_Power"] = 0
		end

		eObjective_power = dObjective["eTotalCFix"] + dObjective["eTotalCVar"] + dObjective["eTotalCNSE"] + dObjective["eTotalCStart"] + dObjective["eTotalHeat_Cost"] - dObjective["eTotalHeat_Rev"] + dObjective["eTotalCRsvPen"] + dObjective["eTotalCNetworkExp"] - dObjective["eTotalexportRev"] + dObjective["eTotalImportCost"] + dObjective["eTotalCO2EmissionCost_Power"]

		dExpressions["eObjective"] = eObjective_power

		# dPowerBalance["ePowerBalance"] = dPowerBalance["ePowerBalance"] - inputs["pD"]
		dPowerBalance["ePowerBalance"] = dPowerBalance["ePowerBalance"] - inputs["pD"] - inputs_H2["pD_extra"] - inputs_H2["BEV_demand"]

		dModuleArgs["dPowerBalance_power_sector"] = dPowerBalance["ePowerBalance"]

		dModuleArgs["dExpressions"] = dExpressions
		dModuleArgs["dObjective"] = dObjective
		dModuleArgs["dPowerBalance"] = dPowerBalance
	end


	@objective(EP,Min, dExpressions["eObjective"]/1e+6)
	@constraint(EP, cPowerBalance[t=1:T, z=1:Z], dPowerBalance["ePowerBalance"][t,z] == 0)

	dModuleArgs["dExpressions"] = dExpressions
	dModuleArgs["dObjective"] = dObjective
	dModuleArgs["dPowerBalance"] = dPowerBalance
	return EP, dModuleArgs
end

function solve_model(EP::Model, dModuleArgs::Dict)
	solver_start_time = time()
	## Solve Model
	status = solve(EP)
	## Record solver time
	solver_time = time() - solver_start_time

	# dfStatus = DataFrame(Status = status, Solve = solver_time,
	# Objval = getobjectivevalue(EP), Objbound= getobjbound(EP),FinalMIPGap =(getobjectivevalue(EP) -getobjbound(EP))/getobjectivevalue(EP) )
	dfStatus = DataFrame(Status = status, Solve = solver_time,
	Objval = getobjectivevalue(EP))

	dModuleArgs["dfStatus"] = dfStatus
	return EP, dModuleArgs
end

function generate_Power_results(EP::Model, setup::Dict, inputs::Dict, dModuleArgs::Dict)
	expressions = dModuleArgs["dExpressions"]
	dObjective = dModuleArgs["dObjective"]

	## Extract data frames from input dictionary
	dfGen = inputs["dfGen"]

	G = inputs["G"]     # Number of resources (generators, storage, DR, and DERs)
	T = inputs["T"]     # Number of time steps (hours)
	Z = inputs["Z"]     # Number of zones
	L = inputs["L"]     # Number of transmission lines
	W = inputs["W"]     # Number of subperiods
    SEG = inputs["SEG"] # Number of load curtailment segments

	DISTR_ZONES = findall(inputs["pDistrZones"].==1)	# Set of distribution network zones

	## Extract and store results...
	# Note: Gurobi excludes constants from solver reported objective function value - MIPGap calculated may be erroneous
	dfStatus = dModuleArgs["dfStatus"]
	status = dfStatus[!,:Status][1]

	## Check if solved sucessfully - time out is included
	# if status!=:Optimal
	# 	if status!=:UserLimit # Model failed to solve, so record solver status and exi
    #         if status!=:Suboptimal
	# 			results = Dict([("STATUS", dfStatus)])
	# 			return results
	# 		end
	# 		# Model reached timelimit but failed to find a feasible solution
	# 	elseif isnan(getobjectivevalue(EP))==true
	# 			# Model failed to solve, so record solver status and exit
	# 		results = Dict([("STATUS", dfStatus)])
	# 		return results
	# 	end
	# end

	if setup["CO2_price_option"] == 1
		eTotalCO2EmissionCost_Power = getvalue(dObjective["eTotalCO2EmissionCost_Power"])
	else
		eTotalCO2EmissionCost_Power = dObjective["eTotalCO2EmissionCost_Power"]
	end


	## Cost results
	if setup["HeatMarket"]==1
		dfCost = DataFrame(Costs = ["cTotal", "cFix", "cVar","cCO2Emission", "cNSE", "cStart", "cUnmetRsv", "cNetworkExp", "cHeatRev","cHeatCost"])
		dfCost[!,Symbol("Total")] = [getvalue(expressions["eObjective"]), getvalue(dObjective["eTotalCFix"]), getvalue(dObjective["eTotalCVar"]),
		eTotalCO2EmissionCost_Power, getvalue(dObjective["eTotalCNSE"]),
								 0, 0, 0, (if dObjective["eTotalHeat_Rev"] == 0 0 else -getvalue(dObjective["eTotalHeat_Rev"]) end), (if dObjective["eTotalHeat_Cost"] == 0 0 else -getvalue(dObjective["eTotalHeat_Cost"]) end)]
	else
		dfCost = DataFrame(Costs = ["cTotal", "cFix", "cVar","cCO2Emission", "cNSE", "cStart", "cUnmetRsv", "cNetworkExp"])
		dfCost[!,Symbol("Total")] = [getvalue(expressions["eObjective"]), getvalue(dObjective["eTotalCFix"]), getvalue(dObjective["eTotalCVar"]),
		eTotalCO2EmissionCost_Power, getvalue(dObjective["eTotalCNSE"]), 0, 0, 0]
	end

	if setup["UCommit"]>=1
		dfCost[!,2][6] = getvalue(dObjective["eTotalCStart"])
	end

	if setup["Reserves"]==1
		dfCost[!,2][7] = getvalue(dObjective["eTotalCRsvPen"])
	end

	if setup["NetworkExpansion"]==1
		dfCost[!,2][8] = getvalue(dObjective["eTotalCNetworkExp"])
	end

	for z in 1:Z
		tempCTotal = 0
		tempCFix = 0
		tempCVar = 0
		tempCStart = 0
		for y in dfGen[dfGen[!,:zone].==z,:][!,:R_ID]
			tempCFix = tempCFix + getvalue(expressions["eCFix"])[y]
			tempCVar = tempCVar + sum(getvalue(expressions["eCVar_in"])[y,:]) + sum(getvalue(expressions["eCVar_out"])[y,:])
			if setup["UCommit"]>=1
				tempCTotal = tempCTotal + getvalue(expressions["eCFix"])[y] + sum(getvalue(expressions["eCVar_in"])[y,:]) + sum(getvalue(expressions["eCVar_out"])[y,:]) + sum(getvalue(expressions["eCStart"])[y,:])
				tempCStart = tempCStart + sum(getvalue(expressions["eCStart"])[y,:])
			else
				tempCTotal = tempCTotal + getvalue(expressions["eCFix"])[y] + sum(getvalue(expressions["eCVar_in"])[y,:]) + sum(getvalue(expressions["eCVar_out"])[y,:])
			end
		end
		if setup["HeatMarket"]==1
			dfCost[!,Symbol("Zone$z")] = [tempCTotal, tempCFix, tempCVar,
			sum(inputs["omega"].*getvalue(dModuleArgs["dExpressions"]["eEmissionsByZone"])[z,:])*inputs["CO2_price"]/1e6,
			sum(getvalue(expressions["eCNSE"])[:,:,z]), tempCStart, "-", "-", "-","-"]
		else
			dfCost[!,Symbol("Zone$z")] = [tempCTotal, tempCFix, tempCVar,
			sum(inputs["omega"].*getvalue(dModuleArgs["dExpressions"]["eEmissionsByZone"])[z,:])*inputs["CO2_price"]/1e6,
			sum(getvalue(expressions["eCNSE"])[:,:,z]), tempCStart, "-", "-"]
		end
	end

	## Extract decision variables
	# Capacity decisions
	capdischarge = zeros(size(inputs["RESOURCES"]))
	capcharge = zeros(size(inputs["RESOURCES"]))
	retcapcharge = zeros(size(inputs["RESOURCES"]))
	aux_i = (dfGen[(dfGen[!,:STOR].==3),:][!,:R_ID])
	for i in aux_i
		capcharge[i] = getvalue(EP[:vCAPCHARGE][i])
		retcapcharge[i] = getvalue(EP[:vRETCAPCHARGE][i])
	end
	if setup["HeatMarket"]==1
		aux_i = (dfGen[(dfGen[!,:HEAT].==2),:][!,:R_ID])
		for i in aux_i
			capdischarge[i] = getvalue(EP[:vCAPDISCHARGE][i])
			capcharge[i] = getvalue(EP[:vCAPCHARGE][i])
		end
	end

	capenergy= zeros(size(inputs["RESOURCES"]))
	retcapenergy = zeros(size(inputs["RESOURCES"]))
	aux_i = (dfGen[(dfGen[!,:STOR].>=2),:][!,:R_ID])
	for i in aux_i
		capenergy[i] = getvalue(EP[:vCAPSTORAGE][i])
		retcapenergy[i] = getvalue(EP[:vRETCAPSTORAGE][i])
	end

	if (setup["UCommit"]>=1)
		dfCap = DataFrame(
			Resource = inputs["RESOURCES"], Zone = dfGen[!,:zone],
			StartCap = dfGen[!,:Existing_Cap_MW],
			RetCap = getvalue(EP[:vRETCAP][:]).*dfGen[!,:Cap_size],
			NewCap = getvalue(EP[:vCAP][:]).*dfGen[!,:Cap_size],
			EndCap = dfGen[!,:Existing_Cap_MW]+getvalue(EP[:vCAP][:]).*dfGen[!,:Cap_size]-getvalue(EP[:vRETCAP][:]).*dfGen[!,:Cap_size],
			StartEnergyCap = dfGen[!,:Existing_Cap_MWh],
			RetEnergyCap = retcapenergy[:],
			NewEnergyCap = capenergy[:],
			EndEnergyCap = dfGen[!,:Existing_Cap_MWh]+capenergy[:]-retcapenergy[:],
			StartChargeCap = dfGen[!,:Existing_Charge_Cap_MW],
			RetChargeCap = retcapcharge[:],
			NewChargeCap = capcharge[:],
			EndChargeCap = dfGen[!,:Existing_Charge_Cap_MW]+capcharge[:]-retcapcharge[:]
		)
	else
		dfCap = DataFrame(
			Resource = inputs["RESOURCES"], Zone = dfGen[!,:zone],
			StartCap = dfGen[!,:Existing_Cap_MW],
			RetCap = getvalue(EP[:vRETCAP][:]),
			NewCap = getvalue(EP[:vCAP][:]),
			EndCap = dfGen[!,:Existing_Cap_MW]+getvalue(EP[:vCAP][:])-getvalue(EP[:vRETCAP][:]),
			StartEnergyCap = dfGen[!,:Existing_Cap_MWh],
			RetEnergyCap = retcapenergy[:],
			NewEnergyCap = capenergy[:],
			EndEnergyCap = dfGen[!,:Existing_Cap_MWh]+capenergy[:]-retcapenergy[:],
			StartChargeCap = dfGen[!,:Existing_Charge_Cap_MW],
			RetChargeCap = retcapcharge[:],
			NewChargeCap = capcharge[:],
			EndChargeCap = dfGen[!,:Existing_Charge_Cap_MW]+capcharge[:]-retcapcharge[:]
		)
	end

	total = DataFrame(
			Resource = "Total", Zone = "n/a",
			StartCap = sum(dfCap[!,:StartCap]), RetCap = sum(dfCap[!,:RetCap]),
			NewCap = sum(dfCap[!,:NewCap]), EndCap = sum(dfCap[!,:EndCap]),
			StartEnergyCap = sum(dfCap[!,:StartEnergyCap]), RetEnergyCap = sum(dfCap[!,:RetEnergyCap]),
			NewEnergyCap = sum(dfCap[!,:NewEnergyCap]), EndEnergyCap = sum(dfCap[!,:EndEnergyCap]),
			StartChargeCap = sum(dfCap[!,:StartChargeCap]), RetChargeCap = sum(dfCap[!,:RetChargeCap]),
			NewChargeCap = sum(dfCap[!,:NewChargeCap]), EndChargeCap = sum(dfCap[!,:EndChargeCap])
		)

	dfCap = vcat(dfCap, total)

	if setup["HeatMarket"]==1
		dfCapheat1 = DataFrame(NewCapDischarge = [capdischarge[:];sum(capdischarge[:])])
		dfCapheat2 = DataFrame(NewCapCharge = [capcharge[:];sum(capcharge[:])])
		dfCap = hcat(dfCap, dfCapheat1,dfCapheat2 )
	end

	# Operational decision variable states

	## Integer Unit Commitment configuration for results
	if (setup["UCommit"]>=1)
		# Commitment state for each resource in each time step
		dfCommit = DataFrame(Resource = inputs["RESOURCES"], Zone = dfGen[!,:zone])
		dfCommit = hcat(dfCommit, convert(DataFrame, getvalue(EP[:vCOMMIT])))
		auxNew_Names=[Symbol("Resource");Symbol("Zone");[Symbol("t$t") for t in 1:T]]
		rename!(dfCommit,auxNew_Names)

		# Startup state for each resource in each time step
		dfStart = DataFrame(Resource = inputs["RESOURCES"], Zone = dfGen[!,:zone], Sum = Array{Union{Missing,Float32}}(undef, G))
		for i in 1:G
			dfStart[!,:Sum][i] = sum(getvalue(EP[:vSTART])[i,:])
		end
		dfStart = hcat(dfStart, convert(DataFrame, getvalue(EP[:vSTART])))
		auxNew_Names=[Symbol("Resource");Symbol("Zone");Symbol("Sum");[Symbol("t$t") for t in 1:T]]
		rename!(dfStart,auxNew_Names)
		total = convert(DataFrame, ["Total" 0 sum(dfStart[!,:Sum]) fill(0.0, (1,T))])
		for t in 1:T
			total[!,t+3] .= sum(dfStart[:,Symbol("t$t")][1:G])
		end
		rename!(total,auxNew_Names)
		dfStart = vcat(dfStart, total)


		# Shutdown state for each resource in each time step
		dfShutdown = DataFrame(Resource = inputs["RESOURCES"], Zone = dfGen[!,:zone], Sum = Array{Union{Missing,Float32}}(undef, G))
		for i in 1:G
			dfShutdown[!,:Sum][i] = sum(getvalue(EP[:vSHUT])[i,:])
		end
		dfShutdown = hcat(dfShutdown, convert(DataFrame, getvalue(EP[:vSHUT])))
		auxNew_Names=[Symbol("Resource");Symbol("Zone");Symbol("Sum");[Symbol("t$t") for t in 1:T]]
		rename!(dfShutdown,auxNew_Names)
		total = convert(DataFrame, ["Total" 0 sum(dfShutdown[!,:Sum]) fill(0.0, (1,T))])
		for t in 1:T
			total[!,t+3] .= sum(dfShutdown[!,Symbol("t$t")][1:G])
		end
		rename!(total,auxNew_Names)
		dfShutdown = vcat(dfShutdown, total)

		if (setup["Reserves"]==1)
			# Regulation up for each resource in each time step
			dfRegUp = DataFrame(Resource = inputs["RESOURCES"], Zone = dfGen[!,:zone], Sum = Array{Union{Missing,Float32}}(undef, G))
			for i in 1:G
				dfRegUp[!,:Sum][i] = sum(getvalue(EP[:vREG_UP])[i,:])
			end
			dfRegUp = hcat(dfRegUp, convert(DataFrame, getvalue(EP[:vREG_UP])))
			auxNew_Names=[Symbol("Resource");Symbol("Zone");Symbol("Sum");[Symbol("t$t") for t in 1:T]]
			rename!(dfRegUp,auxNew_Names)
			total = convert(DataFrame, ["Total" 0 sum(dfRegUp[!,:Sum]) fill(0.0, (1,T))])
			for t in 1:T
				total[!,t+3] .= sum(dfRegUp[!,Symbol("t$t")][1:G])
			end
			rename!(total,auxNew_Names)
			dfRegUp = vcat(dfRegUp, total)

			# Regulation down for each resource in each time step
			dfRegDn = DataFrame(Resource = inputs["RESOURCES"], Zone = dfGen[!,:zone], Sum = Array{Union{Missing,Float32}}(undef, G))
			for i in 1:G
				dfRegDn[!,:Sum][i] = sum(getvalue(EP[:vREG_DN])[i,:])
			end
			dfRegDn = hcat(dfRegDn, convert(DataFrame, getvalue(EP[:vREG_DN])))
			auxNew_Names=[Symbol("Resource");Symbol("Zone");Symbol("Sum");[Symbol("t$t") for t in 1:T]]
			rename!(dfRegDn,auxNew_Names)
			total = convert(DataFrame, ["Total" 0 sum(dfRegDn[!,:Sum]) fill(0.0, (1,T))])
			for t in 1:T
				total[!,t+3] .= sum(dfRegDn[!,Symbol("t$t")][1:G])
			end
			rename!(total,auxNew_Names)
			dfRegDn = vcat(dfRegDn, total)

			#Reserves up for each resource in each time step
			dfRsvUp = DataFrame(Resource = vcat(inputs["RESOURCES"], "unmet"), Zone = vcat(dfGen[!,:zone], "na"), Sum = Array{Union{Missing,Float32}}(undef, G+1))
			# dfRsvUp = DataFrame(Resource = vcat(inputs["RESOURCES"], "unmet"), Zone = vcat(dfGen[!,:zone], "na"), Sum = Array{Union{Missing,Float32}, G+1})

			for i in 1:G
				dfRsvUp[!,:Sum][i] = sum(getvalue(EP[:vRSV_UP])[i,:])
			end
			dfRsvUp[!,:Sum][G+1] = sum(getvalue(EP[:vUNMET_RSV_UP]))
			reg = convert(DataFrame, getvalue(EP[:vRSV_UP]))
			push!(reg, transpose(convert(Array{Union{Missing,Float32}}, getvalue(EP[:vUNMET_RSV_UP]))))
			dfRsvUp = hcat(dfRsvUp, reg)
			auxNew_Names=[Symbol("Resource");Symbol("Zone");Symbol("Sum");[Symbol("t$t") for t in 1:T]]
			rename!(dfRsvUp,auxNew_Names)
			total = convert(DataFrame, ["Total" 0 sum(dfRsvUp[!,:Sum]) fill(0.0, (1,T))])
			for t in 1:T
				total[!,t+3] .= sum(dfRsvUp[!,Symbol("t$t")][1:G])
			end
			rename!(total,auxNew_Names)
			dfRsvUp = vcat(dfRsvUp, total)

			# Reserves down for each resource in each time step
			dfRsvDn = DataFrame(Resource = vcat(inputs["RESOURCES"], "unmet"), Zone = vcat(dfGen[!,:zone], "na"), Sum = Array{Union{Missing,Float32}}(undef, G+1))
			for i in 1:G
				dfRsvDn[!,:Sum][i] = sum(getvalue(EP[:vRSV_DN])[i,:])
			end
			dfRsvDn[!,:Sum][G+1] = sum(getvalue(EP[:vUNMET_RSV_DN]))
			reg = convert(DataFrame, getvalue(EP[:vRSV_DN]))
			push!(reg, transpose(convert(Array{Union{Missing,Float32}}, getvalue(EP[:vUNMET_RSV_DN]))))
			dfRsvDn = hcat(dfRsvDn, reg)
			auxNew_Names=[Symbol("Resource");Symbol("Zone");Symbol("Sum");[Symbol("t$t") for t in 1:T]]
			rename!(dfRsvDn,auxNew_Names)
			total = convert(DataFrame, ["Total" 0 sum(dfRsvDn[!,:Sum]) fill(0.0, (1,T))])
			for t in 1:T
				total[!,t+3] .= sum(dfRsvDn[!,Symbol("t$t")][1:G])
			end
			rename!(total,auxNew_Names)
			dfRsvDn = vcat(dfRsvDn, total)
		end
	end #END unit commitment configuration

	if (setup["HeatMarket"]==1)
		# Power injected by each resource in each time step
		dfPower = DataFrame(Resource = inputs["RESOURCES"], Zone = dfGen[!,:zone], Sum = Array{Union{Missing,Float32}}(undef, G))
		Aux= zeros(Float64,G,T)
		for z in 1:Z
			for y in (dfGen[(dfGen[!,:NACC].==1) .& (dfGen[!,:zone].==z),:][!,:R_ID])
				for x in (dfGen[(dfGen[!,:HEAT].==1) .& (dfGen[!,:zone].==z),:][!,:R_ID])
					Aux[y,:] = dfGen[!,:NACC_Eff][y]*(getvalue(EP[:vHEAT_NG])[y,:] + getvalue(EP[:vHEAT_GENERATION])[x,:])
				end
			end
		end
		for i in 1:G
			dfPower[!,:Sum][i] = sum(getvalue(EP[:vP])[i,:])+sum(Aux[i,:])
		end
		dfPower = hcat(dfPower, convert(DataFrame, getvalue(EP[:vP]) + Aux ) )
		auxNew_Names=[Symbol("Resource");Symbol("Zone");Symbol("Sum");[Symbol("t$t") for t in 1:T]]
		rename!(dfPower,auxNew_Names)
	else
		# Power injected by each resource in each time step
		dfPower = DataFrame(Resource = inputs["RESOURCES"], Zone = dfGen[!,:zone], Sum = Array{Union{Missing,Float32}}(undef, G))
		for i in 1:G
			dfPower[!,:Sum][i] = sum(getvalue(EP[:vP])[i,:])
		end
		dfPower = hcat(dfPower, convert(DataFrame, getvalue(EP[:vP])))
		auxNew_Names=[Symbol("Resource");Symbol("Zone");Symbol("Sum");[Symbol("t$t") for t in 1:T]]
		rename!(dfPower,auxNew_Names)
	end

	total = convert(DataFrame, ["Total" 0 sum(dfPower[!,:Sum]) fill(0.0, (1,T))])
	for t in 1:T
		total[!,t+3] .= sum(dfPower[!,Symbol("t$t")][1:G])
	end
	rename!(total,auxNew_Names)
	dfPower = vcat(dfPower, total)

	# Power withdrawn to charge each resource in each time step
	dfCharge = DataFrame(Resource = inputs["RESOURCES"], Zone = dfGen[!,:zone], Sum = Array{Union{Missing,Float32}}(undef, G))
	for i in 1:G
		dfCharge[!,:Sum][i] = sum(getvalue(EP[:vCHARGE])[i,:])
	end
	dfCharge = hcat(dfCharge, convert(DataFrame, getvalue(EP[:vCHARGE])))
	auxNew_Names=[Symbol("Resource");Symbol("Zone");Symbol("Sum");[Symbol("t$t") for t in 1:T]]
	rename!(dfCharge,auxNew_Names)
	total = convert(DataFrame, ["Total" 0 sum(dfCharge[!,:Sum]) fill(0.0, (1,T))])
	for t in 1:T
		total[!,t+3] .= sum(dfCharge[!,Symbol("t$t")][1:G])
	end
	rename!(total,auxNew_Names)
	dfCharge = vcat(dfCharge, total)


	# Variables related to storage carry over from 1 representative period to the other
	# only valid when operations wrapping =2 and LDS = 1

	if (setup["OperationWrapping"]==2 && setup["LDS"]==1)
		# Initial level of storage in each modeled period
		NPeriods = size(inputs["Period_Map"])[1]
		dfStorageInit = DataFrame(Resource = inputs["RESOURCES"], Zone = dfGen[!,:zone])
		dfStorageInit = hcat(dfStorageInit, convert(DataFrame, getvalue(EP[:vSOCw])))
		auxNew_Names=[Symbol("Resource");Symbol("Zone");[Symbol("n$t") for t in 1:NPeriods]]
		rename!(dfStorageInit,auxNew_Names)

		#Excess inventory of storage period built up during representative period w
		dfdStorage = DataFrame(Resource = inputs["RESOURCES"], Zone = dfGen[!,:zone])
		dfdStorage = hcat(dfdStorage, convert(DataFrame, getvalue(EP[:vdSOC])))
		auxNew_Names=[Symbol("Resource");Symbol("Zone");[Symbol("w$t") for t in 1:W]]
		rename!(dfdStorage,auxNew_Names)
	end

	# Storage level (state of charge) of each resource in each time step
	dfStorage = DataFrame(Resource = inputs["RESOURCES"], Zone = dfGen[!,:zone])
	dfStorage = hcat(dfStorage, convert(DataFrame, getvalue(EP[:vS])))
	auxNew_Names=[Symbol("Resource");Symbol("Zone");[Symbol("t$t") for t in 1:T]]
	rename!(dfStorage,auxNew_Names)

	# Variables related to storage carry over from 1 representative period to the other
	# only valid when operations wrapping =2 and LDS = 1

	if (setup["OperationWrapping"]==2 && setup["LDS"]==1)
		# Initial level of storage in each modeled period
		NPeriods = size(inputs["Period_Map"])[1]
		dfStorageInit = DataFrame(Resource = inputs["RESOURCES"], Zone = dfGen[!,:zone])
		dfStorageInit = hcat(dfStorageInit, convert(DataFrame, getvalue(EP[:vSOCw])))
		auxNew_Names=[Symbol("Resource");Symbol("Zone");[Symbol("n$t") for t in 1:NPeriods]]
		rename!(dfStorageInit,auxNew_Names)

		#Excess inventory of storage period built up during representative period w
		dfdStorage = DataFrame(Resource = inputs["RESOURCES"], Zone = dfGen[!,:zone])
		dfdStorage = hcat(dfdStorage, convert(DataFrame, getvalue(EP[:vdSOC])))
		auxNew_Names=[Symbol("Resource");Symbol("Zone");[Symbol("w$t") for t in 1:W]]
		rename!(dfdStorage,auxNew_Names)
	end


	# Non-served energy/demand curtailment by segment in each time step
	global dfNse = DataFrame()
	for z in 1:Z
		global dfNse
		dfTemp = DataFrame(Segment=zeros(SEG), Zone=zeros(SEG), Sum = Array{Union{Missing,Float32}}(undef, SEG))
		dfTemp[!,:Segment] = (1:SEG)
		dfTemp[!,:Zone] = fill(z,(SEG))
		for i in 1:SEG
			dfTemp[!,:Sum][i] = sum(getvalue(EP[:vNSE])[i,:,z])
		end
		dfTemp = hcat(dfTemp, convert(DataFrame, getvalue(EP[:vNSE])[:,:,z]))
		if z == 1
			dfNse = dfTemp
		else
			dfNse = vcat(dfNse,dfTemp)
		end
	end
	auxNew_Names=[Symbol("Segment");Symbol("Zone");Symbol("Sum");[Symbol("t$t") for t in 1:T]]
	rename!(dfNse,auxNew_Names)
	total = convert(DataFrame, ["Total" 0 sum(dfNse[!,:Sum]) fill(0.0, (1,T))])
	for t in 1:T
		total[!,t+3] .= sum(dfNse[!,Symbol("t$t")][1:Z])
	end
	rename!(total,auxNew_Names)
	dfNse = vcat(dfNse, total)

	# Transmission related values

	# Power flows on transmission lines at each time step
	dfFlow = DataFrame(Line = 1:L, Sum = Array{Union{Missing,Float32}}(undef, L))
	for i in 1:L
		dfFlow[!,:Sum][i] = sum(getvalue(EP[:vFLOW])[i,:])
	end
	dfFlow = hcat(dfFlow, convert(DataFrame, getvalue(EP[:vFLOW])))
	auxNew_Names=[Symbol("Line");Symbol("Sum");[Symbol("t$t") for t in 1:T]]
	rename!(dfFlow,auxNew_Names)
	total = convert(DataFrame, ["Total" sum(dfFlow[!,:Sum]) fill(0.0, (1,T))])
	for t in 1:T
		total[!,t+2] .= sum(dfFlow[!,Symbol("t$t")][1:L])
	end
	rename!(total,auxNew_Names)
	dfFlow = vcat(dfFlow, total)

	# Power losses for transmission between zones at each time step
	dfTLosses = DataFrame(Line = 1:L, Sum = Array{Union{Missing,Float32}}(undef, L))
	for i in 1:L
		dfTLosses[!,:Sum][i] = sum(getvalue(EP[:vTLOSS])[i,:])
	end
	dfTLosses = hcat(dfTLosses, convert(DataFrame, getvalue(EP[:vTLOSS])))
	auxNew_Names=[Symbol("Line");Symbol("Sum");[Symbol("t$t") for t in 1:T]]
	rename!(dfTLosses,auxNew_Names)
	total = convert(DataFrame, ["Total" sum(dfTLosses[!,:Sum]) fill(0.0, (1,T))])
	for t in 1:T
		total[!,t+2] .= sum(dfTLosses[!,Symbol("t$t")][1:L])
	end
	rename!(total,auxNew_Names)
	dfTLosses = vcat(dfTLosses, total)

	# Power losses for within-zone flows for each zone in each time step
	dfWLosses = DataFrame(Zone = 1:Z, Sum = Array{Union{Missing,Float32}}(undef, Z))
	for i in 1:Z
		dfWLosses[!,:Sum][i] = sum(getvalue(EP[:vDLOSS])[i,:])
	end
	dfWLosses = hcat(dfWLosses, convert(DataFrame, getvalue(EP[:vDLOSS])))
	auxNew_Names=[Symbol("Zone");Symbol("Sum");[Symbol("t$t") for t in 1:T]]
	rename!(dfWLosses,auxNew_Names)
	total = convert(DataFrame, ["Total" sum(dfWLosses[!,:Sum]) fill(0.0, (1,T))])
	for t in 1:T
		total[!,t+2] .= sum(dfWLosses[!,Symbol("t$t")][1:Z])
	end
	rename!(total,auxNew_Names)
	dfWLosses = vcat(dfWLosses, total)

	# Total power injection by zone
	dfInject = DataFrame(Zone = Array{Union{Missing,Float64}}(undef, length(DISTR_ZONES)*2), Voltage_Level = Array{Union{Missing,String}}(undef, 2*length(DISTR_ZONES)),Sum = Array{Union{Missing,Float64}}(undef, 2*length(DISTR_ZONES)))
	if length(DISTR_ZONES) != 0
		aux_index_1 = 0
		for i in DISTR_ZONES
			aux_index_1 = aux_index_1+1
			dfInject[!,:Zone][aux_index_1] =  i
			dfInject[!,:Voltage_Level][aux_index_1] =  "LV"
			dfInject[!,:Sum][aux_index_1] = sum(getvalue(EP[:eInjectLV][i,:]))
		end
		aux_index_1 = length(DISTR_ZONES)
		for i in DISTR_ZONES
			aux_index_1 = aux_index_1+1
			dfInject[!,:Zone][aux_index_1] =  i
			dfInject[!,:Voltage_Level][aux_index_1] =  "MV"
			dfInject[!,:Sum][aux_index_1] = sum(getvalue(EP[:eInjectMV][i,:]))
		end
		dfInject_aux1 = DataFrame()
		aux_index_1 = 0
		for i in DISTR_ZONES
			dfTemp = convert(DataFrame, transpose(getvalue(EP[:eInjectLV][i,:]))  )
			if aux_index_1 == 0
				dfInject_aux1 = dfTemp
			else
				dfInject_aux1 = vcat(dfInject_aux1,dfTemp)
			end
			aux_index_1 = aux_index_1+1
		end
		dfInject_aux2 = DataFrame()
		aux_index_1 = 0
		for i in DISTR_ZONES
			dfTemp = convert(DataFrame, transpose(getvalue(EP[:eInjectMV][i,:]))  )
			if aux_index_1 == 0
				dfInject_aux2 = dfTemp
			else
				dfInject_aux2 = vcat(dfInject_aux2,dfTemp)
			end
			aux_index_1 = aux_index_1+1
		end
		dfInject_aux = vcat(dfInject_aux1 ,dfInject_aux2 )

		dfInject = hcat(dfInject, dfInject_aux)
		auxNew_Names=[Symbol("Zone");Symbol("Voltage_Level");Symbol("Sum");[Symbol("t$t") for t in 1:T]]
		rename!(dfInject,auxNew_Names)
	end

	# Total power withdrawal by zone
	dfWithdraw = DataFrame(Zone = Array{Union{Missing,Float64}}(undef, length(DISTR_ZONES)*2), Voltage_Level = Array{Union{Missing,String}}(undef, 2*length(DISTR_ZONES)),Sum = Array{Union{Missing,Float64}}(undef, 2*length(DISTR_ZONES)))
	if length(DISTR_ZONES) != 0
		aux_index_1 = 0
		for i in DISTR_ZONES
			aux_index_1 = aux_index_1+1
			dfWithdraw[!,:Zone][aux_index_1] =  i
			dfWithdraw[!,:Voltage_Level][aux_index_1] =  "LV"
			dfWithdraw[!,:Sum][aux_index_1] = sum(getvalue(EP[:eWithdrawLV][i,:]))
		end
		aux_index_1 = length(DISTR_ZONES)
		for i in DISTR_ZONES
			aux_index_1 = aux_index_1+1
			dfWithdraw[!,:Zone][aux_index_1] =  i
			dfWithdraw[!,:Voltage_Level][aux_index_1] =  "MV"
			dfWithdraw[!,:Sum][aux_index_1] = sum(getvalue(EP[:eWithdrawMV][i,:]))
		end
		dfInject_aux1 = DataFrame()
		aux_index_1 = 0
		for i in DISTR_ZONES
			dfTemp = convert(DataFrame, transpose(getvalue(EP[:eWithdrawLV][i,:]))  )
			if aux_index_1 == 0
				dfInject_aux1 = dfTemp
			else
				dfInject_aux1 = vcat(dfInject_aux1,dfTemp)
			end
			aux_index_1 = aux_index_1+1
		end
		dfInject_aux2 = DataFrame()
		aux_index_1 = 0
		for i in DISTR_ZONES
			dfTemp = convert(DataFrame, transpose(getvalue(EP[:eWithdrawMV][i,:]))  )
			if aux_index_1 == 0
				dfInject_aux2 = dfTemp
			else
				dfInject_aux2 = vcat(dfInject_aux2,dfTemp)
			end
			aux_index_1 = aux_index_1+1
		end
		dfInject_aux = vcat(dfInject_aux1 ,dfInject_aux2 )

		dfWithdraw = hcat(dfWithdraw, dfInject_aux)
		auxNew_Names=[Symbol("Zone");Symbol("Voltage_Level");Symbol("Sum");[Symbol("t$t") for t in 1:T]]
		rename!(dfWithdraw,auxNew_Names)
	end

	## CHP system configuration for results
	if (setup["HeatMarket"]==1)
		# Heat related values
		# Heat sold on the heat market each time step
		aux_y = (dfGen[(dfGen[!,:HEAT].==1),:][!,:R_ID])
		dfHeat_market = DataFrame(Resource = inputs["RESOURCES"][aux_y], Zone = dfGen[!,:zone][aux_y], Sum = Array{Union{Missing,Float32}}(undef, size(aux_y)))
		aux_index = 0
		for i in (dfGen[(dfGen[!,:HEAT].==1),:][!,:R_ID])
			aux_index = aux_index +1
			dfHeat_market[!,:Sum][aux_index] = sum(getvalue(EP[:vHEAT_MARKET])[i,:])
		end
		if size(aux_y)[1] >= 1
			dfHeat_market = hcat(dfHeat_market, convert(DataFrame, getvalue(EP[:vHEAT_MARKET])))
			auxNew_Names=[Symbol("Resource");Symbol("Zone");Symbol("Sum");[Symbol("t$t") for t in 1:T]]
			rename!(dfHeat_market,auxNew_Names)
			total = convert(DataFrame, ["Total" 0 sum(dfHeat_market[!,:Sum]) fill(0.0, (1,T))])
			for t in 1:T
				total[t+3] = sum(dfHeat_market[Symbol("t$t")][1:size(aux_y)[1]])
			end
			rename!(total,auxNew_Names)
			dfHeat_market = vcat(dfHeat_market, total)
		end

		# Heat injected to NACC  market each time step
		aux_y = (dfGen[(dfGen[!,:HEAT].==1),:][!,:R_ID])
		dfHeat_gen = DataFrame(Resource = inputs["RESOURCES"][aux_y], Zone = dfGen[!,:zone][aux_y], Sum = Array{Union{Missing,Float32}}(undef, size(aux_y)))
		aux_index = 0
		for i in (dfGen[(dfGen[!,:HEAT].==1),:][!,:R_ID])
			aux_index = aux_index +1
			dfHeat_gen[!,:Sum][aux_index] = sum(getvalue(EP[:vHEAT_GENERATION])[i,:])
		end
		if size(aux_y)[1] >= 1
			dfHeat_gen = hcat(dfHeat_gen, convert(DataFrame, getvalue(EP[:vHEAT_GENERATION])))
			auxNew_Names=[Symbol("Resource");Symbol("Zone");Symbol("Sum");[Symbol("t$t") for t in 1:T]]
			rename!(dfHeat_gen,auxNew_Names)
			total = convert(DataFrame, ["Total" 0 sum(dfHeat_gen[!,:Sum]) fill(0.0, (1,T))])
			for t in 1:T
				total[t+3] = sum(dfHeat_gen[Symbol("t$t")][1:size(aux_y)[1]])
			end
			rename!(total,auxNew_Names)
			dfHeat_gen = vcat(dfHeat_gen, total)
		end

		# Heat injected to NACC  from NG each time step
		aux_y = (dfGen[(dfGen[!,:NACC].==1),:][!,:R_ID])
		dfHeat_NG = DataFrame(Resource = inputs["RESOURCES"][aux_y], Zone = dfGen[!,:zone][aux_y], Sum = Array{Union{Missing,Float32}}(undef, size(aux_y)))
		aux_index = 0
		for i in (dfGen[(dfGen[!,:NACC].==1),:][!,:R_ID])
			aux_index = aux_index +1
			dfHeat_NG[!,:Sum][aux_index] = sum(getvalue(EP[:vHEAT_NG])[i,:])
		end
		if size(aux_y)[1] >= 1
			dfHeat_NG = hcat(dfHeat_NG, convert(DataFrame, getvalue(EP[:vHEAT_NG])))
			auxNew_Names=[Symbol("Resource");Symbol("Zone");Symbol("Sum");[Symbol("t$t") for t in 1:T]]
			rename!(dfHeat_NG,auxNew_Names)
			total = convert(DataFrame, ["Total" 0 sum(dfHeat_NG[!,:Sum]) fill(0.0, (1,T))])
			for t in 1:T
				total[t+3] = sum(dfHeat_NG[Symbol("t$t")][1:size(aux_y)[1]])
			end
			rename!(total,auxNew_Names)
			dfHeat_NG = vcat(dfHeat_NG, total)
		end

		# power going from heat source to heat storage
		aux_y = (dfGen[(dfGen[!,:NACC].==2),:][!,:R_ID])
		dfHeat_to_Storage = DataFrame(Resource = inputs["RESOURCES"][aux_y], Zone = dfGen[!,:zone][aux_y], Sum = Array{Union{Missing,Float32}}(undef, size(aux_y)))
		aux_index = 0
		for i in (dfGen[(dfGen[!,:NACC].==2),:][!,:R_ID])
			aux_index = aux_index +1
			dfHeat_to_Storage[!,:Sum][aux_index] = sum(getvalue(EP[:vH1])[i,:])
		end
		if size(aux_y)[1] >= 1
			aux_array = zeros(size(aux_y)[1],T)
			for i in 1:size(aux_y)[1]
				aux_array[i,:]=getvalue(EP[:vH1])[aux_y[i],:]
			end
			dfHeat_to_Storage = hcat(dfHeat_to_Storage, convert(DataFrame, aux_array))
			auxNew_Names=[Symbol("Resource");Symbol("Zone");Symbol("Sum");[Symbol("t$t") for t in 1:T]]
			rename!(dfHeat_to_Storage,auxNew_Names)
			total = convert(DataFrame, ["Total" 0 sum(dfHeat_to_Storage[!,:Sum]) fill(0.0, (1,T))])
			for t in 1:T
				total[t+3] = sum(dfHeat_to_Storage[Symbol("t$t")][1:size(aux_y)[1]])
			end
			rename!(total,auxNew_Names)
			dfHeat_to_Storage = vcat(dfHeat_to_Storage, total)
		end

		# power going from heat storage to generation
		aux_y = (dfGen[(dfGen[!,:HEAT].==2),:][!,:R_ID])
		dfHeat_to_turbine= DataFrame(Resource = inputs["RESOURCES"][aux_y], Zone = dfGen[!,:zone][aux_y], Sum = Array{Union{Missing,Float32}}(undef, size(aux_y)))
		aux_index = 0
		for i in (dfGen[(dfGen[!,:HEAT].==2),:][!,:R_ID])
			aux_index = aux_index +1
			dfHeat_to_turbine[!,:Sum][aux_index] = sum(getvalue(EP[:vH2])[i,:])
		end
		if size(aux_y)[1] >= 1
			aux_array = zeros(size(aux_y)[1],T)
			for i in 1:size(aux_y)[1]
				aux_array[i,:]=getvalue(EP[:vH2])[aux_y[i],:]
			end
			dfHeat_to_turbine = hcat(dfHeat_to_turbine, convert(DataFrame, aux_array))
			auxNew_Names=[Symbol("Resource");Symbol("Zone");Symbol("Sum");[Symbol("t$t") for t in 1:T]]
			rename!(dfHeat_to_turbine,auxNew_Names)
			total = convert(DataFrame, ["Total" 0 sum(dfHeat_to_turbine[!,:Sum]) fill(0.0, (1,T))])
			for t in 1:T
				total[t+3] = sum(dfHeat_to_turbine[Symbol("t$t")][1:size(aux_y)[1]])
			end
			rename!(total,auxNew_Names)
			dfHeat_to_turbine = vcat(dfHeat_to_turbine, total)
		end

		# power going reactor or generator to generation
		aux_y = (dfGen[(dfGen[!,:NACC].==2),:][!,:R_ID])
		dfcoreheat_to_turbine= DataFrame(Resource = inputs["RESOURCES"][aux_y], Zone = dfGen[!,:zone][aux_y], Sum = Array{Union{Missing,Float32}}(undef, size(aux_y)))
		aux_index = 0
		for i in (dfGen[(dfGen[!,:NACC].==2),:][!,:R_ID])
			aux_index = aux_index +1
			dfcoreheat_to_turbine[!,:Sum][aux_index] = sum(getvalue(EP[:vP2P])[i,:])
		end
		if size(aux_y)[1] >= 1
			aux_array = zeros(size(aux_y)[1],T)
			for i in 1:size(aux_y)[1]
				aux_array[i,:]=getvalue(EP[:vP2P])[aux_y[i],:]
			end
			dfcoreheat_to_turbine = hcat(dfcoreheat_to_turbine, convert(DataFrame, aux_array))
			auxNew_Names=[Symbol("Resource");Symbol("Zone");Symbol("Sum");[Symbol("t$t") for t in 1:T]]
			rename!(dfcoreheat_to_turbine,auxNew_Names)
			total = convert(DataFrame, ["Total" 0 sum(dfcoreheat_to_turbine[!,:Sum]) fill(0.0, (1,T))])
			for t in 1:T
				total[t+3] = sum(dfcoreheat_to_turbine[Symbol("t$t")][1:size(aux_y)[1]])
			end
			rename!(total,auxNew_Names)
			dfcoreheat_to_turbine = vcat(dfcoreheat_to_turbine, total)
		end
	end #END CHP system configuration

	if (setup["NetworkExpansion"]==1)
		# Transmission network reinforcements
		dfTransCap = DataFrame(
		Line = 1:L, New_Trans_Capacity = convert(Array{Union{Missing,Float32}}, getvalue(EP[:vNEW_TRANS_CAP])),
		Cost_Trans_Capacity = convert(Array{Union{Missing,Float32}}, getvalue(EP[:vNEW_TRANS_CAP]).*inputs["pC_Line_Reinforcement"]),
		)
		# Distribition network reinforcements
		dfDistrCap = DataFrame(
		Zone = DISTR_ZONES,   Peak_Withdraw = Array{Union{Missing,Float32}}(undef, length(DISTR_ZONES)),
		Existing_Withdraw_Capacity = Array{Union{Missing,Float32}}(undef, length(DISTR_ZONES)),
		New_Withdraw_Capacity = Array{Union{Missing,Float32}}(undef, length(DISTR_ZONES)),
		Cost_Withdraw_Capacity = Array{Union{Missing,Float32}}(undef, length(DISTR_ZONES)),
		Withdraw_Margin_Gained =  Array{Union{Missing,Float32}}(undef, length(DISTR_ZONES)),
		Total_Withdraw_Margin = Array{Union{Missing,Float32}}(undef, length(DISTR_ZONES)),
		Peak_Inject =  Array{Union{Missing,Float32}}(undef, length(DISTR_ZONES)),
		Existing_Inject_Capacity =  Array{Union{Missing,Float32}}(undef, length(DISTR_ZONES)),
		New_Inject_Capacity = Array{Union{Missing,Float32}}(undef, length(DISTR_ZONES)),
		Cost_Inject_Capacity =  Array{Union{Missing,Float32}}(undef, length(DISTR_ZONES)),
		Total_Inject_Margin = Array{Union{Missing,Float32}}(undef, length(DISTR_ZONES)),
		)
		aux_index=0
		for z in DISTR_ZONES
			aux_index=aux_index+1
			dfDistrCap[!,:Peak_Withdraw][aux_index]=maximum(getvalue(expressions["ePeakWithdraw"][z,:]))
			dfDistrCap[!,:Existing_Withdraw_Capacity][aux_index]=inputs["pMax_Withdraw"][z]
			dfDistrCap[!,:New_Withdraw_Capacity][aux_index]=getvalue(EP[:vNEW_WITHDRAW_CAP][z])
			dfDistrCap[!,:Cost_Withdraw_Capacity][aux_index]=getvalue(EP[:vNEW_WITHDRAW_CAP][z])*(inputs["pC_MV_Reinforcement"][z]*inputs["pDistr_Share_in_MV"][z] + inputs["pC_LV_Reinforcement"][z]*inputs["pDistr_Share_in_LV"][z])
			dfDistrCap[!,:Withdraw_Margin_Gained][aux_index]=maximum(getvalue(EP[:vWITHDRAW_MARGIN_GAINED][z,:]))
			dfDistrCap[!,:Peak_Inject][aux_index]=maximum(getvalue(expressions["eInjectLV"][z,:])+getvalue(expressions["eInjectMV"][z,:]))
			dfDistrCap[!,:Existing_Inject_Capacity][aux_index]=inputs["pMax_Inject"][z]
			dfDistrCap[!,:New_Inject_Capacity][aux_index]=getvalue(EP[:vNEW_INJECT_CAP][z])
			dfDistrCap[!,:Cost_Inject_Capacity][aux_index]=getvalue(EP[:vNEW_INJECT_CAP][z])*(inputs["pC_MV_Reinforcement"][z]*inputs["pDistr_Share_in_MV"][z] + inputs["pC_LV_Reinforcement"][z]*inputs["pDistr_Share_in_LV"][z])
		end
		dfDistrCap[!,:Total_Withdraw_Margin] = dfDistrCap[!,:Existing_Withdraw_Capacity]+dfDistrCap[!,:New_Withdraw_Capacity]+dfDistrCap[!,:Withdraw_Margin_Gained]
		dfDistrCap[!,:Total_Inject_Margin] = dfDistrCap[!,:Existing_Inject_Capacity]+dfDistrCap[!,:New_Inject_Capacity]
	end

	if (setup["UCommit"]==0 || (setup["UCommit"]==2 && (setup["Reserves"]==0 || (setup["Reserves"]>0 && inputs["pDynamic_Contingency"]==0)))) # fully linear model
		## Extract dual variables of constraints
		# price: Dual variable of hourly power balance constraint = hourly price
		dfPrice = DataFrame(Zone = 1:Z)
		dfPrice = hcat(dfPrice, convert(DataFrame, transpose(getdual(EP[:cPowerBalance]).*1000)))
		auxNew_Names=[Symbol("Zone");[Symbol("t$t") for t in 1:T]]
		rename!(dfPrice,auxNew_Names)

		# reliability: Dual variable of maximum NSE constraint = shadow value of reliability constraint
		dfReliability = DataFrame(Zone = 1:Z)
		dfReliability = hcat(dfReliability, convert(DataFrame, transpose(getdual(EP[:cMaxNSE]))))
		auxNew_Names=[Symbol("Zone");[Symbol("t$t") for t in 1:T]]
		rename!(dfReliability,auxNew_Names)

		# CO2 emissions by zone
		dfEmissions = DataFrame(Zone = 1:Z, CO2_Price = convert(Array{Union{Missing,Float32}}, zeros(Z)), Sum = Array{Union{Missing,Float32}}(undef, Z))
		if setup["CO2Cap"]==1
			# Dual variable of CO2 constraint = shadow price of CO2
			dfEmissions[!,:CO2_Price] = convert(Array{Union{Missing,Float32}}, getdual(EP[:cCO2Emissions_zonal]))
		elseif setup["CO2Cap"]==2
			# Dual variable of CO2 constraint = shadow price of CO2
			dfEmissions[!,:CO2_Price] .= getdual(EP[:cCO2Emissions_systemwide])
		end

		for i in 1:Z
			dfEmissions[!,:Sum][i] = sum(inputs["omega"].*getvalue(dModuleArgs["dExpressions"]["eEmissionsByZone"])[i,:])
		end
		dfEmissions = hcat(dfEmissions, convert(DataFrame, getvalue(expressions["eEmissionsByZone"])))
		auxNew_Names=[Symbol("Zone");Symbol("CO2_Price");Symbol("Sum");[Symbol("t$t") for t in 1:T]]
		rename!(dfEmissions,auxNew_Names)


		# # Dual of storage level (state of charge) balance of each resource in each time step
		dfStorageDual = DataFrame(Resource = inputs["RESOURCES"], Zone = dfGen[!,:zone])
		# # Define an empty
		# x1 = Array{Float64}(undef, G, T)
		#
		# for y in 1:G
		# 	 if y in dfGen[(dfGen[!,:STOR].==1) .| (dfGen[!,:STOR].==2) .| (dfGen[!,:STOR].==3),:][!,:R_ID]
		# ### # script to write duals
		#  		x1[y,:] =getdual(expressions["SoCBal"][:,y])
		# 	else
		# 		x1[y,:] = zeros(T,1) # Empty values for the resource with no ability to store energy
		# 	end
		# end
		#
		# dfStorageDual=hcat(dfStorageDual, convert(DataFrame,x1))


		# auxNew_Names=[Symbol("Resource");Symbol("Zone");[Symbol("t$t") for t in 1:T]]
		# rename!(dfStorageDual,auxNew_Names)
		# rename!(dfStorageDual,[Symbol("Resource");Symbol("Zone");[Symbol("t$t") for t in 1:T]])

	elseif (setup["Dual_MIP"]==1)
		# function to fix integers and linearize problem
		fix_integers(EP)
		# re-solve statement for LP solution
		solve(EP)
		## Extract dual variables of constraints
		# price: Dual variable of hourly power balance constraint = hourly price
		dfPrice = DataFrame(Zone = 1:Z)
		dfPrice = hcat(dfPrice, convert(DataFrame, transpose(getdual(EP[:cPowerBalance]).*1000)))
		auxNew_Names=[Symbol("Zone");[Symbol("t$t") for t in 1:T]]
		rename!(dfPrice,auxNew_Names)

		# reliability: Dual variable of maximum NSE constraint = shadow value of reliability constraint
		dfReliability = DataFrame(Zone = 1:Z)
		dfReliability = hcat(dfReliability, convert(DataFrame, transpose(getdual(EP[:cMaxNSE]))))
		auxNew_Names=[Symbol("Zone");[Symbol("t$t") for t in 1:T]]
		rename!(dfReliability,auxNew_Names)

		# CO2 emissions by zone
		dfEmissions = DataFrame(Zone = 1:Z, CO2_Price = convert(Array{Union{Missing,Float32}}, zeros(Z)), Sum = Array{Union{Missing,Float32}}(undef, Z))
		if setup["CO2Cap"]==1
			# Dual variable of CO2 constraint = shadow price of CO2
			dfEmissions[!,:CO2_Price] = convert(Array{Union{Missing,Float32}}, getdual(EP[:cCO2Emissions_zonal]))
		elseif setup["CO2Cap"]==2
			# Dual variable of CO2 constraint = shadow price of CO2
			dfEmissions[!,:CO2_Price][1:Z] = getdual(EP[:cCO2Emissions_systemwide])
		end

		for i in 1:Z
			dfEmissions[!,:Sum][i] = sum(inputs["omega"].*getvalue(dModuleArgs["dExpressions"]["eEmissionsByZone"])[i,:])
		end
		dfEmissions = hcat(dfEmissions, convert(DataFrame, getvalue(expressions["eEmissionsByZone"])))
		auxNew_Names=[Symbol("Zone");Symbol("CO2_Price");Symbol("Sum");[Symbol("t$t") for t in 1:T]]
		rename!(dfEmissions,auxNew_Names)

		# # Dual of storage level (state of charge) balance of each resource in each time step
		dfStorageDual = DataFrame(Resource = inputs["RESOURCES"], Zone = dfGen[!,:zone])
		# # Define an empty
		# x1 = Array{Float64}(undef, G, T)
		#
		# for y in 1:G
		# 	 if y in dfGen[(dfGen[!,:STOR].==1) .| (dfGen[!,:STOR].==2) .| (dfGen[!,:STOR].==3),:][!,:R_ID]
		# ### # script to write duals
		#  		x1[y,:] =getdual(expressions["SoCBal"][:,y])
		# 	else
		# 		x1[y,:] = zeros(T,1) # Empty values for the resource with no ability to store energy
		# 	end
		# end
		#
		# dfStorageDual=hcat(dfStorageDual, convert(DataFrame,x1 ))

		# auxNew_Names=[Symbol("Resource");Symbol("Zone");[Symbol("t$t") for t in 1:T]]
		# rename!(dfStorageDual,auxNew_Names)
		# rename!(dfStorageDual,[Symbol("Resource");Symbol("Zone");[Symbol("t$t") for t in 1:T]])

	else
		# CO2 emissions by zone
		dfEmissions = DataFrame(Zone = 1:Z, CO2_Price = convert(Array{Union{Missing,Float32}}, zeros(Z)), Sum = Array{Union{Missing,Float32}}(undef, Z))
		for i in 1:Z
			dfEmissions[!,:Sum][i] = sum(inputs["omega"].*getvalue(dModuleArgs["dExpressions"]["eEmissionsByZone"])[i,:])
		end
		dfEmissions = hcat(dfEmissions, convert(DataFrame, getvalue(expressions["eEmissionsByZone"])))
		auxNew_Names=[Symbol("Zone");Symbol("CO2_Price");Symbol("Sum");[Symbol("t$t") for t in 1:T]]
		rename!(dfEmissions,auxNew_Names)

	end

	total = convert(DataFrame, ["Total" 0 sum(dfEmissions[!,:Sum]) fill(0.0, (1,T))])
	for t in 1:T
		total[!,t+3] .= sum(dfEmissions[!,Symbol("t$t")][1:Z])
	end
	rename!(total,auxNew_Names)
	dfEmissions = vcat(dfEmissions, total)

	# Compile final results Dict depending on configuration setup

	# Common elements in all setups
	results = Dict([
	("STATUS", dfStatus)
	("COSTS", dfCost)
	("CAP", dfCap)
	("POWER", dfPower)
	("CHARGE", dfCharge)
	("STORAGE", dfStorage)
	("FLOWS", dfFlow)
	("TLOSSES", dfTLosses)
	("WLOSSES", dfWLosses)
	("INJECT", dfInject)
	("WITHDRAW", dfWithdraw)
	("NSE", dfNse)
	("EMISSIONS",dfEmissions)
	])
	if (setup["UCommit"]==0 || setup["UCommit"]==2)
		# Add additional elements for fully linear case
		results["PRICE"] = dfPrice
		results["RELIABILITY"] = dfReliability
		results["STORAGEDUAL"] = dfStorageDual
	elseif (setup["UCommit"]==1 && setup["Dual_MIP"]==1)
		results["PRICE"] = dfPrice
		results["RELIABILITY"] = dfReliability
		results["STORAGEDUAL"] = dfStorageDual
	end
	if setup["UCommit"]>=1
		# Add additional elements for Unit Commitment set-up
		results["COMMIT"] = dfCommit
		results["START"] = dfStart
		results["SHUTDOWN"] = dfShutdown
		if (setup["Reserves"]==1)
			results["REG_UP"] = dfRegUp
			results["REG_DN"] = dfRegDn
			results["RSV_UP"] = dfRsvUp
			results["RSV_DN"] = dfRsvDn
		end
	end
	if setup["HeatMarket"]==1
		# Add additional elements for Heat Market set-up
		results["HEAT_MARKET"] = dfHeat_market
		results["HEAT_GEN"] = dfHeat_gen
		results["HEAT_NG"] = dfHeat_NG
		results["HEAT_TO_STORAGE"] = dfHeat_to_Storage
		results["HEAT_TO_TOPPING"] = dfHeat_to_turbine
		results["HEAT_TO_BASEGEN"] = dfcoreheat_to_turbine
	end
	if (setup["NetworkExpansion"]==1)
		results["TRANS_CAP"] = dfTransCap
		results["DISTR_CAP"] = dfDistrCap
	end

	# output additional variables related inter-period energy transfer via storage
	if (setup["OperationWrapping"]==2 && setup["LDS"]==1)
		results["dfStorageInit"] =dfStorageInit
		results["dfdStorage"] = dfdStorage
	end

	## Return results object
	return results
end # END solve_model()

function generate_H2_results(EP::Model, dModuleArgs::Dict)
	dfStatus = dModuleArgs["dfStatus"]
	status = dfStatus[!,:Status][1]

	# ## Check if solved sucessfully - time out is included
	# if status!=:Optimal
	# 	if status!=:UserLimit # Model failed to solve, so record solver status and exi
	# 		if status!=:Suboptimal
	# 			results = Dict([("STATUS", dfStatus)])
	# 			return results
	# 		end
	# 		# Model reached timelimit but failed to find a feasible solution
	# 	elseif isnan(getobjectivevalue(EP))==true
	# 			# Model failed to solve, so record solver status and exit
	# 		results = Dict([("STATUS", dfStatus)])
	# 		return results
	# 	end
	# end

	inputs_H2 = dModuleArgs["inputs_H2"]
	setup = dModuleArgs["setup"]
	dExpressions_H2 = dModuleArgs["dExpressions_H2"]
	dObjective_H2 = dModuleArgs["dObjective_H2"]
	dBalance_H2 = dModuleArgs["dBalance_H2"]

	# dExpressions_H2["Total_Compression_Cost"] = dExpressions_H2["CAPEX_Compression_Truck"] + dExpressions_H2["CAPEX_Compression_Pipe"] + dExpressions_H2["CAPEX_Compression_Stor"]
	# dExpressions_H2["Total_Cost"] = dExpressions_H2["eObjective"]


	if setup["single_hydrogen_sector"] == 0
		key_Cost  = ["Total_H2_Cost","CAPEX_H2G", "CAPEX_EL", "CAPEX_SMR", "CAPEX_PfG","CAPEX_Stor","CAPEX_UG","CAPEX_AG", "CAPEX_Pipe", "CAPEX_Truck","Cost_gas", "Cost_start","OPEX_PfG","OPEX_Truck", "Total_Compression_Cost","Total_CO2_Emission_Cost_H2","Cost_unmet_H2"]
		if setup["battery_truck_option"] == 1
			push!(key_Cost,"CAPEX_Bat_Truck")
		end
	else
		key_Cost  = ["Total_H2_Cost","CAPEX_H2G", "CAPEX_EL", "CAPEX_SMR", "CAPEX_PfG","CAPEX_Stor","CAPEX_UG","CAPEX_AG", "CAPEX_Pipe", "CAPEX_Truck","Cost_electricity","Cost_gas", "Cost_start","OPEX_PfG","OPEX_Truck", "Total_Compression_Cost","Total_CO2_Emission_Cost_H2","Cost_unmet_H2"]
		if setup["battery_truck_option"] == 1
			push!(key_Cost,"CAPEX_Bat_Truck")
		end
	end
	value_Total_Cost = []
	value_Unit_Cost = []
	for key in key_Cost
		push!(value_Total_Cost, getvalue(dExpressions_H2[key])/1e+6)
		push!(value_Unit_Cost, (getvalue(dExpressions_H2[key]))/(getvalue(dExpressions_H2["Total_H2_Demand"]*1e+3)))
	end
	dfH2Cost = DataFrame(Costs = key_Cost)
	dfH2Cost[!,Symbol("Total (million \$)")] = value_Total_Cost
	dfH2Cost[!,Symbol("Unit (\$/kg)")] = value_Unit_Cost

	key_Emission  = ["H2_Carbon_Emission","Gen_carbon_emission", "Truck_carbon_emission"]
	value_Emission = []
	for key in key_Emission
		push!(value_Emission, getvalue(dExpressions_H2[key]))
	end
	dfEmissions = DataFrame(Emission = key_Emission)
	dfEmissions[!,Symbol("Total (tonne)")] = value_Emission


	key_PlanVar  = ["vN","vN_int","vM","vH2StorRate", "vH2StorCap","vNPipe","vNTruck","vH2TruckCompressionCap"]
	if setup["battery_truck_option"] == 1
		push!(key_PlanVar,"vN_Bat_Truck")
	end
	value_PlanVar = Dict()
	for key in key_PlanVar
		value_PlanVar[key] = getvalue(EP[Symbol(key)])
	end
	dH2Plan = value_PlanVar

	key_GenVar  = ["vN","vN_int", "vn","vn_start", "vn_int","vn_start_int","vH2Gen","vP2G","vGas"]
	value_GenVar = Dict()
	for key in key_GenVar
		value_GenVar[key] = getvalue(EP[Symbol(key)])
	end
	dH2Gen = value_GenVar

	key_H2PowerGenVar  = ["vM", "vm", "vPfG","vH2P"]
	value_H2PowerGenVar = Dict()
	for key in key_H2PowerGenVar
		value_H2PowerGenVar[key] = getvalue(EP[Symbol(key)])
	end
	dH2Power = value_H2PowerGenVar

	key_StorageVar  = ["vH2StorCap","vH2StorDis", "vH2StorCha", "vH2StorEnergy", "vH2StorRate"]
	if setup["LDS"]==1
		key_StorageVar  = ["vH2StorCap","vH2StorDis", "vH2StorCha", "vH2StorEnergy", "vH2StorRate","vH2SOCw","vH2dSOC"]
	end

	value_StorageVar = Dict()
	for key in key_StorageVar
		value_StorageVar[key] = getvalue(EP[Symbol(key)])
	end
	dH2Storage = value_StorageVar

	key_PipeVar  = ["vNPipe","vH2PipeFlow", "vH2PipeLevel"]
	value_PipeVar = Dict()
	for key in key_PipeVar
		value_PipeVar[key] = getvalue(EP[Symbol(key)])
	end
	dH2Pipe = value_PipeVar

	key_TruckVar  = ["vNTruck","vH2TruckFlow", "vNavail_full","vNtravel_full","vNarrive_full","vNdepart_full","vNavail_empty","vNtravel_empty","vNarrive_empty","vNdepart_empty","vNcharged","vNdischarged","vN_full","vN_empty","vH2TruckCompressionCap"]
	if setup["LDS"]==1
		key_TruckVar  = ["vNTruck","vH2TruckFlow", "vNavail_full","vNtravel_full","vNarrive_full","vNdepart_full","vNavail_empty","vNtravel_empty","vNarrive_empty","vNdepart_empty","vNcharged","vNdischarged","vN_full","vN_empty","vH2TruckCompressionCap","vH2TruckSOCw","vH2TruckdSOC"]
	end
	value_TruckVar = Dict()
	for key in key_TruckVar
		value_TruckVar[key] = getvalue(EP[Symbol(key)])
	end
	dH2Truck = value_TruckVar

	key_DemandVar  = ["vH2H","vH2I","vH2D","vH2DUnmet","vH2Curtail"]
	value_DemandVar = Dict()
	for key in key_DemandVar
		value_DemandVar[key] = getvalue(EP[Symbol(key)])
	end
	dH2Demand = value_DemandVar
	#
	# dModuleArgs["dExpressions_H2"] = dExpressions_H2

	Price_H2 = getdual(EP[:cH2Balance])

	results_H2 = Dict([
	("STATUS", dfStatus)
	("H2COSTS", dfH2Cost)
	("H2PLAN", dH2Plan)
	("H2GEN", dH2Gen)
	("H2STORAGE", dH2Storage)
	("H2PIPE", dH2Pipe)
	("H2TRUCK", dH2Truck)
	("H2Power",dH2Power)
	("H2DEMAND",dH2Demand)
	("EMISSIONS",dfEmissions)
	("H2Price",Price_H2)
	])
	return results_H2
end

function save_H2_results(inputs_H2::Dict, results_H2::Dict,resultpath_H2)
	if (isdir(resultpath_H2)==false)
	   mkdir(resultpath_H2)
	end
	sep = "/"
	T = inputs_H2["T"]
	weights = inputs_H2["weights"]
	CSV.write(string(resultpath_H2,sep,"weights.csv"),DataFrame(transpose(weights)), writeheader=false)

    ## H2 Planning
    dfPlan = results_H2["H2PLAN"]
    json_string = JSON.json(dfPlan)
    open(string(resultpath_H2,sep,"dfPlan.json"),"w") do f
      JSON.print(f, json_string)
    end

    H2vNCap = sum(inputs_H2["H2GenSize"].*results_H2["H2PLAN"]["vN"],dims = 2)[:,1]
    H2vNintCap = sum(inputs_H2["H2GenSize_int"].*results_H2["H2PLAN"]["vN_int"],dims = 2)[:,1]
    H2PfGCap = sum(inputs_H2["PfGSize"].*results_H2["H2PLAN"]["vM"],dims = 2)[:,1]
    H2StorageCap = sum(results_H2["H2PLAN"]["vH2StorCap"],dims = 2)[:,1]
    H2PipeCap = inputs_H2["H2PipeFlowSize"].*sum(sum(results_H2["H2PLAN"]["vNPipe"],dims = 1),dims = 2)[:,:,1]/2
    H2TruckCap = inputs_H2["TruckCap"].*results_H2["H2PLAN"]["vNTruck"]

   H2Cap = [H2vNCap; H2vNintCap;H2PfGCap;H2StorageCap;H2PipeCap;H2TruckCap]
   dfH2Cap = DataFrame(H2Cap)
   CSV.write(string(resultpath_H2,sep,"H2CAP.csv"),dfH2Cap)

    # for key in keys(dfPlan)
    #     # println(DataFrame(dfPlan[key]))
    #     # CSV.write(string(resultpath_H2,sep,key,".csv"),DataFrame(dfPlan[key]))
    #     writedlm(string(resultpath_H2,sep,key,".csv"), dfPlan[key])
    # end
    H2Plan_truck=(results_H2["H2PLAN"]["vNTruck"])
    writedlm(string(resultpath_H2,sep,"H2Plan_vNTruck.csv"), H2Plan_truck)
    H2PLAN_vN_int=DataFrame(results_H2["H2PLAN"]["vN_int"])
    CSV.write(string(resultpath_H2,sep,"H2PLAN_vN_int.csv"),H2PLAN_vN_int)
    H2PLAN_vH2StorCap=DataFrame(results_H2["H2PLAN"]["vH2StorCap"])
    CSV.write(string(resultpath_H2,sep,"H2PLAN_vH2StorCap.csv"),H2PLAN_vH2StorCap)
    H2PLAN_vM=DataFrame(results_H2["H2PLAN"]["vM"])
    CSV.write(string(resultpath_H2,sep,"H2PLAN_vM.csv"),H2PLAN_vM)
    H2PLAN_vN=DataFrame(results_H2["H2PLAN"]["vN"])
    CSV.write(string(resultpath_H2,sep,"H2PLAN_vN.csv"),H2PLAN_vN)
    H2PLAN_vNPipe=DataFrame((results_H2["H2PLAN"]["vNPipe"][:,:,1]))
    CSV.write(string(resultpath_H2,sep,"H2PLAN_vNPipe.csv"),H2PLAN_vNPipe)

    ## H2 Emissions
    EMISSIONS=results_H2["EMISSIONS"]
    CSV.write(string(resultpath_H2,sep,"EMISSIONS.csv"),EMISSIONS)

    # H2 Costs
    H2COSTS=results_H2["H2COSTS"]
    CSV.write(string(resultpath_H2,sep,"H2COSTS.csv"),H2COSTS)

    ## H2 Demand
    H2DEMAND_H2I=DataFrame(results_H2["H2DEMAND"]["vH2I"])
    CSV.write(string(resultpath_H2,sep,"H2DEMAND_H2I.csv"),H2DEMAND_H2I)

    H2DEMAND_H2D=DataFrame(results_H2["H2DEMAND"]["vH2D"])
    CSV.write(string(resultpath_H2,sep,"H2DEMAND_H2D.csv"),H2DEMAND_H2D)

    H2DEMAND_H2H=DataFrame(results_H2["H2DEMAND"]["vH2H"])
    CSV.write(string(resultpath_H2,sep,"H2DEMAND_H2H.csv"),H2DEMAND_H2H)

	H2DEMAND_H2DUnmet=DataFrame(results_H2["H2DEMAND"]["vH2DUnmet"])
    CSV.write(string(resultpath_H2,sep,"H2DEMAND_H2DUnmet.csv"),H2DEMAND_H2DUnmet)

	H2DEMAND_H2Curtail=DataFrame(results_H2["H2DEMAND"]["vH2Curtail"])
    CSV.write(string(resultpath_H2,sep,"H2DEMAND_H2Curtail.csv"),H2DEMAND_H2Curtail)

    ## H2 Gen
    H2Power_vH2Gen=(results_H2["H2GEN"]["vH2Gen"])
    resource_num = size(H2Power_vH2Gen)[1]
    for i in 1:resource_num
        H2Power_vH2Gen_new=DataFrame(H2Power_vH2Gen[i,:,:])
        CSV.write(string(resultpath_H2,sep,"H2GEN_vH2Gen_",i,".csv"),H2Power_vH2Gen_new)
    end

    H2Power_vGas=(results_H2["H2GEN"]["vGas"])
    for i in 1:resource_num
        H2Power_vGas_new=DataFrame(H2Power_vGas[i,:,:])
        CSV.write(string(resultpath_H2,sep,"H2GEN_vGas_",i,".csv"),H2Power_vGas_new)
    end

    H2Power_vP2G=(results_H2["H2GEN"]["vP2G"])
    for i in 1:resource_num
        H2Power_vP2G_new=DataFrame(H2Power_vP2G[i,:,:])
        CSV.write(string(resultpath_H2,sep,"H2GEN_vP2G_",i,".csv"),H2Power_vP2G_new)
    end

    H2Power_vn=(results_H2["H2GEN"]["vn"])
    for i in 1:resource_num
        H2Power_vn_new=DataFrame(H2Power_vn[i,:,:])
        CSV.write(string(resultpath_H2,sep,"H2GEN_vn_",i,".csv"),H2Power_vn_new)
    end

    ## H2 to Power
    H2Power_vH2P=(results_H2["H2Power"]["vH2P"])
    resource_num = size(H2Power_vH2P)[1]
    for i in 1:resource_num
        H2Power_vH2P_new=DataFrame(H2Power_vH2P[i,:,:])
        CSV.write(string(resultpath_H2,sep,"H2Power_vH2P_",i,".csv"),H2Power_vH2P_new)
    end

    H2Power_vPfG=(results_H2["H2Power"]["vPfG"])
    resource_num = size(H2Power_vPfG)[1]
    for i in 1:resource_num
        H2Power_vPfG_new=DataFrame(H2Power_vPfG[i,:,:])
        CSV.write(string(resultpath_H2,sep,"H2Power_vPfG_",i,".csv"),H2Power_vPfG_new)
    end

    H2Power_vm=(results_H2["H2Power"]["vm"])
    resource_num = size(H2Power_vm)[1]
    for i in 1:resource_num
        H2Power_vm_new=DataFrame(H2Power_vm[i,:,:])
        CSV.write(string(resultpath_H2,sep,"H2Power_vm_",i,".csv"),H2Power_vm_new)
    end

    ## H2 Pipeline
    for p in 1:inputs_H2["PT"]
        vPipeFlow = sum(results_H2["H2PIPE"]["vH2PipeFlow"],dims = 2)[:,1,p,:]
        CSV.write(string(resultpath_H2,sep,"vPipeFlow_",p,".csv"),DataFrame(vPipeFlow))
    end

	for p in 1:inputs_H2["PT"]
        vH2PipeLevel = results_H2["H2PIPE"]["vH2PipeLevel"]
		resultpath_vH2PipeLevel = string(resultpath_H2,sep,"vH2PipeLevel_",p)
		if (isdir(resultpath_vH2PipeLevel)==false)
			 mkdir(resultpath_vH2PipeLevel)
		end
		for z in 1:inputs_H2["Z"]
			CSV.write(string(resultpath_vH2PipeLevel,sep,"vH2PipeLevel_",z,".csv"),DataFrame(vH2PipeLevel[z,:,p,:]))
		end
    end

    ## H2 Storage
	dfStorage = results_H2["H2STORAGE"]
    for k in 1:inputs_H2["K_stor"]

		for key in keys(dfStorage)
			if length(size(dfStorage[key])) == 2
				CSV.write(string(resultpath_H2,sep,key,"_",k,".csv"),DataFrame(dfStorage[key]))
			#     println(dfPlan[key])
			elseif length(size(dfStorage[key])) == 3
				CSV.write(string(resultpath_H2,sep,key,"_",k,".csv"),DataFrame(dfStorage[key][k,:,:]))
			# elseif length(size(dfTruck[key])) == 4
			# 	resultpath_H2Truck = string(resultpath_H2,sep,"truck_travel_",k)
			# 	if (isdir(resultpath_H2Truck)==false)
			# 		 mkdir(resultpath_H2Truck)
			# 	end
			# 	for z in 1:inputs_H2["Z"]
			# 		CSV.write(string(resultpath_H2Truck,sep,key,"_",z,".csv"),DataFrame(dfTruck[key][z,:,k,:]))
			# 	end
			end
		end

        # vH2StorDis = results_H2["H2STORAGE"]["vH2StorDis"][k,:,:]
        # vH2StorCha = results_H2["H2STORAGE"]["vH2StorCha"][k,:,:]
        # vH2StorEnergy = results_H2["H2STORAGE"]["vH2StorEnergy"][k,:,:]
        # CSV.write(string(resultpath_H2,sep,"vH2StorDis_",k,".csv"),DataFrame(vH2StorDis))
        # CSV.write(string(resultpath_H2,sep,"vH2StorCha_",k,".csv"),DataFrame(vH2StorCha))
        # CSV.write(string(resultpath_H2,sep,"vH2StorEnergy_",k,".csv"),DataFrame(vH2StorEnergy))
		#
        # vH2StorRate=DataFrame(results_H2["H2STORAGE"]["vH2StorRate"])
        # CSV.write(string(resultpath_H2,sep,"vH2StorRate_",k,".csv"),DataFrame(vH2StorRate))

    end

    ## H2 Truck
    dfTruck = results_H2["H2TRUCK"]
    for r in 1:inputs_H2["RT"]
	    for key in keys(dfTruck)
	        if length(size(dfTruck[key])) == 2
	            CSV.write(string(resultpath_H2,sep,key,"_",r,".csv"),DataFrame(dfTruck[key]))
	        #     println(dfPlan[key])
	        elseif length(size(dfTruck[key])) == 3
	            CSV.write(string(resultpath_H2,sep,key,"_",r,".csv"),DataFrame(dfTruck[key][:,r,:]))
	        elseif length(size(dfTruck[key])) == 4
	            resultpath_H2Truck = string(resultpath_H2,sep,"truck_travel_",r)
	            if (isdir(resultpath_H2Truck)==false)
	        	     mkdir(resultpath_H2Truck)
	        	end
	            for z in 1:inputs_H2["Z"]
	                CSV.write(string(resultpath_H2Truck,sep,key,"_",z,".csv"),DataFrame(dfTruck[key][z,:,r,:]))
	            end
	        end
	    end
    end

    ## H2 Price
    H2Price=DataFrame(results_H2["H2Price"])
    CSV.write(string(resultpath_H2,sep,"H2Price.csv"),H2Price)
	return
end

function generate_H2_model(setup::Dict,inputs_H2::Dict,OPTIMIZER,modeloutput = nothing)
	presolver_start_time = time()
	## Model Definition
	# Define the Energy Portfolio (HY) model
	# Set solver to use Gurobi
	HY = Model(solver = OPTIMIZER)


	# HY=Model(solver=CbcSolver(logLevel = 2))
   	# HY=Model(solver=AmplNLSolver("bonmin"))

	T = inputs_H2["T"]     # Number of time steps (hours)
	Z = inputs_H2["Z"]     # Number of zones

	@variable(HY, vPower2H2_H[t = 1:T,z=1:Z])

	dModuleArgs = Dict()
	dModuleArgs["inputs_H2"] = inputs_H2
	dModuleArgs["setup"] = setup
	dModuleArgs["dExpressions_H2"] = Dict()
	dModuleArgs["dObjective_H2"] = Dict()
	dModuleArgs["dBalance_H2"] = Dict()

	dExpressions_H2 = dModuleArgs["dExpressions_H2"]
	dObjective_H2 = dModuleArgs["dObjective_H2"]
	dBalance_H2 = dModuleArgs["dBalance_H2"]

	@expression(HY, H2Balance[t=1:T, z=1:Z], 0)
	dBalance_H2["H2Balance"] = H2Balance
	@expression(HY, H2PowerBalance[t=1:T, z=1:Z], 0)
	dBalance_H2["H2PowerBalance"] = H2PowerBalance
	@expression(HY, H2_Carbon_Emission, 0)
	dBalance_H2["H2_Carbon_Emission"] = H2_Carbon_Emission

	discount_rate = inputs_H2["discount_rate"]
	life = inputs_H2["life"]
	@expression(HY, discount_factor, discount_rate/(1-(1+discount_rate)^(-life)))
	dExpressions_H2["discount_factor"] = discount_factor

	dModuleArgs["dExpressions_H2"] = dExpressions_H2
	dModuleArgs["dObjective_H2"] = dObjective_H2
	dModuleArgs["dBalance_H2"] = dBalance_H2

	HY, dModuleArgs =  H2_Gen(HY, dModuleArgs)
	HY, dModuleArgs =  H2_PowerGen(HY, dModuleArgs)
	HY, dModuleArgs =  H2_Demand(HY, dModuleArgs)
	HY, dModuleArgs =  H2_Storage(HY, dModuleArgs)
	HY, dModuleArgs =  H2Pipeline(HY, dModuleArgs)
	if inputs_H2["share_truck_options"] == 1
		HY, dModuleArgs =  H2Truck(HY, dModuleArgs)
	else
		HY, dModuleArgs =  H2Truck_old(HY, dModuleArgs)
	end

	dExpressions_H2 = dModuleArgs["dExpressions_H2"]
	dObjective_H2 = dModuleArgs["dObjective_H2"]
	dBalance_H2 = dModuleArgs["dBalance_H2"]
	inputs_H2 = dModuleArgs["inputs_H2"]

	weights = inputs_H2["weights"]

	if setup["single_hydrogen_sector"] == 0

		# setlowerbound(vPower2H2_P[t,z], inputs_H2["Power2H2"])
		# setupperbound(vPower2H2_P[t,z], inputs_H2["Power2H2"])

		# @constraint(HY, cPower2H2_H[t=1:T, z=1:Z], vPower2H2_H[t,z] == inputs_H2["Power2H2"])
		#
		# @variable(HY, 1e+4>= vPower2H2_H_Slack_pos[t = 1:T,z=1:Z] >=0)
		# @variable(HY, 1e+4>= vPower2H2_H_Slack_neg[t = 1:T,z=1:Z] >=0)
		#
		# @expression(HY, ePower2H2_relaxed[t=1:T, z=1:Z],
		# vPower2H2_H[t,z] + vPower2H2_H_Slack_pos[t,z] - vPower2H2_H_Slack_neg[t,z]
		# )
		# dExpressions_H2["ePower2H2_relaxed"] = ePower2H2_relaxed
		# for z in 1:Z
		#     for t in 1:T
		# 		@constraint(HY, ePower2H2_relaxed[t,z] == dExpressions_H2["eH2PowerConsumption"][t,z])
		# 	end
		# end

		electricity_price = inputs_H2["electricity_price"]
		@expression(HY, Cost_electricity, sum(weights[t] * electricity_price[t,z] * dBalance_H2["H2PowerBalance"][t,z] for z=1:Z, t = 1:T))
		dExpressions_H2["Cost_electricity"] = Cost_electricity
		dObjective_H2["Cost_electricity"] = Cost_electricity
	else
		electricity_price = inputs_H2["electricity_price"]
		@expression(HY, Cost_electricity, sum(weights[t] * electricity_price[t,z] * dBalance_H2["H2PowerBalance"][t,z] for z=1:Z, t = 1:T))
		dExpressions_H2["Cost_electricity"] = Cost_electricity
		dObjective_H2["Cost_electricity"] = Cost_electricity
	end

		# @constraint(HY, cPowerBalance[t=1:T, z=1:Z], dBalance_H2["H2PowerBalance"][t,z] == 0)

	@constraint(HY, cH2Balance[t=1:T, z=1:Z], dBalance_H2["H2Balance"][t,z] == 0)

	Total_H2_Demand = dExpressions_H2["Total_H2_Demand"]
	H2_emission_rate = inputs_H2["H2_emission_rate"]

	dExpressions_H2["H2_Carbon_Emission"] = dBalance_H2["H2_Carbon_Emission"]
	dExpressions_H2["Total_CO2_Emission_Cost_H2"] = dBalance_H2["H2_Carbon_Emission"] * inputs_H2["CO2_price"]

	if setup["CO2_price_option"] == 1
		dObjective_H2["eTotalCO2EmissionCost_H2"] = dExpressions_H2["Total_CO2_Emission_Cost_H2"]
	else
		dObjective_H2["Total_CO2_Emission_Cost_H2"] = 0
		@constraint(EP, cCarbonEmission, dBalance_H2["H2_Carbon_Emission"] <= Total_H2_Demand * H2_emission_rate)
	end

	global eObjective = 0
	for key in keys(dObjective_H2)
	    global eObjective
	    eObjective = eObjective + dObjective_H2[key]
	end
	dExpressions_H2["eObjective"] = eObjective
	dExpressions_H2["Total_H2_Cost"] = dExpressions_H2["eObjective"]
	dExpressions_H2["Total_Compression_Cost"] = dExpressions_H2["CAPEX_Compression_Truck"] + dExpressions_H2["CAPEX_Compression_Pipe"] + dExpressions_H2["CAPEX_Compression_Stor"]

	@objective(HY,Min,eObjective)

	dModuleArgs["dExpressions_H2"] = dExpressions_H2
	dModuleArgs["dObjective_H2"] = dObjective_H2
	dModuleArgs["dBalance_H2"] = dBalance_H2
	return HY, dModuleArgs
end

function solve_H2_model(HY::Model, dModuleArgs::Dict)
	solver_start_time = time()

	## Solve Model
	status = solve(HY)
	## Record solver time
	solver_time = time() - solver_start_time

	dfStatus = DataFrame(Status = status, Solve = solver_time,
		Objval = getobjectivevalue(HY), Objbound= getobjbound(HY),FinalMIPGap =(getobjectivevalue(HY) -getobjbound(HY))/getobjectivevalue(HY) )

	## Check if solved sucessfully - time out is included
	if status!=:Optimal
		if status!=:UserLimit # Model failed to solve, so record solver status and exi
			results = Dict([("STATUS", dfStatus)])
			return results
			# Model reached timelimit but failed to find a feasible solution
		elseif isnan(getobjectivevalue(HY))==true
				# Model failed to solve, so record solver status and exit
			results = Dict([("STATUS", dfStatus)])
			return results
		end
	end

	inputs_H2 = dModuleArgs["inputs_H2"]
	setup = dModuleArgs["setup"]
	dExpressions_H2 = dModuleArgs["dExpressions_H2"]
	dObjective_H2 = dModuleArgs["dObjective_H2"]
	dBalance_H2 = dModuleArgs["dBalance_H2"]

	dExpressions_H2["Total_Compression_Cost"] = dExpressions_H2["CAPEX_Compression_Truck"] + dExpressions_H2["CAPEX_Compression_Pipe"] + dExpressions_H2["CAPEX_Compression_Stor"]
	dExpressions_H2["Total_Cost"] = dExpressions_H2["eObjective"]

	key_Cost  = ["Total_Cost","CAPEX_H2G", "CAPEX_EL", "CAPEX_SMR", "CAPEX_Stor","CAPEX_UG","CAPEX_AG","Cost_gas", "Cost_electricity", "CAPEX_Pipe", "CAPEX_Truck", "OPEX_Truck", "Total_Compression_Cost","Cost_unmet_H2"]
	value_Total_Cost = []
	value_Unit_Cost = []
	for key in key_Cost
		push!(value_Total_Cost, getvalue(dExpressions_H2[key])/1000000)
		push!(value_Unit_Cost, getvalue(dExpressions_H2[key])/getvalue(dExpressions_H2["Total_H2_Demand"])/1000)
	end
	dfH2Cost = DataFrame(Costs = key_Cost)
	dfH2Cost[!,Symbol("Total (million \$)")] = value_Total_Cost
	dfH2Cost[!,Symbol("Unit (\$/kg)")] = value_Unit_Cost

	key_Emission  = ["H2_Carbon_Emission","Gen_carbon_emission", "Truck_carbon_emission"]
	value_Emission = []
	for key in key_Emission
		push!(value_Emission, getvalue(dExpressions_H2[key]))
	end
	dfEmissions = DataFrame(Emission = key_Emission)
	dfEmissions[!,Symbol("Total (tonne)")] = value_Emission


	key_PlanVar  = ["vN","vN_int","vM","vH2StorRate", "vH2StorCap","vNPipe","vNTruck","vH2TruckCompressionCap"]
	value_PlanVar = Dict()
	for key in key_PlanVar
		value_PlanVar[key] = getvalue(HY[Symbol(key)])
	end
	dH2Plan = value_PlanVar

	key_GenVar  = ["vN","vN_int", "vn", "vn_int","vH2Gen","vP2G","vGas"]
	value_GenVar = Dict()
	for key in key_GenVar
		value_GenVar[key] = getvalue(HY[Symbol(key)])
	end
	dH2Gen = value_GenVar

	key_StorageVar  = ["vH2StorCap","vH2StorDis", "vH2StorCha", "vH2StorEnergy", "vH2StorRate"]
	value_StorageVar = Dict()
	for key in key_StorageVar
		value_StorageVar[key] = getvalue(HY[Symbol(key)])
	end
	dH2Storage = value_StorageVar

	key_PipeVar  = ["vNPipe","vH2PipeFlow", "vH2PipeLevel"]
	value_PipeVar = Dict()
	for key in key_PipeVar
		value_PipeVar[key] = getvalue(HY[Symbol(key)])
	end
	dH2Pipe = value_PipeVar

	key_TruckVar  = ["vNTruck","vH2TruckFlow", "vNavail_full","vNtravel_full","vNarrive_full","vNdepart_full","vNavail_empty","vNtravel_empty","vNarrive_empty","vNdepart_empty","vNcharged","vNdischarged","vN_full","vN_empty","vH2TruckCompressionCap"]
	value_TruckVar = Dict()
	for key in key_TruckVar
		value_TruckVar[key] = getvalue(HY[Symbol(key)])
	end
	dH2Truck = value_TruckVar

	key_DemandVar  = ["vH2H","vH2I","vH2D","vH2DUnmet","vH2Curtail"]
	value_DemandVar = Dict()
	for key in key_DemandVar
		value_DemandVar[key] = getvalue(HY[Symbol(key)])
	end
	dH2Demand = value_DemandVar


	results_H2 = Dict([
	("STATUS", dfStatus)
	("H2COSTS", dfH2Cost)
	("H2PLAN", dH2Plan)
	("H2GEN", dH2Gen)
	("H2STORAGE", dH2Storage)
	("H2PIPE", dH2Pipe)
	("H2TRUCK", dH2Truck)
	# ("H2Power",dH2Power)
	("H2DEMAND",dH2Demand)
	("EMISSIONS",dfEmissions)
	])
	return HY, dModuleArgs, results_H2
end

# H2 Generation Module
function H2_Gen(HY::Model, dModuleArgs::Dict)

	println("Hydrogen Generation Module")

	inputs_H2 = dModuleArgs["inputs_H2"]
	setup = dModuleArgs["setup"]
	dExpressions_H2 = dModuleArgs["dExpressions_H2"]
	dObjective_H2 = dModuleArgs["dObjective_H2"]
	dBalance_H2 = dModuleArgs["dBalance_H2"]


	Z_total = inputs_H2["Z_total"]
    PT_total = inputs_H2["PT_total"]
    RT_total = inputs_H2["RT_total"]
    K_prod_total = inputs_H2["K_prod_total"]
    K_stor_total = inputs_H2["K_stor_total"]
    G_H2_total = inputs_H2["G_H2_total"]
    T = inputs_H2["T"]
    Z = inputs_H2["Z"]
	Z_set = inputs_H2["Z_set"]
    PT = inputs_H2["PT"]
    RT = inputs_H2["RT"]
    K_prod = inputs_H2["K_prod"]
    K_prod_int = inputs_H2["K_prod_int"]
    K_stor = inputs_H2["K_stor"]
    G_H2 = inputs_H2["G_H2"]
    W = inputs_H2["W"]
    Tw = inputs_H2["Tw"]

	gen_category = inputs_H2["gen_category"]
	no_integer = inputs_H2["no_integer"]
	H2Gen_integer = inputs_H2["H2Gen_integer"]
	weights = inputs_H2["weights"]
	GasPrice= inputs_H2["GasPrice"]
	GasPrice_int= inputs_H2["GasPrice_int"]
	H2Gen_zone_option = inputs_H2["H2Gen_zone_option"]
	H2Gen_int_zone_option = inputs_H2["H2Gen_int_zone_option"]

	etaP2G = inputs_H2["etaP2G"]
	etaGas = inputs_H2["etaGas"]
	rhoH2Gen_min = inputs_H2["rhoH2Gen_min"]
	rhoH2Gen_max = inputs_H2["rhoH2Gen_max"]
	H2GenSize = inputs_H2["H2GenSize"]
	P2GenUnitCapex = inputs_H2["P2GenUnitCapex"]
	start_up_cost = inputs_H2["start_up_cost"]
	tau_up = inputs_H2["tau_up"]
	tau_down = inputs_H2["tau_down"]
	min_up_ratio = inputs_H2["min_up_ratio"]
	gen_emission_rate = inputs_H2["gen_emission_rate"]

	etaP2G_int = inputs_H2["etaP2G_int"]
	etaGas_int = inputs_H2["etaGas_int"]
	rhoH2Gen_min_int = inputs_H2["rhoH2Gen_min_int"]
	rhoH2Gen_max_int = inputs_H2["rhoH2Gen_max_int"]
	H2GenSize_int = inputs_H2["H2GenSize_int"]
	P2GenUnitCapex_int = inputs_H2["P2GenUnitCapex_int"]
	start_up_cost_int = inputs_H2["start_up_cost_int"]
	tau_up_int = inputs_H2["tau_up_int"]
	tau_down_int = inputs_H2["tau_down_int"]
	min_up_ratio_int = inputs_H2["min_up_ratio_int"]
	gen_emission_rate_int = inputs_H2["gen_emission_rate_int"]

	hydrogen_generation_option = inputs_H2["hydrogen_generation_option"]
	SMR_option_avai = inputs_H2["SMR_option_avai"]
	EL_option_avai = inputs_H2["EL_option_avai"]
	hydrogen_demand_option = inputs_H2["hydrogen_demand_option"]
	EL_cost_factor = inputs_H2["EL_cost_factor"]
	# P2GenUnitCapex[1]=P2GenUnitCapex[1]*EL_cost_factor
	P2GenUnitCapex[1]=P2GenUnitCapex[1]/450*1000*EL_cost_factor

	K_EL = 1

	### Variables ###
    @variable(HY, vH2Gen[k=1:K_prod,z=1:Z, t = 1:T] >= 0 )
    @variable(HY, vP2G[k=1:K_prod,z=1:Z, t = 1:T] >= 0 )
    @variable(HY, vGas[k=1:K_prod,z=1:Z, t = 1:T] >= 0 )

	@variable(HY, vn_start[k=1:K_prod,z=1:Z, t = 1:T] >= 0 )
	@variable(HY, vn_shut[k=1:K_prod,z=1:Z, t = 1:T] >= 0 )

    if no_integer == 1
		if H2Gen_integer == 1
			@variable(HY, vn[k=1:K_prod,z=1:Z, t = 1:T] >= 0, Int )
		else
			@variable(HY, vn[k=1:K_prod,z=1:Z, t = 1:T] >= 0 )
		end
        @variable(HY, vN[k=1:K_prod,z=1:Z] >= 0 )
        # @variable(HY, 0 >= vN_new[k=1:K_prod,z=1:Z] >= 0 )
        @variable(HY, 0 >= vN_existing[k=1:K_prod,z=1:Z] >= 0 )
        @variable(HY, 0 >= vN_retired[k=1:K_prod,z=1:Z] >= 0 )
    else
        @variable(HY, vn[k=1:K_prod,z=1:Z, t = 1:T] >= 0, Int )

        @variable(HY, vN[k=1:K_prod,z=1:Z] >= 0, Int )
        # @variable(HY, 0 >= vN_new[k=1:K_prod,z=1:Z] >= 0, Int )
        @variable(HY, 0 >= vN_existing[k=1:K_prod,z=1:Z] >= 0, Int )
        @variable(HY, 0 >= vN_retired[k=1:K_prod,z=1:Z] >= 0, Int )
    end
	# gen_category = 0
	if gen_category == 1
	    @variable(HY, vH2Gen_int[k=1:K_prod_int,z=1:Z, t = 1:T] >= 0 )
	    @variable(HY, vP2G_int[k=1:K_prod_int,z=1:Z, t = 1:T] >= 0 )
	    @variable(HY, vGas_int[k=1:K_prod_int,z=1:Z, t = 1:T] >= 0 )
		@variable(HY, vn_start_int[k=1:K_prod_int,z=1:Z, t = 1:T] >= 0 )
		@variable(HY, vn_shut_int[k=1:K_prod_int,z=1:Z, t = 1:T] >= 0 )

	    if no_integer == 1
			if H2Gen_integer == 1
				@variable(HY, vn_int[k=1:K_prod_int,z=1:Z, t = 1:T] >= 0, Int )
			else
				@variable(HY, vn_int[k=1:K_prod_int,z=1:Z, t = 1:T] >= 0 )
			end


	        @variable(HY, vN_int[k=1:K_prod_int,z=1:Z] >= 0, Int )
	        # @variable(HY, vN_new_int[k=1:K_prod_int,z=1:Z] >= 0, Int )
	        @variable(HY, 0 >= vN_existing_int[k=1:K_prod_int,z=1:Z] >= 0 )
	        @variable(HY, 0 >= vN_retired_int[k=1:K_prod_int,z=1:Z] >= 0 )
	    else
	        @variable(HY, vn_int[k=1:K_prod_int,z=1:Z, t = 1:T] >= 0, Int )
	        @variable(HY, vN_int[k=1:K_prod_int,z=1:Z] >= 0, Int )
	        # @variable(HY, vN_new_int[k=1:K_prod_int,z=1:Z] >= 0, Int )
	        @variable(HY, 0 >= vN_existing_int[k=1:K_prod_int,z=1:Z] >= 0, Int )
	        @variable(HY, 0 >= vN_retired_int[k=1:K_prod_int,z=1:Z] >= 0, Int )
	    end
	else
	    @variable(HY, 0 >= vH2Gen_int[k=1:K_prod_int,z=1:Z, t = 1:T] >= 0 )
	    @variable(HY, 0 >= vP2G_int[k=1:K_prod_int,z=1:Z, t = 1:T] >= 0 )
	    @variable(HY, 0 >= vGas_int[k=1:K_prod_int,z=1:Z, t = 1:T] >= 0 )
	    @variable(HY, 0 >= vn_int[k=1:K_prod_int,z=1:Z, t = 1:T] >= 0 )
	    @variable(HY, 0 >= vN_int[k=1:K_prod_int,z=1:Z] >= 0 )
		@variable(HY, 0 >= vn_start_int[k=1:K_prod_int,z=1:Z, t = 1:T] >= 0 )
		@variable(HY, 0 >= vn_shut_int[k=1:K_prod_int,z=1:Z, t = 1:T] >= 0 )
	end

	if hydrogen_generation_option == 0
		setupperbound.(vH2Gen, 0)
		setupperbound.(vP2G, 0)
		setupperbound.(vGas, 0)
		setupperbound.(vn, 0)
		setupperbound.(vN, 0)
		setupperbound.(vn_start, 0)
		# setupperbound.(vN_new, 0)
	end
	if SMR_option_avai == 0
		setupperbound.(vGas, 0)

	end
	if EL_option_avai == 0
		setupperbound.(vP2G, 0)

	end
	### End Variables ###

	### Expressions ###
	## Objective Function Expressions ##
	discount_factor = dExpressions_H2["discount_factor"]

	@expression(HY, CAPEX_H2G, ( sum(inputs_H2["discount_factor_P2Gen"][k,z]*vN[k,z]* P2GenUnitCapex[k] for k=1:K_prod,z=1:Z) + inputs_H2["discount_factor_P2Gen_int"][k,z]*sum(vN_int[k,z]* P2GenUnitCapex_int[k] for k=1:K_prod_int,z=1:Z) ) )
	dExpressions_H2["CAPEX_H2G"] = CAPEX_H2G
	dObjective_H2["CAPEX_H2G"] = CAPEX_H2G

	@expression(HY, CAPEX_EL, sum(inputs_H2["discount_factor_P2Gen"][k,z]*vN[k,z]* P2GenUnitCapex[k] for k=1:K_EL,z=1:Z))
	dExpressions_H2["CAPEX_EL"] = CAPEX_EL

	@expression(HY, CAPEX_SMR, CAPEX_H2G - CAPEX_EL)
	dExpressions_H2["CAPEX_SMR"] = CAPEX_SMR

	@expression(HY, Cost_gas, sum(weights[t] *vGas[k,z,t]* GasPrice[k,z,t] for k=1:K_prod,z=1:Z,t=1:T) + sum(weights[t] *vGas_int[k,z,t]* GasPrice_int[k,z,t] for k=1:K_prod_int,z=1:Z,t=1:T))
	dExpressions_H2["Cost_gas"] = Cost_gas
	dObjective_H2["Cost_gas"] = Cost_gas

	@expression(HY, Cost_start, sum(weights[t] *vn_start[k,z,t]* start_up_cost[k,z] for k=1:K_prod,z=1:Z,t=1:T) + sum(weights[t] *vn_start_int[k,z,t]* start_up_cost_int[k,z] for k=1:K_prod_int,z=1:Z,t=1:T))
	dExpressions_H2["Cost_start"] = Cost_start
	dObjective_H2["Cost_start"] = Cost_start


	## End Objective Function Expressions ##

	## Balance Expressions ##
	# H2 Power Consumption balance
	@expression(HY, eH2PowerConsumption[t=1:T, z=1:Z],
	sum(vP2G_int[k,z,t] for k=1:K_prod_int) +
	sum(vP2G[k,z,t] for k=1:K_prod)
	)
	dExpressions_H2["eH2PowerConsumption"] = eH2PowerConsumption
	dBalance_H2["H2PowerBalance"] = dBalance_H2["H2PowerBalance"] + eH2PowerConsumption

	# H2 balance
	@expression(HY, eH2Gen[t=1:T, z=1:Z],
	sum(vH2Gen[k,z,t] for k=1:K_prod) +
	sum(vH2Gen_int[k,z,t] for k=1:K_prod_int)
	)
	dExpressions_H2["eH2Gen"] = eH2Gen
	dBalance_H2["H2Balance"] = dBalance_H2["H2Balance"] + eH2Gen

	# Carbon emission balance
	@expression(HY, Gen_carbon_emission, sum(weights[t] *vH2Gen[k,z,t]* gen_emission_rate[k,z] for k=1:K_prod,z=1:Z,t=1:T) + sum(weights[t] *vH2Gen_int[k,z,t]* gen_emission_rate_int[k,z] for k=1:K_prod_int,z=1:Z,t=1:T))
	dExpressions_H2["Gen_carbon_emission"] = Gen_carbon_emission
	dBalance_H2["H2_Carbon_Emission"] = dBalance_H2["H2_Carbon_Emission"] + Gen_carbon_emission
	## End Balance Expressions ##
	### End Expressions ###

	### Constratints ###
	for w in 1:W
	    for h in 1:Tw
	        t = Tw*(w-1)+h
	        tw_min = Tw*(w-1)+1
	        tw_max = Tw*(w-1)+Tw
			for z in 1:Z
        		for k in 1:K_prod
		            @constraint(HY, vP2G[k,z,t] == vH2Gen[k,z,t] * etaP2G[k,z])

		            @constraint(HY, vGas[k,z,t] == vH2Gen[k,z,t] * etaGas[k,z])

		            # Minimum and maximum outputs
		            @constraint(HY, vH2Gen[k,z,t] >= rhoH2Gen_min[k,z] * H2GenSize[k,z] * vn[k,z,t])
		            @constraint(HY, vH2Gen[k,z,t] <= rhoH2Gen_max[k,z] * H2GenSize[k,z] * vn[k,z,t])
		            @constraint(HY, vn[k,z,t] <= vN[k,z])
		            # @constraint(HY, vN[k,z] == vN_existing[k,z] + vN_new[k,z] - vN_retired[k,z])

		            if H2Gen_zone_option[Z_set[z]] == 0
	                    if k > inputs_H2["H2Gen_avai"]
		                	@constraint(HY, vN[k,z] == 0)
						end
		            end

		            @constraint(HY, sum(vn[k,z,t] for t=1:T) >= vN[k,z]*T*min_up_ratio[k,z] )

					# Commitment
					if h > 1
						@constraint(HY, vn[k,z,t] - vn[k,z,t-1] ==  vn_start[k,z,t] - vn_shut[k,z,t])
					else
						@constraint(HY, vn[k,z,t] - vn[k,z,tw_max] ==  vn_start[k,z,t] - vn_shut[k,z,t])
					end

		            # Minimum up and down times
					if (tau_up[k,z]+1) > h >= 1
						@constraint(HY, vn[k,z,t] >=  sum(vn_start[k,z,tt] for tt = (t-(h-1):t)) + sum(vn_start[k,z,tt] for tt=(tw_max-(tau_up[k,z]-h)):tw_max))
					else
						@constraint(HY, vn[k,z,t] >=  sum(vn_start[k,z,tt] for tt = (t-tau_up[k,z]:t)) )
					end


        		end

		        if gen_category == 1
		            for k in 1:K_prod_int
		                @constraint(HY, vP2G_int[k,z,t] == vH2Gen_int[k,z,t] * etaP2G_int[k,z])

		                @constraint(HY, vGas_int[k,z,t] == vH2Gen_int[k,z,t] * etaGas_int[k,z])

		                # Minimum and maximum outputs
		                @constraint(HY, vH2Gen_int[k,z,t] >= rhoH2Gen_min_int[k,z] * H2GenSize_int[k,z] * vn_int[k,z,t])
		                @constraint(HY, vH2Gen_int[k,z,t] <= rhoH2Gen_max_int[k,z] * H2GenSize_int[k,z] * vn_int[k,z,t])
		                @constraint(HY, vn_int[k,z,t] <= vN_int[k,z])
		                # @constraint(HY, vN[k,z] == vN_existing[k,z] + vN_new[k,z] - vN_retired[k,z])

		                if H2Gen_int_zone_option[Z_set[z]] == 0
		                    @constraint(HY, vN_int[k,z] == 0)
		                end

		                @constraint(HY, sum(vn_int[k,z,t] for t=1:T) >= vN_int[k,z]*T*min_up_ratio_int[k,z] )

						# Commitment
						if h > 1
							@constraint(HY, vn_int[k,z,t] - vn_int[k,z,t-1] ==  vn_start_int[k,z,t] - vn_shut_int[k,z,t])
						else
							@constraint(HY, vn_int[k,z,t] - vn_int[k,z,tw_max] ==  vn_start_int[k,z,t] - vn_shut_int[k,z,t])
						end

		                # Minimum up and down times
		                if t > 1 && t > tau_up_int[k,z]
		                @constraint(HY, vn_int[k,z,t] >= vn_int[k,z,t-1] - vn_int[k,z,t - tau_up_int[k,z]])
		                end
		                # if t > 1 && t > tau_down[k,z]
		                # @constraint(HY, 1 - vn[k,z,t] >= vn[k,z,t - tau_down[k,z]] - vn[k,z,t-1])
		                # end
		            end
				end
	        end
	    end
	end

	### End Constratints ###

	dModuleArgs["dExpressions_H2"] = dExpressions_H2
	dModuleArgs["dObjective_H2"] = dObjective_H2
	dModuleArgs["dBalance_H2"] = dBalance_H2
	return HY, dModuleArgs
end # end H2_Gen module

# H2 to Power Module
function H2_PowerGen(HY::Model, dModuleArgs::Dict)

	println("Hydrogen to Power Module")

	inputs_H2 = dModuleArgs["inputs_H2"]
	setup = dModuleArgs["setup"]
	dExpressions_H2 = dModuleArgs["dExpressions_H2"]
	dObjective_H2 = dModuleArgs["dObjective_H2"]
	dBalance_H2 = dModuleArgs["dBalance_H2"]

	Z_total = inputs_H2["Z_total"]
    PT_total = inputs_H2["PT_total"]
    RT_total = inputs_H2["RT_total"]
    K_prod_total = inputs_H2["K_prod_total"]
    K_stor_total = inputs_H2["K_stor_total"]
    G_H2_total = inputs_H2["G_H2_total"]
    T = inputs_H2["T"]
    Z = inputs_H2["Z"]
    PT = inputs_H2["PT"]
    RT = inputs_H2["RT"]
    K_prod = inputs_H2["K_prod"]
    K_prod_int = inputs_H2["K_prod_int"]
    K_stor = inputs_H2["K_stor"]
    G_H2 = inputs_H2["G_H2"]
    W = inputs_H2["W"]
    Tw = inputs_H2["Tw"]
	weights = inputs_H2["weights"]

	no_integer = inputs_H2["no_integer"]
	etaPfG = inputs_H2["etaPfG"]
    rhoPfG_min = inputs_H2["rhoPfG_min"]
    rhoPfG_max = inputs_H2["rhoPfG_max"]
    PfGSize = inputs_H2["PfGSize"]
    PfGUnitCapex = inputs_H2["PfGUnitCapex"]
    tau_up_g = inputs_H2["tau_up_g"]
    tau_down_g = inputs_H2["tau_down_g"]

	PfGVOM = inputs_H2["PfGVOM"]
	PfGFOM = inputs_H2["PfGFOM"]

	H2_to_power_option = inputs_H2["H2_to_power_option"]
	H2_to_power_cost_factor = inputs_H2["H2_to_power_cost_factor"]
	### Variables ###
    @variable(HY, vPfG[g=1:G_H2,z=1:Z, t = 1:T] >= 0 )
    @variable(HY, vH2P[g=1:G_H2,z=1:Z, t = 1:T] >= 0 )
    if no_integer == 1
        @variable(HY, vm[g=1:G_H2,z=1:Z, t = 1:T] >= 0 )
        @variable(HY, vM[g=1:G_H2,z=1:Z] >= 0 )

        # @variable(HY, vM_new[g=1:G_H2,z=1:Z] >= 0 )
        # @variable(HY, vM_existing[g=1:G_H2,z=1:Z] >= 0, Int )
        # @variable(HY, vM_retired[g=1:G_H2,z=1:Z] >= 0, Int )
        # @variable(HY, 0 >= vM_existing[g=1:G_H2,z=1:Z] >= 0 )
        # @variable(HY, 0 >= vM_retired[g=1:G_H2,z=1:Z] >= 0 )
    else
        @variable(HY, vm[g=1:G_H2,z=1:Z, t = 1:T] >= 0, Int )
        @variable(HY, vM[g=1:G_H2,z=1:Z] >= 0, Int )

        # @variable(HY, vM_new[g=1:G_H2,z=1:Z] >= 0, Int )
        # @variable(HY, vM_existing[g=1:G_H2,z=1:Z] >= 0, Int )
        # @variable(HY, vM_retired[g=1:G_H2,z=1:Z] >= 0, Int )
        # @variable(HY, 0 >= vM_existing[g=1:G_H2,z=1:Z] >= 0, Int )
        # @variable(HY, 0 >= vM_retired[g=1:G_H2,z=1:Z] >= 0, Int )
    end

	if H2_to_power_option == 0
		setupperbound.(vm, 0)
		setupperbound.(vM, 0)
		# setupperbound.(vM_new, 0)
		setupperbound.(vPfG, 0)
		setupperbound.(vH2P, 0)
	end

	### Expressions ###
	## Objective Function Expressions ##

	@expression(HY, CAPEX_PfG, H2_to_power_cost_factor*sum(inputs_H2["discount_factor_PfG"][g,z] * vM[g,z]* PfGUnitCapex[g,z] for g=1:G_H2,z=1:Z))
	dExpressions_H2["CAPEX_PfG"] = CAPEX_PfG
	dObjective_H2["CAPEX_PfG"] = CAPEX_PfG

	@expression(HY, OPEX_PfG, H2_to_power_cost_factor*(sum(PfGFOM[g,z]* PfGSize[g,z]*vM[g,z] for g=1:G_H2,z=1:Z) + sum(weights[t]*PfGVOM[g,z]*vPfG[g,z,t] for g=1:G_H2,z=1:Z,t=1:T )  ))
	dExpressions_H2["OPEX_PfG"] = OPEX_PfG
	dObjective_H2["OPEX_PfG"] = OPEX_PfG

	## End Objective Function Expressions ##

	## Balance Expressions ##
	# H2 Power Generation balance
	@expression(HY, eH2PowerGen[t=1:T, z=1:Z], sum(vPfG[g,z,t] for g=1:G_H2))
	dExpressions_H2["eH2PowerGen"] = eH2PowerGen
	dBalance_H2["H2PowerBalance"] = dBalance_H2["H2PowerBalance"] - eH2PowerGen

	# H2 balance
	@expression(HY, eH2Consumption[t=1:T, z=1:Z],
	sum(vH2P[g,z,t] for g=1:G_H2))
	dExpressions_H2["eH2Consumption"] = eH2Consumption
	dBalance_H2["H2Balance"] = dBalance_H2["H2Balance"] - eH2Consumption

	## End Balance Expressions ##
	### End Expressions ###

	### Constratints ###
	for z in 1:Z
        for t in 1:T
            for g in 1:G_H2
            @constraint(HY, vH2P[g,z,t] == vPfG[g,z,t] / etaPfG[g,z])
            # Minimum and maximum outputs
            @constraint(HY, vPfG[g,z,t] >= rhoPfG_min[g,z] * PfGSize[g,z] * vm[g,z,t])
            @constraint(HY, vPfG[g,z,t] <= rhoPfG_max[g,z] * PfGSize[g,z] * vm[g,z,t])
            @constraint(HY, vm[g,z,t] <= vM[g,z])
            # @constraint(HY, vM[g,z] == vM_existing[g,z] + vM_new[g,z] - vM_retired[g,z])

            # Minimum up and down times
            # if t > 1 && t > tau_up_g[g,z]
            # @constraint(HY, vm[g,z,t] >= vm[g,z,t-1] - vm[g,z,t - tau_up_g[g,z]])
            # end
            # if t > 1 && t > tau_down_g[g,z]
            # @constraint(HY, 1 - vm[g,z,t] >= vm[g,z,t - tau_down_g[g,z]] - vm[g,z,t-1])
            # end

            end
        end
    end
	### End Constratints ###

	dModuleArgs["dExpressions_H2"] = dExpressions_H2
	dModuleArgs["dObjective_H2"] = dObjective_H2
	dModuleArgs["dBalance_H2"] = dBalance_H2
	return HY, dModuleArgs
end # end H2Power module

# H2 Demand Module
function H2_Demand(HY::Model, dModuleArgs::Dict)

	println("Hydrogen Demand Module")
	inputs_H2 = dModuleArgs["inputs_H2"]
	setup = dModuleArgs["setup"]
	dExpressions_H2 = dModuleArgs["dExpressions_H2"]
	dObjective_H2 = dModuleArgs["dObjective_H2"]
	dBalance_H2 = dModuleArgs["dBalance_H2"]


	Z_total = inputs_H2["Z_total"]
    PT_total = inputs_H2["PT_total"]
    RT_total = inputs_H2["RT_total"]
    K_prod_total = inputs_H2["K_prod_total"]
    K_stor_total = inputs_H2["K_stor_total"]
    G_H2_total = inputs_H2["G_H2_total"]

    T = inputs_H2["T"]
    Z = inputs_H2["Z"]
    PT = inputs_H2["PT"]
    RT = inputs_H2["RT"]

    K_prod = inputs_H2["K_prod"]
    K_prod_int = inputs_H2["K_prod_int"]
    K_stor = inputs_H2["K_stor"]

    G_H2 = inputs_H2["G_H2"]
    W = inputs_H2["W"]
    Tw = inputs_H2["Tw"]

	H2T = inputs_H2["H2T"]
	H2H_min = inputs_H2["H2H_min"]
    H2H_max = inputs_H2["H2H_max"]
    H2I_min = inputs_H2["H2I_min"]
    H2I_max = inputs_H2["H2I_max"]
	weights = inputs_H2["weights"]
	H2UnmetPrice = inputs_H2["H2UnmetPrice"]

	### Variables ###
	@variable(HY, vH2H[z=1:Z, t = 1:T] >= 0 )
	@variable(HY, vH2I[z=1:Z, t = 1:T] >= 0 )
	@variable(HY, vH2D[z=1:Z, t = 1:T] >= 0 )
	@variable(HY, vH2DUnmet[z=1:Z, t = 1:T] >= 0 )
	@variable(HY, vH2Curtail[z=1:Z, t = 1:T] >= 0 )
	# @variable(HY, vH2D_Slack[z=1:Z, t = 1:T] >= 0 )

	## Balance Expressions ##
	# H2 balance
	@expression(HY, Total_H2_Demand, sum(weights[t]*vH2D[z,t] for z=1:Z,t = 1:T))
	dExpressions_H2["Total_H2_Demand"] = Total_H2_Demand + 0.0001

	@expression(HY, eH2D[t=1:T, z=1:Z],	vH2D[z,t])
	dExpressions_H2["eH2D"] = eH2D
	dBalance_H2["H2Balance"] = dBalance_H2["H2Balance"] - eH2D

	### Constratints ###
	for z in 1:Z
	    for t in 1:T
	        @constraint(HY, vH2D[z,t] == H2T[z,t] + vH2H[z,t] + vH2I[z,t] - vH2DUnmet[z,t] + vH2Curtail[z,t])
	        @constraint(HY, H2H_min[z,t] <= vH2H[z,t] )
	        @constraint(HY, vH2H[z,t] <= H2H_max[z,t])
	        @constraint(HY, H2I_min[z,t] <= vH2I[z,t])
	        @constraint(HY, vH2I[z,t] <= H2I_max[z,t])
	    end
	end

	### Expressions ###
	## Objective Function Expressions ##

	@expression(HY, Cost_unmet_H2, sum(weights[t] *vH2DUnmet[z,t]* H2UnmetPrice for z=1:Z,t=1:T))
	dExpressions_H2["Cost_unmet_H2"] = Cost_unmet_H2
	dObjective_H2["Cost_unmet_H2"] = Cost_unmet_H2

	## End Objective Function Expressions ##

	### End Constratints ###
	dModuleArgs["dExpressions_H2"] = dExpressions_H2
	dModuleArgs["dObjective_H2"] = dObjective_H2
	dModuleArgs["dBalance_H2"] = dBalance_H2
	return HY, dModuleArgs
end # end H2_Demand module

# H2 Storage Module
function H2_Storage(HY::Model, dModuleArgs::Dict)

	println("Hydrogen Storage Module")

	inputs_H2 = dModuleArgs["inputs_H2"]
	inputs = dModuleArgs["inputs"]

	setup = dModuleArgs["setup"]
	dExpressions_H2 = dModuleArgs["dExpressions_H2"]
	dObjective_H2 = dModuleArgs["dObjective_H2"]
	dBalance_H2 = dModuleArgs["dBalance_H2"]

	Z_total = inputs_H2["Z_total"]
    PT_total = inputs_H2["PT_total"]
    RT_total = inputs_H2["RT_total"]
    K_prod_total = inputs_H2["K_prod_total"]
    K_stor_total = inputs_H2["K_stor_total"]
    G_H2_total = inputs_H2["G_H2_total"]
    T = inputs_H2["T"]
    Z = inputs_H2["Z"]
    PT = inputs_H2["PT"]
    RT = inputs_H2["RT"]
    K_prod = inputs_H2["K_prod"]
    K_prod_int = inputs_H2["K_prod_int"]
    K_stor = inputs_H2["K_stor"]
    G_H2 = inputs_H2["G_H2"]
    W = inputs_H2["W"]
    Tw = inputs_H2["Tw"]

	etaStor = inputs_H2["etaStor"]
    rhoH2Stor_min = inputs_H2["rhoH2Stor_min"]
    H2StorUnitCapex = inputs_H2["H2StorUnitCapex"]
    H2StorCap_min = inputs_H2["H2StorCap_min"]
    H2StorCap_max = inputs_H2["H2StorCap_max"]
    H2StorCap_avai = inputs_H2["H2StorCap_avai"]
	H2StorCompressionEnergy = inputs_H2["H2StorCompressionEnergy"]
	H2StorCompressionUnitCapex = inputs_H2["H2StorCompressionUnitCapex"]

	### Variables ###
	@variable(HY, vH2StorDis[k=1:K_stor,z=1:Z, t = 1:T] >= 0 )
	@variable(HY, vH2StorCha[k=1:K_stor,z=1:Z, t = 1:T] >= 0 )
	@variable(HY, vH2StorEnergy[k=1:K_stor, z=1:Z, t = 1:T] >= 0 )
	@variable(HY, vH2StorCap[k=1:K_stor, z=1:Z] >= 0 )
	@variable(HY, vH2StorRate[k=1:K_stor, z=1:Z] >= 0 )


	if setup["LDS"]==1
		NPeriods = size(inputs["Period_Map"])[1] # Number of modeled periods
		# State of charge of storage at beginning of each modeled period n
		@variable(HY, vH2SOCw[k=1:K_stor, z=1:Z,n=1:NPeriods] >= 0)
		# Build up in storage inventory over each representative period w
		# Build up inventory can be positive or negative
		@variable(HY, vH2dSOC[k=1:K_stor, z=1:Z,w=1:W])
	else
		NPeriods = size(inputs["Period_Map"])[1]
		@variable(HY, 0 >= vH2SOCw[k=1:K_stor, z=1:Z,n=1:NPeriods] >= 0)
		@variable(HY, 0 >= vH2dSOC[k=1:K_stor, z=1:Z,w=1:W] >= 0)
	end

	### Expressions ###
	## Objective Function Expressions ##
	discount_factor = dExpressions_H2["discount_factor"]

	@expression(HY, CAPEX_Stor, sum(inputs_H2["discount_factor_H2Stor"][k,z] * vH2StorCap[k,z] * H2StorUnitCapex[k] for k=1:K_stor,z=1:Z))
	dExpressions_H2["CAPEX_Stor"] = CAPEX_Stor
	dObjective_H2["CAPEX_Stor"] = CAPEX_Stor

	@expression(HY, CAPEX_UG, sum(inputs_H2["discount_factor_H2Stor"][k,z] * vH2StorCap[k,z] * H2StorUnitCapex[k] for k=1,z=1:Z))
	dExpressions_H2["CAPEX_UG"] = CAPEX_UG
	@expression(HY, CAPEX_AG, sum(inputs_H2["discount_factor_H2Stor"][k,z] * vH2StorCap[k,z] * H2StorUnitCapex[k] for k=2,z=1:Z))
	dExpressions_H2["CAPEX_AG"] = CAPEX_AG

	@expression(HY, CAPEX_Compression_Stor, sum(inputs_H2["discount_factor_H2Compression"] * vH2StorRate[k,z] * H2StorCompressionUnitCapex[k] for k=1:K_stor,z=1:Z))
	dExpressions_H2["CAPEX_Compression_Stor"] = CAPEX_Compression_Stor
	dObjective_H2["CAPEX_Compression_Stor"] = CAPEX_Compression_Stor

	## End Objective Function Expressions ##

	## Balance Expressions ##
	# H2 Power Consumption balance
	@expression(HY, eH2StorCompressionPowerConsumption[t=1:T, z=1:Z],
	sum(vH2StorCha[k,z,t] * H2StorCompressionEnergy[k] for k = 1:K_stor ) )
	dExpressions_H2["eH2StorCompressionPowerConsumption"] = eH2StorCompressionPowerConsumption
	dBalance_H2["H2PowerBalance"] = dBalance_H2["H2PowerBalance"] + eH2StorCompressionPowerConsumption

	# H2 balance
	@expression(HY, eH2Stor[t=1:T, z=1:Z],
	sum(vH2StorDis[k,z,t] for k=1:K_stor) -
	sum(vH2StorCha[k,z,t] for k=1:K_stor)
	)
	dExpressions_H2["eH2Stor"] = eH2Stor
	dBalance_H2["H2Balance"] = dBalance_H2["H2Balance"] + eH2Stor

	## End Balance Expressions ##
	### End Expressions ###

	### Constraints ###

	# @constraintref SoCBal[1:T,1:G]
	#


	if setup["LDS"]==1
		dfPeriodMap = inputs["Period_Map"] # Dataframe that maps modeled periods to representative periods
		NPeriods = size(inputs["Period_Map"])[1]
		for n in 1:NPeriods # Modeled periods
			for k in 1:K_stor
				for z in 1:Z
					# Storage at beginning of period w = storage at beginning of period w-1+ storage built up in period w (after n representative periods)
					if n<NPeriods  #Multiply storage build up term from prior period with corresponding weight
						@constraint(HY, vH2SOCw[k,z,n+1] == vH2SOCw[k,z,n] + vH2dSOC[k,z,dfPeriodMap[!,:RepPeriod_index][n]])
					else # Last period is linked to first period
						@constraint(HY,vH2SOCw[k,z,1] == vH2SOCw[k,z,n] + vH2dSOC[k,z,dfPeriodMap[!,:RepPeriod_index][n]])
					end

					# Storage at beginning of each modeled period cannot exceed installed energy capacity
					@constraint(HY, vH2SOCw[k,z,n] <= vH2StorCap[k,z])
				end
			end
		end
	end


	for w in 1:W
	    for h in 1:Tw
	        t = Tw*(w-1)+h
	        tw_min = Tw*(w-1)+1
	        tw_max = Tw*(w-1)+Tw

			for z in 1:Z
	    # for t in 1:T
		        for k in 1:K_stor
		            if H2StorCap_avai[k,z] == 0
		                @constraint(HY, vH2StorCap[k,z] == 0)
		                @constraint(HY, vH2StorRate[k,z] == 0)
		            end

		            @constraint(HY, vH2StorCha[k,z,t] <=  vH2StorRate[k,z])
		            # @constraint(HY, vH2StorDis[k,z,t] <=  vH2StorRate[k,z])

		            if h > 1
			            @constraint(HY, vH2StorEnergy[k,z,t] ==  vH2StorEnergy[k,z,t-1] + vH2StorCha[k,z,t] * etaStor[k,z] - vH2StorDis[k,z,t] / etaStor[k,z])
		            else
		                @constraint(HY, vH2StorEnergy[k,z,t] ==  vH2StorEnergy[k,z,tw_max]
						- vH2dSOC[k,z,w]
						 + vH2StorCha[k,z,t] * etaStor[k,z] - vH2StorDis[k,z,t] / etaStor[k,z])
		                # @constraint(HY, vH2StorEnergy[k,z,t] ==  rhoH2Stor_min[k,z] * vH2StorCap[k,z])
		            end

		            @constraint(HY, rhoH2Stor_min[k,z] * vH2StorCap[k,z] <= vH2StorEnergy[k,z,t]  )
		            @constraint(HY, vH2StorEnergy[k,z,t] <= vH2StorCap[k,z]  )
		            @constraint(HY, H2StorCap_min[k,z]*H2StorCap_avai[k,z] <= vH2StorCap[k,z]  )
		            @constraint(HY, vH2StorCap[k,z] <= H2StorCap_max[k,z]*H2StorCap_avai[k,z] )
		        end
		    end
		end
	end
	### End Constraints ###
	dModuleArgs["dExpressions_H2"] = dExpressions_H2
	dModuleArgs["dObjective_H2"] = dObjective_H2
	dModuleArgs["dBalance_H2"] = dBalance_H2
	return HY, dModuleArgs
end # end H2_Storage module

function H2Pipeline(HY::Model, dModuleArgs::Dict)

	println("Hydrogen Pipeline Module")

	inputs_H2 = dModuleArgs["inputs_H2"]
	setup = dModuleArgs["setup"]
	dExpressions_H2 = dModuleArgs["dExpressions_H2"]
	dObjective_H2 = dModuleArgs["dObjective_H2"]
	dBalance_H2 = dModuleArgs["dBalance_H2"]


	Z_total = inputs_H2["Z_total"]
    PT_total = inputs_H2["PT_total"]
    RT_total = inputs_H2["RT_total"]
    K_prod_total = inputs_H2["K_prod_total"]
    K_stor_total = inputs_H2["K_stor_total"]
    G_H2_total = inputs_H2["G_H2_total"]

    T = inputs_H2["T"]
    Z = inputs_H2["Z"]
    PT = inputs_H2["PT"]
    RT = inputs_H2["RT"]

    K_prod = inputs_H2["K_prod"]
    K_prod_int = inputs_H2["K_prod_int"]
    K_stor = inputs_H2["K_stor"]

    G_H2 = inputs_H2["G_H2"]
    W = inputs_H2["W"]
    Tw = inputs_H2["Tw"]

	H2PipeUnitCapex = inputs_H2["H2PipeUnitCapex"]
    PipeLength = inputs_H2["PipeLength"]
    H2PipeFlowSize = inputs_H2["H2PipeFlowSize"]
    H2PipeCap = inputs_H2["H2PipeCap"]
    rhoH2PipeCap_min = inputs_H2["rhoH2PipeCap_min"]

	Number_online_compression = inputs_H2["Number_online_compression"]
	H2PipeCompressionUnitCapex = inputs_H2["H2PipeCompressionUnitCapex"]
	H2PipeCompressionUnitOpex = inputs_H2["H2PipeCompressionUnitOpex"]
	H2PipeCompressionEnergy = inputs_H2["H2PipeCompressionEnergy"]
	H2PipeCompressionOnlineUnitCapex = inputs_H2["H2PipeCompressionOnlineUnitCapex"]
	H2PipeCompressionOnlineUnitOpex = inputs_H2["H2PipeCompressionOnlineUnitOpex"]
	H2PipeCompressionOnlineEnergy = inputs_H2["H2PipeCompressionOnlineEnergy"]

	Transmission_cost_factor = inputs_H2["Transmission_cost_factor"]

	Pipe_cost_factor = inputs_H2["Pipe_cost_factor"]

	### Variables ###
	if inputs_H2["Pipe_option_avai"] == 0
		@variable(HY, 0 >= vNPipe[zz=1:Z,z=1:Z, p=1:PT] >= 0 )
    	@variable(HY, 0 >= vH2PipeLevel[zz=1:Z,z=1:Z, p=1:PT, t = 1:T] >= 0 )
    	@variable(HY, vH2PipeFlow[zz=1:Z,z=1:Z, p=1:PT, t = 1:T] )
    	@variable(HY, 0 >= vH2PipeFlow_pos[zz=1:Z,z=1:Z, p=1:PT, t = 1:T] >= 0 )
    	@variable(HY, 0 >= vH2PipeFlow_neg[zz=1:Z,z=1:Z, p=1:PT, t = 1:T] >= 0)
	else
		if inputs_H2["Pipe_integer"] == 1
	    	@variable(HY, vNPipe[zz=1:Z,z=1:Z, p=1:PT] >= 0, int )
		else
			@variable(HY, vNPipe[zz=1:Z,z=1:Z, p=1:PT] >= 0 )
		end
	    @variable(HY, vH2PipeLevel[zz=1:Z,z=1:Z, p=1:PT, t = 1:T] >= 0 )
	    @variable(HY, vH2PipeFlow[zz=1:Z,z=1:Z, p=1:PT, t = 1:T] )
	    @variable(HY, vH2PipeFlow_pos[zz=1:Z,z=1:Z, p=1:PT, t = 1:T] >= 0 )
	    @variable(HY, vH2PipeFlow_neg[zz=1:Z,z=1:Z, p=1:PT, t = 1:T] >= 0)
	end

	### Expressions ###
	## Objective Function Expressions ##
	discount_factor = dExpressions_H2["discount_factor"]

	@expression(HY, CAPEX_Pipe, Transmission_cost_factor * Pipe_cost_factor * sum(inputs_H2["discount_factor_H2Pipe"][p] * 0.5 * vNPipe[zz,z,p] * PipeLength[zz,z] * H2PipeUnitCapex[p] for zz = 1:Z,z=1:Z,p = 1:PT if zz != z))
	dExpressions_H2["CAPEX_Pipe"] = CAPEX_Pipe
	dObjective_H2["CAPEX_Pipe"] = CAPEX_Pipe

	@expression(HY, CAPEX_Compression_Pipe, Transmission_cost_factor * inputs_H2["discount_factor_H2Compression"] *sum(0.5 * vNPipe[zz,z,p] * H2PipeFlowSize[p] * (H2PipeCompressionUnitCapex[p] + Number_online_compression[zz,z] * H2PipeCompressionOnlineUnitCapex[p]) for zz = 1:Z,z=1:Z,p = 1:PT if zz != z))
	dExpressions_H2["CAPEX_Compression_Pipe"] = CAPEX_Compression_Pipe
	dObjective_H2["CAPEX_Compression_Pipe"] = CAPEX_Compression_Pipe

	## End Objective Function Expressions ##

	## Balance Expressions ##
	# H2 Power Consumption balance
	@expression(HY, eH2PipeCompressionPowerConsumption[t=1:T, z=1:Z],
	sum( (vH2PipeFlow_neg[z,zz,p,t] * (H2PipeCompressionEnergy[p] + Number_online_compression[z,zz] * H2PipeCompressionOnlineEnergy[p]) ) for zz = 1:Z,p = 1:PT if zz != z))
	dExpressions_H2["eH2PipeCompressionPowerConsumption"] = eH2PipeCompressionPowerConsumption
	dBalance_H2["H2PowerBalance"] = dBalance_H2["H2PowerBalance"] + eH2PipeCompressionPowerConsumption

	# H2 balance
	@expression(HY, PipeFlow[t=1:T,z=1:Z], sum(vH2PipeFlow[z,zz,p,t] for zz = 1:Z, p = 1:PT if zz != z))
	dExpressions_H2["PipeFlow"] = PipeFlow
	dBalance_H2["H2Balance"] = dBalance_H2["H2Balance"] + PipeFlow

	## End Balance Expressions ##
	### End Expressions ###

	### Constraints ###
	for w in 1:W
	    for h in 1:Tw
	        t = Tw*(w-1)+h
	        tw_min = Tw*(w-1)+1
	        tw_max = Tw*(w-1)+Tw

			for z in 1:Z
			    # for t in 1:T
		        for zz in 1:Z
		            for p in 1:PT
		                @constraint(HY, vNPipe[zz,z,p] * H2PipeFlowSize[p] >= vH2PipeFlow[zz,z,p,t])
		                @constraint(HY, vH2PipeFlow[zz,z,p,t] >= -vNPipe[zz,z,p] * H2PipeFlowSize[p])

		                @constraint(HY, vNPipe[zz,z,p] * H2PipeFlowSize[p] >= vH2PipeFlow_pos[zz,z,p,t])
		                @constraint(HY, vNPipe[zz,z,p] * H2PipeFlowSize[p] >= vH2PipeFlow_neg[zz,z,p,t])
		                @constraint(HY, vH2PipeFlow[zz,z,p,t] == vH2PipeFlow_pos[zz,z,p,t] - vH2PipeFlow_neg[zz,z,p,t] )

		                if z==zz
		                @constraint(HY, vH2PipeFlow[zz,z,p,t] == 0 )
		                end

		                if h > 1
		                @constraint(HY, vH2PipeLevel[zz,z,p,t] == vH2PipeLevel[zz,z,p,t-1] - vH2PipeFlow[zz,z,p,t] - vH2PipeFlow[z,zz,p,t])
		                else
		                    @constraint(HY, vH2PipeLevel[zz,z,p,t] == vH2PipeLevel[zz,z,p,tw_max] - vH2PipeFlow[zz,z,p,t] - vH2PipeFlow[z,zz,p,t])
		                end
		                @constraint(HY, vH2PipeLevel[zz,z,p,t] >= rhoH2PipeCap_min[zz,z,p] * H2PipeCap[p] * vNPipe[zz,z,p])
		                @constraint(HY, H2PipeCap[p] * vNPipe[zz,z,p] >= vH2PipeLevel[zz,z,p,t])

		                @constraint(HY, vNPipe[zz,z,p] == vNPipe[z,zz,p])
		                @constraint(HY, vH2PipeLevel[zz,z,p,t] == vH2PipeLevel[z,zz,p,t])
		                # @constraint(HY, vH2PipeFlowSize[zz,z,p] == vH2PipeFlowSize[z,zz,p])
		            end
		        end
		    end
		end
	end
	### End Constraints ###
	dModuleArgs["dExpressions_H2"] = dExpressions_H2
	dModuleArgs["dObjective_H2"] = dObjective_H2
	dModuleArgs["dBalance_H2"] = dBalance_H2
	return HY, dModuleArgs
end # end H2Pipeline module

function H2Truck(HY::Model, dModuleArgs::Dict)

	println("Hydrogen Truck Module")

	inputs_H2 = dModuleArgs["inputs_H2"]
	inputs = dModuleArgs["inputs"]

	setup = dModuleArgs["setup"]
	dExpressions_H2 = dModuleArgs["dExpressions_H2"]
	dObjective_H2 = dModuleArgs["dObjective_H2"]
	dBalance_H2 = dModuleArgs["dBalance_H2"]

	Z_total = inputs_H2["Z_total"]
    PT_total = inputs_H2["PT_total"]
    RT_total = inputs_H2["RT_total"]
    K_prod_total = inputs_H2["K_prod_total"]
    K_stor_total = inputs_H2["K_stor_total"]
    G_H2_total = inputs_H2["G_H2_total"]
    T = inputs_H2["T"]
    Z = inputs_H2["Z"]
    PT = inputs_H2["PT"]
    RT = inputs_H2["RT"]
    K_prod = inputs_H2["K_prod"]
    K_prod_int = inputs_H2["K_prod_int"]
    K_stor = inputs_H2["K_stor"]
    G_H2 = inputs_H2["G_H2"]
    W = inputs_H2["W"]
    Tw = inputs_H2["Tw"]

	weights = inputs_H2["weights"]
	RouteLength = inputs_H2["RouteLength"]
	H2TruckUnitCapex = inputs_H2["H2TruckUnitCapex"]
    TruckCap = inputs_H2["TruckCap"]
	Full_weight = inputs_H2["Full_weight_tonne_per_unit"]
	Empty_weight = inputs_H2["Empty_weight_tonne_per_unit"]
    AvgTruckSpeed = inputs_H2["AvgTruckSpeed"]
    H2TruckUnitOpex_full = inputs_H2["H2TruckUnitOpex_full"]
	H2TruckUnitOpex_empty = inputs_H2["H2TruckUnitOpex_empty"]
    TD = inputs_H2["TD"]
    H2TLoss = inputs_H2["H2TLoss"]
    truck_emission_rate = inputs_H2["truck_emission_rate"]
	H2TruckCompressionUnitCapex = inputs_H2["H2TruckCompressionUnitCapex"]
	H2TruckCompressionEnergy = inputs_H2["H2TruckCompressionEnergy"]
	H2TruckCompressionUnitOpex = inputs_H2["H2TruckCompressionUnitOpex"]

	Truck_integer_model = inputs_H2["Truck_integer_model"]
	truck_model_simp = inputs_H2["truck_model_simp"]

	Transmission_cost_factor = inputs_H2["Transmission_cost_factor"]

	### Variables ###
	if setup["LDS"]==1
		NPeriods = size(inputs["Period_Map"])[1]

		# # Number of modeled periods
		# # State of charge of storage at beginning of each modeled period n
		# @variable(HY, vH2TruckSOCw[z=1:Z,r=1:RT, n=1:NPeriods] >= 0)
		# # Build up in storage inventory over each representative period w
		# # Build up inventory can be positive or negative
		# @variable(HY, vH2TruckdSOC[z=1:Z,r=1:RT, w=1:W])
		@variable(HY, 0 >= vH2TruckdSOC[z=1:Z,r=1:RT, w=1:W] >= 0)
		@variable(HY, 0 >= vH2TruckSOCw[z=1:Z,r=1:RT, n=1:NPeriods] >= 0)

	else
		@variable(HY, 0 >= vH2TruckdSOC[z=1:Z,r=1:RT, w=1:W] >= 0)
	end

	@variable(HY, vH2TruckFlow[z=1:Z, r=1:RT, t = 1:T] )
	# @variable(HY, vH2TruckFlow_pos[zz=1:Z,z=1:Z, r=1:RT, t = 1:T] >=0 )
	# @variable(HY, vH2TruckFlow_neg[zz=1:Z,z=1:Z, r=1:RT, t = 1:T] >=0 )
	# @variable(HY, vH2TruckLevel[zz=1:Z,z=1:Z, r=1:RT, t = 1:T] >=0 )


	if inputs_H2["Truck_option_avai"] == 0
		@variable(HY, 0 >= vNTruck[r=1:RT] >= 0 )
        @variable(HY, 0 >= vNTruckRoute[r=1:RT] >= 0 )
        @variable(HY, 0 >= vNavail_full[z=1:Z,r=1:RT,t = 1:T] >= 0  )
        @variable(HY, 0 >= vNtravel_full[zz=1:Z,z=1:Z,r=1:RT,t = 1:T] >= 0 )
        @variable(HY, 0 >= vNarrive_full[zz=1:Z,z=1:Z,r=1:RT,t = 1:T] >= 0 )
        @variable(HY, 0 >= vNdepart_full[zz=1:Z,z=1:Z,r=1:RT,t = 1:T] >= 0 )
        @variable(HY, 0 >= vNavail_empty[z=1:Z,r=1:RT,t = 1:T] >= 0  )
        @variable(HY, 0 >= vNtravel_empty[zz=1:Z,z=1:Z,r=1:RT,t = 1:T] >= 0 )
        @variable(HY, 0 >= vNarrive_empty[zz=1:Z,z=1:Z,r=1:RT,t = 1:T] >= 0 )
        @variable(HY, 0 >= vNdepart_empty[zz=1:Z,z=1:Z,r=1:RT,t = 1:T] >= 0 )
        @variable(HY, 0 >= vNcharged[z=1:Z,r=1:RT,t = 1:T] >= 0 )
        @variable(HY, 0 >= vNdischarged[z=1:Z,r=1:RT,t = 1:T] >= 0 )
        @variable(HY, 0 >= vN_full[r=1:RT,t = 1:T] >= 0 )
        @variable(HY, 0 >= vN_empty[r=1:RT,t = 1:T] >= 0 )
        @variable(HY, 0 >= vH2TruckCompressionCap[z=1:Z,r=1:RT] >= 0 )
	else
		@variable(HY, vNTruck[r=1:RT] >= 0 )
        @variable(HY, vNTruckRoute[r=1:RT] >= 0 )
        @variable(HY, vNavail_full[z=1:Z,r=1:RT,t = 1:T] >= 0  )
        @variable(HY, vNtravel_full[zz=1:Z,z=1:Z,r=1:RT,t = 1:T] >= 0 )
        @variable(HY, vNarrive_full[zz=1:Z,z=1:Z,r=1:RT,t = 1:T] >= 0 )
        @variable(HY, vNdepart_full[zz=1:Z,z=1:Z,r=1:RT,t = 1:T] >= 0 )
        @variable(HY, vNavail_empty[z=1:Z,r=1:RT,t = 1:T] >= 0  )
        @variable(HY, vNtravel_empty[zz=1:Z,z=1:Z,r=1:RT,t = 1:T] >= 0 )
        @variable(HY, vNarrive_empty[zz=1:Z,z=1:Z,r=1:RT,t = 1:T] >= 0 )
        @variable(HY, vNdepart_empty[zz=1:Z,z=1:Z,r=1:RT,t = 1:T] >= 0 )
        @variable(HY, vNcharged[z=1:Z,r=1:RT,t = 1:T] >= 0 )
        @variable(HY, vNdischarged[z=1:Z,r=1:RT,t = 1:T] >= 0 )
        @variable(HY, vN_full[r=1:RT,t = 1:T] >= 0 )
        @variable(HY, vN_empty[r=1:RT,t = 1:T] >= 0 )
        @variable(HY, vH2TruckCompressionCap[z=1:Z,r=1:RT] >= 0 )
	end


	### Expressions ###
	## Objective Function Expressions ##
	discount_factor = dExpressions_H2["discount_factor"]

	@expression(HY, CAPEX_Truck, Transmission_cost_factor * sum(inputs_H2["discount_factor_H2Truck"][r] * vNTruck[r] * H2TruckUnitCapex[r] for r = 1:RT))
	dExpressions_H2["CAPEX_Truck"] = CAPEX_Truck
	dObjective_H2["CAPEX_Truck"] = CAPEX_Truck

    @expression(HY, OPEX_Truck, Transmission_cost_factor * sum(weights[t] *(vNarrive_full[zz,z,r,t]*H2TruckUnitOpex_full[r] + vNarrive_empty[zz,z,r,t]*H2TruckUnitOpex_empty[r]) * RouteLength[zz,z]  for zz = 1:Z,z=1:Z,r = 1:RT,t = 1:T if zz != z))

	dExpressions_H2["OPEX_Truck"] = OPEX_Truck
	dObjective_H2["OPEX_Truck"] = OPEX_Truck

	@expression(HY, CAPEX_Compression_Truck, Transmission_cost_factor * inputs_H2["discount_factor_H2Compression"] *sum(vH2TruckCompressionCap[z,r] * H2TruckCompressionUnitCapex[r] for z=1:Z,r = 1:RT))
	dExpressions_H2["CAPEX_Compression_Truck"] = CAPEX_Compression_Truck
	dObjective_H2["CAPEX_Compression_Truck"] = CAPEX_Compression_Truck

	## End Objective Function Expressions ##

	## Balance Expressions ##
	# H2 Power Consumption balance
	@expression(HY, eH2TruckCompressionPowerConsumption[t=1:T, z=1:Z],
	sum(vNcharged[z,r,t] * TruckCap[r]* H2TruckCompressionEnergy[r] for r = 1:RT)
	)
	dExpressions_H2["eH2TruckCompressionPowerConsumption"] = eH2TruckCompressionPowerConsumption
	dBalance_H2["H2PowerBalance"] = dBalance_H2["H2PowerBalance"] + eH2TruckCompressionPowerConsumption

	# H2 balance
	@expression(HY, TruckFlow[t=1:T,z=1:Z], sum(vH2TruckFlow[z,r,t] for r = 1:RT ))
	dExpressions_H2["TruckFlow"] = TruckFlow
	dBalance_H2["H2Balance"] = dBalance_H2["H2Balance"] + TruckFlow

	# Carbon emission balance
	@expression(HY, Truck_carbon_emission, sum(weights[t] * (vNarrive_full[z,zz,r,t]*Full_weight[r]+vNarrive_empty[z,zz,r,t]*Empty_weight[r]) *RouteLength[z,zz] *truck_emission_rate[r,z] for zz = 1:Z, z = 1:Z, r = 1:RT, t = 1:T if zz != z))
	dExpressions_H2["Truck_carbon_emission"] = Truck_carbon_emission
	dBalance_H2["H2_Carbon_Emission"] = dBalance_H2["H2_Carbon_Emission"] + Truck_carbon_emission
	## End Balance Expressions ##
	### End Expressions ###

	### Constraints ###

	if setup["LDS"]==1
		dfPeriodMap = inputs["Period_Map"] # Dataframe that maps modeled periods to representative periods
		NPeriods = size(inputs["Period_Map"])[1]
		for n in 1:NPeriods # Modeled periods
			for r in 1:RT
				for z in 1:Z
					# Storage at beginning of period w = storage at beginning of period w-1+ storage built up in period w (after n representative periods)
					if n<NPeriods  #Multiply storage build up term from prior period with corresponding weight
						@constraint(HY, vH2TruckSOCw[z,r,n+1] == vH2TruckSOCw[z,r,n] + vH2TruckdSOC[z,r,dfPeriodMap[!,:RepPeriod_index][n]])
					else # Last period is linked to first period
						@constraint(HY,vH2TruckSOCw[z,r,1] == vH2TruckSOCw[z,r,n] + vH2TruckdSOC[z,r,dfPeriodMap[!,:RepPeriod_index][n]])
					end

					# Storage at beginning of each modeled period cannot exceed installed energy capacity
					@constraint(HY, vH2TruckSOCw[z,r,n] <= vNTruck[r])
				end
			end
		end
	end

	for w in 1:W
    	for h in 1:Tw
	        t = Tw*(w-1)+h
	        tw_min = Tw*(w-1)+1
	        tw_max = Tw*(w-1)+Tw
			for z in 1:Z
			    # for t in 1:T
		        for zz in 1:Z
		            for r in 1:RT
						# @constraint(HY, vNtravel[zz,z,r,t] + vNtravel[z,zz,r,t] + vNavail[zz,z,r,t] + vNavail[z,zz,r,t] == vNTruckRoute[z,zz,r])


					    ## Total number
					    @constraint(HY,
					    vN_full[r,t] +
					    vN_empty[r,t] +
					    0
					    ==
					    vNTruck[r]
					    )
					    # The number of total full trucks
					    @constraint(HY, vN_full[r,t] ==
					    sum(vNtravel_full[zz,z,r,t] for zz = 1:Z, z = 1:Z if zz != z) +
					    sum(vNavail_full[z,r,t] for z = 1:Z) +
					    0
					    )
					    ## The number of total empty trucks
						@constraint(HY, vN_empty[r,t] ==
					   sum(vNtravel_empty[zz,z,r,t] for zz = 1:Z, z = 1:Z if zz != z) +
					   sum(vNavail_empty[z,r,t] for z = 1:Z) +
					   0
					   )


					    # @constraint(HY, vNdischarged[z,r,t] <= vNavail_full[z,r,t] )
					    # @constraint(HY, vNcharged[zz,z,r,t] <= vNTruckRoute[zz,z,r] )

					    # @constraint(HY, vNTruckRoute[zz,z,r] == vNTruckRoute[z,zz,r])
						@constraint(HY, vNTruckRoute[r] == vNTruck[r] )
					    if z==zz
						    # @constraint(HY, vNTruckRoute[zz,z,r] == 0 )
						    @constraint(HY, vNtravel_full[zz,z,r,t] == 0 )
						    @constraint(HY, vNtravel_empty[zz,z,r,t] == 0 )
							@constraint(HY, vNarrive_full[zz,z,r,t] == 0 )
						    @constraint(HY, vNdepart_full[zz,z,r,t] == 0 )
						    @constraint(HY, vNarrive_empty[zz,z,r,t] == 0 )
						    @constraint(HY, vNdepart_empty[zz,z,r,t] == 0 )
							# @constraint(HY, vH2TruckFlow[zz,z,r,t] == 0 )
					    else
					    	# @constraint(HY, vNTruckRoute[zz,z,r] >=0 )
					    # @constraint(HY, 5000>= vNTruckRoute[zz,z,r] )
					    end

					    t_arrive = 1
					    t_depart = 1

					    # Change of the number of full available trucks
					    if h > 1
						    @constraint(HY, vNavail_full[z,r,t] - vNavail_full[z,r,t-1] ==
						    vNcharged[z,r,t] - vNdischarged[z,r,t] +
							sum(vNarrive_full[zz,z,r,t-t_arrive] for zz = 1:Z if zz != z) -
							sum(vNdepart_full[z,zz,r,t-t_depart] for zz = 1:Z if zz != z) +
						    0
						    )
					    else
					        @constraint(HY, vNavail_full[z,r,t] - vNavail_full[z,r,tw_max] ==
							-vH2TruckdSOC[z,r,w] +
							vNcharged[z,r,t] - vNdischarged[z,r,t] +
					        sum(vNarrive_full[zz,z,r,tw_max] for zz = 1:Z if zz != z) -
					        sum(vNdepart_full[z,zz,r,tw_max] for zz = 1:Z if zz != z) +
						    0
							)
					    end
					    # Change of the number of empty available trucks
					    if h > 1
						    @constraint(HY, vNavail_empty[z,r,t] - vNavail_empty[z,r,t-1] ==
							-vNcharged[z,r,t] + vNdischarged[z,r,t] +
						    sum(vNarrive_empty[zz,z,r,t-t_arrive] for zz = 1:Z if zz != z) -
						    sum(vNdepart_empty[z,zz,r,t-t_depart] for zz = 1:Z if zz != z)
							)
					    else
					        @constraint(HY, vNavail_empty[z,r,t] - vNavail_empty[z,r,tw_max] ==
							vH2TruckdSOC[z,r,w]
					        -vNcharged[z,r,t] + vNdischarged[z,r,t] +
					        sum(vNarrive_empty[zz,z,r,tw_max] for zz = 1:Z if zz != z) -
					        sum(vNdepart_empty[z,zz,r,tw_max] for zz = 1:Z if zz != z)
							)
					    end

					    # Change of the number of full traveling trucks
					    if h > 1
					    @constraint(HY, vNtravel_full[z,zz,r,t] - vNtravel_full[z,zz,r,t-1] ==
					    vNdepart_full[z,zz,r,t-t_depart] -
					    vNarrive_full[z,zz,r,t-t_arrive]
					    )
					    else
					        @constraint(HY, vNtravel_full[z,zz,r,t] - vNtravel_full[z,zz,r,tw_max] ==
					        vNdepart_full[z,zz,r,tw_max]
					        - vNarrive_full[z,zz,r,tw_max])
					    end

					    # Change of the number of empty traveling trucks
					    if h > 1
						    @constraint(HY, vNtravel_empty[z,zz,r,t] - vNtravel_empty[z,zz,r,t-1] ==
						    vNdepart_empty[z,zz,r,t-t_depart] -
						    vNarrive_empty[z,zz,r,t-t_arrive])
					    else
					        @constraint(HY, vNtravel_empty[z,zz,r,t] - vNtravel_empty[z,zz,r,tw_max] ==
					        vNdepart_empty[z,zz,r,tw_max] -
					        vNarrive_empty[z,zz,r,tw_max])
					    end


					    # Travel delay

					    if t-TD[zz,z]+1 >= tw_min && t-TD[zz,z]+1 <= tw_max && t-TD[zz,z]+1 <= t
					    @constraint(HY, vNtravel_full[zz,z,r,t] >= sum(vNdepart_full[zz,z,r,tt] for tt = (t-TD[zz,z]+1):t) ) # deaprt from zz to z
					    @constraint(HY, vNtravel_empty[zz,z,r,t] >= sum(vNdepart_empty[zz,z,r,tt] for tt = (t-TD[zz,z]+1):t) )
					    end
					    if t+TD[zz,z] >= tw_min && t+TD[zz,z] <= tw_max && t+1 <= t+TD[zz,z]
					    @constraint(HY, vNtravel_full[zz,z,r,t] >= sum(vNarrive_full[zz,z,r,tt] for tt = (t+1):(t+TD[zz,z]) ) ) # arrive to z to zz
					    @constraint(HY, vNtravel_empty[zz,z,r,t] >= sum(vNarrive_empty[zz,z,r,tt] for tt = (t+1):(t+TD[zz,z]) ) )
					    end

					    @constraint(HY, vH2TruckFlow[z,r,t]  == vNdischarged[z,r,t] * TruckCap[r] * (1 - H2TLoss[r]) -
					    vNcharged[z,r,t] * TruckCap[r])

						# @constraint(HY, vH2TruckFlow[z,zz,r,t]  == vH2TruckFlow_pos[z,zz,r,t] -
                        # vH2TruckFlow_neg[z,zz,r,t]  +
                        # 0
                        # )

					    @constraint(HY, vNcharged[z,r,t] * TruckCap[r] <= vH2TruckCompressionCap[z,r]
					    )



					end
				end
			end
		end
	end
	### End Constraints ###

	dModuleArgs["dExpressions_H2"] = dExpressions_H2
	dModuleArgs["dObjective_H2"] = dObjective_H2
	dModuleArgs["dBalance_H2"] = dBalance_H2
	return HY, dModuleArgs
end # end H2Truck module

function H2Truck_old(HY::Model, dModuleArgs::Dict)

	println("Hydrogen Truck Module")

	inputs_H2 = dModuleArgs["inputs_H2"]
	setup = dModuleArgs["setup"]
	dExpressions_H2 = dModuleArgs["dExpressions_H2"]
	dObjective_H2 = dModuleArgs["dObjective_H2"]
	dBalance_H2 = dModuleArgs["dBalance_H2"]

	Z_total = inputs_H2["Z_total"]
    PT_total = inputs_H2["PT_total"]
    RT_total = inputs_H2["RT_total"]
    K_prod_total = inputs_H2["K_prod_total"]
    K_stor_total = inputs_H2["K_stor_total"]
    G_H2_total = inputs_H2["G_H2_total"]
    T = inputs_H2["T"]
    Z = inputs_H2["Z"]
    PT = inputs_H2["PT"]
    RT = inputs_H2["RT"]
    K_prod = inputs_H2["K_prod"]
    K_prod_int = inputs_H2["K_prod_int"]
    K_stor = inputs_H2["K_stor"]
    G_H2 = inputs_H2["G_H2"]
    W = inputs_H2["W"]
    Tw = inputs_H2["Tw"]

	weights = inputs_H2["weights"]
	RouteLength = inputs_H2["RouteLength"]
	H2TruckUnitCapex = inputs_H2["H2TruckUnitCapex"]
    TruckCap = inputs_H2["TruckCap"]
	Full_weight = inputs_H2["Full_weight_tonne_per_unit"]
	Empty_weight = inputs_H2["Empty_weight_tonne_per_unit"]
    AvgTruckSpeed = inputs_H2["AvgTruckSpeed"]
	H2TruckUnitOpex_full = inputs_H2["H2TruckUnitOpex_full"]
	H2TruckUnitOpex_empty = inputs_H2["H2TruckUnitOpex_empty"]
    TD = inputs_H2["TD"]
    H2TLoss = inputs_H2["H2TLoss"]
    truck_emission_rate = inputs_H2["truck_emission_rate"]
	H2TruckCompressionUnitCapex = inputs_H2["H2TruckCompressionUnitCapex"]
	H2TruckCompressionEnergy = inputs_H2["H2TruckCompressionEnergy"]
	H2TruckCompressionUnitOpex = inputs_H2["H2TruckCompressionUnitOpex"]

	Truck_integer_model = inputs_H2["Truck_integer_model"]
	truck_model_simp = inputs_H2["truck_model_simp"]

	Transmission_cost_factor = inputs_H2["Transmission_cost_factor"]

	### Variables ###

	@variable(HY, vH2TruckFlow[zz=1:Z,z=1:Z, r=1:RT, t = 1:T] )
	# @variable(HY, vH2TruckFlow_pos[zz=1:Z,z=1:Z, r=1:RT, t = 1:T] >=0 )
	# @variable(HY, vH2TruckFlow_neg[zz=1:Z,z=1:Z, r=1:RT, t = 1:T] >=0 )
	# @variable(HY, vH2TruckLevel[zz=1:Z,z=1:Z, r=1:RT, t = 1:T] >=0 )

    if inputs_H2["Truck_option_avai"] == 0
		@variable(HY, 0 >= vNTruck[r=1:RT] >= 0 )
        @variable(HY, 0 >= vNTruckRoute[zz=1:Z,z=1:Z,r=1:RT] >= 0 )
        @variable(HY, 0 >= vNavail_full[zz=1:Z,z=1:Z,r=1:RT,t = 1:T] >= 0  )
        @variable(HY, 0 >= vNtravel_full[zz=1:Z,z=1:Z,r=1:RT,t = 1:T] >= 0 )
        @variable(HY, 0 >= vNarrive_full[zz=1:Z,z=1:Z,r=1:RT,t = 1:T] >= 0 )
        @variable(HY, 0 >= vNdepart_full[zz=1:Z,z=1:Z,r=1:RT,t = 1:T] >= 0 )
        @variable(HY, 0 >= vNavail_empty[zz=1:Z,z=1:Z,r=1:RT,t = 1:T] >= 0  )
        @variable(HY, 0 >= vNtravel_empty[zz=1:Z,z=1:Z,r=1:RT,t = 1:T] >= 0 )
        @variable(HY, 0 >= vNarrive_empty[zz=1:Z,z=1:Z,r=1:RT,t = 1:T] >= 0 )
        @variable(HY, 0 >= vNdepart_empty[zz=1:Z,z=1:Z,r=1:RT,t = 1:T] >= 0 )
        @variable(HY, 0 >= vNcharged[zz=1:Z,z=1:Z,r=1:RT,t = 1:T] >= 0 )
        @variable(HY, 0 >= vNdischarged[zz=1:Z,z=1:Z,r=1:RT,t = 1:T] >= 0 )
        @variable(HY, 0 >= vN_full[zz=1:Z,z=1:Z,r=1:RT,t = 1:T] >= 0 )
        @variable(HY, 0 >= vN_empty[zz=1:Z,z=1:Z,r=1:RT,t = 1:T] >= 0 )

        @variable(HY, 0 >= vH2TruckCompressionCap[z=1:Z,r=1:RT] >= 0 )
    else
		@variable(HY, vNTruck[r=1:RT] >= 0 )
        @variable(HY, vNTruckRoute[zz=1:Z,z=1:Z,r=1:RT] >= 0 )
        @variable(HY, vNavail_full[zz=1:Z,z=1:Z,r=1:RT,t = 1:T] >= 0  )
        @variable(HY, vNtravel_full[zz=1:Z,z=1:Z,r=1:RT,t = 1:T] >= 0 )
        @variable(HY, vNarrive_full[zz=1:Z,z=1:Z,r=1:RT,t = 1:T] >= 0 )
        @variable(HY, vNdepart_full[zz=1:Z,z=1:Z,r=1:RT,t = 1:T] >= 0 )
        @variable(HY, vNavail_empty[zz=1:Z,z=1:Z,r=1:RT,t = 1:T] >= 0  )
        @variable(HY, vNtravel_empty[zz=1:Z,z=1:Z,r=1:RT,t = 1:T] >= 0 )
        @variable(HY, vNarrive_empty[zz=1:Z,z=1:Z,r=1:RT,t = 1:T] >= 0 )
        @variable(HY, vNdepart_empty[zz=1:Z,z=1:Z,r=1:RT,t = 1:T] >= 0 )
        @variable(HY, vNcharged[zz=1:Z,z=1:Z,r=1:RT,t = 1:T] >= 0 )
        @variable(HY, vNdischarged[zz=1:Z,z=1:Z,r=1:RT,t = 1:T] >= 0 )
        @variable(HY, vN_full[zz=1:Z,z=1:Z,r=1:RT,t = 1:T] >= 0 )
        @variable(HY, vN_empty[zz=1:Z,z=1:Z,r=1:RT,t = 1:T] >= 0 )

        @variable(HY, vH2TruckCompressionCap[z=1:Z,r=1:RT] >= 0 )
    end

	### Expressions ###
	## Objective Function Expressions ##
	discount_factor = dExpressions_H2["discount_factor"]

	@expression(HY, CAPEX_Truck, Transmission_cost_factor * sum(inputs_H2["discount_factor_H2Truck"][r] * vNTruck[r] * H2TruckUnitCapex[r] for r = 1:RT))
	dExpressions_H2["CAPEX_Truck"] = CAPEX_Truck
	dObjective_H2["CAPEX_Truck"] = CAPEX_Truck

    @expression(HY, OPEX_Truck, Transmission_cost_factor * sum(weights[t] *(vNarrive_full[zz,z,r,t]*H2TruckUnitOpex_full[r] + vNarrive_empty[zz,z,r,t]*H2TruckUnitOpex_empty[r]) * RouteLength[zz,z]  for zz = 1:Z,z=1:Z,r = 1:RT,t = 1:T if zz != z))

	dExpressions_H2["OPEX_Truck"] = OPEX_Truck
	dObjective_H2["OPEX_Truck"] = OPEX_Truck

	@expression(HY, CAPEX_Compression_Truck, Transmission_cost_factor * inputs_H2["discount_factor_H2Compression"] *sum(vH2TruckCompressionCap[z,r] * H2TruckCompressionUnitCapex[r] for z=1:Z,r = 1:RT))
	dExpressions_H2["CAPEX_Compression_Truck"] = CAPEX_Compression_Truck
	dObjective_H2["CAPEX_Compression_Truck"] = CAPEX_Compression_Truck

	## End Objective Function Expressions ##

	## Balance Expressions ##
	# H2 Power Consumption balance
	@expression(HY, eH2TruckCompressionPowerConsumption[t=1:T, z=1:Z],
	sum(vNcharged[z,zz,r,t] * TruckCap[r]* H2TruckCompressionEnergy[r] for zz = 1:Z,r = 1:RT if zz != z))
	dExpressions_H2["eH2TruckCompressionPowerConsumption"] = eH2TruckCompressionPowerConsumption
	dBalance_H2["H2PowerBalance"] = dBalance_H2["H2PowerBalance"] + eH2TruckCompressionPowerConsumption

	# H2 balance
	@expression(HY, TruckFlow[t=1:T,z=1:Z], sum(vH2TruckFlow[z,zz,r,t] for zz = 1:Z, r = 1:RT if zz != z))
	dExpressions_H2["TruckFlow"] = TruckFlow
	dBalance_H2["H2Balance"] = dBalance_H2["H2Balance"] + TruckFlow

	# Carbon emission balance
	@expression(HY, Truck_carbon_emission, sum(weights[t] * (vNarrive_full[z,zz,r,t]*Full_weight[r]+vNarrive_empty[z,zz,r,t]*Empty_weight[r]) *RouteLength[z,zz] *truck_emission_rate[r,z] for zz = 1:Z, z = 1:Z, r = 1:RT, t = 1:T if zz != z))
	dExpressions_H2["Truck_carbon_emission"] = Truck_carbon_emission
	dBalance_H2["H2_Carbon_Emission"] = dBalance_H2["H2_Carbon_Emission"] + Truck_carbon_emission
	## End Balance Expressions ##
	### End Expressions ###

	### Constraints ###
	for w in 1:W
	    for h in 1:Tw
	        t = Tw*(w-1)+h
	        tw_min = Tw*(w-1)+1
	        tw_max = Tw*(w-1)+Tw
	        for z in 1:Z
	            # for t in 1:T
	            for zz in 1:Z
	                for r in 1:RT
                        ## Total number
                        @constraint(HY,
                        vN_full[z,zz,r,t] +
                        vN_empty[z,zz,r,t] +
                        0
                        ==
                        vNTruckRoute[z,zz,r]
                        )

                        # The number of total full trucks
                        @constraint(HY, vN_full[z,zz,r,t] ==
                        vNtravel_full[zz,z,r,t] +
                        vNtravel_full[z,zz,r,t] +
                        vNavail_full[zz,z,r,t] +
                        vNavail_full[z,zz,r,t]
                        )
                        ## The number of total empty trucks
                        @constraint(HY, vN_empty[z,zz,r,t] ==
                        vNtravel_empty[zz,z,r,t] +
                        vNtravel_empty[z,zz,r,t] +
                        vNavail_empty[zz,z,r,t] +
                        vNavail_empty[z,zz,r,t]
                        )

                        if h > 1
                        @constraint(HY, vN_full[z,zz,r,t] - vN_full[z,zz,r,t-1] ==
                        vNcharged[z,zz,r,t] +
                        vNcharged[zz,z,r,t] -
                        vNdischarged[zz,z,r,t] -
                        vNdischarged[z,zz,r,t]
                        )
                        else
                            @constraint(HY, vN_full[z,zz,r,t] - vN_full[z,zz,r,tw_max] ==
                        vNcharged[z,zz,r,t] +
                        vNcharged[zz,z,r,t] -
                        vNdischarged[zz,z,r,t] -
                        vNdischarged[z,zz,r,t]
                            )
                        end
                        ## Change of the number of total empty trucks
                        if h > 1
                        @constraint(HY, vN_empty[z,zz,r,t] - vN_empty[z,zz,r,t-1] ==
                        -vNcharged[z,zz,r,t] -
                        vNcharged[zz,z,r,t] +
                        vNdischarged[z,zz,r,t] +
                        vNdischarged[zz,z,r,t]
                        )
                        else
                            @constraint(HY, vN_empty[z,zz,r,t] - vN_empty[z,zz,r,tw_max] ==
                            -vNcharged[z,zz,r,t] -
                            vNcharged[zz,z,r,t] +
                            vNdischarged[z,zz,r,t] +
                            vNdischarged[zz,z,r,t]
                            )
                        end



                        # Initially empty

                        # @constraint(HY, vNavail_full[zz,z,r,1] == 0)
                        # @constraint(HY, vNtravel_full[zz,z,r,1] == 0)
                        # @constraint(HY, vNavail_empty[1,2,r,1] == 0)

                        ## void cell
                        if h==1
                        # @constraint(HY, vNavail_full[zz,z,r,t] == 0)
                        # @constraint(HY, vNtravel_full[zz,z,r,t] == 0)
                        # @constraint(HY, vNarrive_full[zz,z,r,t] == 0 )
                        # @constraint(HY, vNdepart_full[zz,z,r,t] == 0 )
                        # @constraint(HY, vNarrive_empty[zz,z,r,t] == 0 )

                        # @constraint(HY, vNtravel_empty[zz,z,r,t] == 0)

                        end

                        # @constraint(HY, vNtravel_empty[2,1,r,t] == 0)

                        @constraint(HY, vNdischarged[zz,z,r,t] <= vNTruckRoute[zz,z,r] )
                        @constraint(HY, vNcharged[zz,z,r,t] <= vNTruckRoute[zz,z,r] )

                        @constraint(HY, vNTruckRoute[zz,z,r] == vNTruckRoute[z,zz,r])

						@constraint(HY, vNTruck[r] == sum(0.5*vNTruckRoute[zz,z,r] for zz = 1:Z,z=1:Z,r = 1:RT if zz != z))

                        if z==zz
	                        @constraint(HY, vNTruckRoute[zz,z,r] == 0 )
	                        @constraint(HY, vNtravel_full[zz,z,r,t] == 0 )
	                        @constraint(HY, vNtravel_empty[zz,z,r,t] == 0 )
	                        @constraint(HY, vNavail_full[zz,z,r,t] == 0 )
	                        @constraint(HY, vNavail_empty[zz,z,r,t] == 0 )
	                        @constraint(HY, vNarrive_full[zz,z,r,t] == 0 )
	                        @constraint(HY, vNdepart_full[zz,z,r,t] == 0 )
	                        @constraint(HY, vNarrive_empty[zz,z,r,t] == 0 )
	                        @constraint(HY, vNdepart_empty[zz,z,r,t] == 0 )
	                        @constraint(HY, vH2TruckFlow[zz,z,r,t] == 0 )
	                        else
	                        @constraint(HY, vNTruckRoute[zz,z,r] >=0 )
	                        # @constraint(HY, 5000>= vNTruckRoute[zz,z,r] )
                        end

                        t_arrive = 1
                        t_depart = 1

                        # Change of the number of full available trucks
                        if h > 1
	                        @constraint(HY, vNavail_full[z,zz,r,t] - vNavail_full[z,zz,r,t-1] ==
	                        vNcharged[z,zz,r,t] - vNdischarged[z,zz,r,t] +
	                        vNarrive_full[zz,z,r,t-t_arrive] -
	                        vNdepart_full[z,zz,r,t-t_depart] +
	                        0
	                        )
                        else
                            @constraint(HY, vNavail_full[z,zz,r,t] - vNavail_full[z,zz,r,tw_max] ==
                            vNcharged[z,zz,r,t] - vNdischarged[z,zz,r,t] +
                            vNarrive_full[zz,z,r,tw_max]
                            - vNdepart_full[z,zz,r,tw_max])
                        end
                        # Change of the number of empty available trucks
                        if h > 1
	                        @constraint(HY, vNavail_empty[z,zz,r,t] - vNavail_empty[z,zz,r,t-1] ==
	                        -vNcharged[z,zz,r,t] + vNdischarged[z,zz,r,t] +
	                        vNarrive_empty[zz,z,r,t-t_arrive] -
	                        vNdepart_empty[z,zz,r,t-t_depart])
                        else
                            @constraint(HY, vNavail_empty[z,zz,r,t] - vNavail_empty[z,zz,r,tw_max] ==
                            -vNcharged[z,zz,r,t] + vNdischarged[z,zz,r,t] +
                            vNarrive_empty[zz,z,r,tw_max] -
                            vNdepart_empty[z,zz,r,tw_max])
                        end

                        # Change of the number of full traveling trucks
                        if h > 1
	                        @constraint(HY, vNtravel_full[z,zz,r,t] - vNtravel_full[z,zz,r,t-1] ==
	                        vNdepart_full[z,zz,r,t-t_depart] -
	                        vNarrive_full[z,zz,r,t-t_arrive]
	                        )
                        else
                            @constraint(HY, vNtravel_full[z,zz,r,t] - vNtravel_full[z,zz,r,tw_max] ==
                            vNdepart_full[z,zz,r,tw_max]
                            - vNarrive_full[z,zz,r,tw_max])
                        end

                        # Change of the number of empty traveling trucks
                        if h > 1
	                        @constraint(HY, vNtravel_empty[z,zz,r,t] - vNtravel_empty[z,zz,r,t-1] ==
	                        vNdepart_empty[z,zz,r,t-t_depart] -
	                        vNarrive_empty[z,zz,r,t-t_arrive])
                        else
                            @constraint(HY, vNtravel_empty[z,zz,r,t] - vNtravel_empty[z,zz,r,tw_max] ==
                            vNdepart_empty[z,zz,r,tw_max] -
                            vNarrive_empty[z,zz,r,tw_max])
                        end


                        # Travel delay
                        if t-TD[zz,z]+1 >= tw_min && t-TD[zz,z]+1 <= tw_max && t-TD[zz,z]+1 <= t
	                        @constraint(HY, vNtravel_full[zz,z,r,t] >= sum(vNdepart_full[zz,z,r,tt] for tt = (t-TD[zz,z]+1):t) ) # deaprt from zz to zz
	                        @constraint(HY, vNtravel_empty[zz,z,r,t] >= sum(vNdepart_empty[zz,z,r,tt] for tt = (t-TD[zz,z]+1):t) )
                        end
                        if t+TD[zz,z] >= tw_min && t+TD[zz,z] <= tw_max && t+1 <= t+TD[zz,z]
	                        @constraint(HY, vNtravel_full[zz,z,r,t] >= sum(vNarrive_full[zz,z,r,tt] for tt = (t+1):(t+TD[zz,z]) ) ) # arrive to z to zz
	                        @constraint(HY, vNtravel_empty[zz,z,r,t] >= sum(vNarrive_empty[zz,z,r,tt] for tt = (t+1):(t+TD[zz,z]) ) )
                        end

						# Hydrogen flow by truck
                        @constraint(HY, vH2TruckFlow[z,zz,r,t]  == vNdischarged[z,zz,r,t] * TruckCap[r] * (1 - H2TLoss[r]) -
                        vNcharged[z,zz,r,t] * TruckCap[r]  +
                        0
                        )

                        # @constraint(HY, vH2TruckFlow[z,zz,r,t]  == vH2TruckFlow_pos[z,zz,r,t] -
                        # vH2TruckFlow_neg[z,zz,r,t]  +
                        # 0
                        # )

                        @constraint(HY, sum(vNcharged[z,zz,r,t] for zz = 1:Z if zz != z) * TruckCap[r] <= vH2TruckCompressionCap[z,r]
                        )

	                end
	            end
	        end
	    end
	end
	### End Constraints ###

	dModuleArgs["dExpressions_H2"] = dExpressions_H2
	dModuleArgs["dObjective_H2"] = dObjective_H2
	dModuleArgs["dBalance_H2"] = dBalance_H2
	return HY, dModuleArgs
end # end H2Truck module



function save_results(HY::Model, dModuleArgs::Dict)


end # end save_results module
# dfStatus = DataFrame(Status = status, Presolve = presolver_time, Solve = solver_time,
#     Objval = getobjectivevalue(HY), Objbound= getobjbound(HY),FinalMIPGap =(getobjectivevalue(HY) -getobjbound(HY))/getobjectivevalue(HY) )


# @save "data.jld2" inputs outputs

end # End Module
