module GenX_Modular

# Configurable Capacity Expansion Model with Optimal Dispatch
# Created by Nestor Sepulveda .&  Jesse Jenkins, 2016

#export Module functions
export load_power_inputs, generate_power_system_model, solve_power_system_model, write_power_outputs, economic_dispatch, investment, ucommit, reserves, transmission, dispatchable, nondispatchable, storage, hydro, dr, heat, thermal, hybrid, trading, co2, rps, ces, specificshare

## LOADING PACKAGES

# Define the packages
using JuMP # used for mathematical programming
using DataFrames #This package allows put together data into a matrix
# using Gurobi #Gurobi solver
using MathProgBase #for fix_integers
using CSV
using StatsBase
using LinearAlgebra

#=
working_path = "/Users/jdj2/GenX"
inpath = "$working_path/Input_Short"
outpath = "$working_path/Test_Results"
# # #
mysetup=Dict()
mysetup["Trans_Loss_Segments"]=0 # Number of segments to use in piecewise linear approximation of losses; 1 = linear, >2 = piecewise quadratic
mysetup["Distr_Loss_Segments"]=0 # Number of segments to use in piecewise linear approximation of losses; 1 = linear, >2 = piecewise quadratic
mysetup["Distr_Margin_LV_Segments"]=0 # Number of segments to use in piecewise linear approximation of quadratic term in LV distribution margin equation
mysetup["CO2Cap"]=0 # 0 = no cap; 1 = zonal caps; 2 = overall cap;
mysetup["UCommit"]=2 # 0 = no UC, 1 integer clestering, 2 linearized clustering
mysetup["Specific_Share"] =0 # 1 = active otherwise no active
mysetup["RPS"]=2 #minimum qualifying renewables penetration; 0 = no active; 1 = active per node; 2 = active as system
mysetup["CES"]=2 #minimum qualifying clean energy resource penetration; 0 = no active; 1 = active per node; 2 = active as system
mysetup["NetworkExpansion"]=0 # 1 = active otherwise no active
mysetup["HeatMarket"]=0 # 1 = active otherwise no active
mysetup["Externalities"]=0 # 1 = active otherwise no active
mysetup["pTolerance"] = 0.001 # Wrap-up tolerance
mysetup["Reserves"] = 0 # 0 = no reserves, 1 regulation (primary) and operating (secondary) reserves
mysetup["MIPGap"] = 0.01 # 1.0% MIP Gap
mysetup["MacOrWindows"] = "Mac"  # Set to either "Mac" (also works for Linux) or "Windows" to ensure use of proper file directory separator "\" or "/
mysetup["Dual_MIP"] = 0 # get dual of MIP
mysetup["Feasib_Tol"] = 1e-6 # All constraints must be satisfied to a tolerance of
mysetup["Pre_Solve"] = 2 # Controls the presolve level. A value of -1 corresponds to an automatic setting. Other options are off (0), conservative (1), or aggressive (2).
mysetup["Crossover"] = -1 # Determines the crossover strategy used to transform the interior solution produced by barrier into a basic solution.
mysetup["Method"] = -1 # Algorithm used to solve continuous models or the root node of a MIP model. Options are: -1=automatic, 0=primal simplex, 1=dual simplex, 2=barrier, 3=concurrent, 4=deterministic concurrent.
mysetup["NumericFocus"] = 0 #The NumericFocus parameter controls the degree to which the code attempts to detect and manage numerical issues. The default setting (0)
mysetup["Optimal_Tol"] = 1e-6 #Reduced costs must all be smaller than OptimalityTol in the improving direction in order for a model to be declared optimal.
mysetup["TimeLimit"] = Inf # Time limit (in seconds) for Gurobi solver. Will terminate if optimal solution not found before this specified limit. Default is Infinity.
#NOTE: new setup option for time wrapping
mysetup["OperationWrapping"]=0 # 0 = Full Year First-Last Hour Year; 1 = Full Year First-Last Hour Week; 2 = Sample Weeks First-Last Hour Week;
mysetup["StorageLosses"]=0 # 0 VRE and CO2 constraint DO NOT account for energy lost; 1 constraint DO  account for energy lost
mysetup["PrintModel"]=0 # 1 for including the model equation as an output; 0 for the model equation not being included as an output
mysetup["Trading"]=0 # Jesse Trading function

 path = inpath
 setup = mysetup
 inputs = load_inputs(setup,path)
 myresults = solve_model(setup,inputs)
 write_outputs(setup,outpath,myresults,inputs)
=#

################################################################################
## Function load_inputs(setup, path)
##
## inputs:
##  setup - dict object containing setup parameters
##  path - [string] path to working directory
##
## description: Loads various data inputs from multiple input .csv files in path
## directory and stores variables in a Dict (dictionary) object for use in
## model() function
##
## returns: Dict (dictionary) object containing all data inputs
################################################################################
function load_power_inputs(setup::Dict,path::AbstractString)

	## Use appropriate directory separator depending on Mac or Windows config
	if setup["MacOrWindows"]=="Mac"
		sep = "/"
	else
		sep = "\U005c"
	end

	## Read input files
	println("About to Read In CSV Files")
	# Load related inputs
	load_in = DataFrame(CSV.File(string(path,sep,"Load_data.csv"), header=true), copycols=true)
	println("Load_data.csv Successfully Read!")
	# Generator related inputs
	gen_in = DataFrame(CSV.File(string(path,sep,"Generators_data.csv"), header=true), copycols=true)
	println("Generators_data.csv Successfully Read!")
	# Fuel related inputs
	fuels_in = DataFrame(CSV.File(string(path,sep,"Fuels_data.csv"), header=true), copycols=true)
	println("Fuels_data.csv Successfully Read!")
	# Hourly capacity factor for variable
	gen_var = DataFrame(CSV.File(string(path,sep,"Generators_variability.csv"), header=true), copycols=true)
	println("Generators_variability.csv Successfully Read!")
	# Network zones inputs and Network topology inputs
	network_var = DataFrame(CSV.File(string(path,sep,"Network.csv"), header=true), copycols=true)
	println("Network.csv Successfully Read!")

	println("CSV Files Successfully Read In")

	if setup["Reserves"]==1
		# Reserve related inputs
		res_in = DataFrame(CSV.File(string(path,sep,"Reserves.csv"), header=true), copycols=true)
	end

	if setup["Trading"]==1
		# Reserve related inputs
		trade_in = CSV.read(string(path,sep,"Trading.csv"), copycols=true, header=true)
	end
	if setup["Fuel_Price_Var"]==1
		# Hourly scaling factor for fuel price
		fuels_var = DataFrame!(CSV.File(string(path,sep,"Fuels_variability.csv"), header=true), copycols=true)
	end
	## Declare Dict (dictionary) object used to store parameters
	inputs = Dict()
	# inputs["Power2H2"] = 0

	## Create sets and set indices
	# Number of lines in the network
	inputs["L"]=size(collect(skipmissing(network_var[!,:Network_lines])),1)
	# Number of zones in the network
	inputs["Z"]=size(collect(skipmissing(network_var[!,:Network_zones])),1)
	# Number of resources
	inputs["G"] = size(collect(skipmissing(gen_in[!,:R_ID])),1)
	# Names of resources
	inputs["RESOURCES"] = collect(skipmissing(gen_in[!,:Resource][1:inputs["G"]]))
	# Number of time steps (hours)
	inputs["T"] = size(collect(skipmissing(load_in[!,:Time_index])),1)
	# Number of demand curtailment/lost load segments
	inputs["SEG"]=size(collect(skipmissing(load_in[!,:Demand_segment])),1)

	## (from Period_Map.csv)
	# read in mapping of modeled periods to representative periods
	if (setup["OperationWrapping"]==2 && setup["LDS"]==1)
		inputs["Period_Map"] = DataFrame(CSV.File(string(path,sep,"Period_Map.csv"), header=true), copycols=true)
		println("Perio_Map.csv Successfully Read!")
	end


	# inputs["omega"] = zeros(Float64, T)
	# inputs["W"] = 1
	# inputs["H"] = 1
	if setup["OperationWrapping"]==0
		# Total number of subtime periods
		inputs["H"] = collect(skipmissing(load_in[!,:Hours_per_period]))[1]
		inputs["W"] = 1
		inputs["T"] = inputs["W"] * inputs["H"]
		# Simple scaling factor for number of hours
		inputs["omega"] = zeros(Float64, inputs["T"])
		inputs["omega"][:] .= 1
	elseif setup["OperationWrapping"]==1
		inputs["H"] = collect(skipmissing(load_in[!,:Hours_per_period]))[1]
		inputs["W"] = div.(T,inputs["H"])[1]  # Total number of subtime periods
		# Simple scaling factor for number of hours
		inputs["T"] = inputs["W"] * inputs["H"]
		inputs["omega"] = zeros(Float64, inputs["T"])
		inputs["omega"][:] .= 1
	elseif setup["OperationWrapping"]==2
		# Weights for each subperiod
		inputs["Weights"] = collect(skipmissing(load_in[!,:Sub_Weights])) # Weights each period

		# Total number of subtime periods
		inputs["W"] = convert(Int8, collect(skipmissing(load_in[!,:Subperiods]))[1])
		inputs["H"] = convert(Int16, collect(skipmissing(load_in[!,:Hours_per_period]))[1])

		# inputs["T"] = size(collect(skipmissing(load_in[!,:Time_index])),1)
		inputs["T"] = inputs["W"] * inputs["H"]
		inputs["omega"] = zeros(Float64, inputs["T"])
		# creating hourly weights from weekly weights
		for w in 1:inputs["W"]
			for h in 1:inputs["H"]
				t = inputs["H"]*(w-1)+h
				inputs["omega"][t] = inputs["Weights"][w]/inputs["H"]*8736/sum(inputs["Weights"][1:inputs["W"]])
			end
		end
	end



	## Set indices for internal use
	G = inputs["G"]   # Number of resources (generators, storage, DR, and DERs)
	T = inputs["T"]   # Total number of time steps (hours)
	Z = inputs["Z"]   # Number of zones
	L = inputs["L"]   # Number of lines
	SEG = inputs["SEG"]  # Number of lines

	## (from Load_data.csv)
	println("from Load_data.csv")
	# Demand in MW for each zone
	start = findall(s -> s == "Load_MW_z1", names(load_in))[1]
	if setup["ParameterScale"] ==1  # Parameter scaling turned on
		# Max value of non-served energy
		inputs["Voll"] = collect(skipmissing(load_in[!,:Voll])) /1e+3 # convert from $/MWh $ million/GWh
		# Demand in MW
		inputs["pD"] =convert(Matrix, load_in[1:inputs["T"],start:start-1+inputs["Z"]])*setup["power_load_factor"]/1e+3  # convert to GW

	else # No scaling
		# Max value of non-served energy
		inputs["Voll"] = collect(skipmissing(load_in[!,:Voll]))
		# Demand in MW
		inputs["pD"] =convert(Matrix, load_in[1:inputs["T"],start:start-1+inputs["Z"]])*setup["power_load_factor"]
	end

	# Cost of non-served energy/demand curtailment (for each segment)
	inputs["pC_D_Curtail"] = zeros(SEG)
	inputs["pMax_D_Curtail"] = zeros(SEG)
	for s in 1:SEG
		inputs["pC_D_Curtail"][s] = collect(skipmissing(load_in[!,:Cost_of_demand_curtailment_perMW]))[s]*inputs["Voll"][1]
		# Maximum hourly demand curtailable as % of the max demand (for each segment)
		inputs["pMax_D_Curtail"][s] = collect(skipmissing(load_in[!,:Max_demand_curtailment]))[s]
	end

	## (from Fuels_data.csv)
	println("from Fuels_data.csv")
	# Fuel costs .&  CO2 emissions rate for each fuel type (stored in dictionary objects)
	fuels = fuels_in[!,:Fuel] # fuel type indexes
	if setup["ParameterScale"] ==1  # Parameter scaling turned on
		# scaled to convert objective function to $million with power in GW
		costs = convert(Array{Float64}, fuels_in[!,:Cost_per_MMBtu] )/1e+3
	else # No scaling
		costs = convert(Array{Float64}, fuels_in[!,:Cost_per_MMBtu] )   # $/MMBtu
	end

	CO2_content = convert(Array{Float64}, fuels_in[!,:CO2_content_tons_perMMBtu] )    # tons CO2/MMBtu
	fuel_costs = Dict{AbstractString,Float64}()
	fuel_CO2 = Dict{AbstractString,Float64}()
	for i = 1:length(fuels)
		fuel_costs[fuels[i]] = costs[i]
		fuel_CO2[fuels[i]] = CO2_content[i]
	end

	## (from Generators_data.csv)'
	println("from Generators_data.csv")
	# Number of generators/resources eligible per zone
	inputs["pGensPerZone"]=countmap(collect(skipmissing(gen_in[!,:zone])))
	# Store DataFrame of generators/resources input data for use in model
	inputs["dfGen"]=gen_in
	if setup["ParameterScale"] ==1  # Parameter scaling turned on - adjust values of subset of parameter values
		inputs["dfGen"][!,:Existing_Charge_Cap_MW] = gen_in[!,:Existing_Charge_Cap_MW]/1e+3 # Convert to GW
		inputs["dfGen"][!,:Existing_Cap_MWh] = gen_in[!,:Existing_Cap_MWh]/1e+3 # Convert to GWh
		inputs["dfGen"][!,:Existing_Cap_MW] = gen_in[!,:Existing_Cap_MW]/1e+3 # Convert to GW

		# Cap_size scale only capacities for those technologies with capacity >1
		# Step 1: convert vector to float
		inputs["dfGen"][!,:Cap_size] =convert(Array{Float64}, gen_in[!,:Cap_size] )
		for g in 1:G  # Scale only those capacities for which cap_size > 1
			if inputs["dfGen"][!,:Cap_size][g]>1.0
				inputs["dfGen"][!,:Cap_size][g] = gen_in[!,:Cap_size][g]/1e+3 # Convert to GW
				# println(gen_in[!,:Cap_size][g])
			end
		end
		inputs["dfGen"][!,:Min_Cap_MW] = gen_in[!,:Min_Cap_MW]/1e+3 # Convert to GW
		inputs["dfGen"][!,:Max_Cap_MW] = gen_in[!,:Max_Cap_MW]/1e+3 # Convert to GW
		# Investment cost terms
		inputs["dfGen"][!,:Inv_cost_per_MWyr] = gen_in[!,:Inv_cost_per_MWyr]/1e+3 # Convert to $ million/GW/yr with objective function in millions
		inputs["dfGen"][!,:Inv_cost_per_MWhyr] = gen_in[!,:Inv_cost_per_MWhyr]/1e+3 # Convert to $ million/GWh/yr  with objective function in millions
		inputs["dfGen"][!,:Inv_cost_charge_per_MWyr] = gen_in[!,:Inv_cost_charge_per_MWyr]/1e+3 # Convert to $ million/GWh/yr  with objective function in millions
		# Fixed O&M cost terms
		inputs["dfGen"][!,:Fixed_OM_cost_per_MWyr] = gen_in[!,:Fixed_OM_cost_per_MWyr]/1e+3 # Convert to $ million/GW/yr with objective function in millions
		inputs["dfGen"][!,:Fixed_OM_cost_per_MWhyr] = gen_in[!,:Fixed_OM_cost_per_MWhyr]/1e+3 # Convert to $ million/GW/yr with objective function in millions
		inputs["dfGen"][!,:Fixed_OM_cost_charge_per_MWyr] = gen_in[!,:Fixed_OM_cost_charge_per_MWyr]/1e+3 # Convert to $ million/GW/yr with objective function in millions
		# Variable O&M cost terms
		inputs["dfGen"][!,:Var_OM_cost_per_MWh] = gen_in[!,:Var_OM_cost_per_MWh]/1e+3 # Convert to $ million/GWh with objective function in millions
		inputs["dfGen"][!,:Var_OM_cost_per_MWh_in] = gen_in[!,:Var_OM_cost_per_MWh_in]/1e+3 # Convert to $ million/GWh with objective function in millions
		inputs["dfGen"][!,:Externality_cost_MWh] = gen_in[!,:Externality_cost_MWh]/1e+3 # Convert to $ million/GWh with objective function in millions

		inputs["dfGen"][!,:Start_cost_per_MW] = gen_in[!,:Start_cost_per_MW]/1e+3 # Convert to $ million/GW with objective function in millions
		inputs["dfGen"][!,:Reg_Cost] = gen_in[!,:Reg_Cost]/1e+3 # Convert to $ million/GW with objective function in millions
		inputs["dfGen"][!,:Rsv_Cost] = gen_in[!,:Rsv_Cost]/1e+3 # Convert to $ million/GW with objective function in millions

	end

	if setup["UCommit"]>=1
		# Fuel consumed on start-up (million BTUs per MW per start) if unit commitment is modelled
		Start_Fuel = convert(Array{Float64}, collect(skipmissing(gen_in[!,:Start_fuel_MMBTU_per_MW])) )
	# Fixed cost per start-up ($ per MW per start) if unit commitment is modelled
		Start_Cost = convert(Array{Float64}, collect(skipmissing(gen_in[!,:Start_cost_per_MW])) )
		inputs["dfGen"][!,:C_Start] = zeros(Float64, G)
		inputs["dfGen"][!,:CO2_per_Start] = zeros(Float64, G)
	end

	# Heat rate of all resources (million BTUs/MWh)
	Heat_Rate = convert(Array{Float64}, collect(skipmissing(gen_in[!,:Heat_rate_MMBTU_per_MWh])) )
	# Fuel used by each resource
	Fuel_Type = collect(skipmissing(gen_in[!,:Fuel]))
	# Maximum fuel cost in $ per MWh and CO2 emissions in tons per MWh
	inputs["dfGen"][!,:C_Fuel_per_MWh] = zeros(Float64, G)
	inputs["dfGen"][!,:CO2_per_MWh] = zeros(Float64, G)
	for g in 1:G
		inputs["dfGen"][!,:C_Fuel_per_MWh][g] = fuel_costs[Fuel_Type[g]]*Heat_Rate[g]
		inputs["dfGen"][!,:CO2_per_MWh][g] = fuel_CO2[Fuel_Type[g]]*Heat_Rate[g] *1e+3
		if setup["UCommit"]>=1
			# Start-up cost is sum of fixed cost per start plus cost of fuel consumed on startup.
			# CO2 from fuel consumption during startup also calculated
			inputs["dfGen"][!,:C_Start][g] = inputs["dfGen"][!,:Cap_size][g]*(fuel_costs[Fuel_Type[g]]*Start_Fuel[g] + Start_Cost[g])
			inputs["dfGen"][!,:CO2_per_Start][g] = inputs["dfGen"][!,:Cap_size][g]*(fuel_CO2[Fuel_Type[g]]*Start_Fuel[g])
		end
	end

	## (from Generators_variability.csv)
	println("from Generators_variability.csv")
	# Maximum power output and variability of each energy resource
	inputs["pP_Max"] = transpose(convert(Matrix{Float64}, gen_var[1:inputs["T"],2:(inputs["G"]+1)]))

	## (from Network.csv)
	println("from Network.csv")
	# Topology of the network
	start = findall(s -> s == "z1", names(network_var))[1]
	inputs["pNet_Map"] = convert(Matrix{Float64}, network_var[1:inputs["L"],start:start+inputs["Z"]-1])
	# Max zone/bus angle
	## DC OPF not currently implemented; code segment commented out below; commenting this out as well to avoid having to include column in input file
	# inputs["pTheta_Max"] = convert(Array{Float64}, collect(skipmissing(network_var[!,:Thetha_max])))

	if setup["CO2Cap"] >= 1
		# Max zone/bus emissions in tons
		inputs["pMaxCO2"] = convert(Array{Float64}, collect(skipmissing(network_var[!,:CO_2_Max_Mtons]))*10^6)
		# Max zone/bus emissions in ton/MWh
		inputs["pMaxCO2Rate"] = convert(Array{Float64}, collect(skipmissing(network_var[!,:CO_2_Max_ton_MWh])))
		# CO2 emissions rate applied per MWh of imports to zone for purposes of emissions cap
		inputs["pCO2ImportsRate"] = convert(Array{Float64}, collect(skipmissing(network_var[!,:CO_2_Import_ton_MWh])))
		# CO2 emissions rate credited per MWh of exports to zone for purposes of emissions cap
		inputs["pCO2ExportsRate"] = convert(Array{Float64}, collect(skipmissing(network_var[!,:CO_2_Export_ton_MWh])))

	end

	# Min zone/bus qualifying renewables share (Renewable Portfolio Standard)
	if setup["RPS"] >= 1
		inputs["pRPS"] = convert(Array{Float64}, collect(skipmissing(network_var[!,:RPS])))
	end
	# Min zone/bus qualifying clean energy share (Clean Energy Standard)
	if setup["CES"] >= 1
		inputs["pCES"] = convert(Array{Float64}, collect(skipmissing(network_var[!,:CES])))
	end

	# Transmission capacity of the network (in MW)
	if setup["ParameterScale"] ==1  # Parameter scaling turned on - adjust values of subset of parameter values
		inputs["pTrans_Max"] = convert(Array{Float64}, collect(skipmissing(network_var[!,:Line_Max_Flow_MW])))/1e+3  # convert to GW
	else # no scaling
		inputs["pTrans_Max"] = convert(Array{Float64}, collect(skipmissing(network_var[!,:Line_Max_Flow_MW])))
	end

	if setup["Trans_Loss_Segments"] == 1
		# Line percentage Loss - valid for case when modeling losses as a fixed percent of absolute value of power flows
		inputs["pPercent_Loss"] = convert(Array{Float64}, collect(skipmissing(network_var[!,:Line_Loss_Percentage])))
	elseif setup["Trans_Loss_Segments"] >= 2
		# Transmission line voltage (in kV)
		inputs["kV"] = convert(Array{Float64}, collect(skipmissing(network_var[!,:Line_Voltage_kV])))
		# Transmission line resistance (in Ohms) - Used in transmission losses
		inputs["Ohms"] = convert(Array{Float64}, collect(skipmissing(network_var[!,:Line_Resistance_ohms])))
		# Transmission line resistance (in Ohms) - [used in DC power flow. Appears to be unused input so commented out]
		# inputs["line_R"] = convert(Array{Float64}, collect(skipmissing(network_var[!,:Line_R_ohms])))
		# Transmission line reactance (in Ohms) - [used in DC power flow. Appears to be unused input so commented out]
		# inputs["line_X"] = convert(Array{Float64}, collect(skipmissing(network_var[!,:Line_X_ohms])))
	end


	# Maximum possible flow after reinforcement for use in linear segments of piecewise approximation
	inputs["pTrans_Max_Possible"] = zeros(Float64, inputs["L"])
	if setup["NetworkExpansion"]==1
		if setup["ParameterScale"] ==1  # Parameter scaling turned on - adjust values of subset of parameter values
			# Read between zone network reinforcement costs per peak MW of capacity added
			inputs["pC_Line_Reinforcement"] = convert(Array{Float64}, collect(skipmissing(network_var[!,:Line_Reinforcement_Cost_per_MW_yr])))/1e+3 # convert to million $/GW/yr with objective function in millions
			# Maximum reinforcement allowed in MW
			#NOTE: values <0 indicate no expansion possible
			inputs["pMax_Line_Reinforcement"] = convert(Array{Float64}, collect(skipmissing(network_var[!,:Line_Max_Reinforcement_MW])))/1e+3 # convert to GW
		else
			# Read between zone network reinforcement costs per peak MW of capacity added
			inputs["pC_Line_Reinforcement"] = convert(Array{Float64}, collect(skipmissing(network_var[!,:Line_Reinforcement_Cost_per_MW_yr])))
			# Maximum reinforcement allowed in MW
			#NOTE: values <0 indicate no expansion possible
			inputs["pMax_Line_Reinforcement"] = convert(Array{Float64}, collect(skipmissing(network_var[!,:Line_Max_Reinforcement_MW])))
		end
		for l in 1:inputs["L"]
			if inputs["pMax_Line_Reinforcement"][l] > 0
				inputs["pTrans_Max_Possible"][l] = inputs["pTrans_Max"][l] + inputs["pMax_Line_Reinforcement"][l]
			else
				inputs["pTrans_Max_Possible"][l] = inputs["pTrans_Max"][l]
			end
		end
	else
		inputs["pTrans_Max_Possible"] = inputs["pTrans_Max"]
	end
	# Transmission line (between zone) loss coefficient (resistance/voltage^2)
	inputs["pTrans_Loss_Coef"] = zeros(Float64, inputs["L"])
	for l in 1:inputs["L"]
		#for cases with only one segment
		if setup["Trans_Loss_Segments"] == 1
			inputs["pTrans_Loss_Coef"][l] = inputs["pPercent_Loss"][l]/inputs["pTrans_Max_Possible"][l]
		elseif setup["Trans_Loss_Segments"] >= 2
			# If zones are connected, loss coefficient is R/V^2 where R is resistance in Ohms and V is voltage in Volts
			if setup["ParameterScale"] ==1  # Parameter scaling turned on - adjust values of subset of parameter values
					inputs["pTrans_Loss_Coef"][l] = (inputs["Ohms"][l]/10^6)/(inputs["kV"][l]/10^3)^2 *1e+3 # 1/GW
			else
				inputs["pTrans_Loss_Coef"][l] = (inputs["Ohms"][l]/10^6)/(inputs["kV"][l]/10^3)^2 # 1/MW
			end

		end
	end


	## Distribution network zone parameters
	# NOTE: the scaling of distribution parameters as per parameter units noted in Table 4.3 of JDJ thesis. This has not been tested.
	###  None of the loss coefficients or distribution margin coefficient terms have been scaled since all variables and constraints are also written in GWs
	# Record subset of distribution voltage zones
	inputs["pDistrZones"] = convert(Array{Float64}, collect(skipmissing(network_var[!,:DistrZones])))
	if sum(inputs["pDistrZones"]) > 0
	# There are distribution zones, read relevant inputs...
		# Within zone (distribution) loss coefficients
		inputs["pDistr_Loss_LV_Net_Quad"] = convert(Array{Float64}, collect(skipmissing(network_var[!,:DistrLossFact_LV_Net_Quad])))
		inputs["pDistr_Loss_MV_Net_Linear"] = convert(Array{Float64}, collect(skipmissing(network_var[!,:DistrLossFact_MV_Net_Linear])))
		inputs["pDistr_Loss_LV_Total_Linear"] = convert(Array{Float64}, collect(skipmissing(network_var[!,:DistrLossFact_LV_Total_Linear])))
		# Share of Distr zone demand that is in medium and low voltage
		inputs["pDistr_Share_in_MV"]= convert(Array{Float64}, collect(skipmissing(network_var[!,:Share_in_MV])))
		inputs["pDistr_Share_in_LV"]= 1 .- inputs["pDistr_Share_in_MV"]
		## Parameters used for within zone losses:
		if setup["NetworkExpansion"]==1

			inputs["pDistr_Max_Inject_Possible"] = zeros(Float64, Z)  # maximum power injections within a zone (set below)
			inputs["pDistr_Max_Withdraw_Possible"] = zeros(Float64, Z)  # maximum power withdrawals within a zone (set below)

			if setup["ParameterScale"] ==1  # Parameter scaling turned on - adjust values of subset of parameter values
				# Read within zone network reinforcement cost per peak MW of capacity added due to MV and LV demand
				inputs["pC_MV_Reinforcement"] = convert(Array{Float64}, collect(skipmissing(network_var[!,:Distr_MV_Reinforcement_Cost_per_MW_yr])))/1e+3 # convert to million $/GW/yr with objective in millions
				inputs["pC_LV_Reinforcement"] = convert(Array{Float64}, collect(skipmissing(network_var[!,:Distr_LV_Reinforcement_Cost_per_MW_yr])))/1e+3 # convert to million $/GW/yr with objective in millions
				# Maximum net withdrawals and maximum net injection by zone (distribution network capacity)
				inputs["pMax_Inject"] = convert(Array{Float64}, collect(skipmissing(network_var[!,:Distr_Max_Inject])))/1e+3 # Convert to GW
				inputs["pMax_Withdraw"] = convert(Array{Float64}, collect(skipmissing(network_var[!,:Distr_Max_Withdraw])))/1e+3 # Convert to GW
				# Maximum zonal network capacity reinforcement allowed in MW (for both injection and withdrawal)
				# (note: negative values indicate zone cannot be reinforced)
				inputs["pMax_Inject_Reinforcement"] = convert(Array{Float64}, collect(skipmissing(network_var[!,:Distr_Inject_Max_Reinforcement_MW])))/1e+3 # Convert to GW
				inputs["pMax_Withdraw_Reinforcement"] = convert(Array{Float64}, collect(skipmissing(network_var[!,:Distr_Withdraw_Max_Reinforcement_MW])))/1e+3 # Convert to GW

			else
				# Read within zone network reinforcement cost per peak MW of capacity added due to MV and LV demand
				inputs["pC_MV_Reinforcement"] = convert(Array{Float64}, collect(skipmissing(network_var[!,:Distr_MV_Reinforcement_Cost_per_MW_yr])))
				inputs["pC_LV_Reinforcement"] = convert(Array{Float64}, collect(skipmissing(network_var[!,:Distr_LV_Reinforcement_Cost_per_MW_yr])))
				# Maximum net withdrawals and maximum net injection by zone (distribution network capacity)
				inputs["pMax_Inject"] = convert(Array{Float64}, collect(skipmissing(network_var[!,:Distr_Max_Inject])))
				inputs["pMax_Withdraw"] = convert(Array{Float64}, collect(skipmissing(network_var[!,:Distr_Max_Withdraw])))
				# Maximum zonal network capacity reinforcement allowed in MW (for both injection and withdrawal)
				# (note: negative values indicate zone cannot be reinforced)
				inputs["pMax_Inject_Reinforcement"] = convert(Array{Float64}, collect(skipmissing(network_var[!,:Distr_Inject_Max_Reinforcement_MW])))
				inputs["pMax_Withdraw_Reinforcement"] = convert(Array{Float64}, collect(skipmissing(network_var[!,:Distr_Withdraw_Max_Reinforcement_MW])))

			end
			# Length of segments used in segment-wise approximation of quadratic term for net demand reductions required by LV network margin needed
			inputs["pDistr_Margin_LV_Segment_Length"] = inputs["pMax_Withdraw_Reinforcement"].*inputs["pDistr_Share_in_LV"]/setup["Distr_Margin_LV_Segments"]
			# Read set of hours that may contain peak withdrawal and injection periods
			# Specify "all" in first entry to include all hours, otherwise list specific hours in these column
			inputs["pPeak_Inject_Hrs"] = collect(skipmissing(network_var[!,:Peak_Injection_Hours]))
			inputs["pPeak_Withdraw_Hrs"] = collect(skipmissing(network_var[!,:Peak_Withdrawal_Hours]))
			if (occursin(r"(?i)all", convert(AbstractString,inputs["pPeak_Inject_Hrs"][1])))
				inputs["pPeak_Inject_Hrs"] = 1:T
			end
			if (occursin(r"(?i)all", convert(AbstractString, inputs["pPeak_Withdraw_Hrs"][1])))
				inputs["pPeak_Withdraw_Hrs"] = 1:T
			end
			# Maximum possible injection/withdrawal after reinforcement for use in linear segments of piecewise approximation
			for z in 1:Z
				if inputs["pMax_Inject_Reinforcement"][z] > 0
					inputs["pDistr_Max_Inject_Possible"][z] = inputs["pMax_Inject"][z] + inputs["pMax_Inject_Reinforcement"][z]
				else
					inputs["pDistr_Max_Inject_Possible"][z] = inputs["pMax_Inject"][z]
				end
				if inputs["pMax_Withdraw_Reinforcement"][z] > 0
					inputs["pDistr_Max_Withdraw_Possible"][z] = inputs["pMax_Withdraw"][z] + inputs["pMax_Withdraw_Reinforcement"][z]
				else
					inputs["pDistr_Max_Withdraw_Possible"][z] = inputs["pMax_Withdraw"][z]
				end
			end
			# Network margin gained per MW of net demand reduction parameters
			# LV and MV coefficients
			inputs["pDistr_Margin_LV_Linear"]= convert(Array{Float64}, collect(skipmissing(network_var[!,:DistrMarginFact_LV_Linear])))
			inputs["pDistr_Margin_LV_Quad"]= convert(Array{Float64}, collect(skipmissing(network_var[!,:DistrMarginFact_LV_Quad])))
			inputs["pDistr_Margin_MV_Linear"]= convert(Array{Float64}, collect(skipmissing(network_var[!,:DistrMarginFact_MV_Linear])))
			# Maximum contribution of MV to margin
			inputs["pDistr_Margin_MV_Max"]= convert(Array{Float64}, collect(skipmissing(network_var[!,:DistrMargin_MV_Max])))
			# `Discount rate` to apply to MV reductions to account for marginal losses in distribution
			inputs["pDistr_Margin_MV_Discount"]= convert(Array{Float64}, collect(skipmissing(network_var[!,:DistrMargin_MV_DiscountFact])))
		else
			for z in 1:Z
				maxExport = dot(abs.(inputs["pNet_Map"][:,z]), inputs["pTrans_Max"][:])
				inputs["pDistr_Max_Inject_Possible"][z] = maximum(inputs["pD"][:,z]) + maxExport# max for segments is largest demand and largest possible export from zone
				inputs["pDistr_Max_Withdraw_Possible"][z] = maximum(inputs["pD"][:,z]) + maxExport# max for segments is largest demand and largest possible export from zone
			end
		end
	end


	## Heat data depending on setup
	if setup["HeatMarket"] == 1
		# heat market related inputs
		heat_in = DataFrame(CSV.File(string(path,sep,"Heat_data.csv"), header=true))
		if setup["ParameterScale"] ==1  # Parameter scaling turned on - adjust values of subset of parameter values
			# Heat demand in MWh for each zone
			inputs["pHeat_Demand"] =convert(Array, heat_in[1:inputs["T"],3:2+inputs["Z"]])/1e+3 # convert MWh to GWh
			# Cost of heat for heat market
			inputs["pHeat_Price"] = collect(skipmissing(heat_in[!,:Heat_Price_dollar_per_MWh]))/1e+3  # Convert $/MWh to million $/GWh with objective function in millions
			# Cost of heat for heat market
		else
			# Heat demand in MWh for each zone
			inputs["pHeat_Demand"] =convert(Array, heat_in[1:inputs["T"],3:2+inputs["Z"]])
			# Cost of heat for heat market
			inputs["pHeat_Price"] = collect(skipmissing(heat_in[!,:Heat_Price_dollar_per_MWh]))
		end
		# CO2 emissions associated with heat generation for heat market - assumed to come from natural gas
		inputs["pHeat_CO2"] =  fuel_CO2["natural_gas"]/0.29307107
	end #END heat data inputs

	##Reserve inputs
	if setup["Reserves"]==1
		inputs["pReg_Req_Load"] = convert(Float64, res_in[!,:Reg_Req_Percent_Load][1] )
		inputs["pReg_Req_VRE"] = convert(Float64, res_in[!,:Reg_Req_Percent_VRE][1] )
		inputs["pRsv_Up_Req_Load"] = convert(Float64, res_in[!,:Rsv_Up_Req_Percent_Load][1] )
		inputs["pRsv_Up_Req_VRE"] = convert(Float64, res_in[!,:Rsv_Up_Req_Percent_VRE][1] )
		inputs["pRsv_Dn_Req_Load"] = convert(Float64, res_in[!,:Rsv_Dn_Req_Percent_Load][1] )
		inputs["pRsv_Dn_Req_VRE"] = convert(Float64, res_in[!,:Rsv_Dn_Req_Percent_VRE][1] )
		if setup["ParameterScale"] ==1  # Parameter scaling turned on - adjust values of subset of parameter values
			inputs["pC_Rsv_Penalty"] = convert(Float64, res_in[!,:Unmet_Rsv_Penalty_dollar_per_MW][1] )/1e+3 # convert to million $/GW with objective function in millions
		else
			inputs["pC_Rsv_Penalty"] = convert(Float64, res_in[!,:Unmet_Rsv_Penalty_dollar_per_MW][1] )
		end
		inputs["pDynamic_Contingency"] = convert(Int8, res_in[!,:Dynamic_Contingency][1] )
		inputs["pNetwork_Contingency"] = convert(Int8, res_in[!,:Network_Contingency][1] )
		# Set BigM value used for dynamic contingencies cases to be largest possible cluster size
		if inputs["pDynamic_Contingency"] > 0
			inputs["pContingency_BigM"] = inputs["dfGen"][!,:Max_Cap_MW].*inputs["dfGen"][!,:Cap_size]
			for y in inputs["dfGen"][(inputs["dfGen"][!,:Commit].==1),:][!,:R_ID]
				if inputs["pContingency_BigM"][y] < 0
					# NOTE: this effectively acts as a maximum cluster size when not otherwise specified, adjust accordingly
					inputs["pContingency_BigM"][y] = 5000*inputs["dfGen"][!,:Cap_size][y]
				end
			end
		end
	end


	## Print confirmation
	print("Loaded inputs from $path$sep")

	return inputs

end # END data_input()
################################################################################


################################################################################
## function fix_integers()
##
## inputs: jump_model - a model object containing that has been previously solved.
##
## description: fixes the iteger variables ones the model has been solved in order
## to calculate approximations of dual variables
##
## returns: no result since it modifies an existing-solved model in the memory.
## solve() must be run again to solve and getdual veriables
##
################################################################################
function fix_integers(jump_model)
	N = MathProgBase.numvar(jump_model)
	for i in 1:N
		v = Variable(jump_model,i)
		if getcategory(v) != :Cont
			setcategory(v,:Cont)
			setlowerbound(v, getvalue(v))
			setupperbound(v, getvalue(v))
		end
	end
end
################################################################################



################################################################################
## function solve_model()
##
## inputs: in - a Dict() object containing all input parameters (see
## data_inputs() function).
##
## description: Sets up and solves constrained optimization model of electricity
## generation capacity expansion .&  operation problem and extracts solution variables
## for later processing
##
## returns: results - Dict() object with a set of DataFrames containing key results
##
################################################################################
function generate_power_system_model(setup::Dict,inputs::Dict,OPTIMIZER,modeloutput = nothing)
	## Start pre-solve timer
	println("Entered Method")
	presolver_start_time = time()
	println("Tock!")

	## Model Definition
	# Define the Energy Portfolio (EP) model
	# Set solver to use Gurobi
	EP=Model(solver=OPTIMIZER)

	# EP=Model(solver=CbcSolver(logLevel = 2))
   	# EP=Model(solver=AmplNLSolver("bonmin"))

	T = inputs["T"]     # Number of time steps (hours)
	Z = inputs["Z"]     # Number of zones

	@variable(EP, vPower2H2_P[t = 1:T,z=1:Z])

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
	EP, dModuleArgs = economic_dispatch(EP, dModuleArgs)

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

	if setup["single_power_sector"] == 1
		@constraint(EP, cPower2H2_P[t=1:T, z=1:Z], vPower2H2_P[t,z] == inputs["Power2H2"])

		@variable(EP, 1e+4>= vPower2H2_P_Slack_pos[t = 1:T,z=1:Z] >=0)
		@variable(EP, 1e+4>= vPower2H2_P_Slack_neg[t = 1:T,z=1:Z] >=0)

		@expression(EP, ePower2H2_relaxed_P[t=1:T, z=1:Z],
		vPower2H2_P[t,z] + vPower2H2_P_Slack_pos[t,z] - vPower2H2_P_Slack_neg[t,z]
		)
		dExpressions["ePower2H2_relaxed_P"] = ePower2H2_relaxed_P

		@expression(EP, H2PowerPenalty, sum(inputs["omega"][t] * 1e+5 * (vPower2H2_P_Slack_pos[t,z] + vPower2H2_P_Slack_neg[t,z]) for z=1:Z, t = 1:T))
		dExpressions["H2PowerPenalty"] = H2PowerPenalty
		dObjective["H2PowerPenalty"] = H2PowerPenalty
		# for t = 1:T
		# 	for z=1:Z
				# setlowerbound(vPower2H2_P[t,z], inputs["Power2H2"])
				# setupperbound(vPower2H2_P[t,z], inputs["Power2H2"])
		# 	end
		# end

		## Define the objective function
		@objective(EP,Min,dObjective["eTotalCFix"] + dObjective["eTotalCVar"] + dObjective["eTotalCNSE"] + dObjective["eTotalCStart"] + dObjective["eTotalHeat_Cost"] - dObjective["eTotalHeat_Rev"] + dObjective["eTotalCRsvPen"] + dObjective["eTotalCNetworkExp"] - dObjective["eTotalexportRev"] + dObjective["eTotalImportCost"] +
		dObjective["H2PowerPenalty"]
		)

		## Power balance constraints
		# demand = generation + storage discharge - storage charge - demand deferral + deferred demand satisfaction - demand curtailment (NSE)
		#          + incoming power flows - outgoing power flows - flow losses - charge of heat storage + generation from NACC
		@constraint(EP, cPowerBalance[t=1:T, z=1:Z], dPowerBalance["ePowerBalance"][t,z] == inputs["pD"][t,z] + ePower2H2_relaxed_P[t,z]/1e+3
		)

		## Record pre-solver time
		presolver_time = time() - presolver_start_time
	    #### Question - What do we do with this time now that we've split this function into 2?
		# if setup["PrintModel"] == 1
		# 	if modeloutput == nothing
		# 		filepath = joinpath(pwd(), "YourModel.lp")
		# 		JuMP.writeLP(EP, filepath; genericnames=false)
		# 	else
		# 		filepath = joinpath(modeloutput, "YourModel.lp")
		# 		JuMP.writeLP(EP, filepath; genericnames=false)
		# 	end
		# 	println("Model Printed")
	    # end
	else
		if setup["CO2_price_option"] == 1
			dObjective["eTotalCO2EmissionCost_Power"] = dExpressions["eCO2Emissions_power"] * inputs["CO2_price"]/1e+6
		else
			dObjective["eTotalCO2EmissionCost"] = 0
		end

		eObjective = dObjective["eTotalCFix"] + dObjective["eTotalCVar"] + dObjective["eTotalCNSE"] + dObjective["eTotalCStart"] + dObjective["eTotalHeat_Cost"] - dObjective["eTotalHeat_Rev"] + dObjective["eTotalCRsvPen"] + dObjective["eTotalCNetworkExp"] - dObjective["eTotalexportRev"] + dObjective["eTotalImportCost"] + dObjective["eTotalCO2EmissionCost_Power"]

		dExpressions["eObjective"] = eObjective

		dPowerBalance["ePowerBalance"] = dPowerBalance["ePowerBalance"] - inputs["pD"]

		@objective(EP,Min, dExpressions["eObjective"])
		@constraint(EP, cPowerBalance[t=1:T, z=1:Z], dPowerBalance["ePowerBalance"][t,z] == 0)
	end

	dModuleArgs["dExpressions"] = dExpressions
	dModuleArgs["dObjective"] = dObjective
	dModuleArgs["dPowerBalance"] = dPowerBalance
    return EP, dModuleArgs
end

    #### This is where we want to split this funtion into two

################################################################################
## function solve_model()
##
## inputs: in - a JuMP model representing the energy optimization problem
##
##
## description: Solves and extracts solution variables for later processing
##
##
## returns: results - Dict() object with a set of DataFrames containing key results
##
################################################################################
function solve_power_system_model(EP::Model, setup::Dict, inputs::Dict, dModuleArgs::Dict)

	expressions = dModuleArgs["dExpressions"]
	dObjective = dModuleArgs["dObjective"]

	## Start solve timer
	solver_start_time = time()

	## Solve Model
	status = solve(EP)
	## Record solver time
	solver_time = time() - solver_start_time

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
	dfStatus = DataFrame(Status = status, Solve = solver_time,
		Objval = getobjectivevalue(EP), Objbound= getobjbound(EP),FinalMIPGap =(getobjectivevalue(EP) -getobjbound(EP))/getobjectivevalue(EP) )

	## Check if solved sucessfully - time out is included
	if status!=:Optimal
		if status!=:UserLimit # Model failed to solve, so record solver status and exi
			results = Dict([("STATUS", dfStatus)])
			return results
			# Model reached timelimit but failed to find a feasible solution
		elseif isnan(getobjectivevalue(EP))==true
				# Model failed to solve, so record solver status and exit
			results = Dict([("STATUS", dfStatus)])
			return results
		end
	end


	## Cost results
	if setup["HeatMarket"]==1
		dfCost = DataFrame(Costs = ["cTotal", "cFix", "cVar", "cNSE", "cStart", "cUnmetRsv", "cNetworkExp", "cHeatRev","cHeatCost"])
		dfCost[!,Symbol("Total")] = [getobjectivevalue(EP), getvalue(dObjective["eTotalCFix"]), getvalue(dObjective["eTotalCVar"]), getvalue(dObjective["eTotalCNSE"]),
								 0, 0, 0, (if dObjective["eTotalHeat_Rev"] == 0 0 else -getvalue(dObjective["eTotalHeat_Rev"]) end), (if dObjective["eTotalHeat_Cost"] == 0 0 else -getvalue(dObjective["eTotalHeat_Cost"]) end)]
	else
		dfCost = DataFrame(Costs = ["cTotal", "cFix", "cVar", "cNSE", "cStart", "cUnmetRsv", "cNetworkExp"])
		dfCost[!,Symbol("Total")] = [getobjectivevalue(EP), getvalue(dObjective["eTotalCFix"]), getvalue(dObjective["eTotalCVar"]), getvalue(dObjective["eTotalCNSE"]), 0, 0, 0]
	end

	if setup["UCommit"]>=1
		dfCost[!,2][5] = getvalue(dObjective["eTotalCStart"])
	end

	if setup["Reserves"]==1
		dfCost[!,2][6] = getvalue(dObjective["eTotalCRsvPen"])
	end

	if setup["NetworkExpansion"]==1
		dfCost[!,2][7] = getvalue(dObjective["eTotalCNetworkExp"])
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
			dfCost[!,Symbol("Zone$z")] = [tempCTotal, tempCFix, tempCVar, sum(getvalue(expressions["eCNSE"])[:,:,z]), tempCStart, "-", "-", "-","-"]
		else
			dfCost[!,Symbol("Zone$z")] = [tempCTotal, tempCFix, tempCVar, sum(getvalue(expressions["eCNSE"])[:,:,z]), tempCStart, "-", "-"]
		end
	end

	## Extract decision variables
	# Capacity decisions
	capdischarge = zeros(size(inputs["RESOURCES"]))
	capcharge = zeros(size(inputs["RESOURCES"]))
	retcapcharge = zeros(size(inputs["RESOURCES"]))
	aux_i = (dfGen[(dfGen[!,:STOR].==3) .| (dfGen[!,:STOR].==5) ,:][!,:R_ID])
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

	# Storage level (state of charge) of each resource in each time step
	dfStorage = DataFrame(Resource = inputs["RESOURCES"], Zone = dfGen[!,:zone])
	dfStorage = hcat(dfStorage, convert(DataFrame, getvalue(EP[:vS])))
	auxNew_Names=[Symbol("Resource");Symbol("Zone");[Symbol("t$t") for t in 1:T]]
	rename!(dfStorage,auxNew_Names)

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
	dfNse = DataFrame()
	for z in 1:Z
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
		# Dividing dual variable for each hour with corresponding hourly weight to retrieve marginal cost of generation
		dfPrice = hcat(dfPrice, convert(DataFrame, transpose(getdual(EP[:cPowerBalance])./inputs["omega"])))
		# dfPrice = hcat(dfPrice, convert(DataFrame, transpose(getdual(EP[:cPowerBalance]))))
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
			dfEmissions[!,:Sum][i] = sum(getvalue(expressions["eEmissionsByZone"])[i,:])
		end
		dfEmissions = hcat(dfEmissions, convert(DataFrame, getvalue(expressions["eEmissionsByZone"])))
		auxNew_Names=[Symbol("Zone");Symbol("CO2_Price");Symbol("Sum");[Symbol("t$t") for t in 1:T]]
		rename!(dfEmissions,auxNew_Names)

		# # Dual of storage level (state of charge) balance of each resource in each time step
		dfStorageDual = DataFrame(Resource = inputs["RESOURCES"], Zone = dfGen[!,:zone])
		# Define an empty
		x1 = Array{Float64}(undef, G, T)

		for y in 1:G
			 if y in dfGen[(dfGen[!,:STOR].==1) .| (dfGen[!,:STOR].==2) .| (dfGen[!,:STOR].==3),:][!,:R_ID]
		### # script to write duals
		 		x1[y,:] =getdual(expressions["SoCBal"][:,y])
			else
				x1[y,:] = zeros(T,1) # Empty values for the resource with no ability to store energy
			end
		end

		dfStorageDual=hcat(dfStorageDual, convert(DataFrame,x1))
		# auxNew_Names=[Symbol("Resource");Symbol("Zone");[Symbol("t$t") for t in 1:T]]
		# rename!(dfStorageDual,auxNew_Names)
		rename!(dfStorageDual,[Symbol("Resource");Symbol("Zone");[Symbol("t$t") for t in 1:T]])

	elseif (setup["Dual_MIP"]==1)
		# function to fix integers and linearize problem
		fix_integers(EP)
		# re-solve statement for LP solution
		println("Solving LP solution for duals")
		solve(EP)
		## Extract dual variables of constraints
		# price: Dual variable of hourly power balance constraint = hourly price
		dfPrice = DataFrame(Zone = 1:Z)
		dfPrice = hcat(dfPrice, convert(DataFrame, transpose(getdual(EP[:cPowerBalance])./inputs["omega"])))
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
			dfEmissions[!,:Sum][i] = sum(getvalue(expressions["eEmissionsByZone"])[i,:])
		end
		dfEmissions = hcat(dfEmissions, convert(DataFrame, getvalue(expressions["eEmissionsByZone"])))
		auxNew_Names=[Symbol("Zone");Symbol("CO2_Price");Symbol("Sum");[Symbol("t$t") for t in 1:T]]
		rename!(dfEmissions,auxNew_Names)

		# # Dual of storage level (state of charge) balance of each resource in each time step
		dfStorageDual = DataFrame(Resource = inputs["RESOURCES"], Zone = dfGen[!,:zone])
		# Define an empty
		x1 = Array{Float64}(undef, G, T)

		for y in 1:G
			 if y in dfGen[(dfGen[!,:STOR].==1) .| (dfGen[!,:STOR].==2) .| (dfGen[!,:STOR].==3),:][!,:R_ID]
		### # script to write duals
		 		x1[y,:] =getdual(expressions["SoCBal"][:,y])
			else
				x1[y,:] = zeros(T,1) # Empty values for the resource with no ability to store energy
			end
		end

		dfStorageDual=hcat(dfStorageDual, convert(DataFrame,x1 ))
		# auxNew_Names=[Symbol("Resource");Symbol("Zone");[Symbol("t$t") for t in 1:T]]
		# rename!(dfStorageDual,auxNew_Names)
		rename!(dfStorageDual,[Symbol("Resource");Symbol("Zone");[Symbol("t$t") for t in 1:T]])

	else
		# CO2 emissions by zone
		dfEmissions = DataFrame(Zone = 1:Z, CO2_Price = convert(Array{Union{Missing,Float32}}, zeros(Z)), Sum = Array{Union{Missing,Float32}}(undef, Z))
		for i in 1:Z
			dfEmissions[!,:Sum][i] = sum(getvalue(expressions["eEmissionsByZone"])[i,:])
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
################################################################################

################################################################################
## function dftranspose(df)
##
## inputs: df - [DataFrame] a DataFrame object to be transposed
## results t - [DataFrame] a transposed version of the df DataFrame.
##		   withhead - [Boolean] if True, first column of df will become column
##		   names for t. Otherwise, first column will first row and column names
##		   will be generic (e.g. x1:xN)
##
## Note this function is necessary because no stock function to transpose
## DataFrames appears to exist.
################################################################################
function dftranspose(df::DataFrame, withhead::Bool)
	# Extract old column names (as Symbols)
	oldnames_temp = names(df)
	# Convert to String vector and save for use as new row names
	oldnames = Vector{Union{Nothing,String}}(nothing, length(oldnames_temp))
	for r in 1:length(oldnames_temp)
		oldnames[r] = String(oldnames_temp[r])
	end
	if(withhead)
		# Extract first row of data frame (Resources names) (as Strings) and save as new column names
		newnames = string.(df[:,1])
		startcol=2
	else
		startcol=1
	end
	# Collect each of the old columns and tranpose to new rows
	t = DataFrame(permutedims(df[:,startcol]))
	for c in (startcol+1):ncol(df)
		t = vcat(t,DataFrame(permutedims(df[:,c])))
	end
	# Set new column names
	if(withhead)
		t = DataFrame(t,Symbol(newnames[c]))
	end
	# Add new row names vector to data frame
	t = hcat(DataFrame(Row=oldnames[startcol:length(oldnames)]), t)
	# Return transposed data frame
	return t
end # End dftranpose()
################################################################################

################################################################################
## function output(path, results)
##
## inputs: path - [string] path to working directory
## results - results dict data structure containing results to write to files.
##           See model() function for description of results dict.
##
## description: Writes results to multiple .csv output files in path directory
##
## returns: n/a
################################################################################
function write_power_outputs(setup::Dict, path::AbstractString, results::Dict, inputs::Dict)

	## Use appropriate directory separator depending on Mac or Windows config
	if setup["MacOrWindows"]=="Mac"
		sep = "/"
	else
		sep = "\U005c"
	end

	# If output directory does not exist, create it
	if !(isdir(path))
		mkdir(path)
	end

	CSV.write(string(path,sep,"status.csv"),results["STATUS"])


	# if results["STATUS"][!,:Status][1]!=:Optimal
	# 	if results["STATUS"][!,:Status][1]!=:UserLimit
	# 	# Model failed to solve succesfully, so exit here
	# 		println("Wrote outputs to $path$sep - Model not solved to optimality")
	# 		return
	# 	elseif isnan(results["STATUS"][!,Objval][1])==true # Model reached timelimit but failed to find a feasible solution
	# 		println("Wrote outputs to $path$sep - Model not solved to optimality")
	# 		return
	# 	end
	# end


	Z = inputs["Z"]   # Number of zones
	T = inputs["T"]   # Number of time periods
	CSV.write(string(path,sep,"costs.csv"), results["COSTS"])
	CSV.write(string(path,sep,"capacity.csv"), results["CAP"])
	CSV.write(string(path,sep,"power.csv"), dftranspose(results["POWER"], false), writeheader=false)
	CSV.write(string(path,sep,"charge.csv"), dftranspose(results["CHARGE"], false), writeheader=false)
	CSV.write(string(path,sep,"storage.csv"), dftranspose(results["STORAGE"], false), writeheader=false)
	CSV.write(string(path,sep,"nse.csv"),  dftranspose(results["NSE"], false), writeheader=false)
	if Z > 1
		CSV.write(string(path,sep,"flow.csv"), dftranspose(results["FLOWS"], false), writeheader=false)
		CSV.write(string(path,sep,"tlosses.csv"), dftranspose(results["TLOSSES"], false), writeheader=false)
	end
	CSV.write(string(path,sep,"wlosses.csv"), dftranspose(results["WLOSSES"], false), writeheader=false)
	CSV.write(string(path,sep,"injections.csv"), dftranspose(results["INJECT"], false), writeheader=false)
	CSV.write(string(path,sep,"withdrawals.csv"), dftranspose(results["WITHDRAW"], false), writeheader=false)
	CSV.write(string(path,sep,"emissions.csv"), dftranspose(results["EMISSIONS"], false), writeheader=false)

	if setup["HeatMarket"]==1
		## CHP system configuration final output

		CSV.write(string(path,sep,"heat_market.csv"), dftranspose(results["HEAT_MARKET"], false), writeheader=false)
		CSV.write(string(path,sep,"heat_gen.csv"), dftranspose(results["HEAT_GEN"], false), writeheader=false)
		CSV.write(string(path,sep,"heat_ng.csv"), dftranspose(results["HEAT_NG"], false), writeheader=false)

		CSV.write(string(path,sep,"heat_to_storage.csv"), dftranspose(results["HEAT_TO_STORAGE"], false), writeheader=false)
		CSV.write(string(path,sep,"heat_to_topping.csv"), dftranspose(results["HEAT_TO_TOPPING"], false), writeheader=false)
		CSV.write(string(path,sep,"heat_to_basegen.csv"), dftranspose(results["HEAT_TO_BASEGEN"], false), writeheader=false)
	end #END CHP system configuration


	if (setup["UCommit"]==0)
		## Linear configuration final output

		CSV.write(string(path,sep,"prices.csv"), dftranspose(results["PRICE"], false), writeheader=false)
		CSV.write(string(path,sep,"reliabilty.csv"), dftranspose( results["RELIABILITY"], false), writeheader=false)
		CSV.write(string(path,sep,"storagebal_duals.csv"), dftranspose(results["STORAGEDUAL"], false), writeheader=false)
	elseif (setup["UCommit"]==2)

		## Linear configuration final output
		CSV.write(string(path,sep,"prices.csv"), dftranspose(results["PRICE"], false), writeheader=false)
		CSV.write(string(path,sep,"reliabilty.csv"), dftranspose(results["RELIABILITY"], false), writeheader=false)
		CSV.write(string(path,sep,"storagebal_duals.csv"), dftranspose(results["STORAGEDUAL"], false), writeheader=false)
		## Unit Commitment configuration final output

		CSV.write(string(path,sep,"commit.csv"), dftranspose(results["COMMIT"], false), writeheader=false)
		CSV.write(string(path,sep,"start.csv"), dftranspose(results["START"], false), writeheader=false)
		CSV.write(string(path,sep,"shutdown.csv"), dftranspose(results["SHUTDOWN"], false), writeheader=false)

		if setup["Reserves"]==1
			CSV.write(string(path,sep,"reg_up.csv"), dftranspose(results["REG_UP"], false), writeheader=false)
			CSV.write(string(path,sep,"reg_dn.csv"), dftranspose(results["REG_DN"], false), writeheader=false)
			CSV.write(string(path,sep,"rsv_up.csv"), dftranspose(results["RSV_UP"], false), writeheader=false)
			CSV.write(string(path,sep,"rsv_dn.csv"), dftranspose(results["RSV_DN"], false), writeheader=false)
		end
	else
		## Unit Commitment configuration final output
		CSV.write(string(path,sep,"commit.csv"), dftranspose(results["COMMIT"], false), writeheader=false)
		CSV.write(string(path,sep,"start.csv"), dftranspose(results["START"], false), writeheader=false)
		CSV.write(string(path,sep,"shutdown.csv"),dftranspose( results["SHUTDOWN"], false), writeheader=false)
		if setup["Reserves"]==1
			CSV.write(string(path,sep,"reg_up.csv"), dftranspose(results["REG_UP"], false), writeheader=false)
			CSV.write(string(path,sep,"reg_dn.csv"), dftranspose(results["REG_DN"], false), writeheader=false)
			CSV.write(string(path,sep,"rsv_up.csv"), dftranspose(results["RSV_UP"], false), writeheader=false)
			CSV.write(string(path,sep,"rsv_dn.csv"), dftranspose(results["RSV_DN"], false), writeheader=false)
		end

		if setup["Dual_MIP"]==1 # IF user requested duals of LP version of the problem - provde them
			println("exporting prices")
			CSV.write(string(path,sep,"prices.csv"), dftranspose(results["PRICE"], false), writeheader=false)
			CSV.write(string(path,sep,"storagebal_duals.csv"), dftranspose(results["STORAGEDUAL"], false), writeheader=false)
		end

	end #END unit commitment configuration

	if setup["NetworkExpansion"]==1
		CSV.write(string(path,sep,"network_expansion.csv"), results["TRANS_CAP"])
		if sum(inputs["pDistrZones"]) > 0
			CSV.write(string(path,sep,"distribution_expansion.csv"), results["DISTR_CAP"])
		end
	end

	if (setup["OperationWrapping"]==2 && setup["LDS"]==1)
		CSV.write(string(path,sep,"StorageInit.csv"), dftranspose(results["dfStorageInit"], false), writeheader=false)
		CSV.write(string(path,sep,"dStorage.csv"), dftranspose(results["dfdStorage"], false), writeheader=false)
	end

	## Print confirmation
	println("Wrote outputs to $path$sep")

end # END output()
################################################################################

# Begin Infrastructure Functions

function economic_dispatch(EP::Model, dModuleArgs::Dict)

	println("Economic Dispatch Module")

	inputs = dModuleArgs["inputs"]
	setup = dModuleArgs["setup"]
	dExpressions = dModuleArgs["dExpressions"]
	dObjective = dModuleArgs["dObjective"]
	dPowerBalance = dModuleArgs["dPowerBalance"]

	dfGen = inputs["dfGen"]

	G = inputs["G"]     # Number of resources (generators, storage, DR, and DERs)
	T = inputs["T"]     # Number of time steps (hours)
	Z = inputs["Z"]     # Number of zones
	SEG = inputs["SEG"] # Number of load curtailment segments

	### Variables ###

	# Energy injected into the grid by resource "y" at hour "t" [MWh]
	@variable(EP, vP[y=1:G,t=1:T] >=0);
	# Energy withdrawn from grid by resource "y" at hour "t" [MWh] on zone "z"
	@variable(EP, vCHARGE[y=1:G,t=1:T] >= 0);
	# Non-served energy/curtailed demand in the segment "s" at hour "t" [MWh] on zone "z"
	@variable(EP, vNSE[s=1:SEG,t=1:T,z=1:Z] >= 0);

	# Storage level of resource "y" at hour "t" [MWh] on zone "z" - unbounded
	@variable(EP, vS[y=1:G,t=1:T]);

	### End Variables ###

	### Expressions ###

	#Variable costs of "generation" for resource "y" during hour "t" = variable O&M plus fuel cost
	if setup["Fuel_Price_Var"]==1
		@expression(EP, eCVar_out[y=1:G,t=1:T], (inputs["omega"][t]*(dfGen[!,:Var_OM_cost_per_MWh][y]+dfGen[!,:C_Fuel_per_MWh][y]*inputs["dfFuelsVar"][!,Symbol(dfGen[!,:Fuel][y])][t])*vP[y,t]))
	else
		@expression(EP, eCVar_out[y=1:G,t=1:T], (inputs["omega"][t]*(dfGen[!,:Var_OM_cost_per_MWh][y]+dfGen[!,:C_Fuel_per_MWh][y])*vP[y,t]))
	end
	dExpressions["eCVar_out"] = eCVar_out

	#NOTE: are we going to assume some cost for charging and discharging at the heat storage
	#Variable costs of "charging" for technologies "y" during hour "t" in zone "z"
	@expression(EP, eCVar_in[y=1:G,t=1:T], inputs["omega"][t]*dfGen[!,:Var_OM_cost_per_MWh_in][y]*vCHARGE[y,t])
	dExpressions["eCVar_in"] = eCVar_in

	# energy losses related to technologies (increase in effective demand)
	@expression(EP, eELOSS[y=1:G], 0)
	for y in (dfGen[(dfGen[!,:STOR].>=1),:][!,:R_ID])
		eELOSS[y] = sum(inputs["omega"][t]*vCHARGE[y,t] for t in 1:T) - sum(inputs["omega"][t]*vP[y,t] for t in 1:T)
	end
	dExpressions["eELOSS"] = eELOSS

	# Cost of non-served energy/curtailed demand at hour "t" in zone "z"
	@expression(EP, eCNSE[s=1:SEG,t=1:T,z=1:Z], (inputs["omega"][t]*inputs["pC_D_Curtail"][s]*vNSE[s,t,z]))
	dExpressions["eCNSE"] = eCNSE

	## Objective Function Expressions ##

	@expression(EP, eTotalCNSE, sum(eCNSE[s,t,z] for s in 1:SEG,t in 1:T,z in 1:Z) )
    dObjective["eTotalCNSE"] = eTotalCNSE

	@expression(EP, eTotalCVar, sum(eCVar_out[y,t] + eCVar_in[y,t] for y in 1:G, t in 1:T) )
	dObjective["eTotalCVar"] = eTotalCVar

	## End Objective Function Expressions ##

	## Power Balance Expressions ##

	# Expression for "baseline" power balance constraint?
	@expression(EP, ePowerBalance[t=1:T, z=1:Z], 0)
	dPowerBalance["ePowerBalance"] = ePowerBalance

	@expression(EP, ePowerBalanceNse[t=1:T, z=1:Z],
		sum(vNSE[s,t,z] for s=1:SEG))

	# Add to power balance expression
	dPowerBalance["ePowerBalance"] = dPowerBalance["ePowerBalance"] + ePowerBalanceNse

	## End Power Balance Expressions ##

	### End Expressions ###

	### Constratints ###

	## Non-served energy (demand curtailment) constraints
	# Demand curtailed in each segment of curtailable demands cannot exceed maximum allowable share of demand
	@constraint(EP, cNSEPerSeg[s=1:SEG, t=1:T, z=1:Z], vNSE[s,t,z] <= inputs["pMax_D_Curtail"][s]*inputs["pD"][t,z])
	# Total demand curtailed in each time step (hourly) cannot exceed total demand
	@constraint(EP, cMaxNSE[t=1:T, z=1:Z], sum(vNSE[s,t,z] for s=1:SEG) <= inputs["pD"][t,z])

	### End Constraints ###

	dModuleArgs["dExpressions"] = dExpressions
	dModuleArgs["dObjective"] = dObjective
	dModuleArgs["dPowerBalance"] = dPowerBalance

	return EP, dModuleArgs
end

function investment(EP::Model, dModuleArgs::Dict)

	println("Investment Module")

	inputs = dModuleArgs["inputs"]
	setup = dModuleArgs["setup"]
	dExpressions = dModuleArgs["dExpressions"]
	dObjective = dModuleArgs["dObjective"]
	dPowerBalance = dModuleArgs["dPowerBalance"]

	dfGen = inputs["dfGen"]

	G = inputs["G"]     # Number of resources (generators, storage, DR, and DERs)

	### Variables ###

	# Retired capacity of resource "y" from existing capacity [MW]
    @variable(EP, vRETCAP[y=1:G] >= 0);
    # Installed capacity of resource "y" capacity [MW]
	@variable(EP, vCAP[y=1:G] >= 0);

	# Energy storage reservoir capacity (MWh capacity) built/retired for storage with variable power to energy ratio (STOR=2 or STOR=3)
	@variable(EP, vCAPSTORAGE[y=(dfGen[(dfGen[!,:STOR].>=2),:][!,:R_ID])] >= 0)
	@variable(EP, vRETCAPSTORAGE[y=(dfGen[(dfGen[!,:STOR].>=2),:][!,:R_ID])] >= 0)
	# Storage capacity built and retired for storage resources with independent charge and discharge power capacities (STOR=3)
	# Note that discharge capacity built/retired for such resources is designated by vCAP and vRETCAP variables
	# Retired storage capacity tracked via power capacity
	@variable(EP, vCAPCHARGE[y=(dfGen[(dfGen[!,:STOR].==3) .| (dfGen[!,:STOR].==5) ,:][!,:R_ID])] >= 0)
	@variable(EP, vRETCAPCHARGE[y=(dfGen[(dfGen[!,:STOR].==3) .| (dfGen[!,:STOR].==5) ,:][!,:R_ID])] >= 0)

	### End Variables ###

	### Expressions ###

	# Fixed costs for resource "y" = annuitized investment cost plus fixed O&M costs
	@expression(EP, eCFix[y=1:G], ( dfGen[!,:Inv_cost_per_MWyr][y]*dfGen[!,:Cap_size][y]*vCAP[y])
		+ dfGen[!,:Fixed_OM_cost_per_MWyr][y]*(dfGen[!,:Existing_Cap_MW][y] + dfGen[!,:Cap_size][y]*vCAP[y] - dfGen[!,:Cap_size][y]*vRETCAP[y])
		)

	# Special cases for fixed cost expressions follow...
	# Costs for heat storage resources
	## NOTE: currently only greenfield capacity costs considered for heat storage; need to update
	## Also appears that vCAP variable applies to storage capacity (energy) not power, as it does for all other resources
	## Rewrite for consistency.
	#NOTE: just greenfield fix cost update later
	for y in (dfGen[(dfGen[!,:HEAT].==2),:][!,:R_ID])
		eCFix[y] = ( dfGen[!,:Inv_cost_per_MWhyr][y]*dfGen[!,:Cap_size][y]*vCAP[y]
			+ (dfGen[!,:Inv_cost_per_MWyr][y] + dfGen[!,:Fixed_OM_cost_per_MWyr][y])*(vCAPCHARGE[y] + vCAPDISCHARGE[y])
		)
	end
	# Costs for storage resources with independent power/energy capacities
	# STOR = 2 corresponds to battery with variable power and energy ratios but symmetric charge/discharg power ratings
	for y in (dfGen[(dfGen[!,:STOR].==2),:][!,:R_ID])
		eCFix[y] = ( dfGen[!,:Inv_cost_per_MWyr][y]*dfGen[!,:Cap_size][y]*vCAP[y]
			+ dfGen[!,:Fixed_OM_cost_per_MWyr][y]*(dfGen[!,:Existing_Cap_MW][y] + dfGen[!,:Cap_size][y]*vCAP[y] - dfGen[!,:Cap_size][y]*vRETCAP[y])
			+ dfGen[!,:Inv_cost_per_MWhyr][y]*vCAPSTORAGE[y]
			+ dfGen[!,:Fixed_OM_cost_per_MWhyr][y]*(dfGen[!,:Existing_Cap_MWh][y] + vCAPSTORAGE[y] - vRETCAPSTORAGE[y]) )
	end


	# Costs for storage resources with independent charge and discharge power capacities and energy storage capacity
	for y in (dfGen[(dfGen[!,:STOR].==3),:][!,:R_ID])
		eCFix[y] = ( dfGen[!,:Inv_cost_per_MWyr][y]*dfGen[!,:Cap_size][y]*vCAP[y]
			+ dfGen[!,:Fixed_OM_cost_per_MWyr][y]*(dfGen[!,:Existing_Cap_MW][y] + dfGen[!,:Cap_size][y]*vCAP[y] - dfGen[!,:Cap_size][y]*vRETCAP[y])
			+ dfGen[!,:Inv_cost_per_MWhyr][y]*vCAPSTORAGE[y]
			+ dfGen[!,:Fixed_OM_cost_per_MWhyr][y]*(dfGen[!,:Existing_Cap_MWh][y] + vCAPSTORAGE[y] - vRETCAPSTORAGE[y]) +
			+ dfGen[!,:Inv_cost_charge_per_MWyr][y]*dfGen[!,:Cap_size][y]*vCAPCHARGE[y]
			+ dfGen[!,:Fixed_OM_cost_charge_per_MWyr][y]*(dfGen[!,:Existing_Charge_Cap_MW][y] + dfGen[!,:Cap_size][y]*vCAPCHARGE[y] - dfGen[!,:Cap_size][y]*vRETCAPCHARGE[y])
		)
	end

	dExpressions["eCFix"] = eCFix

	## Objective Function Expressions ##

	@expression(EP, eTotalCFix, sum(eCFix[y] for y in 1:G))
	dObjective["eTotalCFix"] = eTotalCFix

	## End Objective Function Expressions ##

	## Power Balance Expressions ##
	## End Power Balance Expressions ##

	### End Expressions ###

	### Constratints ###

	## Constraints on retirements and capacity additions
	#Cannot retire more capacity than existing capacity
	@constraint(EP, cMaxRet[y=1:G], dfGen[!,:Cap_size][y]*vRETCAP[y] <= dfGen[!,:Existing_Cap_MW][y])
	@constraint(EP, cMaxRetStorage[y=(dfGen[(dfGen[!,:STOR].>=2),:][!,:R_ID])], vRETCAPSTORAGE[y] <= dfGen[!,:Existing_Cap_MWh][y])
 	@constraint(EP, cMaxRetCharge[y=(dfGen[(dfGen[!,:STOR].==3),:][!,:R_ID])], vRETCAPCHARGE[y] <= dfGen[!,:Existing_Charge_Cap_MW][y])

  	#Constraints on new built capacity
	for y in 1:G
		# Can only install resources eligible for new capacity
		if (dfGen[!,:New_Build][y] == 0 || dfGen[!,:Max_Cap_MW][y] == 0 || dfGen[!,:New_Build][y] == -1 )
			setupperbound(vCAP[y], 0.0)
			setlowerbound(vCAP[y], 0.0)
			if (y in (dfGen[(dfGen[!,:STOR].>=2) ,:][!,:R_ID]))
				setupperbound(vCAPSTORAGE[y], 0.0)
				setlowerbound(vCAPSTORAGE[y], 0.0)
			end
 			if (y in (dfGen[(dfGen[!,:STOR].==3) ,:][!,:R_ID]))
 				setupperbound(vCAPCHARGE[y], 0.0)
 				setlowerbound(vCAPCHARGE[y], 0.0)
 			end
		end
		# If no existing capacity eliminate retirement variable or
		# if dfGen[!,:New_Build]==-1, no retirement is allowed
		if (dfGen[!,:New_Build][y] == -1 || dfGen[!,:Existing_Cap_MW][y] == 0)
			setupperbound(vRETCAP[y], 0.0)
			setlowerbound(vRETCAP[y], 0.0)
			if (y in (dfGen[(dfGen[!,:STOR].>=2) ,:][!,:R_ID]))
				setupperbound(vRETCAPSTORAGE[y], 0.0)
				setlowerbound(vRETCAPSTORAGE[y], 0.0)
			end
 			if (y in (dfGen[(dfGen[!,:STOR].==3) ,:][!,:R_ID]))
 				setupperbound(vRETCAPCHARGE[y], 0.0)
 				setlowerbound(vRETCAPCHARGE[y], 0.0)
 			end
		end

		# Constraint on maximum capacity (if applicable) [set input to -1 if no constraint on maximum capacity]
		if dfGen[!,:Max_Cap_MW][y] > 0
			@constraint(EP, dfGen[!,:Cap_size][y]*(dfGen[!,:Existing_Cap_MW][y] + vCAP[y] - vRETCAP[y]) <= dfGen[!,:Max_Cap_MW][y])
			if (y in (dfGen[(dfGen[!,:STOR].==3) ,:][!,:R_ID]))
				# For storage with independent charge/discharge power capacities, max charge capacity is constrainted to same value as max discharge capacity
				@constraint(EP, dfGen[!,:Cap_size][y]*(dfGen[!,:Existing_Cap_MW][y] + vCAPCHARGE[y] - vRETCAPCHARGE[y]) <= dfGen[!,:Max_Cap_MW][y])
			end
		end

		# Constraint on minimum capacity (if applicable) [set input to -1 if no constraint on minimum capacity]
		if dfGen[!,:Min_Cap_MW][y] > 0
			@constraint(EP, dfGen[!,:Cap_size][y]*(dfGen[!,:Existing_Cap_MW][y] + vCAP[y] - vRETCAP[y]) >= dfGen[!,:Min_Cap_MW][y])
			if (y in (dfGen[(dfGen[!,:STOR].==3) ,:][!,:R_ID]))
				# For storage with independent charge/discharge power capacities, min charge capacity is constrainted to same value as min discharge capacity
				@constraint(EP, dfGen[!,:Cap_size][y]*(dfGen[!,:Existing_Cap_MW][y] + vCAPCHARGE[y] - vRETCAPCHARGE[y]) >= dfGen[!,:Min_Cap_MW][y])
			end
		end
	end

	### End Constraints ###

	dModuleArgs["dExpressions"] = dExpressions
	dModuleArgs["dObjective"] = dObjective
	dModuleArgs["dPowerBalance"] = dPowerBalance

	return EP, dModuleArgs
end

function ucommit(EP::Model, dModuleArgs::Dict)

	println("Unit Commitment Module")

	inputs = dModuleArgs["inputs"]
	setup = dModuleArgs["setup"]
	dExpressions = dModuleArgs["dExpressions"]
	dObjective = dModuleArgs["dObjective"]
	dPowerBalance = dModuleArgs["dPowerBalance"]

	dfGen = inputs["dfGen"]

	G = inputs["G"]     # Number of resources (generators, storage, DR, and DERs)
	T = inputs["T"]     # Number of time steps (hours)

	### Variables ###

	## Integer Unit Commitment configuration for variables
	if (setup["UCommit"]>=1)
		## Decision variables for unit commitment
		# commitment state variable
		@variable(EP, vCOMMIT[y=1:G,t=1:T] >= 0)
		# startup event variable
		@variable(EP, vSTART[y=1:G,t=1:T] >= 0)
		# shutdown event variable
		@variable(EP, vSHUT[y=1:G,t=1:T] >= 0)
	else
		dfGen[!,:Cap_size].=1
	end

	### End Variables ###

	### Expressions ###

	## Objective Function Expressions ##

	if (setup["UCommit"]>=1)
		# Startup costs of "generation" for resource "y" during hour "t"
		@expression(EP, eCStart[y=1:G,t=1:T],(inputs["omega"][t]*dfGen[!,:C_Start][y]*vSTART[y,t]))
		dExpressions["eCStart"] = eCStart

		@expression(EP, eTotalCStart, sum(eCStart[y,t] for y=1:G, t=1:T))
	else
		@expression(EP, eTotalCStart, 0)
	end
	dObjective["eTotalCStart"] = eTotalCStart

	## End Objective Function Expressions ##

	## Power Balance Expressions ##
	## End Power Balance Expressions ##

	### End Expressions ###

	### Constratints ###
	if (setup["UCommit"]>=1)
		## Declaration of integer/binary variables
		for t in 1:T
			if (setup["UCommit"]==1) # Integer UC constraints
					for y in dfGen[(dfGen[!,:Commit].==1),:][!,:R_ID]
					setcategory(vCOMMIT[y,t], :Int )
					setcategory(vSTART[y,t], :Int )
					setcategory(vSHUT[y,t], :Int )
					setcategory(EP[:vRETCAP][y], :Int )
					setcategory(EP[:vCAP][y], :Int )
					end
			else # Linear relaxation of UC constraints
				for y in dfGen[(dfGen[!,:Commit].==1),:][!,:R_ID]
					setcategory(vCOMMIT[y,t], :Cont )
					setcategory(vSTART[y,t], :Cont )
					setcategory(vSHUT[y,t], :Cont )
					setcategory(EP[:vRETCAP][y], :Cont )
					setcategory(EP[:vCAP][y], :Cont )
					end
			end
				for y in dfGen[(dfGen[!,:Commit].==0),:][!,:R_ID]
					setlowerbound(vCOMMIT[y,t], 0.0)
					setlowerbound(vSTART[y,t], 0.0)
					setlowerbound(vSHUT[y,t], 0.0)
					setupperbound(vCOMMIT[y,t], 0.0)
					setupperbound(vSTART[y,t], 0.0)
					setupperbound(vSHUT[y,t], 0.0)
				end
		end #END unit commitment configuration
	end

	### End Constraints ###

	dModuleArgs["dExpressions"] = dExpressions
	dModuleArgs["dObjective"] = dObjective
	dModuleArgs["dPowerBalance"] = dPowerBalance

	return EP, dModuleArgs
end

function reserves(EP::Model, dModuleArgs::Dict)

	println("Reserves Module")

	inputs = dModuleArgs["inputs"]
	setup = dModuleArgs["setup"]
	dExpressions = dModuleArgs["dExpressions"]
	dObjective = dModuleArgs["dObjective"]
	dPowerBalance = dModuleArgs["dPowerBalance"]

	dfGen = inputs["dfGen"]

	G = inputs["G"]     # Number of resources (generators, storage, DR, and DERs)
	T = inputs["T"]     # Number of time steps (hours)
	Z = inputs["Z"]     # Number of zones

	### Variables ###

	## Integer Unit Commitment configuration for variables
	if (setup["UCommit"]>=1)
		## Decision variables for reserves
		if (setup["Reserves"]==1)
			@variable(EP, vREG_UP[y=1:G,t=1:T] >= 0) # Contribution to regulation (primary reserves) up
			@variable(EP, vREG_DN[y=1:G,t=1:T] >= 0) # Contribution to regulation (primary reserves) down
			@variable(EP, vRSV_UP[y=1:G,t=1:T] >= 0) # Contribution to operating reserves (secondary reserves) up
			@variable(EP, vRSV_DN[y=1:G,t=1:T] >= 0) # Contribution to operating reserves (secondary reserves) down

				# Storage techs have two pairs of auxilary variables to reflect contributions to regulation and reserves
				# when charging and discharging (primary variable becomes equal to sum of these auxilary variables)
			@variable(EP, vREG_UP_discharge[y=dfGen[(dfGen[!,:STOR].>=1),:][!,:R_ID],t=1:T] >= 0) # Contribution to regulation (primary reserves) up (mirrored variable used for storage devices)
			@variable(EP, vREG_DN_discharge[y=dfGen[(dfGen[!,:STOR].>=1),:][!,:R_ID],t=1:T] >= 0) # Contribution to regulation (primary reserves) down (mirrored variable used for storage devices)
			@variable(EP, vRSV_UP_discharge[y=dfGen[(dfGen[!,:STOR].>=1),:][!,:R_ID],t=1:T] >= 0) # Contribution to operating reserves (secondary reserves) up (mirrored variable used for storage devices)
			@variable(EP, vRSV_DN_discharge[y=dfGen[(dfGen[!,:STOR].>=1),:][!,:R_ID],t=1:T] >= 0) # Contribution to operating reserves (secondary reserves) down (mirrored variable used for storage devices)
			@variable(EP, vREG_UP_charge[y=dfGen[(dfGen[!,:STOR].>=1),:][!,:R_ID],t=1:T] >= 0) # Contribution to regulation (primary reserves) up (mirrored variable used for storage devices)
			@variable(EP, vREG_DN_charge[y=dfGen[(dfGen[!,:STOR].>=1),:][!,:R_ID],t=1:T] >= 0) # Contribution to regulation (primary reserves) down (mirrored variable used for storage devices)
			@variable(EP, vRSV_UP_charge[y=dfGen[(dfGen[!,:STOR].>=1),:][!,:R_ID],t=1:T] >= 0) # Contribution to operating reserves (secondary reserves) up (mirrored variable used for storage devices)
			@variable(EP, vRSV_DN_charge[y=dfGen[(dfGen[!,:STOR].>=1),:][!,:R_ID],t=1:T] >= 0) # Contribution to operating reserves (secondary reserves) down (mirrored variable used for storage devices)


			@variable(EP, vUNMET_RSV_UP[t=1:T] >= 0) # Unmet operating reserves up
			@variable(EP, vUNMET_RSV_DN[t=1:T] >= 0) # Unmet operating reserves down
			if (inputs["pDynamic_Contingency"] == 2)
				if (inputs["pNetwork_Contingency"] == 1)
					# Contingency = largest committed thermal unit in each time period or largest transmission line
					@variable(EP, vLARGEST_CONTINGENCY[t=1:T] >= maximum(inputs["pTrans_Max"]))
				else
					# Contingency = largest committed thermal unit in each time period
					@variable(EP, vLARGEST_CONTINGENCY[t=1:T] >= 0)
				end
				# Auxiliary variable that is 0 if vCOMMIT = 0, 1 otherwise
				@variable(EP, vCONTINGENCY_AUX[y=1:G, t=1:T], Bin)
			elseif (inputs["pDynamic_Contingency"] == 1)
				if (inputs["pNetwork_Contingency"] == 1)
					# Contingency = largest installed thermal unit or largest transmission line
					@variable(EP, vLARGEST_CONTINGENCY[1] >= maximum(inputs["pTrans_Max"]))
				else
					# Contingency = largest installed thermal unit
					@variable(EP, vLARGEST_CONTINGENCY[1] >= 0)
				end
				# Auxiliary variable that is 0 if vCAP = 0, 1 otherwise
				@variable(EP, vCONTINGENCY_AUX[y=1:G], Bin)
			end
			# NOTE: If Dynamic_Contingency == 0, then contingency is fixed parameter equal to capacity
			# of largest possible thermal unit or largest transmission line, whichever is greater (if Network_Contingency=1)
			# or the largest possible thermal unit (if Network_Contingency=0)
		end
	end

	### End Variables ###

	### Expressions ###
	## Total system reserve expressions (if active)
	if (setup["Reserves"]==1)
		# Regulation requirements as a percentage of load and scheduled variable renewable energy production in each hour
		# Reg up and down requirements are symmetric
		@expression(EP, eRegReq[t=1:T], inputs["pReg_Req_Load"]*sum(inputs["pD"][t,z] for z=1:Z) +
			inputs["pReg_Req_VRE"]*sum(EP[:vCAP][y]*inputs["pP_Max"][y,t] for y=dfGen[(dfGen[!,:DISP].==1) .| (dfGen[!,:NDISP].==1),:][!,:R_ID]) )
		# Operating reserve up requirements as a percentage of load and scheduled variable renewable energy production in each hour
		# and the largest single contingency
		if (inputs["pDynamic_Contingency"]==2)
			@expression(EP, eRsvUpReq[t=1:T], inputs["pRsv_Up_Req_Load"]*sum(inputs["pD"][t,z], z=1:Z) +
					inputs["pRsv_Up_Req_VRE"]*sum(EP[:vCAP][y]*inputs["pP_Max"][y,t], y=dfGen[(dfGen[!,:DISP].==1) .| (dfGen[!,:NDISP].==1),:][!,:R_ID]) +
					vLARGEST_CONTINGENCY[t] )
		elseif (inputs["pDynamic_Contingency"]==1)
			@expression(EP, eRsvUpReq[t=1:T], inputs["pRsv_Up_Req_Load"]*sum(inputs["pD"][t,z] for z=1:Z) +
				inputs["pRsv_Up_Req_VRE"]*sum(EP[:vCAP][y]*inputs["pP_Max"][y,t] for y=dfGen[(dfGen[!,:DISP].==1) .| (dfGen[!,:NDISP].==1),:][!,:R_ID]) +
				vLARGEST_CONTINGENCY[1] )
		else
			if (inputs["pNetwork_Contingency"]==1)
				# Largest contingency defined fixed as max of either largest thermal generator eligible for construction or largest transmission line
				@expression(EP, eRsvUpReq[t=1:T], inputs["pRsv_Up_Req_Load"]*sum(inputs["pD"][t,z] for z=1:Z) +
					inputs["pRsv_Up_Req_VRE"]*sum(EP[:vCAP][y]*inputs["pP_Max"][y,t] for y=dfGen[(dfGen[!,:DISP].==1) .| (dfGen[!,:NDISP].==1),:][!,:R_ID]) +
					max(maximum(inputs["dfGen"][!,:Cap_size]), maximum(inputs["pTrans_Max"])) )
			else
				# Largest contingency defined fixed as largest thermal generator eligible for construction
				@expression(EP, eRsvUpReq[t=1:T], inputs["pRsv_Up_Req_Load"]*sum(inputs["pD"][t,z] for z=1:Z) +
					inputs["pRsv_Up_Req_VRE"]*sum(EP[:vCAP][y]*inputs["pP_Max"][y,t] for y=dfGen[(dfGen[!,:DISP].==1) .| (dfGen[!,:NDISP].==1),:][!,:R_ID]) +
					maximum(inputs["dfGen"][!,:Cap_size]) )
			end
		end
		# Operating reserve down requirements as a percentage of load and scheduled variable renewable energy production in each hour
		@expression(EP, eRsvDnReq[t=1:T], inputs["pRsv_Dn_Req_Load"]*sum(inputs["pD"][t,z] for z=1:Z) +
			inputs["pRsv_Dn_Req_VRE"]*sum(EP[:vCAP][y]*inputs["pP_Max"][y,t] for y=dfGen[(dfGen[!,:DISP].==1) .| (dfGen[!,:NDISP].==1),:][!,:R_ID]))
	end ## End reserve expressions block

	## Objective Function Expressions ##

	## Integer Unit Commitment configuration for expressions
	if (setup["UCommit"]>=1 && setup["Reserves"]==1)
		# Penalty for unmet operating reserves
		@expression(EP, eCRsvPen[t=1:T], (inputs["omega"][t]*inputs["pC_Rsv_Penalty"]*(vUNMET_RSV_UP[t]+vUNMET_RSV_DN[t])))
			dExpressions["eCRsvPen"] = eCRsvPen
	end #END unit commitment configuration

	# If reserves used, include cost of unmet reserves
	if (setup["Reserves"]==1)
		@expression(EP, eTotalCRsvPen, sum(eCRsvPen[t] for t=1:T) +
			sum(dfGen[!,:Reg_Cost][y]*vRSV_UP[y,t] + dfGen[!,:Reg_Cost][y]*vRSV_DN[y,t] for y=1:G, t=1:T) +
			sum(dfGen[!,:Rsv_Cost][y]*vREG_UP[y,t] + dfGen[!,:Rsv_Cost][y]*vREG_DN[y,t] for y=1:G, t=1:T) )
	else
		@expression(EP, eTotalCRsvPen, 0)
	end
	dObjective["eTotalCRsvPen"] = eTotalCRsvPen

	## End Objective Function Expressions ##

	## Power Balance Expressions ##
	## End Power Balance Expressions ##

	### End Expressions ###

	### Constratints ###
	### End Constraints ###

	### Expressions & Constraints ###

	## Total system reserve constraints (if active)
	if (setup["Reserves"]==1)
		# Regulation requirements as a percentage of load and scheduled variable renewable energy production in each hour
		# Reg up and down requirements are symmetric
		@constraint(EP, cRegUp[t=1:T], sum(vREG_UP[y,t] for y=1:G) >= eRegReq[t])
		@constraint(EP, cRegDn[t=1:T], sum(vREG_DN[y,t] for y=1:G) >= eRegReq[t])

		# Operating reserve up requirements as a percentage of load and scheduled variable renewable energy production in each hour
		# and the largest single contingency
		if (inputs["pDynamic_Contingency"]==2)
			# Largest contingency defined for each hour as max of either largest committed generator or largest transmission line
			@constraint(EP, cContingency[y=dfGen[(dfGen[!,:Commit].==1),:][!,:R_ID], t=1:T], vLARGEST_CONTINGENCY[t] >=
				inputs["dfGen"][!,:Cap_size][y]*vCONTINGENCY_AUX[y,t] )

			# Ensure vCONTINGENCY_AUX = 0 if vCOMMIT = 0
			@constraint(EP, cContAux[y=dfGen[(dfGen[!,:Commit].==1),:][!,:R_ID], t=1:T], vCONTINGENCY_AUX[y,t] <= EP[:vCOMMIT][y,t])

			# Ensure vCONTINGENCY_AUX = 1 if vCOMMIT > 0
			@constraint(EP, cContAux2[y=dfGen[(dfGen[!,:Commit].==1),:][!,:R_ID], t=1:T], EP[:vCOMMIT][y, t] <= inputs["pContingency_BigM"][y]*vCONTINGENCY_AUX[y,t])
		elseif (inputs["pDynamic_Contingency"]==1)
			# Largest contingency defined as max of either largest installed generator or largest transmission line
			@constraint(EP, cContingency[y=dfGen[(dfGen[!,:Commit].==1),:][!,:R_ID]], vLARGEST_CONTINGENCY[1] >=
				inputs["dfGen"][!,:Cap_size][y]*vCONTINGENCY_AUX[y] )

			# Ensure vCONTINGENCY_AUX = 0 if vCAP = 0
			@constraint(EP, cContAux1[y=dfGen[(dfGen[!,:Commit].==1),:][!,:R_ID]], vCONTINGENCY_AUX[y] <= EP[:vCAP][y])

			# Ensure vCONTINGENCY_AUX = 1 if vCAP > 0
			@constraint(EP, cContAux2[y=dfGen[(dfGen[!,:Commit].==1),:][!,:R_ID]], EP[:vCAP][y] <= inputs["pContingency_BigM"][y]*vCONTINGENCY_AUX[y])
		end
		@constraint(EP, cRsvUp[t=1:T], sum(vRSV_UP[y,t] for y=1:G) + vUNMET_RSV_UP[t] >= eRsvUpReq[t])

		@constraint(EP, cRsvDn[t=1:T], sum(vRSV_DN[y,t] for y=1:G) + vUNMET_RSV_DN[t] >= eRsvDnReq[t] )
	end ## End reserve constraints block

	dModuleArgs["dExpressions"] = dExpressions
	dModuleArgs["dObjective"] = dObjective
	dModuleArgs["dPowerBalance"] = dPowerBalance

	return EP, dModuleArgs
end

function transmission(EP::Model, dModuleArgs::Dict)

	println("Transmission Module")

	inputs = dModuleArgs["inputs"]
	setup = dModuleArgs["setup"]
	dExpressions = dModuleArgs["dExpressions"]
	dObjective = dModuleArgs["dObjective"]
	dPowerBalance = dModuleArgs["dPowerBalance"]

	dfGen = inputs["dfGen"]

	T = inputs["T"]     # Number of time steps (hours)
	Z = inputs["Z"]     # Number of zones
	L = inputs["L"]     # Number of transmission lines
	SEG = inputs["SEG"] # Number of load curtailment segments

	## sets and indices for transmission and distribution losses and expansion
	TRANS_LOSS_SEGS = setup["Trans_Loss_Segments"] # Number of segments used in piecewise linear approximations quadratic loss functions
	DISTR_LOSS_SEGS = setup["Distr_Loss_Segments"] # Number of segments used in piecewise linear approximations quadratic loss functions
	DISTR_MARGIN_LV_SEGS = setup["Distr_Margin_LV_Segments"] # Number of segments used in piecewise linear approximations quadratic loss functions
	LOSS_LINES = findall(inputs["pTrans_Loss_Coef"].!=0) # Lines for which loss coefficients apply (are non-zero);
	NOLOSS_LINES = findall(inputs["pTrans_Loss_Coef"].==0) # Other lines will not have loss-related constraints and losses are set automatically to zero
	DISTR_ZONES = findall(inputs["pDistrZones"].==1)	# Set of distribution network zones
	if sum(inputs["pDistrZones"]) > 0
		LOSS_ZONES = findall(inputs["pDistr_Loss_LV_Net_Quad"].!=0) # Zones for which within-zone loss coefficients apply (are non-zero); other zones will not have loss-related constraints
		NOLOSS_ZONES = findall(inputs["pDistr_Loss_LV_Net_Quad"].==0) # Other zones will not have loss-related constraints and losses are set automatically to zero
	end
	if setup["NetworkExpansion"] == 1
		# Network lines and zones that are expandable have non-negative maximum reinforcement inputs
		EXPANSION_LINES = findall(inputs["pMax_Line_Reinforcement"].>=0)
		NO_EXPANSION_LINES = findall(inputs["pMax_Line_Reinforcement"].<0)
		if sum(inputs["pDistrZones"]) > 0
			EXPANSION_INJECT = findall(inputs["pMax_Inject_Reinforcement"].>=0)
			NO_EXPANSION_INJECT = intersect(findall(inputs["pMax_Inject_Reinforcement"].<0), DISTR_ZONES)
			EXPANSION_WITHDRAW =  intersect(findall(inputs["pMax_Withdraw_Reinforcement"].>=0), DISTR_ZONES)
			NO_EXPANSION_WITHDRAW =  intersect(findall(inputs["pMax_Withdraw_Reinforcement"].<0), DISTR_ZONES)
		end
	end

	### Variables ###

	# Power flow on each transmission line "l" at hour "t"
	@variable(EP, vFLOW[l=1:L,t=1:T]);

	if (setup["NetworkExpansion"]==1)
		# Transmission network capacity reinforcements per line
		@variable(EP, vNEW_TRANS_CAP[l=1:L] >= 0)
		# Distribution (within-zone) network capacity reinforcements per zone
		@variable(EP, vNEW_INJECT_CAP[z=DISTR_ZONES] >= 0)
		@variable(EP, vNEW_WITHDRAW_CAP[z=DISTR_ZONES] >= 0)
		if (DISTR_LOSS_SEGS > 0)
			# Distribution network `margin' gained via net demand reductions
			# (currently no option for generation curtailment to contribute to injection network margin
			# except in a 1:1 manner as curtailable DG during coincident peak injection periods reduces the
			# peak injection variable
			@variable(EP, vWITHDRAW_MARGIN_GAINED[z=DISTR_ZONES, t=inputs["pPeak_Inject_Hrs"]] >= 0)
			@variable(EP, vMV_NET_REDUCTION_AUX[z=DISTR_ZONES, t=inputs["pPeak_Inject_Hrs"]] >= 0) # Auxiliary variable for MV net demand reduction
			@variable(EP, vDISTR_MARGIN_LV_AUX[z=DISTR_ZONES,t=inputs["pPeak_Inject_Hrs"], s=1:DISTR_MARGIN_LV_SEGS] >= 0)
		end
	end

	if (TRANS_LOSS_SEGS == 0)
		# zone or bus angles at hour "t" in rad
		# @variable(EP, vTHETHA[z=1:Z,t=1:T]);
  	elseif (TRANS_LOSS_SEGS==1)  #loss is a constant times absolute value of power flow
		# Positive and negative flow variables
		@variable(EP, vTAUX_NEG[l=LOSS_LINES,t=1:T] >= 0)
		@variable(EP, vTAUX_POS[l=LOSS_LINES,t=1:T] >= 0)

		if (setup["UCommit"]==1)
			# Single binary variable to ensure positive or negative flows only
			@variable(EP, vTAUX_POS_ON[l=LOSS_LINES,t=1:T],Bin)
		end
	else
		# Auxiliary variables for linear piecewise interpolation of quadratic losses
		@variable(EP, vTAUX_NEG[l=LOSS_LINES,s=0:TRANS_LOSS_SEGS,t=1:T] >= 0)
		@variable(EP, vTAUX_POS[l=LOSS_LINES,s=0:TRANS_LOSS_SEGS,t=1:T] >= 0)
		if (setup["UCommit"]==1)
			# Binary auxilary variables for each segment >1 to ensure segments fill in order
			@variable(EP, vTAUX_POS_ON[l=LOSS_LINES,s=1:TRANS_LOSS_SEGS,t=1:T], Bin)
			@variable(EP, vTAUX_NEG_ON[l=LOSS_LINES,s=1:TRANS_LOSS_SEGS,t=1:T], Bin)
		end
    end

	# Transmission losses on each transmission line "l" at hour "t"
	@variable(EP, vTLOSS[l=1:L,t=1:T] >= 0)

	if (DISTR_LOSS_SEGS > 0 && sum(inputs["pDistrZones"]) > 0)
		# Electrical losses for flows within zone "i" in hour "t"
		@variable(EP, vDLOSS[z=1:Z,t=1:T] >= 0)
		# Auxiliary variables for linear piecewise interpolation of quadratic losses
		@variable(EP, vWAUX_NEG[z=LOSS_ZONES,s=0:DISTR_LOSS_SEGS,t=1:T] >= 0)
		@variable(EP, vWAUX_POS[z=LOSS_ZONES,s=0:DISTR_LOSS_SEGS,t=1:T] >= 0)
		if (setup["UCommit"]==1)
			# Binary auxilary variables for each segment >1 to ensure segments fill in order
			@variable(EP, vWAUX_POS_ON[z=LOSS_ZONES,s=1:DISTR_LOSS_SEGS,t=1:T], Bin)
			@variable(EP, vWAUX_NEG_ON[z=LOSS_ZONES,s=1:DISTR_LOSS_SEGS,t=1:T], Bin)
		end
	else # Either no distribution zones at all or not modeling distribution losses, so set vDLOSS to zero
		@variable(EP, vDLOSS[z=1:Z,t=1:T] == 0)
	end

	### End Variables ###

	### Expressions ###

	## Transmission power flow and loss related expressions:
	# Net power flow outgoing from zone "z" at hour "t" in MW
    @expression(EP, eNet_Export_Flows[z=1:Z,t=1:T], sum(inputs["pNet_Map"][l,z] * vFLOW[l,t] for l=1:L))
	dExpressions["eNet_Export_Flows"] = eNet_Export_Flows

	# Losses from power flows into or out of zone "z" in MW
    @expression(EP, eLosses_By_Zone[z=1:Z,t=1:T], sum(abs(inputs["pNet_Map"][l,z]) * vTLOSS[l,t] for l=1:L))
	dExpressions["eLosses_By_Zone"] = eLosses_By_Zone

	## Distribution network zone related expressions
	if (setup["NetworkExpansion"]==1)
		# Total power withdrawal and injection expressions used for distribution loss calculations
		# Total injection expressions also currently used for determining injection related distribution network capacity needed
		# Total power withdrawn from grid in node z during time t in MW (equal to sum of demand + storage charging + deferred demand satisfied - demand deffered - nonserved energy)
		# Withdrawal in LV (low voltage level)
		@expression(EP, eWithdrawLV[z=DISTR_ZONES, t=1:T],  inputs["pD"][t,z]*inputs["pDistr_Share_in_LV"][z]
															+ sum(EP[:vCHARGE][y,t] for y=dfGen[(dfGen[!,:STOR].>=1) .& (dfGen[!,:zone].==z) .& (dfGen[!,:voltage_level].=="LV"),:][!,:R_ID])
															+ sum(EP[:vP][y,t] for y=dfGen[(dfGen[!,:DR].==1) .& (dfGen[!,:zone].==z) .& (dfGen[!,:voltage_level].=="LV"),:][!,:R_ID])
															+ sum(EP[:vCHARGE][y,t] for y=dfGen[(dfGen[!,:HEAT].==1) .& (dfGen[!,:zone].==z) .& (dfGen[!,:voltage_level].=="LV"),:][!,:R_ID])
															- sum(EP[:vCHARGE][y,t] for y=dfGen[(dfGen[!,:DR].==1) .& (dfGen[!,:zone].==z) .& (dfGen[!,:voltage_level].=="LV"),:][!,:R_ID])
															- sum(EP[:vNSE][s,t,z] for s=1:SEG)
															)
															# Note: assumes any price responsive demand curtailment is in LV as this has greatest impact on losses per MW and there is not currently any differentiation
                                                            # between curtailable load segments in each voltage level.
        dExpressions["eWithdrawLV"] = eWithdrawLV

		# Withdrawal in MV (medium voltage level)
		@expression(EP, eWithdrawMV[z=DISTR_ZONES, t=1:T],  inputs["pD"][t,z]*inputs["pDistr_Share_in_MV"][z]
															+ sum(EP[:vCHARGE][y,t] for y=dfGen[(dfGen[!,:STOR].>=1) .& (dfGen[!,:zone].==z) .& (dfGen[!,:voltage_level].=="MV"),:][!,:R_ID])
															+ sum(EP[:vP][y,t] for y=dfGen[(dfGen[!,:DR].==1) .& (dfGen[!,:zone].==z) .& (dfGen[!,:voltage_level].=="MV"),:][!,:R_ID])
															+ sum(EP[:vCHARGE][y,t] for y=dfGen[(dfGen[!,:HEAT].==1) .& (dfGen[!,:zone].==z) .& (dfGen[!,:voltage_level].=="MV"),:][!,:R_ID])
															- sum(EP[:vCHARGE][y,t] for y=dfGen[(dfGen[!,:DR].==1) .& (dfGen[!,:zone].==z) .& (dfGen[!,:voltage_level].=="MV"),:][!,:R_ID])
															)
															# Note: assumes any price responsive demand curtailment is in LV as this has greatest impact on losses per MW and there is not currently any differentiation
															# between curtailable load segments in each voltage level.

        dExpressions["eWithdrawMV"] = eWithdrawMV
		# Total power injected into grid in node z during time t in MW (equal to sum of power generated by generators + discharged by storage)
		# Injection in LV (low voltage level)
		@expression(EP, eInjectLV[z=DISTR_ZONES, t=1:T], sum(EP[:vP][y,t] for y=dfGen[(dfGen[!,:THERM].==1) .& (dfGen[!,:zone].==z) .& (dfGen[!,:voltage_level].=="LV"),:][!,:R_ID])
															+ sum(EP[:vP][y,t] for y=dfGen[(dfGen[!,:DISP].==1) .& (dfGen[!,:zone].==z) .& (dfGen[!,:voltage_level].=="LV"),:][!,:R_ID])
															+ sum(EP[:vP][y,t] for y=dfGen[(dfGen[!,:NDISP].==1) .&  (dfGen[!,:zone].==z) .& (dfGen[!,:voltage_level].=="LV"),:][!,:R_ID])
															+ sum(EP[:vP][y,t] for y=dfGen[(dfGen[!,:STOR].>=1) .&  (dfGen[!,:zone].==z) .& (dfGen[!,:voltage_level].=="LV"),:][!,:R_ID])
															+ sum(EP[:vP][y,t] for y=dfGen[(dfGen[!,:NACC].==1) .&  (dfGen[!,:zone].==z) .& (dfGen[!,:voltage_level].=="LV"),:][!,:R_ID]) )

        dExpressions["eInjectLV"] = eInjectLV
        # Injection in MV (medium voltage level)
		@expression(EP, eInjectMV[z=DISTR_ZONES, t=1:T], sum(EP[:vP][y,t] for y=dfGen[(dfGen[!,:THERM].==1) .&  (dfGen[!,:zone].==z) .&  (dfGen[!,:voltage_level].=="MV"),:][!,:R_ID])
															+ sum(EP[:vP][y,t] for y=dfGen[(dfGen[!,:DISP].==1) .&  (dfGen[!,:zone].==z) .&  (dfGen[!,:voltage_level].=="MV"),:][!,:R_ID])
															+ sum(EP[:vP][y,t] for y=dfGen[(dfGen[!,:NDISP].==1) .&  (dfGen[!,:zone].==z) .&  (dfGen[!,:voltage_level].=="MV"),:][!,:R_ID])
															+ sum(EP[:vP][y,t] for y=dfGen[(dfGen[!,:STOR].>=1) .&  (dfGen[!,:zone].==z) .&  (dfGen[!,:voltage_level].=="MV"),:][!,:R_ID])
															+ sum(EP[:vP][y,t] for y=dfGen[(dfGen[!,:NACC].==1) .&  (dfGen[!,:zone].==z) .&  (dfGen[!,:voltage_level].=="MV"),:][!,:R_ID]) )
        dExpressions["eInjectMV"] = eInjectMV

		if sum(inputs["pDistrZones"]) > 0
			# Peak withdrawal expression for determining total withdrawal capacity needed either via network reinforcement or network margin gained from net demand reductions
			# Peak withdrawal includes demand plus all deferred demand satisfied and all storage charging during that period.
			# Note that demand deferred and price responsive demand curtailment appear in the eNetReduction expressions, which is why this differs from the eWithdraw expression above
			# used for losses calculations. (Demand reduction has a greater than 1:1 effect on network margins needed)
			@expression(EP, ePeakWithdraw[z=DISTR_ZONES, t=inputs["pPeak_Withdraw_Hrs"]],
																inputs["pD"][t,z]
																+ sum(EP[:vCHARGE][y,t] for y=dfGen[(dfGen[!,:STOR].>=1) .&  (dfGen[!,:zone].==z),:][!,:R_ID])
																+ sum(EP[:vP][y,t] for y=dfGen[(dfGen[!,:DR].==1) .&  (dfGen[!,:zone].==z),:][!,:R_ID])
																+ sum(EP[:vCHARGE][y,t] for y=dfGen[(dfGen[!,:HEAT].==1) .&  (dfGen[!,:zone].==z),:][!,:R_ID])
																)
            dExpressions["ePeakWithdraw"] = ePeakWithdraw
			# Net demand reduction expressions for each distribution voltage level
			# Net reductions are equal to all demand deferral and price responsive demand (NSE) plus all injections from DERs in the voltage level and zone.
			# Note that all price responsive demand is assumed to contribute to LV net reductions at this point, as there is no way to disaggregate
			# demand response by voltage level at the moment, except via the fixed shares.
			@expression(EP, eNetReductionMV[z=DISTR_ZONES, t=inputs["pPeak_Withdraw_Hrs"]],
																sum(EP[:vCHARGE][y,t] for y=dfGen[(dfGen[!,:DR].==1) .&  (dfGen[!,:zone].==z) .&  (dfGen[!,:voltage_level].=="MV"),:][!,:R_ID])
																+ eInjectMV[z,t])
			@expression(EP, eNetReductionLV[z=DISTR_ZONES, t=inputs["pPeak_Withdraw_Hrs"]],
																sum(EP[:vCHARGE][y,t] for y=dfGen[(dfGen[!,:DR].==1) .&  (dfGen[!,:zone].==z) .&  (dfGen[!,:voltage_level].=="LV"),:][!,:R_ID])
																+ sum(EP[:vNSE][s,t,z] for s=1:SEG)
                                                                + eInjectLV[z,t])
            dExpressions["eNetReductionMV"] = eNetReductionMV
        	dExpressions["eNetReductionLV"] = eNetReductionLV
				# These two expressions break out the total net demand reduction across the zone (in either MV or LV) that is required by network margin gained in each voltage level
				# Functions based on Jenkins (2018) PhD thesis work
			@expression(EP, eReductionRequiredByMV[z=DISTR_ZONES, t=inputs["pPeak_Withdraw_Hrs"]],
															vWITHDRAW_MARGIN_GAINED[z,t]*inputs["pDistr_Margin_MV_Linear"][z]*inputs["pDistr_Share_in_MV"][z])

			@expression(EP, eReductionRequiredByLV[z=DISTR_ZONES, t=inputs["pPeak_Withdraw_Hrs"]],
																vWITHDRAW_MARGIN_GAINED[z,t]*inputs["pDistr_Margin_LV_Linear"][z]*inputs["pDistr_Share_in_LV"][z]
																+ sum((2*s-1)*vDISTR_MARGIN_LV_AUX[z,t,s] for
																	s=1:DISTR_MARGIN_LV_SEGS)*inputs["pDistr_Margin_LV_Quad"][z])
            dExpressions["eReductionRequiredByMV"] = eReductionRequiredByMV
            dExpressions["eReductionRequiredByLV"] = eReductionRequiredByLV
            end # End block of distribution related expressions
	end

	## Objective Function Expressions ##

	if (setup["NetworkExpansion"]==1)
		if sum(inputs["pDistrZones"]) > 0
			@expression(EP, eTotalCNetworkExp, sum( vNEW_TRANS_CAP[l]*inputs["pC_Line_Reinforcement"][l] for l=EXPANSION_LINES) +
											sum( vNEW_INJECT_CAP[z]*inputs["pDistr_Share_in_MV"][z]*inputs["pC_MV_Reinforcement"][z] +
												 vNEW_INJECT_CAP[z]*inputs["pDistr_Share_in_LV"][z]*inputs["pC_LV_Reinforcement"][z] for z=EXPANSION_INJECT) +
											sum( vNEW_WITHDRAW_CAP[z]*inputs["pDistr_Share_in_MV"][z]*inputs["pC_MV_Reinforcement"][z] +
												 vNEW_WITHDRAW_CAP[z]*inputs["pDistr_Share_in_LV"][z]*inputs["pC_LV_Reinforcement"][z] for z=EXPANSION_WITHDRAW)
												 )
		else
			@expression(EP, eTotalCNetworkExp, sum( vNEW_TRANS_CAP[l]*inputs["pC_Line_Reinforcement"][l] for l=EXPANSION_LINES) )
		end
	else
		@expression(EP, eTotalCNetworkExp, 0)
    end

    dObjective["eTotalCNetworkExp"] = eTotalCNetworkExp

	## End Objective Function Expressions ##

	## Power Balance Expressions ##

	@expression(EP, ePowerBalanceDLoss[t=1:T, z=1:Z],
		-vDLOSS[z,t])
	@expression(EP, ePowerBalanceNetExportFlows[t=1:T, z=1:Z],
		-eNet_Export_Flows[z,t])
	@expression(EP, ePowerBalanceLossesByZone[t=1:T, z=1:Z],
		-(1/2)*eLosses_By_Zone[z,t])

	dPowerBalance["ePowerBalance"] = dPowerBalance["ePowerBalance"] + ePowerBalanceLossesByZone
	dPowerBalance["ePowerBalance"] = dPowerBalance["ePowerBalance"] + ePowerBalanceNetExportFlows
	dPowerBalance["ePowerBalance"] = dPowerBalance["ePowerBalance"] + ePowerBalanceDLoss

	## End Power Balance Expressions ##

	### End Expressions ###

	### Constraints ###

  	## Power flow and transmission (between zone) loss related constraints
	# setting transmission variables to zero if just one zone
	if Z == 1
		for t in 1:T
			for l in 1:L
				setupperbound(vFLOW[l,t], 0.0)
				setlowerbound(vFLOW[l,t], 0.0)
				setupperbound(vTLOSS[l,t], 0.0)
				setlowerbound(vTLOSS[l,t], 0.0)
 				# setupperbound(vNEW_TRANS_CAP[l], 0.0)
          		# setlowerbound(vNEW_TRANS_CAP[l], 0.0)
			end
			# setupperbound(vTHETHA[Z,t], 0.0)
			# setlowerbound(vTHETHA[Z,t], 0.0)
		end
		# Multizonal transmission and distribution zone/line related constraints if more than one zone
	else

		# Maximum power flows, power flow on each transmission line cannot exceed maximum capacity of the line at any hour "t"
		# If network expansion is used, power injections and withdrawals by zone also constrained
		if (setup["NetworkExpansion"]==1)
			# Transmission network related power flow and capacity constraints
			# Allow expansion of transmission capacity for lines eligible for reinforcement
			@constraint(EP, cMaxFlow_out[l=1:L, t=1:T], vFLOW[l,t] <= inputs["pTrans_Max"][l] + vNEW_TRANS_CAP[l])
			@constraint(EP, cMaxFlow_in[l=1:L, t=1:T], vFLOW[l,t] >= -inputs["pTrans_Max"][l] - vNEW_TRANS_CAP[l])
			# Constrain maximum line capacity reinforcement for lines eligible for expansion
			@constraint(EP, cMaxLineReinforcement[l=EXPANSION_LINES], vNEW_TRANS_CAP[l] <= inputs["pMax_Line_Reinforcement"][l])
			# Constrain lines ineligible for expansion (solver will drop variable from model)
			for l in NO_EXPANSION_LINES
				setupperbound(vNEW_TRANS_CAP[l], 0.0)
				setlowerbound(vNEW_TRANS_CAP[l], 0.0)
			end

			if sum(inputs["pDistrZones"]) > 0
				# Distribution network zone related withdrawal/injection and capacity constraints
				# Note constraints applied only to hours in the relevant peak hours set to limit # of constraints added
				# Peak net injection cannot exceed max network capacity
				# (note that any withdrawals during these periods offset injections one for one at moment, pending experiments to confirm if injection
				# margins gained by net injection reductions mirror withdrawal margins gained by net demand reductions)
				@constraint(EP, cMaxInject[z=DISTR_ZONES, t=inputs["pPeak_Inject_Hrs"]], ( eInjectMV[z,t] + eInjectLV[z,t] ) <=
																			inputs["pMax_Inject"][z] + vNEW_INJECT_CAP[z] + eWithdrawMV[z,t] + eWithdrawMV[z,t])
				# Peak withdrawals cannot exceed max network capacity plus network margin gained via targeted net demand reductions as per Jenkins (2018) PhD thesis
				@constraint(EP, cMaxWithdraw[z=DISTR_ZONES, t=inputs["pPeak_Withdraw_Hrs"]],  ePeakWithdraw[z,t] <=
																			inputs["pMax_Withdraw"][z] + vNEW_WITHDRAW_CAP[z] + vWITHDRAW_MARGIN_GAINED[z,t])
				# Net demand reductions must exceed network margin gain needed in each time period as per Jenkins (2018) PhD thesis
				@constraint(EP, cNetReductions[z=DISTR_ZONES, t=inputs["pPeak_Withdraw_Hrs"]],
													(vMV_NET_REDUCTION_AUX[z,t]*inputs["pDistr_Margin_MV_Discount"][z] + eNetReductionLV[z,t]) >=
													eReductionRequiredByMV[z,t] + eReductionRequiredByLV[z,t])

				# Constraints on MV net reduction auxiliary variable (must be <= max contribution from MV or net demand reductions in MV, whichever is lower)
				@constraint(EP, cMVMarginMax1[z=DISTR_ZONES, t=inputs["pPeak_Withdraw_Hrs"]], vMV_NET_REDUCTION_AUX[z,t] <= eNetReductionMV[z,t])
				@constraint(EP, cMVMarginMax2[z=DISTR_ZONES, t=inputs["pPeak_Withdraw_Hrs"]], vMV_NET_REDUCTION_AUX[z,t] <=
																				inputs["pDistr_Margin_MV_Max"][z]*eReductionRequiredByLV[z,t] +
																				eReductionRequiredByMV[z,t])
				# Constraints on auxiliary variable used to approximate quadratic term in net reduction required by margin needed in LV
					# Segment length constraint
				@constraint(EP, cMaxDistrSegLength[z=DISTR_ZONES, t=inputs["pPeak_Withdraw_Hrs"], s=1:DISTR_MARGIN_LV_SEGS],
																				vDISTR_MARGIN_LV_AUX[z,t,s] <= inputs["pDistr_Margin_LV_Segment_Length"][z] )
					# Sum of each auxiliary variable must equal total LV margin gained (used to approximate quadratic term in eNetReductionLV)
				@constraint(EP, cSumDistrSegs[z=DISTR_ZONES, t=inputs["pPeak_Withdraw_Hrs"]],
																				sum(vDISTR_MARGIN_LV_AUX[z,t,s] for s=1:DISTR_MARGIN_LV_SEGS) ==
																				vWITHDRAW_MARGIN_GAINED[z,t]*inputs["pDistr_Share_in_LV"][z] )

				# Constrain maximum distribution zone capacity reinforcement for zones eligible for expansion
				@constraint(EP, cMaxInjectReinforcement[z=EXPANSION_INJECT], vNEW_INJECT_CAP[z] <= inputs["pMax_Inject_Reinforcement"][z])
				@constraint(EP, cMaxWithdrawReinforcement[z=EXPANSION_WITHDRAW], vNEW_WITHDRAW_CAP[z] <= inputs["pMax_Withdraw_Reinforcement"][z])
				# Constrain distribution zones ineligible for expansion (solver will drop variable from model)
				for z in NO_EXPANSION_INJECT
					setupperbound(vNEW_INJECT_CAP[z], 0.0)
					setlowerbound(vNEW_INJECT_CAP[z], 0.0)
				end
				for z in NO_EXPANSION_WITHDRAW
					setupperbound(vNEW_WITHDRAW_CAP[z], 0.0)
					setlowerbound(vNEW_WITHDRAW_CAP[z], 0.0)
				end
			end # End block of distribution zone related constraints
		else
			# Capacity of lines is fixed and zones do not have constraints on injection and withdrawls
			@constraint(EP, cMaxFlow_out[l=1:L, t=1:T], vFLOW[l,t] <= inputs["pTrans_Max"][l])
			@constraint(EP, cMaxFlow_in[l=1:L, t=1:T], vFLOW[l,t] >= -inputs["pTrans_Max"][l])
		end
		#END network expansion contraints

		# If no loss segments, then set vTLOSS variables to 0 to eliminate them from model
		if TRANS_LOSS_SEGS == 0
			# DC power flow related constraints; currently not implemented
			# Maximum zone/bus angle at any hour "t"
			#@constraint(EP, cThetha_up[l=1:L, t=1:T], sum(inputs["pNet_Map"][l,z]*vTHETHA[z,t],z=1:Z) <= inputs["pTheta_Max"][l])
			#@constraint(EP, cThetha_down[l=1:L, t=1:T], sum(inputs["pNet_Map"][l,z]*vTHETHA[z,t],z=1:Z) >= -inputs["pTheta_Max"][l])
			# Maximum power flow as function of thetha at any hour "t"
			#@constraint(EP, cMaxFlow_thetha[l=1:L, t=1:T],  inputs["line_X"][1]*vFLOW[l,t] == (10^3*inputs["kV"][1])^2*sum(inputs["pNet_Map"][l,z]*vTHETHA[z,t],z=1:Z))
			# reference zone
			#@constraint(EP, cReferenceZone[z=1:Z,t=1:T],vTHETHA[z=1,t] == 0)
			for t in 1:T
				for l in 1:L
					setupperbound(vTLOSS[l,t], 0.0)
					setlowerbound(vTLOSS[l,t], 0.0)
				end
				#=
				for z in 1:Z
					setupperbound(vTHETHA[z,t], 0.0)
					setlowerbound(vTHETHA[z,t], 0.0)
				end
				=#
			end
		end

		# Transmission loss related constraints - linear losses as a function of absolute value
		if TRANS_LOSS_SEGS == 1
		    # Losses are alpha times abosolute values
		    @constraint(EP, cTLoss[l=LOSS_LINES, t=1:T], vTLOSS[l,t] ==
		                inputs["pPercent_Loss"][l]*(vTAUX_POS[l,t] + vTAUX_NEG[l,t]))

		    # power flow is sum of positive and negative components
		    @constraint(EP, cTAuxSum[l=LOSS_LINES, t=1:T], vTAUX_POS[l,t] - vTAUX_NEG[l,t]  == vFLOW[l,t])

			if (setup["UCommit"]==0 || setup["UCommit"]==2)
		    # positive and negative power flows are bounded
		    	@constraint(EP, cTAuxPosUB[l=LOSS_LINES, t=1:T], vTAUX_POS[l,t] <= inputs["pTrans_Max_Possible"][l])
		    	@constraint(EP, cTAuxNegUB[l=LOSS_LINES, t=1:T], vTAUX_NEG[l,t] <= inputs["pTrans_Max_Possible"][l])

		    # when modeling unit commitment, we need additional binary variables to avoid phantom losses
			else # setup["UCommit"]==1
		        @constraint(EP, cTAuxPosUB[l=LOSS_LINES, t=1:T], vTAUX_POS[l,t] <=
		                            inputs["pTrans_Max_Possible"][l]*vTAUX_POS_ON[l,t])
		        # Either negative or positive flows are activated, not both
		        @constraint(EP, cTAuxNegUB[l=LOSS_LINES, t=1:T], vTAUX_NEG[l,t] <=
		                            inputs["pTrans_Max_Possible"][l]*(1-vTAUX_POS_ON[l,t]))
		    end

		end # End if(TRANS_LOSS_SEGS == 1) block

		# When number of segments is greater than 1
		if (TRANS_LOSS_SEGS > 1)
				## Between zone transmission loss constraints
					# Losses are expressed as a piecewise approximation of a quadratic function of power flows across each line
					# Eq 1: Total losses are function of loss coefficient times the sum of auxilary segment variables across all segments of piecewise approximation
					# (Includes both positive domain and negative domain segments)
			@constraint(EP, cTLoss[l=LOSS_LINES, t=1:T], vTLOSS[l,t] ==
								(inputs["pTrans_Loss_Coef"][l]*sum((2*s-1)*(inputs["pTrans_Max_Possible"][l]/TRANS_LOSS_SEGS)*vTAUX_POS[l,s,t] for s=1:TRANS_LOSS_SEGS)) +
								(inputs["pTrans_Loss_Coef"][l]*sum((2*s-1)*(inputs["pTrans_Max_Possible"][l]/TRANS_LOSS_SEGS)*vTAUX_NEG[l,s,t] for s=1:TRANS_LOSS_SEGS)) )
				# Eq 2: Sum of auxilary segment variables (s >= 1) minus the "zero" segment (which allows values to go negative)
				# from both positive and negative domains must total the actual power flow across the line
			@constraint(EP, cTAuxSumPos[l=LOSS_LINES, t=1:T], sum(vTAUX_POS[l,s,t] for s=1:TRANS_LOSS_SEGS) - vTAUX_POS[l,0,t]  == vFLOW[l,t])
			@constraint(EP, cTAuxSumNeg[l=LOSS_LINES, t=1:T], sum(vTAUX_NEG[l,s,t] for s=1:TRANS_LOSS_SEGS) - vTAUX_NEG[l,0,t]  == -vFLOW[l,t])
			if (setup["UCommit"]==0 || setup["UCommit"]==2)
				# Eq 3: Each auxilary segment variables (s >= 1) must be less than the maximum power flow in the zone / number of segments
				@constraint(EP, cTAuxMaxPos[l=LOSS_LINES, s=1:TRANS_LOSS_SEGS, t=1:T], vTAUX_POS[l,s,t] <= (inputs["pTrans_Max_Possible"][l]/TRANS_LOSS_SEGS))
				@constraint(EP, cTAuxMaxNeg[l=LOSS_LINES, s=1:TRANS_LOSS_SEGS, t=1:T], vTAUX_NEG[l,s,t] <= (inputs["pTrans_Max_Possible"][l]/TRANS_LOSS_SEGS))
			else # Constraints that can be ommitted if problem is convex (i.e. if not using MILP unit commitment constraints)
					# Eqs 3-4: Ensure that auxilary segment variables do not exceed maximum value per segment and that they
					# "fill" in order: i.e. one segment cannot be non-zero unless prior segment is at it's maximum value
					# (These constraints are necessary to prevents phantom losses in MILP problems)
				@constraint(EP, cTAuxOrderPos1[l=LOSS_LINES, s=1:TRANS_LOSS_SEGS, t=1:T], vTAUX_POS[l,s,t] <=  (inputs["pTrans_Max_Possible"][l]/TRANS_LOSS_SEGS)*vTAUX_POS_ON[l,s,t])
				@constraint(EP, cTAuxOrderNeg1[l=LOSS_LINES, s=1:TRANS_LOSS_SEGS, t=1:T], vTAUX_NEG[l,s,t] <=  (inputs["pTrans_Max_Possible"][l]/TRANS_LOSS_SEGS)*vTAUX_NEG_ON[l,s,t])
				@constraint(EP, cTAuxOrderPos2[l=LOSS_LINES, s=1:(TRANS_LOSS_SEGS-1), t=1:T], vTAUX_POS[l,s,t] >=  (inputs["pTrans_Max_Possible"][l]/TRANS_LOSS_SEGS)*vTAUX_POS_ON[l,s+1,t])
				@constraint(EP, cTAuxOrderNeg2[l=LOSS_LINES, s=1:(TRANS_LOSS_SEGS-1), t=1:T], vTAUX_NEG[l,s,t] >=  (inputs["pTrans_Max_Possible"][l]/TRANS_LOSS_SEGS)*vTAUX_NEG_ON[l,s+1,t])
				# Eq 5: Binary constraints to deal with absolute value of vFLOW.
						# If flow is positive, vTAUX_POS segment 0 must be zero; If flow is negative, vTAUX_POS segment 0 must be positive
						# (and takes on value of the full negative flow), forcing all vTAUX_POS other segments (s>=1) to be zero
				@constraint(EP, cTAuxSegmentZeroPos[l=LOSS_LINES, t=1:T], vTAUX_POS[l,0,t] <= inputs["pTrans_Max_Possible"][l]*(1-vTAUX_POS_ON[l,1,t]))
						# If flow is negative, vTAUX_NEG segment 0 must be zero; If flow is positive, vTAUX_NEG segment 0 must be positive
						# (and takes on value of the full positive flow), forcing all other vTAUX_NEG segments (s>=1) to be zero
				@constraint(EP, cTAuxSegmentZeroNeg[l=LOSS_LINES, t=1:T], vTAUX_NEG[l,0,t] <= inputs["pTrans_Max_Possible"][l]*(1-vTAUX_NEG_ON[l,1,t]))
			end
			# Lines without a loss coefficient have no losses in all time periods (and above constraints are not applied)
			for t in 1:T
				for l in NOLOSS_LINES
					setupperbound(vTLOSS[l,t], 0.0)
					setlowerbound(vTLOSS[l,t], 0.0)
				end
			end
		end # End if(TRANS_LOSS_SEGS > 0) block

		if sum(inputs["pDistrZones"]) > 0
		## Distribution zone related block:
			# If no distribution loss segments, then set vDLOSS variables to 0 to eliminate them from model
			if DISTR_LOSS_SEGS == 0
				for t in 1:T
					for z in 1:Z
						setupperbound(vDLOSS[z,t], 0.0)
						setlowerbound(vDLOSS[z,t], 0.0)
					end
				end
			end

			# If distribution loss segments, then set distribution loss related constraints
			if (DISTR_LOSS_SEGS > 0)
				## Distribution zone loss constraints:
				# Losses are expressed as a polynomial function with a piecewise approximation of a quadratic term (withdrawal_lv - injection_lv)^2
				# plus linear terms for (withdrawal_mv - injection_mv) and (Withdrawal_lv + Injection_lv) within each zone

				# Eq 1: Total losses are function of loss coefficient times the sum of auxilary segment variables across all segments of piecewise approximation of (net LV withdrawals)^2
				# (Includes both positive domain and negative domain segments) plus linear terms for net withdrawals in MV and total withdrawal + injection in LV
				@constraint(EP, cWLoss[z=LOSS_ZONES, t=1:T], vDLOSS[z,t] ==
						(inputs["pDistr_Loss_LV_Net_Quad"][z]*sum((2*s-1)*(inputs["pDistr_Max_Withdraw_Possible"][z]/DISTR_LOSS_SEGS)*vWAUX_POS[z,s,t] for s=1:DISTR_LOSS_SEGS)) +
						(inputs["pDistr_Loss_LV_Net_Quad"][z]*sum((2*s-1)*(inputs["pDistr_Max_Inject_Possible"][z]/DISTR_LOSS_SEGS)*vWAUX_NEG[z,s,t] for s=1:DISTR_LOSS_SEGS)) +
						inputs["pDistr_Loss_MV_Net_Linear"][z]*(eWithdrawMV[z,t]-eInjectMV[z,t]) +
						inputs["pDistr_Loss_LV_Total_Linear"][z]*(eWithdrawLV[z,t]+eInjectLV[z,t])
						)
				# Eq 2: Sum of auxilary segment variables (s >= 1) minus the "zero" segment (which allows values to go negative)
				# from both positive and negative domains must total the actual power flow across the line
				@constraint(EP, cWAuxSumPos[z=LOSS_ZONES, t=1:T], sum(vWAUX_POS[z,s,t] for s=1:DISTR_LOSS_SEGS) - vWAUX_POS[z,0,t]  == (eWithdrawLV[z,t] - eInjectLV[z,t]))
				@constraint(EP, cWAuxSumNeg[z=LOSS_ZONES, t=1:T], sum(vWAUX_NEG[z,s,t] for s=1:DISTR_LOSS_SEGS) - vWAUX_NEG[z,0,t]  == -(eWithdrawLV[z,t] - eInjectLV[z,t]))
				if (setup["UCommit"]==0 || setup["UCommit"]==2)
					# Eq 3: Each auxilary segment variables (s >= 1) must be less than the maximum power flow in the zone / number of segments
					@constraint(EP, cWAuxMaxPos[z=LOSS_ZONES, s=1:DISTR_LOSS_SEGS, t=1:T], vWAUX_POS[z,s,t] <= inputs["pDistr_Max_Withdraw_Possible"][z]/DISTR_LOSS_SEGS)
					@constraint(EP, cWAuxMaxNeg[z=LOSS_ZONES, s=1:DISTR_LOSS_SEGS, t=1:T], vWAUX_NEG[z,s,t] <= inputs["pDistr_Max_Inject_Possible"][z]/DISTR_LOSS_SEGS)
				else # Constraints that can be ommitted if problem is convex (i.e. if not using MILP unit commitment constraints)
					# Eqs 3-4: Ensure that auxilary segment variables do not exceed maximum value per segment and that they
						# "fill" in order: i.e. one segment cannot be non-zero unless prior segment is at it's maximum value
						# (These constraints are necessary to prevents phantom losses in MILP problems)
					@constraint(EP, cWAuxOrderPos1[z=LOSS_ZONES, s=1:DISTR_LOSS_SEGS, t=1:T], vWAUX_POS[z,s,t] <=  (inputs["pDistr_Max_Withdraw_Possible"][z]/DISTR_LOSS_SEGS)*vWAUX_POS_ON[z,s,t])
					@constraint(EP, cWAuxOrderNeg1[z=LOSS_ZONES, s=1:DISTR_LOSS_SEGS, t=1:T], vWAUX_NEG[z,s,t] <=  (inputs["pDistr_Max_Inject_Possible"][z]/DISTR_LOSS_SEGS)*vWAUX_NEG_ON[z,s,t])
					@constraint(EP, cWAuxOrderPos2[z=LOSS_ZONES, s=1:(DISTR_LOSS_SEGS-1), t=1:T], vWAUX_POS[z,s,t] >=  (inputs["pDistr_Max_Withdraw_Possible"][z]/DISTR_LOSS_SEGS)*vWAUX_POS_ON[z,s+1,t])
					@constraint(EP, cWAuxOrderNeg2[z=LOSS_ZONES, s=1:(DISTR_LOSS_SEGS-1), t=1:T], vWAUX_NEG[z,s,t] >=  (inputs["pDistr_Max_Inject_Possible"][z]/DISTR_LOSS_SEGS)*vWAUX_NEG_ON[z,s+1,t])
					# Eq 5: Binary constraints to deal with absolute value of (withdrawal - injection).
							# If (withdrawal - injection) is positive, vWAUX_POS segment 0 must be zero; If (withdrawal - injection) is negative, vWAUX_POS segment 0 must be positive
							# (and takes on value of the full net injection), forcing all vWAUX_POS other segments (s>=1) to be zero
					@constraint(EP, cWAuxSegmentZeroPos[z=LOSS_ZONES, t=1:T], vWAUX_POS[z,0,t] <= inputs["pDistr_Max_Inject_Possible"][z]*(1-vWAUX_POS_ON[z,1,t])) #inputs["pDistr_MaxP"][z]*(1-vWAUX_POS_ON[z,1,t]))
							# If (withdrawal - injection) is negative, vWAUX_NEG segment 0 must be zero; If (withdrawal - injection) is positive, vWAUX_NEG segment 0 must be positive
							# (and takes on value of the full net withdrawal), forcing all other vWAUX_NEG segments (s>=1) to be zero
					@constraint(EP, cWAuxSegmentZeroNeg[z=LOSS_ZONES, t=1:T], vWAUX_NEG[z,0,t] <= inputs["pDistr_Max_Withdraw_Possible"][z]*(1-vWAUX_NEG_ON[z,1,t])) #inputs["pDistr_MaxP"][z]*(1-vWAUX_NEG_ON[z,1,t]))
				end
				# Zones without a loss coefficient have no losses in all time periods (and above constraints are not applied)
				for t in 1:T
					for z in NOLOSS_ZONES
						setupperbound(vDLOSS[z,t], 0.0)
						setlowerbound(vDLOSS[z,t], 0.0)
					end
				end
			end # End if(DISTR_LOSS_SEGS > 0) block

		end # End Distribution zone related block

	end # End if(Z>1) multizonal constraints block

	### End Constraints ###

	dModuleArgs["dExpressions"] = dExpressions
	dModuleArgs["dObjective"] = dObjective
	dModuleArgs["dPowerBalance"] = dPowerBalance

	return EP, dModuleArgs
end

# End Infrastructure Functions

# Begin Technologies Functions
# Constraints specific to dispatchable renewables technologies
function dispatchable(EP::Model, dModuleArgs::Dict)

	println("Dispatchable Resources Module")

	inputs = dModuleArgs["inputs"]
	setup = dModuleArgs["setup"]
	dExpressions = dModuleArgs["dExpressions"]
	dObjective = dModuleArgs["dObjective"]
	dPowerBalance = dModuleArgs["dPowerBalance"]

	dfGen = inputs["dfGen"]

	T = inputs["T"]     # Number of time steps (hours)
	Z = inputs["Z"]     # Number of zones
	W = inputs["W"]     # Number of subperiods

	Tw = div.(T,W) #total number of hours per subperiod

	### Variables ###

	## Decision variables hybrid solar dc-coupled storage behind the inverter
	# Power from PV to inverter

	# @variable(EP, vPPVI[y=(dfGen[(dfGen[!,:DISP].==2) ,:][!,:R_ID]),t=1:T] >= 0)

	### End Variables ###

	### Expressions ###

	## Objective Function Expressions ##
	## End Objective Function Expressions ##

	## Power Balance Expressions ##

	@expression(EP, ePowerBalanceDisp[t=1:T, z=1:Z],
	sum(EP[:vP][y,t] for y=dfGen[(dfGen[!,:DISP].>=1) .&  (dfGen[!,:zone].==z),:][!,:R_ID]))

	dPowerBalance["ePowerBalance"] = dPowerBalance["ePowerBalance"] + ePowerBalanceDisp

	## End Power Balance Expressions ##

	### End Expressions ###

	### Constratints ###

	for w in 1:W
		for h in 1:Tw
			t = Tw*(w-1)+h
			## Dispatchable (curtailable) variable renewable generators
			### Option of modeling VRE generators with multiple availability profiles and capacity limits -  Num_VRE_bins in Generators_data.csv  >1
			## Default value of Num_VRE_bins ==1
			for y in dfGen[(dfGen[!,:DISP].==1),:][!,:R_ID]

				# Define the set of indices making up the different VRE bins
				VRE_BINS = dfGen[(dfGen[!,:R_ID].>=y) .&  (dfGen[!,:R_ID].<=y+dfGen[!,:Num_VRE_bins][y]-1),:][!,:R_ID]

				if (dfGen[!,:Num_VRE_bins][y] >=1) # For resource for which we are modeling hourly power output
					if (setup["Reserves"]==1)
						if (dfGen[!,:Reg_Up][y] == 0 && dfGen[!,:Reg_Dn][y] == 0 && dfGen[!,:Rsv_Up][y] ==0 && dfGen[!,:Rsv_Dn][y] == 0)
							setupperbound(EP[:vREG_UP][y,t], 0.0)
							setlowerbound(EP[:vREG_UP][y,t], 0.0)
							setupperbound(EP[:vREG_DN][y,t], 0.0)
							setlowerbound(EP[:vREG_DN][y,t], 0.0)
							setupperbound(EP[:vRSV_DN][y,t], 0.0)
							setlowerbound(EP[:vRSV_DN][y,t], 0.0)
							setupperbound(EP[:vRSV_UP][y,t], 0.0)
							setlowerbound(EP[:vRSV_UP][y,t], 0.0)
						end
						# Constraints on contribution to regulation and reserves
						# For VRE, reserve contributions must be less than the specified percentage of the capacity factor for the hour
						@constraint(EP, EP[:vREG_UP][y,t] <=
									dfGen[!,:Reg_Up][y]*sum(inputs["pP_Max"][yy,t]*(dfGen[!,:Existing_Cap_MW][yy] + EP[:vCAP][yy] - EP[:vRETCAP][yy])
										for yy =VRE_BINS))
						@constraint(EP, EP[:vREG_DN][y,t] <=
									dfGen[!,:Reg_Dn][y]*sum(inputs["pP_Max"][yy,t]*(dfGen[!,:Existing_Cap_MW][yy] + EP[:vCAP][yy] - EP[:vRETCAP][yy])
										for yy =VRE_BINS))
						@constraint(EP, EP[:vRSV_DN][y,t] <=
									dfGen[!,:Rsv_Up][y]*sum(inputs["pP_Max"][yy,t]*(dfGen[!,:Existing_Cap_MW][yy] + EP[:vCAP][yy] - EP[:vRETCAP][yy])
										for yy =VRE_BINS))
						@constraint(EP, EP[:vRSV_UP][y,t] <=
									dfGen[!,:Rsv_Dn][y]*sum(inputs["pP_Max"][yy,t]*(dfGen[!,:Existing_Cap_MW][yy] + EP[:vCAP][yy] - EP[:vRETCAP][yy])
										for yy =VRE_BINS))
						# Power generated and reserve contributions down per hour must be greater than zero
						@constraint(EP, EP[:vP][y,t] - EP[:vREG_DN][y,t] - EP[:vRSV_DN][y,t] >= 0)
						# Power generated and reserve contributions up per hour by renewable generators must be less than
						# hourly capacity factor. Note: inequality constraint allows curtailment of output
						# below maximum level.
						@constraint(EP, EP[:vP][y,t] + EP[:vREG_UP][y,t] + EP[:vREG_DN][y,t] <=
															sum(inputs["pP_Max"][yy,t]*(dfGen[!,:Existing_Cap_MW][yy] + EP[:vCAP][yy] - EP[:vRETCAP][yy])
						 									for yy =VRE_BINS))
					else
						# Maximum power generated per hour by renewable generators must be less than
						# sum of product of hourly capacity factor for each bin times its the bin installed capacity
						# Note: inequality constraint allows curtailment of output below maximum level.
						@constraint(EP, EP[:vP][y,t] <= sum(inputs["pP_Max"][yy,t]*(dfGen[!,:Existing_Cap_MW][yy] + EP[:vCAP][yy] - EP[:vRETCAP][yy])
						 									for yy =VRE_BINS))
					end

					# Set power variables for all bins that are not being modeled for hourly output to be zero
				else #  (dfGen[!,:Num_VRE_bins][y] == 0)
					setupperbound(EP[:vP][y,t], 0.0)
					setlowerbound(EP[:vP][y,t], 0.0)
					if (setup["Reserves"]==1)
						setupperbound(EP[:vREG_UP][y,t], 0.0)
						setlowerbound(EP[:vREG_UP][y,t], 0.0)
						setupperbound(EP[:vREG_DN][y,t], 0.0)
						setlowerbound(EP[:vREG_DN][y,t], 0.0)
						setupperbound(EP[:vRSV_DN][y,t], 0.0)
						setlowerbound(EP[:vRSV_DN][y,t], 0.0)
						setupperbound(EP[:vRSV_UP][y,t], 0.0)
						setlowerbound(EP[:vRSV_UP][y,t], 0.0)
					end
					# println("Zero Power VRE resources")
					# println(y)
				end # END for "if" loop on Power_Var ==0 or 1

				# Dispatchable technologies cannot store energy
				setupperbound(EP[:vCHARGE][y,t], 0.0)
				setlowerbound(EP[:vCHARGE][y,t], 0.0)
			end #END "for y"

		end #END "for h"
	end #END "for w"

	### End Constraints ###

	dModuleArgs["dExpressions"] = dExpressions
	dModuleArgs["dObjective"] = dObjective
	dModuleArgs["dPowerBalance"] = dPowerBalance

	return EP, dModuleArgs
end

function nondispatchable(EP::Model, dModuleArgs::Dict)

	println("Non-Dispatchable Resources Module")

	inputs = dModuleArgs["inputs"]
	setup = dModuleArgs["setup"]
	dExpressions = dModuleArgs["dExpressions"]
	dObjective = dModuleArgs["dObjective"]
	dPowerBalance = dModuleArgs["dPowerBalance"]

	dfGen = inputs["dfGen"]

	T = inputs["T"]     # Number of time steps (hours)
	Z = inputs["Z"]     # Number of zones
	W = inputs["W"]     # Number of subperiods

	Tw = div.(T,W) #total number of hours per subperiod

	### Variables ###
	### End Variables ###

	### Expressions ###

	## Objective Function Expressions ##
	## End Objective Function Expressions ##

	## Power Balance Expressions ##

	@expression(EP, ePowerBalanceNdisp[t=1:T, z=1:Z],
		sum(EP[:vP][y,t] for y=dfGen[(dfGen[!,:NDISP].==1) .&  (dfGen[!,:zone].==z),:][!,:R_ID]))

	dPowerBalance["ePowerBalance"] = dPowerBalance["ePowerBalance"] + ePowerBalanceNdisp

	## End Power Balance Expressions ##

	### End Expressions ###

	### Constratints ###

	for w in 1:W
		for h in 1:Tw
			t = Tw*(w-1)+h
			## Non-dispatchable (must-take)  renewable generators
			for y in dfGen[(dfGen[!,:NDISP].==1),:][!,:R_ID]
				# for y in dfGen[dfGen[!,:NDISP].==1),:][!,:R_ID]
					# Power generated per hour by renewable generators exactly equals
					# hourly capacity factor times installed capacity.
					@constraint(EP, EP[:vP][y,t] == inputs["pP_Max"][y,t]*(dfGen[!,:Existing_Cap_MW][y] + EP[:vCAP][y] - EP[:vRETCAP][y]))
					# Non-dispatchable renewable technologies cannot store energy
					setupperbound(EP[:vCHARGE][y,t], 0.0)
					setlowerbound(EP[:vCHARGE][y,t], 0.0)
					# No constraint on hourly ramp rate
					if (setup["Reserves"]==1)
						# Non-dispatchable renewables can't contribute to reserves
						setupperbound(EP[:vREG_UP][y,t], 0.0)
						setlowerbound(EP[:vREG_UP][y,t], 0.0)
						setupperbound(EP[:vREG_DN][y,t], 0.0)
						setlowerbound(EP[:vREG_DN][y,t], 0.0)
						setupperbound(EP[:vRSV_DN][y,t], 0.0)
						setlowerbound(EP[:vRSV_DN][y,t], 0.0)
						setupperbound(EP[:vRSV_UP][y,t], 0.0)
						setlowerbound(EP[:vRSV_UP][y,t], 0.0)
					end
			end #END "for y"
		end #END "for h"
	end #END "for w"

	### End Constraints ###

	dModuleArgs["dExpressions"] = dExpressions
	dModuleArgs["dObjective"] = dObjective
	dModuleArgs["dPowerBalance"] = dPowerBalance

	return EP, dModuleArgs
end

function storage(EP::Model, dModuleArgs::Dict)

	println("Storage Resources Module")

	inputs = dModuleArgs["inputs"]
	setup = dModuleArgs["setup"]
	dExpressions = dModuleArgs["dExpressions"]
	dObjective = dModuleArgs["dObjective"]
	dPowerBalance = dModuleArgs["dPowerBalance"]

	dfGen = inputs["dfGen"]

	G = inputs["G"]     # Number of resources (generators, storage, DR, and DERs)
	T = inputs["T"]     # Number of time steps (hours)
	Z = inputs["Z"]     # Number of zones
	W = inputs["W"]     # Number of subperiods

	Tw = div.(T,W) #total number of hours per subperiod

	###  New Variables related to storage ###

	# Bound storage level of resource "y" at hour "t" [MWh] on zone "z" for storage technologies
	for y in (dfGen[(dfGen[!,:STOR].>=1) ,:][!,:R_ID])
		for t in 1:T
			setlowerbound(EP[:vS][y,t], 0.0)
		end
	end
	# Variable to define direct conversion of power to gas (via electrolyzer that is part of storage unit)
	@variable(EP, vPower2Gas[y=1:G,t=1:T]>=0)
	# Variable to define storage discharge to supply gas demand (E.g. H2 storage)
	@variable(EP, vStorage2Gas[y=1:G,t=1:T]>=0)

	# Variables to define inter-period energy transferred between modeled periods
	# valid only when modeling representative periods and long-duration storage (STOR =3 or STOR=5)
	if (setup["OperationWrapping"]==2 && setup["LDS"]==1)
		NPeriods = size(inputs["Period_Map"])[1] # Number of modeled periods
		# State of charge of storage at beginning of each modeled period n
		# Apply only to (dfGen[(dfGen[!,:STOR].==2).| (dfGen[!,:STOR].==3) .| (dfGen[!,:STOR].==5) ,:][!,:R_ID])
		@variable(EP, vSOCw[y=1:G,n=1:NPeriods] >= 0)

		# Build up in storage inventory over each representative period w
		# Build up inventory can be positive or negative
		# Apply only to (dfGen[(dfGen[!,:STOR].==2).| (dfGen[!,:STOR].==3) .| (dfGen[!,:STOR].==5) ,:][!,:R_ID])
		@variable(EP, vdSOC[y=1:G,w=1:W])
	end
	### End Variables ###

	### Expressions ###

	## Objective Function Expressions ##
	## End Objective Function Expressions ##

	## Power Balance Expressions ##
	# Term to represent net dispatch from storage in any period
	@expression(EP, ePowerBalanceStor[t=1:T, z=1:Z],
		sum(EP[:vP][y,t]-EP[:vCHARGE][y,t] for y=dfGen[(dfGen[!,:STOR].>=1) .&  (dfGen[!,:zone].==z),:][!,:R_ID]))

	dPowerBalance["ePowerBalance"] = dPowerBalance["ePowerBalance"] + ePowerBalanceStor

	## End Power Balance Expressions ##

	##  RE Gas Demand Balance Expressions ##
	# Total Renewable gas suppply = supply from charging device + supply from storage discharging
	@expression(EP, eREGasSupply[t=1:T, z=1:Z],
		sum(dfGen[!,:Eff_up][y]*vPower2Gas[y,t] + vStorage2Gas[y,t]  for y=dfGen[(dfGen[!,:STOR].==5) .&  (dfGen[!,:zone].==z),:][!,:R_ID]))

	## End RE Gas Demand Balance Expressions ##

	### End Expressions ###

	### Constraints ###

	#### Define constraint name for storage inventory balance
	@constraintref SoCBal[1:T,1:G]

	if (setup["OperationWrapping"]==2 && setup["LDS"]==1)
		dfPeriodMap = inputs["Period_Map"] # Dataframe that maps modeled periods to representative periods
		NPeriods = size(inputs["Period_Map"])[1]
		# inter-period linking constraints connecting storage state of charge across Modeled periods
		# NOTE: These constraints are adapted from  1) Kotzur et al. 2018 - https://www.sciencedirect.com/science/article/pii/S0306261918300242
		# and 2) Zhang et al. 2019 - https://www.sciencedirect.com/science/article/pii/S0098135418306525
		for r  in 1:NPeriods # Modeled periods
			for y in dfGen[(dfGen[!,:STOR].==2) .| (dfGen[!,:STOR].==3) .| (dfGen[!,:STOR].==5),:][!,:R_ID]
				# Storage at beginning of period w = storage at beginning of period w-1+ storage built up in period w (after n representative periods)
				if r<NPeriods  #Multiply storage build up term from prior period with corresponding weight
					@constraint(EP, vSOCw[y,r+1] == vSOCw[y,r] + vdSOC[y,dfPeriodMap[!,:RepPeriod_index][r]])
				else # Last period is linked to first period
					@constraint(EP,vSOCw[y,1] == vSOCw[y,r] + vdSOC[y,dfPeriodMap[!,:RepPeriod_index][r]])
				end

				# Storage at beginning of each modeled period cannot exceed installed energy capacity
				@constraint(EP, vSOCw[y,r] <=dfGen[!,:Existing_Cap_MWh][y] + EP[:vCAPSTORAGE][y] - EP[:vRETCAPSTORAGE][y])
			end
		end
	end

	for w in 1:W  # Representative periods
		tw_min = Tw*(w-1)+1 # Beginning hour in period w
		tw_max = Tw*(w-1)+Tw # Ending  hour in period w

		for h in 1:Tw  # Looping of each hour of each week
			t = Tw*(w-1)+h  # time index for each model period

			# Electrochemical and mechanical energy storage
			## Constraints applied to all energy storage technologies
			for y in dfGen[(dfGen[!,:STOR].==1) .| (dfGen[!,:STOR].==2) .| (dfGen[!,:STOR].==3),:][!,:R_ID]

				if (setup["Reserves"]==1)
					if (dfGen[!,:Reg_Up][y] == 0 && dfGen[!,:Reg_Dn][y] == 0 && dfGen[!,:Rsv_Up][y] ==0 && dfGen[!,:Rsv_Dn][y] == 0)
						setupperbound(EP[:vREG_UP][y,t], 0.0)
						setlowerbound(EP[:vREG_UP][y,t], 0.0)
						setupperbound(EP[:vREG_DN][y,t], 0.0)
						setlowerbound(EP[:vREG_DN][y,t], 0.0)
						setupperbound(EP[:vRSV_DN][y,t], 0.0)
						setlowerbound(EP[:vRSV_DN][y,t], 0.0)
						setupperbound(EP[:vRSV_UP][y,t], 0.0)
						setlowerbound(EP[:vRSV_UP][y,t], 0.0)
					end
					# Maximum storage contribution to reserves is a specified fraction of installed capacity
					@constraint(EP, EP[:vREG_UP][y,t] <= dfGen[!,:Reg_Up][y]*(dfGen[!,:Existing_Cap_MW][y] + EP[:vCAP][y] - EP[:vRETCAP][y]))
					@constraint(EP, EP[:vREG_DN][y,t] <= dfGen[!,:Reg_Dn][y]*(dfGen[!,:Existing_Cap_MW][y] + EP[:vCAP][y] - EP[:vRETCAP][y]))
					@constraint(EP, EP[:vRSV_UP][y,t] <= dfGen[!,:Rsv_Up][y]*(dfGen[!,:Existing_Cap_MW][y] + EP[:vCAP][y] - EP[:vRETCAP][y]))
					@constraint(EP, EP[:vRSV_DN][y,t] <= dfGen[!,:Rsv_Dn][y]*(dfGen[!,:Existing_Cap_MW][y] + EP[:vCAP][y] - EP[:vRETCAP][y]))
					# Actual contribution to regulation and reserves is sum of auxilary variables for portions contributed during
					# charging and discharging
					@constraint(EP, EP[:vREG_UP][y,t] == EP[:vREG_UP_charge][y,t] + EP[:vREG_UP_discharge][y,t])
					@constraint(EP, EP[:vREG_DN][y,t] == EP[:vREG_DN_charge][y,t] + EP[:vREG_DN_discharge][y,t])
					@constraint(EP, EP[:vRSV_UP][y,t] == EP[:vRSV_UP_charge][y,t] + EP[:vRSV_UP_discharge][y,t])
					@constraint(EP, EP[:vRSV_DN][y,t] == EP[:vRSV_DN_charge][y,t] + EP[:vRSV_DN_discharge][y,t])
					# Maximum charging rate plus contribution to reserves up must be greater than zero
				    @constraint(EP, EP[:vCHARGE][y,t] - EP[:vREG_UP_charge][y,t] - EP[:vRSV_UP_charge][y,t] >= 0)
					# Maximum discharging rate and contribution to reserves down must be greater than zero
				    @constraint(EP, EP[:vP][y,t] - EP[:vREG_DN_discharge][y,t] - EP[:vRSV_DN_discharge][y,t] >= 0)
					# Maximum discharging rate and contribution to reserves up must be less than power rating OR available stored energy, whichever is less
				    @constraint(EP, EP[:vP][y,t] + EP[:vREG_UP_discharge][y,t] + EP[:vRSV_UP_discharge][y,t] <= (dfGen[!,:Existing_Cap_MW][y] + EP[:vCAP][y] - EP[:vRETCAP][y]))
				else
					# Maximum discharging rate must be less than power rating OR available stored energy, whichever is less
					@constraint(EP, EP[:vP][y,t] <= (dfGen[!,:Existing_Cap_MW][y] + EP[:vCAP][y] - EP[:vRETCAP][y]))
				end
			end

			## Energy storage technologies with fixed power/energy ratio
			for y in dfGen[(dfGen[!,:STOR].==1),:][!,:R_ID]

				# State of charge constraint (equals previous state + charge - discharge)
				if h == 1
					# Links last time step with first time step, ensuring position in hour 1 is within eligible change from final hour position
					#  Wrap-up Tolerance
						 SoCBal[t,y] = @constraint(EP, EP[:vS][y,t] == (1 - setup["pTolerance"])*EP[:vS][y,tw_max]-(1/dfGen[!,:Eff_down][y]*EP[:vP][y,t])+(dfGen[!,:Eff_up][y]*EP[:vCHARGE][y,t])-(dfGen[!,:Self_disch][y]*EP[:vS][y,tw_max]))
						 # energy stored for the next hour
						# SoCBal[y,t] @constraint(EP, EP[:vS][y,t] <= (1 + setup["pTolerance"])*EP[:vS][y,tw_max]-(1/dfGen[!,:Eff_down][y]*EP[:vP][y,t])+(dfGen[!,:Eff_up][y]*EP[:vCHARGE][y,t])-(dfGen[!,:Self_disch][y]*EP[:vS][y,tw_max])) #energy stored for the next hour
				else
						 SoCBal[t,y] =@constraint(EP, EP[:vS][y,t] == EP[:vS][y,t-1]-(1/dfGen[!,:Eff_down][y]*EP[:vP][y,t])+(dfGen[!,:Eff_up][y]*EP[:vCHARGE][y,t])-(dfGen[!,:Self_disch][y]*EP[:vS][y,t-1]))
						 #energy stored for the next hour
				end
					# Maximum energy stored must be less than energy capacity
					@constraint(EP, EP[:vS][y,t] <= (1/dfGen[!,:Ratio_power_to_energy][y])*(dfGen[!,:Existing_Cap_MW][y] + EP[:vCAP][y] - EP[:vRETCAP][y]))
					if (setup["Reserves"]==1)
						 # Storage units charging can charge faster to provide reserves down and charge slower to provide reserves up
						 # Maximum charging rate plus contribution to reserves down must be less than power rating OR available storage capacity, whichever is less
						@constraint(EP, EP[:vCHARGE][y,t] + EP[:vREG_DN_charge][y,t] + EP[:vRSV_DN_charge][y,t] <=
							(dfGen[!,:Existing_Cap_MW][y] + EP[:vCAP][y] - EP[:vRETCAP][y]))
					 	 # Max simultaneous charge and discharge rates cannot be greater than capacity
						@constraint(EP, (EP[:vP][y,t] + EP[:vREG_UP_discharge][y,t] + EP[:vRSV_UP_discharge][y,t])
						 +(EP[:vCHARGE][y,t] + EP[:vREG_DN_charge][y,t] + EP[:vRSV_DN_charge][y,t]) <= (dfGen[!,:Existing_Cap_MW][y] + EP[:vCAP][y] - EP[:vRETCAP][y]))
					else
						# Maximum charging rate must be less than power rating OR available storage capacity, whichever is less
						@constraint(EP, EP[:vCHARGE][y,t] <= (dfGen[!,:Existing_Cap_MW][y] + EP[:vCAP][y] - EP[:vRETCAP][y]))
						# Max simultaneous charge and discharge cannot be greater than capacity
						@constraint(EP, EP[:vP][y,t]+EP[:vCHARGE][y,t] <= (dfGen[!,:Existing_Cap_MW][y] + EP[:vCAP][y] - EP[:vRETCAP][y]))
					end
			end #END "for y"

			## Energy storage technologies with independent power and energy storage capacities (but symmetrical charge/discharge capacity)
			for y in dfGen[(dfGen[!,:STOR].==2),:][!,:R_ID]

				# State of charge constraint (equals previous state + charge - discharge)
				#NOTE: for longer-duration storage we might want to connect the weeks
				if h == 1
					# Links last time step with first time step, ensuring position in hour 1 is within eligible change from final hour position
					#  Wrap-up Tolerance
					# Modified initial state of storage for long-duration storage - initialize wth value carried over from last period
					if (setup["OperationWrapping"]==2 && setup["LDS"]==1)
						SoCBal[t,y] = @constraint(EP, EP[:vS][y,t] == (1-dfGen[!,:Self_disch][y])*(EP[:vS][y,tw_max] -vdSOC[y,w])-(1/dfGen[!,:Eff_down][y]*EP[:vP][y,t])+(dfGen[!,:Eff_up][y]*EP[:vCHARGE][y,t]))
					else
						 SoCBal[t,y] = @constraint(EP, EP[:vS][y,t] == (1 - setup["pTolerance"])*EP[:vS][y,tw_max]-(1/dfGen[!,:Eff_down][y]*EP[:vP][y,t])+(dfGen[!,:Eff_up][y]*EP[:vCHARGE][y,t])-(dfGen[!,:Self_disch][y]*EP[:vS][y,tw_max]))
					end
				else
					SoCBal[t,y] =@constraint(EP, EP[:vS][y,t] == EP[:vS][y,t-1]-(1/dfGen[!,:Eff_down][y]*EP[:vP][y,t])+(dfGen[!,:Eff_up][y]*EP[:vCHARGE][y,t])-(dfGen[!,:Self_disch][y]*EP[:vS][y,t-1]))
						 #energy stored for the next hour
				end
					# Maximum energy stored must be less than energy capacity
					@constraint(EP, EP[:vS][y,t] <= dfGen[!,:Existing_Cap_MWh][y] + EP[:vCAPSTORAGE][y] - EP[:vRETCAPSTORAGE][y])
					# Max and min constraints on energy storage capacity built (as proportion to power capacity)
					@constraint(EP, dfGen[!,:Min_Duration][y] * (dfGen[!,:Existing_Cap_MW][y] + EP[:vCAP][y] - EP[:vRETCAP][y]) <=  dfGen[!,:Existing_Cap_MWh][y] + EP[:vCAPSTORAGE][y] - EP[:vRETCAPSTORAGE][y])
					@constraint(EP, dfGen[!,:Max_Duration][y] * (dfGen[!,:Existing_Cap_MW][y] + EP[:vCAP][y] - EP[:vRETCAP][y]) >=  dfGen[!,:Existing_Cap_MWh][y] + EP[:vCAPSTORAGE][y] - EP[:vRETCAPSTORAGE][y])
					if (setup["Reserves"]==1)
						# Storage units charging can charge faster to provide reserves down and charge slower to provide reserves up
						 # Maximum charging rate plus contribution to reserves down must be less than power rating OR available storage capacity, whichever is less
						@constraint(EP, EP[:vCHARGE][y,t] + EP[:vREG_DN_charge][y,t] + EP[:vRSV_DN_charge][y,t] <= dfGen[!,:Existing_Cap_MW][y] + EP[:vCAP][y] - EP[:vRETCAP][y])
						 # Max simultaneous charge and discharge rates cannot be greater than capacity
						@constraint(EP, (EP[:vP][y,t] + EP[:vREG_UP_discharge][y,t] + EP[:vRSV_UP_discharge][y,t])
						 + (EP[:vCHARGE][y,t] + EP[:vREG_DN_charge][y,t] + EP[:vRSV_DN_charge][y,t]) <= (dfGen[!,:Existing_Cap_MW][y] + EP[:vCAP][y] - EP[:vRETCAP][y]))
					else
						# Maximum charging rate must be less than power rating OR available storage capacity, whichever is less
						@constraint(EP, EP[:vCHARGE][y,t] <= (dfGen[!,:Existing_Cap_MW][y] + EP[:vCAP][y] - EP[:vRETCAP][y]))
						# Max simultaneous charge and discharge cannot be greater than capacity
						@constraint(EP, EP[:vP][y,t]+EP[:vCHARGE][y,t] <= (dfGen[!,:Existing_Cap_MW][y] + EP[:vCAP][y] - EP[:vRETCAP][y]))
					end
			end #END "for y"

			## Energy storage technologies with independent charge and discharge power capacities and independent energy storage capacity
			for y in dfGen[(dfGen[!,:STOR].==3),:][!,:R_ID]

				# State of charge constraint (equals previous state + charge - discharge)
				#NOTE: for longer-duration storage we might want to connect the weeks
				if h == 1
					# Links last time step with first time step, ensuring position in hour 1 is within eligible change from final hour position
					#  Wrap-up Tolerance
					# Modified initial state of storage for long-duration storage - initialize wth value carried over from last representative period
					if (setup["OperationWrapping"]==2 && setup["LDS"]==1)
						SoCBal[t,y] = @constraint(EP, EP[:vS][y,t] == (1-dfGen[!,:Self_disch][y])*(EP[:vS][y,tw_max] -vdSOC[y,w])-(1/dfGen[!,:Eff_down][y]*EP[:vP][y,t])+(dfGen[!,:Eff_up][y]*EP[:vCHARGE][y,t]))
					else
						 SoCBal[t,y] = @constraint(EP, EP[:vS][y,t] == (1 - setup["pTolerance"])*EP[:vS][y,tw_max]-(1/dfGen[!,:Eff_down][y]*EP[:vP][y,t])+(dfGen[!,:Eff_up][y]*EP[:vCHARGE][y,t])-(dfGen[!,:Self_disch][y]*EP[:vS][y,tw_max]))
					end
				else
					SoCBal[t,y] =@constraint(EP, EP[:vS][y,t] == EP[:vS][y,t-1]-(1/dfGen[!,:Eff_down][y]*EP[:vP][y,t])+(dfGen[!,:Eff_up][y]*EP[:vCHARGE][y,t])-(dfGen[!,:Self_disch][y]*EP[:vS][y,t-1]))
						 #energy stored for the next hour
				end
					# Maximum energy stored must be less than energy capacity
					@constraint(EP, EP[:vS][y,t] <= (dfGen[!,:Existing_Cap_MWh][y] + EP[:vCAPSTORAGE][y] - EP[:vRETCAPSTORAGE][y]))
					# Max and min constraints on energy storage capacity built (as proportion to discharge power capacity)
					@constraint(EP, dfGen[!,:Min_Duration][y] * (dfGen[!,:Existing_Cap_MW][y] + EP[:vCAP][y] - EP[:vRETCAP][y]) <=  (dfGen[!,:Existing_Cap_MWh][y] + EP[:vCAPSTORAGE][y] - EP[:vRETCAPSTORAGE][y]))
					@constraint(EP, dfGen[!,:Max_Duration][y] * (dfGen[!,:Existing_Cap_MW][y] + EP[:vCAP][y] - EP[:vRETCAP][y]) >=  (dfGen[!,:Existing_Cap_MWh][y] + EP[:vCAPSTORAGE][y] - EP[:vRETCAPSTORAGE][y]))
					if (setup["Reserves"]==1)
						# Storage units charging can charge faster to provide reserves down and charge slower to provide reserves up
						 # Maximum charging rate plus contribution to reserves down must be less than charge power rating OR available storage capacity, whichever is less
						@constraint(EP, EP[:vCHARGE][y,t] + EP[:vREG_DN_charge][y,t] + EP[:vRSV_DN_charge][y,t] <= (dfGen[!,:Existing_Charge_Cap_MW][y] + EP[:vCAPCHARGE][y] - EP[:vRETCAPCHARGE][y]))
					else
						# Maximum charging rate must be less than charge power rating OR available storage capacity, whichever is less
						@constraint(EP, EP[:vCHARGE][y,t] <= (dfGen[!,:Existing_Charge_Cap_MW][y] + EP[:vCAPCHARGE][y] - EP[:vRETCAPCHARGE][y]))
					end
			end #END "for y"

			#  Storage technologies exporting product (e.g. H2) along with providing electrity storage AND independent charge, discharge power capacity and energy storage  capacity
			# STOR TYPE = 5
			# NOTE: Implemented August 7, 2020
			for y in dfGen[(dfGen[!,:STOR].==5) ,:][!,:R_ID]
				# Maximum energy stored must be less than energy capacity
				@constraint(EP, EP[:vS][y,t] <= (dfGen[!,:Existing_Cap_MWh][y] + EP[:vCAPSTORAGE][y] - EP[:vRETCAPSTORAGE][y]))
				# Max and min constraints on energy storage capacity built (as proportion to discharge power capacity)
				@constraint(EP, dfGen[!,:Min_Duration][y] * (dfGen[!,:Existing_Cap_MW][y] + EP[:vCAP][y] - EP[:vRETCAP][y]) <=  (dfGen[!,:Existing_Cap_MWh][y] + EP[:vCAPSTORAGE][y] - EP[:vRETCAPSTORAGE][y]))
				@constraint(EP, dfGen[!,:Max_Duration][y] * (dfGen[!,:Existing_Cap_MW][y] + EP[:vCAP][y] - EP[:vRETCAP][y]) >=  (dfGen[!,:Existing_Cap_MWh][y] + EP[:vCAPSTORAGE][y] - EP[:vRETCAPSTORAGE][y]))

				# State of charge constraint (equals previous state + charge to storage (excluding direct supply to end use) - discharge to power gen - discharge to direct end use)
				#NOTE: for longer-duration storage we might want to connect the weeks
				if h == 1
					# Links last time step with first time step, ensuring position in hour 1 is within eligible change from final hour position
					#  Wrap-up Tolerance
					# Modified initial state of storage for long-duration storage - initialize wth value carried over from last representative period
					if (setup["OperationWrapping"]==2 && setup["LDS"]==1)
						SoCBal[t,y] = @constraint(EP, EP[:vS][y,t] == (1-dfGen[!,:Self_disch][y])*(EP[:vS][y,tw_max] -vdSOC[y,w])-(1/dfGen[!,:Eff_down][y]*EP[:vP][y,t]) -(vStorage2Gas[y,t])+
						(dfGen[!,:Eff_up][y]*(EP[:vCHARGE][y,t] -vPower2Gas[y,t])))
					else
						SoCBal[t,y] = @constraint(EP, EP[:vS][y,t] == (1 - setup["pTolerance"])*EP[:vS][y,tw_max]-(1/dfGen[!,:Eff_down][y]*EP[:vP][y,t]) -(vStorage2Gas[y,t])+
						(dfGen[!,:Eff_up][y]*(EP[:vCHARGE][y,t] -vPower2Gas[y,t]))-(dfGen[!,:Self_disch][y]*EP[:vS][y,tw_max]))
					end
				else
						 SoCBal[t,y] =@constraint(EP, EP[:vS][y,t] == EP[:vS][y,t-1]-(1/dfGen[!,:Eff_down][y]*EP[:vP][y,t]) -(vStorage2Gas[y,t])+
						 					(dfGen[!,:Eff_up][y]*(EP[:vCHARGE][y,t] -vPower2Gas[y,t]))-(dfGen[!,:Self_disch][y]*EP[:vS][y,t-1]))
						 #energy stored for the next hour
				end
				if (setup["Reserves"]==1)
					if (dfGen[!,:Reg_Up][y] == 0 && dfGen[!,:Reg_Dn][y] == 0 && dfGen[!,:Rsv_Up][y] ==0 && dfGen[!,:Rsv_Dn][y] == 0)
						setupperbound(EP[:vREG_UP][y,t], 0.0)
						setlowerbound(EP[:vREG_UP][y,t], 0.0)
						setupperbound(EP[:vREG_DN][y,t], 0.0)
						setlowerbound(EP[:vREG_DN][y,t], 0.0)
						setupperbound(EP[:vRSV_DN][y,t], 0.0)
						setlowerbound(EP[:vRSV_DN][y,t], 0.0)
						setupperbound(EP[:vRSV_UP][y,t], 0.0)
						setlowerbound(EP[:vRSV_UP][y,t], 0.0)
					end
					# Maximum storage contribution to reserves is a specified fraction of installed capacity
					@constraint(EP, EP[:vREG_UP][y,t] <= dfGen[!,:Reg_Up][y]*(dfGen[!,:Existing_Cap_MW][y] + EP[:vCAP][y] - EP[:vRETCAP][y]))
					@constraint(EP, EP[:vREG_DN][y,t] <= dfGen[!,:Reg_Dn][y]*(dfGen[!,:Existing_Cap_MW][y] + EP[:vCAP][y] - EP[:vRETCAP][y]))
					@constraint(EP, EP[:vRSV_UP][y,t] <= dfGen[!,:Rsv_Up][y]*(dfGen[!,:Existing_Cap_MW][y] + EP[:vCAP][y] - EP[:vRETCAP][y]))
					@constraint(EP, EP[:vRSV_DN][y,t] <= dfGen[!,:Rsv_Dn][y]*(dfGen[!,:Existing_Cap_MW][y] + EP[:vCAP][y] - EP[:vRETCAP][y]))
					# Actual contribution to regulation and reserves is sum of auxilary variables for portions contributed during
					# charging and discharging
					@constraint(EP, EP[:vREG_UP][y,t] == EP[:vREG_UP_charge][y,t] + EP[:vREG_UP_discharge][y,t])
					@constraint(EP, EP[:vREG_DN][y,t] == EP[:vREG_DN_charge][y,t] + EP[:vREG_DN_discharge][y,t])
					@constraint(EP, EP[:vRSV_UP][y,t] == EP[:vRSV_UP_charge][y,t] + EP[:vRSV_UP_discharge][y,t])
					@constraint(EP, EP[:vRSV_DN][y,t] == EP[:vRSV_DN_charge][y,t] + EP[:vRSV_DN_discharge][y,t])
					# Maximum charging rate plus contribution to reserves up must be greater than zero
					@constraint(EP, EP[:vCHARGE][y,t] - EP[:vREG_UP_charge][y,t] - EP[:vRSV_UP_charge][y,t] >= 0)
					# Maximum discharging rate and contribution to reserves down must be greater than zero
					@constraint(EP, EP[:vP][y,t] - EP[:vREG_DN_discharge][y,t] - EP[:vRSV_DN_discharge][y,t] >= 0)
					# Maximum discharging rate and contribution to reserves up must be less than power rating OR available stored energy, whichever is less
					@constraint(EP, EP[:vP][y,t] + EP[:vREG_UP_discharge][y,t] + EP[:vRSV_UP_discharge][y,t] <= (dfGen[!,:Existing_Cap_MW][y] + EP[:vCAP][y] - EP[:vRETCAP][y]))

					# Storage units charging can charge faster to provide reserves down and charge slower to provide reserves up
					 # Maximum charging rate plus contribution to reserves down must be less than charge power rating OR available storage capacity, whichever is less
					@constraint(EP, EP[:vCHARGE][y,t] + EP[:vREG_DN_charge][y,t] + EP[:vRSV_DN_charge][y,t] <= (dfGen[!,:Existing_Charge_Cap_MW][y] + EP[:vCAPCHARGE][y] - EP[:vRETCAPCHARGE][y]))

				else
					# Maximum charging rate must be less than charge power rating OR available storage capacity, whichever is less
					@constraint(EP, EP[:vCHARGE][y,t] <= (dfGen[!,:Existing_Charge_Cap_MW][y] + EP[:vCAPCHARGE][y] - EP[:vRETCAPCHARGE][y]))

					# Maximum discharging rate must be less than power rating OR available stored energy, whichever is less
					@constraint(EP, EP[:vP][y,t] <= (dfGen[!,:Existing_Cap_MW][y] + EP[:vCAP][y] - EP[:vRETCAP][y]))
				end
			end #END "for y"

		end
	end

	### End Constraints ###
	dExpressions["SoCBal"] = SoCBal # Passing constraint names for SoC Balance constraint
	dExpressions["eREGasSupply"] =eREGasSupply # Passing expressions for Renewable GAs supply



	dModuleArgs["dExpressions"] = dExpressions
	dModuleArgs["dObjective"] = dObjective
	dModuleArgs["dPowerBalance"] = dPowerBalance


	return EP, dModuleArgs
end

function hydro(EP::Model, dModuleArgs::Dict)

	println("Hydro Resources Module")

	inputs = dModuleArgs["inputs"]
	setup = dModuleArgs["setup"]
	dExpressions = dModuleArgs["dExpressions"]
	dObjective = dModuleArgs["dObjective"]
	dPowerBalance = dModuleArgs["dPowerBalance"]

	dfGen = inputs["dfGen"]

	G = inputs["G"]     # Number of resources (generators, storage, DR, and DERs)
	T = inputs["T"]     # Number of time steps (hours)
	Z = inputs["Z"]     # Number of zones
	W = inputs["W"]     # Number of subperiods

	Tw = div.(T,W) #total number of hours per subperiod

	### Variables ###

	# State of charge variable is non-negative for hydro
	for y in (dfGen[(dfGen[!,:HYDRO].==1) ,:][!,:R_ID])
		for t in 1:T
			setlowerbound(EP[:vS][y,t], 0.0)
		end
	end

	# Hydro reservoir overflow
	@variable(EP, vSPILL[y=(dfGen[(dfGen[!,:HYDRO].==1) ,:][!,:R_ID]),t=1:T] >= 0)

	### End Variables ###

	### Expressions ###

	## Objective Function Expressions ##
	## End Objective Function Expressions ##

	## Power Balance Expressions ##
	@expression(EP, ePowerBalanceHydro[t=1:T, z=1:Z],
		sum(EP[:vP][y,t] for y=dfGen[(dfGen[!,:HYDRO].>=1) .&  (dfGen[!,:zone].==z),:][!,:R_ID]))

	dPowerBalance["ePowerBalance"] = dPowerBalance["ePowerBalance"] + ePowerBalanceHydro

	## End Power Balance Expressions ##

	### End Expressions ###

	### Constratints ###

	for w in 1:W
		for h in 1:Tw
			t = Tw*(w-1)+h
			tw_min = Tw*(w-1)+1
			tw_max = Tw*(w-1)+Tw
			## Hydro reservoir resources - with reservoir modeling
			for y in dfGen[(dfGen[!,:HYDRO].==1),:][!,:R_ID]
				# State of charge constraint (equals previous state + charge - discharge)
				if h == 1
					# Links last time step with first time step, ensuring position in hour 1 is within eligible change from final hour position
					#  Wrap-up Tolerance
					@constraint(EP, EP[:vS][y,t] >= (1 - setup["pTolerance"])*EP[:vS][y,tw_max]-(1/dfGen[!,:Eff_down][y]*EP[:vP][y,t]) - vSPILL[y,t] + inputs["pP_Max"][y,t]*(dfGen[!,:Existing_Cap_MW][y] + EP[:vCAP][y] - EP[:vRETCAP][y]))
					@constraint(EP, EP[:vS][y,t] <= (1 + setup["pTolerance"])*EP[:vS][y,tw_max]-(1/dfGen[!,:Eff_down][y]*EP[:vP][y,t]) - vSPILL[y,t] + inputs["pP_Max"][y,t]*(dfGen[!,:Existing_Cap_MW][y] + EP[:vCAP][y] - EP[:vRETCAP][y]))
					#Energy stored at the first hour as fraction of installed capacity less first period generation and spill + inflows
					@constraint(EP, EP[:vS][y,t] == (dfGen[!,:Hydro_level][y])*(1/dfGen[!,:Ratio_power_to_energy][y])*(dfGen[!,:Existing_Cap_MW][y] + EP[:vCAP][y] - EP[:vRETCAP][y]))
												- (1/dfGen[!,:Eff_down][y]*EP[:vP][y,t]) - vSPILL[y,t] +  inputs["pP_Max"][y,t]*(dfGen[!,:Existing_Cap_MW][y] + EP[:vCAP][y] - EP[:vRETCAP][y])
					# Maximum ramp up and down between consecutive hours
					@constraint(EP, EP[:vP][y,t] - EP[:vP][y,tw_max] <= dfGen[!,:Ramp_Up_percentage][y]*(dfGen[!,:Existing_Cap_MW][y] + EP[:vCAP][y] - EP[:vRETCAP][y]))
					@constraint(EP, EP[:vP][y,tw_max] - EP[:vP][y,t] <= dfGen[!,:Ramp_Dn_percentage][y]*(dfGen[!,:Existing_Cap_MW][y] + EP[:vCAP][y] - EP[:vRETCAP][y]))
					# Links last time step with first time step, ensuring position in hour 1 is within eligible ramp of final hour position
						# Superseded version which included spill rates in ramp constraints to limit fluctuations in streamflow from hour to hour;
						# revisit if appropriate but can produce infeasibilities as currently formulated so removed for now
						# @constraint(EP, EP[:vP][y,t] + vSPILL[y,t] - EP[:vP][y,tw_max] - vSPILL[y,tw_max] <= dfGen[!,:Ramp_Up_percentage][y]*(dfGen[!,:Existing_Cap_MW][y] + EP[:vCAP][y] - EP[:vRETCAP][y]))
						# @constraint(EP, EP[:vP][y,tw_max] + vSPILL[y,tw_max] - EP[:vP][y,t] - vSPILL[y,t] <= dfGen[!,:Ramp_Dn_percentage][y]*(dfGen[!,:Existing_Cap_MW][y] + EP[:vCAP][y] - EP[:vRETCAP][y]))
				else
					# Energy stored for the next hour (includes inflows)
					@constraint(EP, EP[:vS][y,t] == EP[:vS][y,t-1] - (1/dfGen[!,:Eff_down][y]*EP[:vP][y,t]) - vSPILL[y,t] +  inputs["pP_Max"][y,t]*(dfGen[!,:Existing_Cap_MW][y] + EP[:vCAP][y] - EP[:vRETCAP][y]))
					# Maximum ramp up and down between consecutive hours (applies only to power generation, not total streamflow inclusive of water spilled)
					@constraint(EP, EP[:vP][y,t] - EP[:vP][y,t-1] <= dfGen[:Ramp_Up_percentage][y]*(dfGen[:Existing_Cap_MW][y] + EP[:vCAP][y] - EP[:vRETCAP][y]))
					@constraint(EP, EP[:vP][y,t-1] - EP[:vP][y,t] <= dfGen[:Ramp_Dn_percentage][y]*(dfGen[:Existing_Cap_MW][y] + EP[:vCAP][y] - EP[:vRETCAP][y]))
						# Superseded version which included spill rates in ramp constraints to limit fluctuations in streamflow from hour to hour;
						# revisit if appropriate but can produce infeasibilities as currently formulated so removed for now
						# Maximum ramp up and down between consecutive hours (applies to total streamflow, including power generation & water spilled)
						# @constraint(EP, EP[:vP][y,t] + vSPILL[y,t] - EP[:vP][y,t-1] - vSPILL[y,t-1] <= dfGen[:Ramp_Up_percentage][y]*(dfGen[:Existing_Cap_MW][y] + EP[:vCAP][y] - EP[:vRETCAP][y]))
						# @constraint(EP, EP[:vP][y,t-1] + vSPILL[y,t-1] - EP[:vP][y,t] + vSPILL[y,t] <= dfGen[:Ramp_Dn_percentage][y]*(dfGen[:Existing_Cap_MW][y] + EP[:vCAP][y] - EP[:vRETCAP][y]))
				end
					# Maximum energy stored must be less than energy capacity
					@constraint(EP, EP[:vS][y,t] <= (1/dfGen[!,:Ratio_power_to_energy][y])*(dfGen[!,:Existing_Cap_MW][y] + EP[:vCAP][y] - EP[:vRETCAP][y]))
					# Minimum running requirements (does not currently include spilled water)
					@constraint(EP, EP[:vP][y,t] >= dfGen[!,:Min_power][y]*(dfGen[!,:Existing_Cap_MW][y] + EP[:vCAP][y] - EP[:vRETCAP][y]))
				if (setup["Reserves"]==1)
					if (dfGen[!,:Reg_Up][y] == 0 && dfGen[!,:Reg_Dn][y] == 0 && dfGen[!,:Rsv_Up][y] ==0 && dfGen[!,:Rsv_Dn][y] == 0)
						setupperbound(EP[:vREG_UP][y,t], 0.0)
						setlowerbound(EP[:vREG_UP][y,t], 0.0)
						setupperbound(EP[:vREG_DN][y,t], 0.0)
						setlowerbound(EP[:vREG_DN][y,t], 0.0)
						setupperbound(EP[:vRSV_DN][y,t], 0.0)
						setlowerbound(EP[:vRSV_DN][y,t], 0.0)
						setupperbound(EP[:vRSV_UP][y,t], 0.0)
						setlowerbound(EP[:vRSV_UP][y,t], 0.0)
					end
					# Maximum storage contribution to reserves is a specified fraction of installed capacity
					@constraint(EP, EP[:vREG_UP][y,t] <= dfGen[!,:Reg_Up][y]*(dfGen[!,:Existing_Cap_MW][y] + EP[:vCAP][y] - EP[:vRETCAP][y]))
					@constraint(EP, EP[:vREG_DN][y,t] <= dfGen[!,:Reg_Dn][y]*(dfGen[!,:Existing_Cap_MW][y] + EP[:vCAP][y] - EP[:vRETCAP][y]))
					@constraint(EP, EP[:vRSV_UP][y,t] <= dfGen[!,:Rsv_Up][y]*(dfGen[!,:Existing_Cap_MW][y] + EP[:vCAP][y] - EP[:vRETCAP][y]))
					@constraint(EP, EP[:vRSV_DN][y,t] <= dfGen[!,:Rsv_Dn][y]*(dfGen[!,:Existing_Cap_MW][y] + EP[:vCAP][y] - EP[:vRETCAP][y]))
					# Maximum discharging rate and contribution to reserves up must be less than power rating OR available stored energy, whichever is less
					@constraint(EP, EP[:vP][y,t] + EP[:vREG_UP][y,t] + EP[:vRSV_UP][y,t] <= (dfGen[!,:Existing_Cap_MW][y] + EP[:vCAP][y] - EP[:vRETCAP][y]))
						# Maximum discharging rate and contribution to reserves down must be greater than zero
					@constraint(EP, EP[:vP][y,t] - EP[:vREG_DN][y,t] - EP[:vRSV_DN][y,t] >= 0)
				else
					# Maximum discharging rate must be less than power rating OR available stored energy, whichever is less
					@constraint(EP, EP[:vP][y,t] <= (dfGen[!,:Existing_Cap_MW][y] + EP[:vCAP][y] - EP[:vRETCAP][y]))
					@constraint(EP, EP[:vP][y,t] <= EP[:vS][y,t])
				end
				# hydro resources cannot charge
				setupperbound(EP[:vCHARGE][y,t], 0.0)
				setlowerbound(EP[:vCHARGE][y,t], 0.0)
			end #END "for y"

			## Hydro power with reservoir resources - without reservoir modeling, only subperiod flow limits (outside time loop!)
			# Hydro = 2
			for y in dfGen[(dfGen[!,:HYDRO].==2),:][!,:R_ID]
				# Minimum running requirements (does not currently include spilled water)
				@constraint(EP, EP[:vP][y,t] >= dfGen[!,:Min_power][y]*(dfGen[!,:Existing_Cap_MW][y] + EP[:vCAP][y] - EP[:vRETCAP][y]))

				if (setup["Reserves"]==1)
					if (dfGen[!,:Reg_Up][y] == 0 && dfGen[!,:Reg_Dn][y] == 0 && dfGen[!,:Rsv_Up][y] ==0 && dfGen[!,:Rsv_Dn][y] == 0)
						setupperbound(EP[:vREG_UP][y,t], 0.0)
						setlowerbound(EP[:vREG_UP][y,t], 0.0)
						setupperbound(EP[:vREG_DN][y,t], 0.0)
						setlowerbound(EP[:vREG_DN][y,t], 0.0)
						setupperbound(EP[:vRSV_DN][y,t], 0.0)
						setlowerbound(EP[:vRSV_DN][y,t], 0.0)
						setupperbound(EP[:vRSV_UP][y,t], 0.0)
						setlowerbound(EP[:vRSV_UP][y,t], 0.0)
					end
					# Maximum storage contribution to reserves is a specified fraction of installed capacity
					@constraint(EP, EP[:vREG_UP][y,t] <= dfGen[!,:Reg_Up][y]*(dfGen[!,:Existing_Cap_MW][y] + EP[:vCAP][y] - EP[:vRETCAP][y]))
					@constraint(EP, EP[:vREG_DN][y,t] <= dfGen[!,:Reg_Dn][y]*(dfGen[!,:Existing_Cap_MW][y] + EP[:vCAP][y] - EP[:vRETCAP][y]))
					@constraint(EP, EP[:vRSV_UP][y,t] <= dfGen[!,:Rsv_Up][y]*(dfGen[!,:Existing_Cap_MW][y] + EP[:vCAP][y] - EP[:vRETCAP][y]))
					@constraint(EP, EP[:vRSV_DN][y,t] <= dfGen[!,:Rsv_Dn][y]*(dfGen[!,:Existing_Cap_MW][y] + EP[:vCAP][y] - EP[:vRETCAP][y]))

					# Maximum discharging rate and contribution to reserves up must be less than power rating
					@constraint(EP, EP[:vP][y,t] + EP[:vREG_UP][y,t] + EP[:vRSV_UP][y,t] <= (dfGen[!,:Existing_Cap_MW][y] + EP[:vCAP][y] - EP[:vRETCAP][y]))

						# Maximum discharging rate and contribution to reserves down must be greater than zero
					@constraint(EP, EP[:vP][y,t] - EP[:vREG_DN][y,t] - EP[:vRSV_DN][y,t] >= 0)
				else
					# Maximum discharging rate must be less than power rating OR available stored energy, whichever is less
					@constraint(EP, EP[:vP][y,t] <= (dfGen[!,:Existing_Cap_MW][y] + EP[:vCAP][y] - EP[:vRETCAP][y]))

				end
				# hydro resources cannot store energy
				setupperbound(EP[:vCHARGE][y,t], 0.0)
				setlowerbound(EP[:vCHARGE][y,t], 0.0)
			end #END "for y"

			## Non-dispatchable (must-take) Run-of-river hydropower plants (HYDRO=3)
			for y in dfGen[(dfGen[!,:HYDRO].==3),:][!,:R_ID]
				# Power generated per hour by renewable generators exactly equals
				# hourly capacity factor times installed capacity.
				@constraint(EP, EP[:vP][y,t] == inputs["pP_Max"][y,t]*(dfGen[!,:Existing_Cap_MW][y] + EP[:vCAP][y] - EP[:vRETCAP][y]))
				# Non-dispatchable renewable technologies cannot store energy
				setupperbound(EP[:vCHARGE][y,t], 0.0)
				setlowerbound(EP[:vCHARGE][y,t], 0.0)
				# No constraint on hourly ramp rate
				if (setup["Reserves"]==1)
					# Non-dispatchable renewables can't contribute to reserves
					setupperbound(EP[:vREG_UP][y,t], 0.0)
					setlowerbound(EP[:vREG_UP][y,t], 0.0)
					setupperbound(EP[:vREG_DN][y,t], 0.0)
					setlowerbound(EP[:vREG_DN][y,t], 0.0)
					setupperbound(EP[:vRSV_DN][y,t], 0.0)
					setlowerbound(EP[:vRSV_DN][y,t], 0.0)
					setupperbound(EP[:vRSV_UP][y,t], 0.0)
					setlowerbound(EP[:vRSV_UP][y,t], 0.0)
				end
			end #END "for y"
		end
	end

	#### SPECIAL CONSTRAINT for Reservoir Hydropower modeled without reservoir capacity limit - energy generation bounds for each day or each supperiod
	# NOTE: this constraint will only work in the case of modeling representative weeks
	if setup["OperationWrapping"]==2
		# Case when energy generation must adhere to bounds for each subperiod (Hydro_Res_flexibility_days=0)
			for y in dfGen[(dfGen[!,:HYDRO].==2) .& (dfGen[!,:Hydro_Res_flexibility_days].==0),:][!,:R_ID]
				for w in 1:W
					tw_min = Tw*(w-1)+1 # starting index of each subperiod
					tw_max = Tw*(w-1)+Tw # ending index of each subperiod
				# limits on energy generation for each sub period - total number of constraint = total number of subperiods
					@constraint(EP, sum(EP[:vP][y,e] for e=tw_min:tw_max) <= sum(inputs["pP_Max"][y,e] for e=tw_min:tw_max)*(dfGen[!,:Existing_Cap_MW][y] + EP[:vCAP][y] - EP[:vRETCAP][y]))
				end # end of subperiod loop
			end # End of y loop

		# Case when energy generation must adhere to bounds for each day (Hydro_Res_flexibility_days=1)
			for y in dfGen[(dfGen[!,:HYDRO].==2) .& (dfGen[!,:Hydro_Res_flexibility_days].==1),:][!,:R_ID]
				NumDays = div.(Tw,24)  # Number of days per subperiod
				for w in 1:W
					for d in 1:NumDays #total number of days per subperiod
						tw_min = Tw*(w-1)+ 24*(d-1)+1 # starting index of each day
						tw_max = Tw*(w-1)+ 24*(d-1)+24 # ending index of each day
				# limits on energy generation for each day - total number of constraint = total number of subperiods*days per subperiod
						@constraint(EP, sum(EP[:vP][y,e] for e=tw_min:tw_max) <= sum(inputs["pP_Max"][y,e] for e=tw_min:tw_max)*(dfGen[!,:Existing_Cap_MW][y] + EP[:vCAP][y] - EP[:vRETCAP][y]))
					end  # for daily loop
				end # end of subperiod loop
			end # End of y loop
		end # End if condition

	### End Constraints ###

	dModuleArgs["dExpressions"] = dExpressions
	dModuleArgs["dObjective"] = dObjective
	dModuleArgs["dPowerBalance"] = dPowerBalance

	return EP, dModuleArgs
end

function dr(EP::Model, dModuleArgs::Dict)

	println("Demand Response Resources Module")

	inputs = dModuleArgs["inputs"]
	setup = dModuleArgs["setup"]
	dExpressions = dModuleArgs["dExpressions"]
	dObjective = dModuleArgs["dObjective"]
	dPowerBalance = dModuleArgs["dPowerBalance"]

	dfGen = inputs["dfGen"]

	T = inputs["T"]     # Number of time steps (hours)
	Z = inputs["Z"]     # Number of zones
	W = inputs["W"]     # Number of subperiods

	Tw = div.(T,W) #total number of hours per subperiod

	### Variables ###
	### End Variables ###

	### Expressions ###

	## Objective Function Expressions ##
	## End Objective Function Expressions ##

	## Power Balance Expressions ##
	@expression(EP, ePowerBalanceDr[t=1:T, z=1:Z],
		sum(-EP[:vP][y,t]+EP[:vCHARGE][y,t] for y=dfGen[(dfGen[!,:DR].>=1) .&  (dfGen[!,:zone].==z),:][!,:R_ID]))

	dPowerBalance["ePowerBalance"] = dPowerBalance["ePowerBalance"] + ePowerBalanceDr

	## End Power Balance Expressions ##

	### End Expressions ###

	### Constratints ###

	for w in 1:W
		for h in 1:Tw
			t = Tw*(w-1)+h
			tw_min = Tw*(w-1)+1
			tw_max = Tw*(w-1)+Tw

			## Flexible demand side management available during all hours and can be either delayed or advanced (virtual storage-shiftable demand) - DR ==1
			for z in 1:Z
				# NOTE: DSM operates by zone since capacity is now related to zone demand
				for y in dfGen[(dfGen[!,:DR].==1) .&  (dfGen[!,:zone].==z),:][!,:R_ID]
					# State of "charge" constraint (equals previous state + charge - discharge)
					if h <= 1
						# Links last time step with first time step, ensuring position in hour 1 is within eligible change from final hour position
						# DSM_Energy_Eff corresponds to energy loss due to time shifting
						@constraint(EP, EP[:vS][y,t] == EP[:vS][y,tw_max]-dfGen[!,:DSM_energy_eff][y]*(EP[:vP][y,t])+(EP[:vCHARGE][y,t]))
					else
						@constraint(EP, EP[:vS][y,t] == EP[:vS][y,t-1]-dfGen[!,:DSM_energy_eff][y]*(EP[:vP][y,t])+(EP[:vCHARGE][y,t]))
					end # END if
					# Maximum energy "stored" or deferred for later hours
					# NOTE: the maximum capacity is given by the hours demand can be deferred and maximum amount of demand thta can be deferred
					#@constraint(EP, EP[:vS][y,t] <= (dfGen[!,:Existing_Cap_MW][y] + EP[:vCAP][y] - EP[:vRETCAP][y])/dfGen[!,:Ratio_power_to_energy][y])
					# Maximum charging rate
					# NOTE: the maximum amount that can be shifted is given by a percentage of the hourly demand
					@constraint(EP, EP[:vCHARGE][y,t] <= dfGen[!,:Ratio_power_to_energy][y]*inputs["pD"][t,z])
					# Maximum discharging rate
					# NOTE: no maximum discharge rate unless constrained by other factors like transmission, etc.
					#@constraint(EP, EP[:vP][y,t] <= (dfGen[!,:Existing_Cap_MW][y] + EP[:vCAP][y] - EP[:vRETCAP][y]))
					# Require deferred demands to be satisfied within the specified time delay
					if (Tw-dfGen[!,:Max_DSM_delay][y]) < h < Tw
						# Constraint wraps around to first hours of time series
							@constraint(EP, sum(EP[:vP][y,e] for e=(t+1):tw_max)+sum(EP[:vP][y,e] for e=tw_min:(tw_min-1+dfGen[!,:Max_DSM_delay][y]-(Tw-h))) >= EP[:vS][y,t])
					elseif h == Tw
						# Constraint wraps around to first hours of time series
							@constraint(EP, sum(EP[:vP][y,e] for e=tw_min:(tw_min-1+dfGen[!,:Max_DSM_delay][y])) >= EP[:vS][y,t])
					else
						# Constraint looks back over last n hours, where n = inputs["Max_DSM_delay"][y]
							@constraint(EP, sum(EP[:vP][y,e] for e=(t+1):(t+dfGen[!,:Max_DSM_delay][y])) >= EP[:vS][y,t])
					end #END if

					# Require deferred demands to be satisfied within the specified time advanced
					if (Tw-dfGen[!,:Max_DSM_advance][y]) < h < Tw
						# Constraint wraps around to first hours of time series
							@constraint(EP, sum(EP[:vCHARGE][y,e] for e=(t+1):tw_max)+sum(EP[:vCHARGE][y,e] for e=tw_min:(tw_min-1+dfGen[!,:Max_DSM_advance][y]-(Tw-h))) >= -EP[:vS][y,t])
					elseif h == Tw
						# Constraint wraps around to first hours of time series
							@constraint(EP, sum(EP[:vCHARGE][y,e] for e=tw_min:(tw_min-1+dfGen[!,:Max_DSM_advance][y])) >= -EP[:vS][y,t])
					else
						# Constraint looks for over next n hours, where n = inputs["Max_DSM_advance"][y]
							@constraint(EP, sum(EP[:vCHARGE][y,e] for e=(t+1):(t+dfGen[!,:Max_DSM_advance][y])) >= -EP[:vS][y,t])
					end #END if

					if (setup["Reserves"]==1)
						# NOTE: DMS can't contribute to reserves
						setupperbound(EP[:vREG_UP][y,t], 0.0)
						setlowerbound(EP[:vREG_UP][y,t], 0.0)
						setupperbound(EP[:vREG_DN][y,t], 0.0)
						setlowerbound(EP[:vREG_DN][y,t], 0.0)
						setupperbound(EP[:vRSV_DN][y,t], 0.0)
						setlowerbound(EP[:vRSV_DN][y,t], 0.0)
						setupperbound(EP[:vRSV_UP][y,t], 0.0)
						setlowerbound(EP[:vRSV_UP][y,t], 0.0)
					end
				end #END "for y"
			end #END "for z"

			## Flexible demand side management available only during specified hours with time delay or time advance  (virtual storage-shiftable demand) - DR ==2
			for z in 1:Z
				# NOTE: DSM operates by zone since capacity is now related to zone demand
				for y in dfGen[(dfGen[!,:DR].==2) .&  (dfGen[!,:zone].==z),:][!,:R_ID]
					# State of "charge" constraint (equals previous state + charge - discharge)
					if h <= 1
							# Links last time step with first time step, ensuring position in hour 1 is within eligible change from final hour position
							# DSM_Energy_Eff corresponds to energy loss due to time shifting
							@constraint(EP, EP[:vS][y,t] == EP[:vS][y,tw_max]-dfGen[!,:DSM_energy_eff][y]*(EP[:vP][y,t])+(EP[:vCHARGE][y,t]))
					else
							@constraint(EP, EP[:vS][y,t] == EP[:vS][y,t-1]-dfGen[!,:DSM_energy_eff][y]*(EP[:vP][y,t])+(EP[:vCHARGE][y,t]))
					end # END if
					# Maximum energy "stored" or deferred for later hours
					# NOTE: the maximum capacity is given by the hours demand can be deferred and maximum amount of demand thta can be deferred
					#@constraint(EP, EP[:vS][y,t] <= (dfGen[!,:Existing_Cap_MW][y] + EP[:vCAP][y] - EP[:vRETCAP][y])/dfGen[!,:Ratio_power_to_energy][y])
					# Maximum charging rate
					# NOTE: the maximum amount that can be shifted is given by a percentage of the hourly demand
					@constraint(EP, EP[:vCHARGE][y,t] <= inputs["pP_Max"][y,t]*inputs["pD"][t,z])
					# Maximum discharging rate
					# NOTE: no maximum discharge rate unless constrained by other factors like transmission, etc.
					#@constraint(EP, EP[:vP][y,t] <= (dfGen[!,:Existing_Cap_MW][y] + EP[:vCAP][y] - EP[:vRETCAP][y]))
					# Require deferred demands to be satisfied within the specified time limit
					if (Tw-dfGen[!,:Max_DSM_delay][y]) < h < Tw
						# Constraint wraps around to first hours of time series
							@constraint(EP, sum(EP[:vP][y,e] for e=(t+1):tw_max)+sum(EP[:vP][y,e] for e=tw_min:(tw_min-1+dfGen[!,:Max_DSM_delay][y]-(Tw-h))) >= EP[:vS][y,t])
					elseif h == Tw
						# Constraint wraps around to first hours of time series
							@constraint(EP, sum(EP[:vP][y,e] for e=tw_min:(tw_min-1+dfGen[!,:Max_DSM_delay][y]-(Tw-h))) >= EP[:vS][y,t])
					else
						# Constraint looks back over last n hours, where n = inputs["pDMS_Time"][y]
							@constraint(EP, sum(EP[:vP][y,e] for e=(t+1):(t+dfGen[!,:Max_DSM_delay][y])) >= EP[:vS][y,t])
					end #END if

					# Require deferred demands to be satisfied within the specified time advanced
					if (Tw-dfGen[!,:Max_DSM_advance][y]) < h < Tw
						# Constraint wraps around to first hours of time series
							@constraint(EP, sum(EP[:vCHARGE][y,e] for e=(t+1):tw_max)+sum(EP[:vCHARGE][y,e] for e=tw_min:(tw_min-1+dfGen[!,:Max_DSM_advance][y]-(Tw-h))) >= -EP[:vS][y,t])
					elseif h == Tw
						# Constraint wraps around to first hours of time series
							@constraint(EP, sum(EP[:vCHARGE][y,e] for e=tw_min:(tw_min-1+dfGen[!,:Max_DSM_advance][y])) >= -EP[:vS][y,t])
					else
						# Constraint looks for over next n hours, where n = inputs["Max_DSM_advance"][y]
							@constraint(EP, sum(EP[:vCHARGE][y,e] for e=(t+1):(t+dfGen[!,:Max_DSM_advance][y])) >= -EP[:vS][y,t])
					end #END if

					if (setup["Reserves"]==1)
						# NOTE: DMS can't contribute to reserves
						setupperbound(EP[:vREG_UP][y,t], 0.0)
						setlowerbound(EP[:vREG_UP][y,t], 0.0)
						setupperbound(EP[:vREG_DN][y,t], 0.0)
						setlowerbound(EP[:vREG_DN][y,t], 0.0)
						setupperbound(EP[:vRSV_DN][y,t], 0.0)
						setlowerbound(EP[:vRSV_DN][y,t], 0.0)
						setupperbound(EP[:vRSV_UP][y,t], 0.0)
						setlowerbound(EP[:vRSV_UP][y,t], 0.0)
					end
				end #END "for y"
			end #END "for z"
		end
	end

	### End Constraints ###

	dModuleArgs["dExpressions"] = dExpressions
	dModuleArgs["dObjective"] = dObjective
	dModuleArgs["dPowerBalance"] = dPowerBalance

	return EP, dModuleArgs
end

function heat(EP::Model, dModuleArgs::Dict)

	println("Heat Resources Module")

	inputs = dModuleArgs["inputs"]
	setup = dModuleArgs["setup"]
	dExpressions = dModuleArgs["dExpressions"]
	dObjective = dModuleArgs["dObjective"]
	dPowerBalance = dModuleArgs["dPowerBalance"]

	dfGen = inputs["dfGen"]

	T = inputs["T"]     # Number of time steps (hours)
	Z = inputs["Z"]     # Number of zones
	W = inputs["W"]     # Number of subperiods

	Tw = div.(T,W) #total number of hours per subperiod

	### Variables ###

	## Heat-electriciy and NACC interaction configuration for variables
	if (setup["HeatMarket"]==1)
		## Decision variables for heat-electricity interaction
		# Heat going to heat market
		@variable(EP, vHEAT_MARKET[y=(dfGen[(dfGen[:HEAT].==1) ,:][!,:R_ID]),t=1:T] >= 0)
		# Heat from heat storage that is going to generation
		@variable(EP, vHEAT_GENERATION[y=(dfGen[(dfGen[!,:HEAT].==1) ,:][!,:R_ID]),t=1:T] >= 0)
		## NOTE Variables for on-site heat storage
		# Storage charge capacity
		@variable(EP, vCAPCHARGE[y=(dfGen[(dfGen[!,:HEAT].==2) ,:][!,:R_ID])] >= 0)
		# Storage discharge capacity
		@variable(EP, vCAPDISCHARGE[y=(dfGen[(dfGen[!,:HEAT].==2) ,:][!,:R_ID])] >= 0)

		# Heat from supplemental natural gas that is going to generation
		@variable(EP, vHEAT_NG[y=(dfGen[(dfGen[!,:NACC].==1) ,:][!,:R_ID]),t=1:T] >= 0)
		## NOTE Variables for on-site heat storage
		# P2p variable
		@variable(EP, vP2P[y=(dfGen[(dfGen[!,:NACC].==2) ,:][!,:R_ID]),t=1:T] >= 0)
		# Power charging TES
		@variable(EP, vH1[y=(dfGen[(dfGen[!,:NACC].==2) ,:][!,:R_ID]),t=1:T] >= 0)
		# Power discharging TES
		@variable(EP, vH2[y=(dfGen[(dfGen[!,:HEAT].==2) ,:][!,:R_ID]),t=1:T] >= 0)
	end #END heat-electricity configuration

	### End Variables ###

	### Expressions ###

	## Objective Function Expressions ##

	## Heat-electriciy interaction configuration for expressions
	if (setup["HeatMarket"]==1)
		# Revenues from heat market during hour "t"
		@expression(EP, eHeat_Rev[y=(dfGen[(dfGen[!,:HEAT].==1) ,:][!,:R_ID]),t=1:T],inputs["omega"][t]*(inputs["pHeat_Price"][1]*vHEAT_MARKET[y,t]))
		dExpressions["eHeat_Rev"] = eHeat_Rev

		@expression(EP, eTotalHeat_Rev, sum(eHeat_Rev[y,t] for  y in (dfGen[(dfGen[!,:HEAT].==1) ,:][!,:R_ID]), t in 1:T))

		# Cost from heat market during hour "t"
		@expression(EP, eHeat_Cost[y=(dfGen[(dfGen[!,:NACC].==1) ,:][!,:R_ID]),t=1:T],inputs["omega"][t]*(inputs["pHeat_Price"][1]*vHEAT_NG[y,t]))
		dExpressions["eHeat_Cost"] = eHeat_Cost

		@expression(EP, eTotalHeat_Cost, sum(eHeat_Cost[y,t] for  y in (dfGen[(dfGen[!,:NACC].==1) ,:][!,:R_ID]), t in 1:T))
	else
		@expression(EP, eTotalHeat_Rev, 0)
		@expression(EP, eTotalHeat_Cost, 0)
	end

	dObjective["eTotalHeat_Cost"] = eTotalHeat_Cost
	dObjective["eTotalHeat_Rev"] = eTotalHeat_Rev
	## End Objective Function Expressions ##

	## Power Balance Expressions ##

	if setup["HeatMarket"]==1
		# Adding term related to charging heat storage
		@expression(EP, ePowerBalanceHeatCharge[t=1:T, z=1:Z],
			-sum(EP[:vCHARGE][y,t] for y in dfGen[(dfGen[!,:HEAT].==1) .&  (dfGen[!,:zone].==z),:][!,:R_ID]))
		# Terms related to power supply from coupled power plant + heat storage unit
		# 1st term - baseload power (Nuclear with NACC=1) + efficiency of using NG from pre-heated air from nuclear
		# 2nd term - power generated via stored heat (without NG)
		# 3rd term - power generated via NACC of type 2
		# NOTE: not sure what is NACC of type 2 is - consider deleting later on
		@expression(EP, ePowerBalanceHeat[t=1:T, z=1:Z],
			sum(EP[:vP][y,t] + dfGen[!,:NACC_Eff][y]*vHEAT_NG[y,t] for y in dfGen[(dfGen[!,:NACC].==1) .&  (dfGen[!,:zone].==z),:][!,:R_ID])
			+ sum(dfGen[!,:NACC_Eff][y]*vHEAT_GENERATION[x,t] for y in dfGen[(dfGen[!,:NACC].==1) .&  (dfGen[!,:zone].==z),:][!,:R_ID], x in dfGen[(dfGen[!,:HEAT].==1) .&  (dfGen[!,:zone].==z),:][!,:R_ID])
			+ sum(EP[:vP][y,t]  for y in dfGen[(dfGen[!,:NACC].==2) .&  (dfGen[!,:zone].==z),:][!,:R_ID]))

		dPowerBalance["ePowerBalance"] = dPowerBalance["ePowerBalance"] + ePowerBalanceCharge
		dPowerBalance["ePowerBalance"] = dPowerBalance["ePowerBalance"] + ePowerBalanceHeat
	end

	## End Power Balance Expressions ##

	### End Expressions ###

	### Constraints ###
	for w in 1:W
		for h in 1:Tw
			t = Tw*(w-1)+h
			tw_min = Tw*(w-1)+1
			tw_max = Tw*(w-1)+Tw

			## Heat-electricity configuration for operational constraints
			if (setup["UCommit"]==0 && setup["HeatMarket"]==1)
				## NACC generators WHITOUT unit commitment constraints
				for z in 1:Z
						for y in (dfGen[(dfGen[!,:NACC].==1) .&  (dfGen[!,:zone].==z),:][!,:R_ID])
						for x in (dfGen[(dfGen[!,:HEAT].==1) .&  (dfGen[!,:zone].==z),:][!,:R_ID])
							if (setup["Reserves"]==1)
								if (dfGen[!,:Reg_Up][y] == 0 && dfGen[!,:Reg_Dn][y] == 0 && dfGen[!,:Rsv_Up][y] ==0 && dfGen[!,:Rsv_Dn][y] == 0)
									setupperbound(EP[:vREG_UP][y,t], 0.0)
									setlowerbound(EP[:vREG_UP][y,t], 0.0)
									setupperbound(EP[:vREG_DN][y,t], 0.0)
									setlowerbound(EP[:vREG_DN][y,t], 0.0)
									setupperbound(EP[:vRSV_DN][y,t], 0.0)
									setlowerbound(EP[:vRSV_DN][y,t], 0.0)
									setupperbound(EP[:vRSV_UP][y,t], 0.0)
									setlowerbound(EP[:vRSV_UP][y,t], 0.0)
								end
								# Constraints on contribution to regulation and reserves
								@constraint(EP, EP[:vREG_UP][y,t] <= dfGen[!,:Reg_Up][y]*dfGen[!,:NACC_Peak_to_Base][y]*(dfGen[!,:Existing_Cap_MW][y] + EP[:vCAP][y] - EP[:vRETCAP][y]))
								@constraint(EP, EP[:vREG_DN][y,t] <= dfGen[!,:Reg_Dn][y]*dfGen[!,:NACC_Peak_to_Base][y]*(dfGen[!,:Existing_Cap_MW][y] + EP[:vCAP][y] - EP[:vRETCAP][y]))
								@constraint(EP, EP[:vRSV_DN][y,t] <= dfGen[!,:Rsv_Up][y]*dfGen[!,:NACC_Peak_to_Base][y]*(dfGen[!,:Existing_Cap_MW][y] + EP[:vCAP][y] - EP[:vRETCAP][y]))
								@constraint(EP, EP[:vRSV_UP][y,t] <= dfGen[!,:Rsv_Dn][y]*dfGen[!,:NACC_Peak_to_Base][y]*(dfGen[!,:Existing_Cap_MW][y] + EP[:vCAP][y] - EP[:vRETCAP][y]))
								# Minimum and maximum stable power generated per technology "y" at hour "t" Min_power for base capacity
								@constraint(EP, EP[:vP][y,t] == dfGen[!,:Min_power][y]*(dfGen[!,:Existing_Cap_MW][y] + EP[:vCAP][y] - EP[:vRETCAP][y]))
								# Minimum stable power generated per technology "y" at hour "t" Min_power
								@constraint(EP, dfGen[!,:NACC_Eff][y]*(vHEAT_NG[y,t] + vHEAT_GENERATION[x,t]) - EP[:vREG_DN][y,t] - EP[:vRSV_DN][y,t] >= 0)
								# Maximum power generated per peak capacity technology "y" at hour "t"
								@constraint(EP, dfGen[!,:NACC_Eff][y]*(vHEAT_NG[y,t] + vHEAT_GENERATION[x,t]) + EP[:vREG_UP][y,t] + EP[:vRSV_UP][y,t] <= dfGen[!,:NACC_Peak_to_Base][y]*(dfGen[!,:Existing_Cap_MW][y] + EP[:vCAP][y] - EP[:vRETCAP][y]))
							else
								# Minimum stable power generated per technology "y" at hour "t" Min_power
								@constraint(EP, EP[:vP][y,t] == dfGen[!,:Min_power][y]*(dfGen[!,:Existing_Cap_MW][y] + EP[:vCAP][y] - EP[:vRETCAP][y]))
								@constraint(EP, dfGen[!,:NACC_Eff][y]*(vHEAT_NG[y,t] + vHEAT_GENERATION[x,t]) <= dfGen[!,:NACC_Peak_to_Base][y]*(dfGen[!,:Existing_Cap_MW][y] + EP[:vCAP][y] - EP[:vRETCAP][y]))
							end
						end #END "x" dim for heat storage
						end #END "y" dim for NACC technology
						## Heat storage technologies
					for y in (dfGen[(dfGen[!,:HEAT].==1),:][!,:R_ID])
						# State of charge constraint (equals previous state + charge - discharge)
						if h <= 1
							# Links last time step with first time step, ensuring position in hour 1 is within eligible change from final hour position
							#  Wrap-up Tolerance
							@constraint(EP, EP[:vS][y,t] == (1 - setup["pTolerance"])*EP[:vS][y,tw_max]-(1/dfGen[!,:Eff_down][y]*EP[:vP][y,t])+(dfGen[!,:Eff_up][y]*EP[:vCHARGE][y,t])-(dfGen[!,:Self_disch][y]*EP[:vS][y,tw_max])) #energy stored for the next hour
						else
							@constraint(EP, EP[:vS][y,t] == EP[:vS][y,t-1]-(1/dfGen[!,:Eff_down][y]*EP[:vP][y,t])+(dfGen[!,:Eff_up][y]*EP[:vCHARGE][y,t])-(dfGen[!,:Self_disch][y]*EP[:vS][y,t]))  #energy stored for the next hour
						end #END if
						# Maximum energy stored must be less than energy capacity
						@constraint(EP, EP[:vS][y,t] <= (1/dfGen[!,:Ratio_power_to_energy][y])*(dfGen[!,:Existing_Cap_MW][y] + EP[:vCAP][y] - EP[:vRETCAP][y]))
						if (setup["Reserves"]==1)
							if (dfGen[!,:Reg_Up][y] == 0 && dfGen[!,:Reg_Dn][y] == 0 && dfGen[!,:Rsv_Up][y] ==0 && dfGen[!,:Rsv_Dn][y] == 0)
								setupperbound(EP[:vREG_UP][y,t], 0.0)
								setlowerbound(EP[:vREG_UP][y,t], 0.0)
								setupperbound(EP[:vREG_DN][y,t], 0.0)
								setlowerbound(EP[:vREG_DN][y,t], 0.0)
								setupperbound(EP[:vRSV_DN][y,t], 0.0)
								setlowerbound(EP[:vRSV_DN][y,t], 0.0)
								setupperbound(EP[:vRSV_UP][y,t], 0.0)
								setlowerbound(EP[:vRSV_UP][y,t], 0.0)
							end
							# Maximum storage contribution to reserves is a specified fraction of installed capacity
							@constraint(EP, EP[:vREG_UP][y,t] <= dfGen[!,:Reg_Up][y]*(dfGen[!,:Existing_Cap_MW][y] + EP[:vCAP][y] - EP[:vRETCAP][y]))
							@constraint(EP, EP[:vREG_DN][y,t] <= dfGen[!,:Reg_Dn][y]*(dfGen[!,:Existing_Cap_MW][y] + EP[:vCAP][y] - EP[:vRETCAP][y]))
							@constraint(EP, EP[:vRSV_UP][y,t] <= dfGen[!,:Rsv_Up][y]*(dfGen[!,:Existing_Cap_MW][y] + EP[:vCAP][y] - EP[:vRETCAP][y]))
							@constraint(EP, EP[:vRSV_DN][y,t] <= dfGen[!,:Rsv_Dn][y]*(dfGen[!,:Existing_Cap_MW][y] + EP[:vCAP][y] - EP[:vRETCAP][y]))
							# Storage units charging can charge faster to provide reserves down and charge slower to provide reserves up
							# Maximum charging rate plus contribution to reserves down must be less than power rating
							@constraint(EP, EP[:vCHARGE][y,t] + EP[:vREG_DN][y,t] + EP[:vRSV_DN][y,t] <= (dfGen[!,:Existing_Cap_MW][y] + EP[:vCAP][y] - EP[:vRETCAP][y]))
							# Maximum charging rate plus contribution to reserves up must be greater than zero
							@constraint(EP, EP[:vCHARGE][y,t] - EP[:vREG_UP][y,t] - EP[:vRSV_UP][y,t] >= 0)

						else
							# Maximum charging rate must be less than power rating
							@constraint(EP, EP[:vCHARGE][y,t] <= (dfGen[!,:Existing_Cap_MW][y] + EP[:vCAP][y] - EP[:vRETCAP][y]))

						end

						# Maximum discharging rate and contribution to reserves up must be less than power rating
						# NOTE: this is heat discharging so it cannot supply reserves to electricity market
						@constraint(EP, EP[:vP][y,t] <= (dfGen[!,:Existing_Cap_MW][y] + EP[:vCAP][y] - EP[:vRETCAP][y]))

						# Heat variables from storage discharge
						@constraint(EP, vHEAT_MARKET[y,t] + vHEAT_GENERATION[y,t] == EP[:vP][y,t])
						@constraint(EP, sum(vHEAT_MARKET[y,t] for y=dfGen[(dfGen[!,:HEAT].==1) .&  (dfGen[!,:zone].==z),:][!,:R_ID]) <= inputs["pHeat_Demand"][t,z])

					end #END "for y"
				end #END "for z"

			# setting variables in case no CHP system is allowed
			elseif (setup["UCommit"]==0 && setup["HeatMarket"]==0)
				for y in dfGen[(dfGen[!,:NACC].==1) .| (dfGen[!,:HEAT].==1) ,:][!,:R_ID]
					setupperbound(EP[:vP][y,t], 0.0)
					setlowerbound(EP[:vP][y,t], 0.0)
					setupperbound(EP[:vCHARGE][y,t], 0.0)
					setlowerbound(EP[:vCHARGE][y,t], 0.0)
					setupperbound(EP[:vS][y,t], 0.0)
					setlowerbound(EP[:vS][y,t], 0.0)
					setupperbound(EP[:vCAP][y], 0.0)
					setlowerbound(EP[:vCAP][y], 0.0)
				end #END "for y"
			end #END CHP system configuration for operational constraints # Loop of operational constraints

			## NACC Stuff
			## Advanced nuclear generators (NACC) WITH unit commitment constraints
			if (setup["HeatMarket"]==1 && setup["UCommit"]>=1 )
				for z in 1:Z
						for y in (dfGen[(dfGen[!,:NACC].==1) .&  (dfGen[!,:zone].==z),:][!,:R_ID])
							for x in (dfGen[(dfGen[!,:HEAT].==1) .&  (dfGen[!,:zone].==z),:][!,:R_ID])
								if (setup["Reserves"]==1)
									if (dfGen[!,:Reg_Up][y] == 0 && dfGen[!,:Reg_Dn][y] == 0 && dfGen[!,:Rsv_Up][y] ==0 && dfGen[!,:Rsv_Dn][y] == 0)
										setupperbound(EP[:vREG_UP][y,t], 0.0)
										setlowerbound(EP[:vREG_UP][y,t], 0.0)
										setupperbound(EP[:vREG_DN][y,t], 0.0)
										setlowerbound(EP[:vREG_DN][y,t], 0.0)
										setupperbound(EP[:vRSV_DN][y,t], 0.0)
										setlowerbound(EP[:vRSV_DN][y,t], 0.0)
										setupperbound(EP[:vRSV_UP][y,t], 0.0)
										setlowerbound(EP[:vRSV_UP][y,t], 0.0)
									end
									# Constraints on contribution to regulation and reserves
									@constraint(EP, EP[:vREG_UP][y,t] <= dfGen[!,:Reg_Up][y]*dfGen[!,:NACC_Eff][y]*dfGen[!,:NACC_Peak_to_Base][y]*dfGen[!,:Cap_size][y]*EP[:vCOMMIT][y,t])
									@constraint(EP, EP[:vREG_DN][y,t] <= dfGen[!,:Reg_Dn][y]*dfGen[!,:NACC_Eff][y]*dfGen[!,:NACC_Peak_to_Base][y]*dfGen[!,:Cap_size][y]*EP[:vCOMMIT][y,t])
									@constraint(EP, EP[:vRSV_DN][y,t] <= dfGen[!,:Rsv_Up][y]*dfGen[!,:NACC_Eff][y]*dfGen[!,:NACC_Peak_to_Base][y]*dfGen[!,:Cap_size][y]*EP[:vCOMMIT][y,t])
									@constraint(EP, EP[:vRSV_UP][y,t] <= dfGen[!,:Rsv_Dn][y]*dfGen[!,:NACC_Eff][y]*dfGen[!,:NACC_Peak_to_Base][y]*dfGen[!,:Cap_size][y]*EP[:vCOMMIT][y,t])
									# Minimum and maximum stable power generated per technology "y" at hour "t" Min_power for base capacity
									# If committed, NACC operates at fixed output that cannot be adjusted; flexibility comes from peaker capability
									@constraint(EP, EP[:vP][y,t] == dfGen[!,:Cap_size][y]*EP[:vCOMMIT][y,t])
									# Minimum stable power generated per technology "y" at hour "t" Min_power
									@constraint(EP, dfGen[!,:NACC_Eff][y]*(vHEAT_NG[y,t] + vHEAT_GENERATION[x,t]) - EP[:vREG_DN][y,t] - EP[:vRSV_DN][y,t] >= 0)
									# Maximum power generated per peak capacity technology "y" at hour "t"
									@constraint(EP, dfGen[!,:NACC_Eff][y]*(vHEAT_NG[y,t] + vHEAT_GENERATION[x,t]) + EP[:vREG_UP][y,t] + EP[:vRSV_UP][y,t] <= dfGen[!,:NACC_Peak_to_Base][y]*(dfGen[!,:Cap_size][y]*EP[:vCOMMIT][y,t]))
								else
									# Minimum stable power generated per technology "y" at hour "t" Min_power
									@constraint(EP, EP[:vP][y,t] == dfGen[!,:Cap_size][y]*EP[:vCOMMIT][y,t])
									@constraint(EP, dfGen[!,:NACC_Eff][y]*(vHEAT_NG[y,t] + vHEAT_GENERATION[x,t]) <= dfGen[!,:NACC_Peak_to_Base][y]*(dfGen[!,:Cap_size][y]*EP[:vCOMMIT][y,t]))
								end
							end #END "x" dim for heat storage
						end #END "y" dim for NACC technology
					## Heat storage technologies
					for y in (dfGen[(dfGen[!,:HEAT].==1),:][!,:R_ID])
							# State of charge constraint (equals previous state + charge - discharge)
							if h <= 1
									# Links last time step with first time step, ensuring position in hour 1 is within eligible change from final hour position
									#  Wrap-up Tolerance
									@constraint(EP, EP[:vS][y,t] == (1 - setup["pTolerance"])*EP[:vS][y,tw_max]-(1/dfGen[!,:Eff_down][y]*EP[:vP][y,t])+(dfGen[!,:Eff_up][y]*EP[:vCHARGE][y,t])-(dfGen[!,:Self_disch][y]*EP[:vS][y,tw_max])) #energy stored for the next hour
							else
									@constraint(EP, EP[:vS][y,t] == EP[:vS][y,t-1]-(1/dfGen[!,:Eff_down][y]*EP[:vP][y,t])+(dfGen[!,:Eff_up][y]*EP[:vCHARGE][y,t])-(dfGen[!,:Self_disch][y]*EP[:vS][y,t]))  #energy stored for the next hour
							end #END if
							# Maximum energy stored must be less than energy capacity
							@constraint(EP, EP[:vS][y,t] <= (1/dfGen[!,:Ratio_power_to_energy][y])*(dfGen[!,:Existing_Cap_MW][y] + EP[:vCAP][y] - EP[:vRETCAP][y]))
							if (setup["Reserves"]==1)
								if (dfGen[!,:Reg_Up][y] == 0 && dfGen[!,:Reg_Dn][y] == 0 && dfGen[!,:Rsv_Up][y] ==0 && dfGen[!,:Rsv_Dn][y] == 0)
									setupperbound(EP[:vREG_UP][y,t], 0.0)
									setlowerbound(EP[:vREG_UP][y,t], 0.0)
									setupperbound(EP[:vREG_DN][y,t], 0.0)
									setlowerbound(EP[:vREG_DN][y,t], 0.0)
									setupperbound(EP[:vRSV_DN][y,t], 0.0)
									setlowerbound(EP[:vRSV_DN][y,t], 0.0)
									setupperbound(EP[:vRSV_UP][y,t], 0.0)
									setlowerbound(EP[:vRSV_UP][y,t], 0.0)
								end
								# Maximum storage contribution to reserves is a specified fraction of installed capacity
								@constraint(EP, EP[:vREG_UP][y,t] <= dfGen[!,:Reg_Up][y]*(dfGen[!,:Existing_Cap_MW][y] + EP[:vCAP][y] - EP[:vRETCAP][y]))
								@constraint(EP, EP[:vREG_DN][y,t] <= dfGen[!,:Reg_Dn][y]*(dfGen[!,:Existing_Cap_MW][y] + EP[:vCAP][y] - EP[:vRETCAP][y]))
								@constraint(EP, EP[:vRSV_UP][y,t] <= dfGen[!,:Rsv_Up][y]*(dfGen[!,:Existing_Cap_MW][y] + EP[:vCAP][y] - EP[:vRETCAP][y]))
								@constraint(EP, EP[:vRSV_DN][y,t] <= dfGen[!,:Rsv_Dn][y]*(dfGen[!,:Existing_Cap_MW][y] + EP[:vCAP][y] - EP[:vRETCAP][y]))
								# Storage units charging can charge faster to provide reserves down and charge slower to provide reserves up
								# Maximum charging rate plus contribution to reserves down must be less than power rating
								@constraint(EP, EP[:vCHARGE][y,t] + EP[:vREG_DN][y,t] + EP[:vRSV_DN][y,t] <= (dfGen[!,:Existing_Cap_MW][y] + EP[:vCAP][y] - EP[:vRETCAP][y]))
								# Maximum charging rate plus contribution to reserves up must be greater than zero
								@constraint(EP, EP[:vCHARGE][y,t] - EP[:vREG_UP][y,t] - EP[:vRSV_UP][y,t] >= 0)
								# Maximum discharge rate (in terms of heat output) cannot exceed installed heat capacity
								@constraint(EP, EP[:vP][y,t] <= (dfGen[!,:Existing_Cap_MW][y] + EP[:vCAP][y] - EP[:vRETCAP][y]))

							else
								# Maximum charging rate must be less than power rating
								@constraint(EP, EP[:vCHARGE][y,t] <= ((1/dfGen[!,:Eff_up][y])*(dfGen[!,:Existing_Cap_MW][y] + EP[:vCAP][y] - EP[:vRETCAP][y])))
							end
							# Maximum discharging rate must be less than power rating
							@constraint(EP, EP[:vP][y,t] <= (dfGen[!,:Eff_down][y])*(dfGen[!,:Existing_Cap_MW][y] + EP[:vCAP][y] - EP[:vRETCAP][y]))

							# Heat variables from storage discharge
							@constraint(EP, vHEAT_MARKET[y,t] + vHEAT_GENERATION[y,t] == EP[:vP][y,t])
							@constraint(EP, sum(vHEAT_MARKET[y,t] for y=dfGen[(dfGen[!,:HEAT].==1) .&  (dfGen[!,:zone].==z),:][!,:R_ID]) <= inputs["pHeat_Demand"][t,z])
					end #END "for y"f
					# Heat variables are zero for non heat storage technologies
					for y in dfGen[(dfGen[!,:HEAT].==0),:][!,:R_ID]
							setupperbound(vHEAT_MARKET[y,t], 0.0)
							setlowerbound(vHEAT_MARKET[y,t], 0.0)
							setupperbound(vHEAT_GENERATION[y,t], 0.0)
							setlowerbound(vHEAT_GENERATION[y,t], 0.0)
					end #END "for y"
				end #END "for z"
				for y in dfGen[(dfGen[!,:NACC].==1) .&  (dfGen[!,:Commit].==1),:][!,:R_ID]
					@constraint(EP, EP[:vCOMMIT][y,t] <= ((dfGen[!,:Existing_Cap_MW][y]/dfGen[!,:Cap_size][y]) + EP[:vCAP][y] - EP[:vRETCAP][y]))
					@constraint(EP, EP[:vSTART][y,t]  <= ((dfGen[!,:Existing_Cap_MW][y]/dfGen[!,:Cap_size][y]) + EP[:vCAP][y] - EP[:vRETCAP][y]))
					@constraint(EP, EP[:vSHUT][y,t]   <= ((dfGen[!,:Existing_Cap_MW][y]/dfGen[!,:Cap_size][y]) + EP[:vCAP][y] - EP[:vRETCAP][y]))
					# Thermal technologies cannot store energy
					setupperbound(EP[:vCHARGE][y,t], 0.0)
					setlowerbound(EP[:vCHARGE][y,t], 0.0)
					# Minimum up and down times
					# up time
						if  1 < h < (dfGen[!,:Up_time][y]+1)
							# Constraint wraps around to first hours of time series
							@constraint(EP, EP[:vCOMMIT][y,t] >= sum(EP[:vSTART][y,e] for e=(t-(h-1):t))+sum(EP[:vSTART][y,e] for e=(tw_max-(dfGen[!,:Up_time][y]-h)):tw_max))
						elseif h == 1
							# Constraint looks back over last n hours, where n = inputs["pDMS_Time"][y]
							@constraint(EP, EP[:vCOMMIT][y,t]  >= EP[:vSTART][y,t] + sum(EP[:vSTART][y,e] for e=(tw_max-(dfGen[!,:Up_time][y] -1 )):tw_max))
						else
							# Constraint looks back over last n hours, where n = inputs["pDMS_Time"][y]
							@constraint(EP, EP[:vCOMMIT][y,t]  >= sum(EP[:vSTART][y,e] for e=(t-dfGen[!,:Up_time][y]):t))
						end
					# down time
						if  1 < h < (dfGen[!,:Down_time][y]+1)
							# Constraint wraps around to first hours of time series
							@constraint(EP, ((dfGen[!,:Existing_Cap_MW][y]/dfGen[!,:Cap_size][y]) + EP[:vCAP][y] - EP[:vRETCAP][y]) - EP[:vCOMMIT][y,t]
													>= sum(EP[:vSHUT][y,e] for e=(t-(h-1):t))+sum(EP[:vSHUT][y,e] for e=(tw_max-(dfGen[!,:Down_time][y]-h)):tw_max))
						elseif h == 1
							# Constraint looks back over last n hours, where n = inputs["pDMS_Time"][y]
							@constraint(EP, ((dfGen[!,:Existing_Cap_MW][y]/dfGen[!,:Cap_size][y]) + EP[:vCAP][y] - EP[:vRETCAP][y]) - EP[:vCOMMIT][y,t]
													>= EP[:vSHUT][y,t] + sum(EP[:vSHUT][y,e] for e=(tw_max-(dfGen[!,:Down_time][y] - 1 )):tw_max))
						else
							# Constraint looks back over last n hours, where n = inputs["pDMS_Time"][y]
							@constraint(EP, ((dfGen[!,:Existing_Cap_MW][y]/dfGen[!,:Cap_size][y]) + EP[:vCAP][y] - EP[:vRETCAP][y]) - EP[:vCOMMIT][y,t]
													>= sum(EP[:vSHUT][y,e] for e=(t-dfGen[!,:Down_time][y]):t))
						end
					# Startup and shutdown events
						if t <= 1
							@constraint(EP, EP[:vCOMMIT][y,t] == EP[:vCOMMIT][y,tw_max] + EP[:vSTART][y,t] - EP[:vSHUT][y,t] )
						else
							@constraint(EP, EP[:vCOMMIT][y,t] == EP[:vCOMMIT][y,t-1] + EP[:vSTART][y,t] - EP[:vSHUT][y,t] )
						end
				end #END "for y"
			# setting variables in case no CHP system is allowed
			else
				for y in dfGen[(dfGen[!,:NACC].==1) .| (dfGen[!,:HEAT].==1) ,:][!,:R_ID]
					setupperbound(EP[:vP][y,t], 0.0)
					setlowerbound(EP[:vP][y,t], 0.0)
					setupperbound(EP[:vCHARGE][y,t], 0.0)
					setlowerbound(EP[:vCHARGE][y,t], 0.0)
					setupperbound(EP[:vS][y,t], 0.0)
					setlowerbound(EP[:vS][y,t], 0.0)
					setupperbound(EP[:vCAP][y], 0.0)
					setlowerbound(EP[:vCAP][y], 0.0)
				end #END "for y"
			end #END "if" CHP system configuration for operational constraints

		end
	end

	dModuleArgs["dExpressions"] = dExpressions
	dModuleArgs["dObjective"] = dObjective
	dModuleArgs["dPowerBalance"] = dPowerBalance

	return EP, dModuleArgs
end

function thermal(EP::Model, dModuleArgs::Dict)

	println("Thermal Resources Module")

	inputs = dModuleArgs["inputs"]
	setup = dModuleArgs["setup"]
	dExpressions = dModuleArgs["dExpressions"]
	dObjective = dModuleArgs["dObjective"]
	dPowerBalance = dModuleArgs["dPowerBalance"]

	dfGen = inputs["dfGen"]

	T = inputs["T"]     # Number of time steps (hours)
	Z = inputs["Z"]     # Number of zones
	W = inputs["W"]     # Number of subperiods

	Tw = div.(T,W) #total number of hours per subperiod

	### Variables ###
	### End Variables ###

	### Expressions ###

	## Objective Function Expressions ##
	## End Objective Function Expressions ##

	## Power Balance Expressions ##
	@expression(EP, ePowerBalanceTherm[t=1:T, z=1:Z],
	sum(EP[:vP][y,t] for y in  dfGen[(dfGen[!,:THERM].==1) .&  (dfGen[!,:zone].==z),:][!,:R_ID]))

	dPowerBalance["ePowerBalance"] = dPowerBalance["ePowerBalance"] + ePowerBalanceTherm

	## End Power Balance Expressions ##

	### End Expressions ###

	### Constratints ###

	for w in 1:W
		for h in 1:Tw
			t = Tw*(w-1)+h
			tw_min = Tw*(w-1)+1
			tw_max = Tw*(w-1)+Tw

			## Integer Unit Commitment configuration for operational constraints
			if (setup["UCommit"]>=1)
				## Thermal generators WITH unit commitment constraints
				for y in dfGen[(dfGen[!,:THERM].==1) .& (dfGen[!,:Commit].==1),:][!,:R_ID]
					@constraint(EP, EP[:vCOMMIT][y,t] <= ((dfGen[!,:Existing_Cap_MW][y]/dfGen[!,:Cap_size][y]) + EP[:vCAP][y] - EP[:vRETCAP][y]))
					@constraint(EP, EP[:vSTART][y,t] <= ((dfGen[!,:Existing_Cap_MW][y]/dfGen[!,:Cap_size][y]) + EP[:vCAP][y] - EP[:vRETCAP][y]))
					@constraint(EP, EP[:vSHUT][y,t] <= ((dfGen[!,:Existing_Cap_MW][y]/dfGen[!,:Cap_size][y]) + EP[:vCAP][y] - EP[:vRETCAP][y]))

					# Maximum ramp up and down between consecutive hours
					if h <= 1
						# Links last time step with first time step, ensuring position in hour 1 is within eligible ramp of final hour position
						# DEVELOPMENT NOTEs: We should make this a configurable option...
							#rampup constraints
							@constraint(EP, EP[:vP][y,t]-EP[:vP][y,tw_max] <= dfGen[!,:Ramp_Up_percentage][y]*dfGen[!,:Cap_size][y]*(EP[:vCOMMIT][y,t]-EP[:vSTART][y,t])
																											 + min(inputs["pP_Max"][y,t],max(dfGen[!,:Min_power][y],dfGen[!,:Ramp_Up_percentage][y]))*dfGen[!,:Cap_size][y]*EP[:vSTART][y,t]
																											 - dfGen[!,:Min_power][y]*dfGen[!,:Cap_size][y]*EP[:vSHUT][y,t])
							#rampdown constraints
							@constraint(EP, EP[:vP][y,tw_max] - EP[:vP][y,t] <= dfGen[!,:Ramp_Dn_percentage][y]*dfGen[!,:Cap_size][y]*(EP[:vCOMMIT][y,t]-EP[:vSTART][y,t])
																											 - dfGen[!,:Min_power][y]*dfGen[!,:Cap_size][y]*EP[:vSTART][y,t]
																											 + min(inputs["pP_Max"][y,t],max(dfGen[!,:Min_power][y],dfGen[!,:Ramp_Dn_percentage][y]))*dfGen[!,:Cap_size][y]*EP[:vSHUT][y,t])
					else
							#rampup constraints
							@constraint(EP, EP[:vP][y,t]-EP[:vP][y,t-1] <= dfGen[!,:Ramp_Up_percentage][y]*dfGen[!,:Cap_size][y]*(EP[:vCOMMIT][y,t]-EP[:vSTART][y,t])
																											 + min(inputs["pP_Max"][y,t],max(dfGen[!,:Min_power][y],dfGen[!,:Ramp_Up_percentage][y]))*dfGen[!,:Cap_size][y]*EP[:vSTART][y,t]
																											 -dfGen[!,:Min_power][y]*dfGen[!,:Cap_size][y]*EP[:vSHUT][y,t])
							#rampdown constraints
							@constraint(EP, EP[:vP][y,t-1] - EP[:vP][y,t] <= dfGen[!,:Ramp_Dn_percentage][y]*dfGen[!,:Cap_size][y]*(EP[:vCOMMIT][y,t]-EP[:vSTART][y,t])
																											 -dfGen[!,:Min_power][y]*dfGen[!,:Cap_size][y]*EP[:vSTART][y,t]
																											 + min(inputs["pP_Max"][y,t],max(dfGen[!,:Min_power][y],dfGen[!,:Ramp_Dn_percentage][y]))*dfGen[!,:Cap_size][y]*EP[:vSHUT][y,t])
					end

					if (setup["Reserves"]==1)
						# Constraints on contribution to regulation and reserves
						if (dfGen[!,:Reg_Up][y] == 0 && dfGen[!,:Reg_Dn][y] == 0 && dfGen[!,:Rsv_Up][y] ==0 && dfGen[!,:Rsv_Dn][y] == 0)
							setupperbound(EP[:vREG_UP][y,t], 0.0)
							setlowerbound(EP[:vREG_UP][y,t], 0.0)
							setupperbound(EP[:vREG_DN][y,t], 0.0)
							setlowerbound(EP[:vREG_DN][y,t], 0.0)
							setupperbound(EP[:vRSV_DN][y,t], 0.0)
							setlowerbound(EP[:vRSV_DN][y,t], 0.0)
							setupperbound(EP[:vRSV_UP][y,t], 0.0)
							setlowerbound(EP[:vRSV_UP][y,t], 0.0)
						end
						@constraint(EP, EP[:vREG_UP][y,t] <= dfGen[!,:Reg_Up][y]*dfGen[!,:Cap_size][y]*EP[:vCOMMIT][y,t])
						@constraint(EP, EP[:vREG_DN][y,t] <= dfGen[!,:Reg_Dn][y]*dfGen[!,:Cap_size][y]*EP[:vCOMMIT][y,t])
						@constraint(EP, EP[:vRSV_DN][y,t] <= dfGen[!,:Rsv_Up][y]*dfGen[!,:Cap_size][y]*EP[:vCOMMIT][y,t])
						@constraint(EP, EP[:vRSV_UP][y,t] <= dfGen[!,:Rsv_Dn][y]*dfGen[!,:Cap_size][y]*EP[:vCOMMIT][y,t])
						# Minimum stable power generated per technology "y" at hour "t" Min_power
						@constraint(EP, EP[:vP][y,t] - EP[:vREG_DN][y,t] - EP[:vRSV_DN][y,t] >= dfGen[!,:Min_power][y]*dfGen[!,:Cap_size][y]*EP[:vCOMMIT][y,t])
						# Maximum power generated per technology "y" at hour "t"
						@constraint(EP, EP[:vP][y,t] + EP[:vREG_UP][y,t] + EP[:vRSV_UP][y,t] <= inputs["pP_Max"][y,t]*dfGen[!,:Cap_size][y]*EP[:vCOMMIT][y,t])
					else
						# Minimum stable power generated per technology "y" at hour "t" Min_power
						@constraint(EP, EP[:vP][y,t] >= dfGen[!,:Min_power][y]*dfGen[!,:Cap_size][y]*EP[:vCOMMIT][y,t])
						# Maximum power generated per technology "y" at hour "t"
						@constraint(EP, EP[:vP][y,t] <= inputs["pP_Max"][y,t]*dfGen[!,:Cap_size][y]*EP[:vCOMMIT][y,t])
					end
					# Thermal technologies cannot store energy
					setupperbound(EP[:vCHARGE][y,t], 0.0)
					setlowerbound(EP[:vCHARGE][y,t], 0.0)

					# Minimum up and down times
					# up time
					if  1 < h < (dfGen[!,:Up_time][y]+1)
						# Constraint wraps around to first hours of time series
						@constraint(EP, EP[:vCOMMIT][y,t] >= sum(EP[:vSTART][y,e] for e=(t-(h-1):t))+sum(EP[:vSTART][y,e] for e=(tw_max-(dfGen[!,:Up_time][y]-h)):tw_max))
					elseif h == 1
						# Constraint looks back over last n hours, where n = inputs["pDMS_Time"][y]
						@constraint(EP, EP[:vCOMMIT][y,t]  >= EP[:vSTART][y,t] + sum(EP[:vSTART][y,e] for e=(tw_max-(dfGen[!,:Up_time][y] -1 )):tw_max))
					else
						# Constraint looks back over last n hours, where n = inputs["pDMS_Time"][y]
						@constraint(EP, EP[:vCOMMIT][y,t]  >= sum(EP[:vSTART][y,e] for e=(t-dfGen[!,:Up_time][y]):t))
					end
					# down time
					if  1 < h < (dfGen[!,:Down_time][y]+1)
						# Constraint wraps around to first hours of time series
						@constraint(EP, ((dfGen[!,:Existing_Cap_MW][y]/dfGen[!,:Cap_size][y]) + EP[:vCAP][y] - EP[:vRETCAP][y]) - EP[:vCOMMIT][y,t] >=
						sum(EP[:vSHUT][y,e] for e=(t-(h-1):t))+sum(EP[:vSHUT][y,e] for e=(tw_max-(dfGen[!,:Down_time][y]-h)):tw_max))
					elseif h == 1
						# Constraint looks back over last n hours, where n = inputs["pDMS_Time"][y]
						@constraint(EP, ((dfGen[!,:Existing_Cap_MW][y]/dfGen[!,:Cap_size][y]) + EP[:vCAP][y] - EP[:vRETCAP][y]) - EP[:vCOMMIT][y,t]  >=
						EP[:vSHUT][y,t] + sum(EP[:vSHUT][y,e] for e=(tw_max-(dfGen[!,:Down_time][y] -1 )):tw_max))
					else
						# Constraint looks back over last n hours, where n = inputs["pDMS_Time"][y]
						@constraint(EP, ((dfGen[!,:Existing_Cap_MW][y]/dfGen[!,:Cap_size][y]) + EP[:vCAP][y] - EP[:vRETCAP][y]) - EP[:vCOMMIT][y,t]  >=
						sum(EP[:vSHUT][y,e] for e=(t-dfGen[!,:Down_time][y]):t))
					end
					# Startup and shutdown events
					if h <= 1
						@constraint(EP, EP[:vCOMMIT][y,t] == EP[:vCOMMIT][y,tw_max] + EP[:vSTART][y,t] - EP[:vSHUT][y,t] )
					else
						@constraint(EP, EP[:vCOMMIT][y,t] == EP[:vCOMMIT][y,t-1] + EP[:vSTART][y,t] - EP[:vSHUT][y,t] )
					end
				end #END "for y"
			end
			if (setup["UCommit"]==0)
				## Thermal generators WITHOUT unit commitment constraints - this includes resources with Commit=0 and Commit =1
				for y in dfGen[(dfGen[!,:THERM].==1),:][!,:R_ID]
					# Maximum ramp up and down between consecutive hours
					if h <= 1
						# Links last time step with first time step, ensuring position in hour 1 is within eligible ramp of final hour position
						# NOTE: We should make wrap-around a configurable option
						@constraint(EP, EP[:vP][y,t]-EP[:vP][y,tw_max] <= dfGen[!,:Ramp_Up_percentage][y]*(dfGen[!,:Existing_Cap_MW][y] + EP[:vCAP][y] - EP[:vRETCAP][y]))
						@constraint(EP, EP[:vP][y,tw_max] - EP[:vP][y,t] <= dfGen[!,:Ramp_Dn_percentage][y]*(dfGen[!,:Existing_Cap_MW][y] + EP[:vCAP][y] - EP[:vRETCAP][y]))
					else
						@constraint(EP, EP[:vP][y,t]-EP[:vP][y,t-1] <= dfGen[!,:Ramp_Up_percentage][y]*(dfGen[!,:Existing_Cap_MW][y] + EP[:vCAP][y] - EP[:vRETCAP][y]))
						@constraint(EP, EP[:vP][y,t-1] - EP[:vP][y,t] <= dfGen[!,:Ramp_Dn_percentage][y]*(dfGen[!,:Existing_Cap_MW][y] + EP[:vCAP][y] - EP[:vRETCAP][y]))
					end
					if (setup["Reserves"]==1)
						if (dfGen[!,:Reg_Up][y] == 0 && dfGen[!,:Reg_Dn][y] == 0 && dfGen[!,:Rsv_Up][y] ==0 && dfGen[!,:Rsv_Dn][y] == 0)
							setupperbound(EP[:vREG_UP][y,t], 0.0)
							setlowerbound(EP[:vREG_UP][y,t], 0.0)
							setupperbound(EP[:vREG_DN][y,t], 0.0)
							setlowerbound(EP[:vREG_DN][y,t], 0.0)
							setupperbound(EP[:vRSV_DN][y,t], 0.0)
							setlowerbound(EP[:vRSV_DN][y,t], 0.0)
							setupperbound(EP[:vRSV_UP][y,t], 0.0)
							setlowerbound(EP[:vRSV_UP][y,t], 0.0)
						end
						# Constraints on contribution to regulation and reserves
						# Thermal units without commitment constraints assumed to be quick-start units, so can contribute a specified fraction of their
						# total installed capacity to reserves
						@constraint(EP, EP[:vREG_UP][y,t] <= dfGen[!,:Reg_Up][y]*(dfGen[!,:Existing_Cap_MW][y] + EP[:vCAP][y] - EP[:vRETCAP][y]))
						@constraint(EP, EP[:vREG_DN][y,t] <= dfGen[!,:Reg_Dn][y]*(dfGen[!,:Existing_Cap_MW][y] + EP[:vCAP][y] - EP[:vRETCAP][y]))
						@constraint(EP, EP[:vRSV_DN][y,t] <= dfGen[!,:Rsv_Up][y]*(dfGen[!,:Existing_Cap_MW][y] + EP[:vCAP][y] - EP[:vRETCAP][y]))
						@constraint(EP, EP[:vRSV_UP][y,t] <= dfGen[!,:Rsv_Dn][y]*(dfGen[!,:Existing_Cap_MW][y] + EP[:vCAP][y] - EP[:vRETCAP][y]))
						# Minimum stable power generated per technology "y" at hour "t" Min_power
						@constraint(EP, EP[:vP][y,t] - EP[:vREG_DN][y,t] - EP[:vRSV_DN][y,t] >= dfGen[!,:Min_power][y]*(dfGen[!,:Existing_Cap_MW][y] + EP[:vCAP][y] - EP[:vRETCAP][y]))
						# Maximum power generated per technology "y" at hour "t"
						@constraint(EP, EP[:vP][y,t] + EP[:vREG_UP][y,t] + EP[:vRSV_UP][y,t] <= inputs["pP_Max"][y,t]*(dfGen[!,:Existing_Cap_MW][y] + EP[:vCAP][y] - EP[:vRETCAP][y]))
					else
						# Minimum stable power generated per technology "y" at hour "t" Min_power
						@constraint(EP, EP[:vP][y,t] >= dfGen[!,:Min_power][y]*(dfGen[!,:Existing_Cap_MW][y] + EP[:vCAP][y] - EP[:vRETCAP][y]))
						# Maximum power generated per technology "y" at hour "t"
						@constraint(EP, EP[:vP][y,t] <= inputs["pP_Max"][y,t]*(dfGen[!,:Existing_Cap_MW][y] + EP[:vCAP][y] - EP[:vRETCAP][y]))
					end
					# Thermal technologies cannot store energy
					setupperbound(EP[:vCHARGE][y,t], 0.0)
					setlowerbound(EP[:vCHARGE][y,t], 0.0)
				end #END "for y"
			elseif (setup["UCommit"]>=1)
				## Thermal generators WHITOUT unit commitment constraints
					for y in dfGen[(dfGen[!,:THERM].==1) .&  (dfGen[!,:Commit].==0),:][!,:R_ID]
						# Maximum ramp up and down between consecutive hours
						if h <= 1
							# Links last time step with first time step, ensuring position in hour 1 is within eligible ramp of final hour position
							# NOTE: We should make wrap-around a configurable option
							@constraint(EP, EP[:vP][y,t]-EP[:vP][y,tw_max] <= dfGen[!,:Ramp_Up_percentage][y]*(dfGen[!,:Existing_Cap_MW][y] + EP[:vCAP][y] - EP[:vRETCAP][y]))
							@constraint(EP, EP[:vP][y,tw_max] - EP[:vP][y,t] <= dfGen[!,:Ramp_Dn_percentage][y]*(dfGen[!,:Existing_Cap_MW][y] + EP[:vCAP][y] - EP[:vRETCAP][y]))
						else
							@constraint(EP, EP[:vP][y,t]-EP[:vP][y,t-1] <= dfGen[!,:Ramp_Up_percentage][y]*(dfGen[!,:Existing_Cap_MW][y] + EP[:vCAP][y] - EP[:vRETCAP][y]))
							@constraint(EP, EP[:vP][y,t-1] - EP[:vP][y,t] <= dfGen[!,:Ramp_Dn_percentage][y]*(dfGen[!,:Existing_Cap_MW][y] + EP[:vCAP][y] - EP[:vRETCAP][y]))
						end
						if (setup["Reserves"]==1)
							if (dfGen[!,:Reg_Up][y] == 0 && dfGen[!,:Reg_Dn][y] == 0 && dfGen[!,:Rsv_Up][y] ==0 && dfGen[!,:Rsv_Dn][y] == 0)
								setupperbound(EP[:vREG_UP][y,t], 0.0)
								setlowerbound(EP[:vREG_UP][y,t], 0.0)
								setupperbound(EP[:vREG_DN][y,t], 0.0)
								setlowerbound(EP[:vREG_DN][y,t], 0.0)
								setupperbound(EP[:vRSV_DN][y,t], 0.0)
								setlowerbound(EP[:vRSV_DN][y,t], 0.0)
								setupperbound(EP[:vRSV_UP][y,t], 0.0)
								setlowerbound(EP[:vRSV_UP][y,t], 0.0)
							else
								# Constraints on contribution to regulation and reserves
								# Thermal units without commitment constraints assumed to be quick-start units, so can contribute a specified fraction of their
								# total installed capacity to reserves
								@constraint(EP, EP[:vREG_UP][y,t] <= dfGen[!,:Reg_Up][y]*(dfGen[!,:Existing_Cap_MW][y] + EP[:vCAP][y] - EP[:vRETCAP][y]))
								@constraint(EP, EP[:vREG_DN][y,t] <= dfGen[!,:Reg_Dn][y]*(dfGen[!,:Existing_Cap_MW][y] + EP[:vCAP][y] - EP[:vRETCAP][y]))
								@constraint(EP, EP[:vRSV_DN][y,t] <= dfGen[!,:Rsv_Up][y]*(dfGen[!,:Existing_Cap_MW][y] + EP[:vCAP][y] - EP[:vRETCAP][y]))
								@constraint(EP, EP[:vRSV_UP][y,t] <= dfGen[!,:Rsv_Dn][y]*(dfGen[!,:Existing_Cap_MW][y] + EP[:vCAP][y] - EP[:vRETCAP][y]))
								# Minimum stable power generated per technology "y" at hour "t" Min_power
								@constraint(EP, EP[:vP][y,t] - EP[:vREG_DN][y,t] - EP[:vRSV_DN][y,t] >= dfGen[!,:Min_power][y]*(dfGen[!,:Existing_Cap_MW][y] + EP[:vCAP][y] - EP[:vRETCAP][y]))
								# Maximum power generated per technology "y" at hour "t"
								@constraint(EP, EP[:vP][y,t] + EP[:vREG_UP][y,t] + EP[:vRSV_UP][y,t] <= inputs["pP_Max"][y,t]*(dfGen[!,:Existing_Cap_MW][y] + EP[:vCAP][y] - EP[:vRETCAP][y]))
							end
						else
							# Minimum stable power generated per technology "y" at hour "t" Min_power
							@constraint(EP, EP[:vP][y,t] >= dfGen[!,:Min_power][y]*(dfGen[!,:Existing_Cap_MW][y] + EP[:vCAP][y] - EP[:vRETCAP][y]))
							# Maximum power generated per technology "y" at hour "t"
							@constraint(EP, EP[:vP][y,t] <= inputs["pP_Max"][y,t]*(dfGen[!,:Existing_Cap_MW][y] + EP[:vCAP][y] - EP[:vRETCAP][y]))
						end
						# Thermal technologies cannot store energy
						setupperbound(EP[:vCHARGE][y,t], 0.0)
						setlowerbound(EP[:vCHARGE][y,t], 0.0)
					end #END "for y"
			end #END configuration
		end
	end

	### End Constraints ###

	dModuleArgs["dExpressions"] = dExpressions
	dModuleArgs["dObjective"] = dObjective
	dModuleArgs["dPowerBalance"] = dPowerBalance

	return EP, dModuleArgs
end

function hybrid(EP::Model, dModuleArgs::Dict)

	println("Hybrid Resources Module")

	inputs = dModuleArgs["inputs"]
	setup = dModuleArgs["setup"]
	dExpressions = dModuleArgs["dExpressions"]
	dObjective = dModuleArgs["dObjective"]
	dPowerBalance = dModuleArgs["dPowerBalance"]

	dfGen = inputs["dfGen"]

	T = inputs["T"]     # Number of time steps (hours)
	Z = inputs["Z"]     # Number of zones
	W = inputs["W"]     # Number of subperiods

	Tw = div.(T,W) #total number of hours per subperiod

	### Variables ###
	## Decision variables hybrid solar dc-coupled storage behind the inverter
	# Power from PV to inverter
	@variable(EP, vPPVI[y=(dfGen[(dfGen[!,:DISP].==2) ,:][!,:R_ID]),t=1:T] >= 0)
	### End Variables ###

	### Expressions ###

	## Objective Function Expressions ##
	## End Objective Function Expressions ##

	## Power Balance Expressions ##
	## End Power Balance Expressions ##

	### End Expressions ###

	### Constratints ###

	for w in 1:W
		for h in 1:Tw
			t = Tw*(w-1)+h
			tw_min = Tw*(w-1)+1
			tw_max = Tw*(w-1)+Tw
			## Hybrid VRE and Storage behind the inverter
			## NOTE: Changed numerical flag for hybrid VRE+storage behind the inverter to :STOR==4
			##  Currently system assumes integrated system can only provide reserves based on PV generation in each hour and not storage capacity of DC system
			## Future work should revisit this formulation to account for potential for dc-coupled storage to enhance reserve provision from this system

			for z in 1:Z
				for y in (dfGen[(dfGen[!,:DISP].==2) .&  (dfGen[!,:zone].==z),:][!,:R_ID])
					# Dispatchable technologies cannot store energy
					setupperbound(EP[:vCHARGE][y,t], 0.0)
					setlowerbound(EP[:vCHARGE][y,t], 0.0)
				   for x in (dfGen[(dfGen[!,:STOR].==4) .&  (dfGen[!,:zone].==z),:][!,:R_ID]) # Storage dc-coupled with PV
					   if (setup["Reserves"]==1)
						   if (dfGen[!,:Reg_Up][y] == 0 && dfGen[!,:Reg_Dn][y] == 0 && dfGen[!,:Rsv_Up][y] ==0 && dfGen[!,:Rsv_Dn][y] == 0)
							   setupperbound(EP[:vREG_UP][y,t], 0.0)
							   setlowerbound(EP[:vREG_UP][y,t], 0.0)
							   setupperbound(EP[:vREG_DN][y,t], 0.0)
							   setlowerbound(EP[:vREG_DN][y,t], 0.0)
							   setupperbound(EP[:vRSV_DN][y,t], 0.0)
							   setlowerbound(EP[:vRSV_DN][y,t], 0.0)
							   setupperbound(EP[:vRSV_UP][y,t], 0.0)
							   setlowerbound(EP[:vRSV_UP][y,t], 0.0)
						   end
						   # # Constraints on contribution to regulation and reserves
						   # # For VRE, reserve contributions must be less than the specified percentage of the capacity factor for the hour
						   @constraint(EP, EP[:vREG_UP][y,t] <= dfGen[!,:Reg_Up][y]*inputs["pP_Max"][y,t]*(dfGen[!,:Existing_Cap_MW][y] + EP[:vCAP][y] - EP[:vRETCAP][y]))
						   @constraint(EP, EP[:vREG_DN][y,t] <= dfGen[!,:Reg_Dn][y]*inputs["pP_Max"][y,t]*(dfGen[!,:Existing_Cap_MW][y] + EP[:vCAP][y] - EP[:vRETCAP][y]))
						   @constraint(EP, EP[:vRSV_DN][y,t] <= dfGen[!,:Rsv_Up][y]*inputs["pP_Max"][y,t]*(dfGen[!,:Existing_Cap_MW][y] + EP[:vCAP][y] - EP[:vRETCAP][y]))
						   @constraint(EP, EP[:vRSV_UP][y,t] <= dfGen[!,:Rsv_Dn][y]*inputs["pP_Max"][y,t]*(dfGen[!,:Existing_Cap_MW][y] + EP[:vCAP][y] - EP[:vRETCAP][y]))
						   # Power generated and reserve contributions down per hour must be greater than zero
						   @constraint(EP, EP[:vPPVI][y,t] + EP[:vCHARGE][x,t] - EP[:vREG_DN][y,t] - EP[:vRSV_DN][y,t] >= 0)
						   # Power flow to inverter + storage charging and reserve contributions up per hour by renewable generators must be less than
						   # hourly capacity factor. Note: inequality constraint allows curtailment of output
						   # below maximum level.
						   @constraint(EP, EP[:vPPVI][y,t] + EP[:vCHARGE][x,t] + EP[:vREG_UP][y,t] + EP[:vREG_DN][y,t] <= inputs["pP_Max"][y,t]*(dfGen[!,:Existing_Cap_MW][y] + EP[:vCAP][y] - EP[:vRETCAP][y]))
					   else
						   # Maximum power generated per hour by renewable generators must be less than
						   # hourly capacity factor. Note: inequality constraint allows curtailment of output
						   # below maximum level.
						   # Sum of power flow from PV to inverter and storage charging cannot exceed availability times installed DC capacity
						   @constraint(EP, EP[:vPPVI][y,t] + EP[:vCHARGE][x,t] <= inputs["pP_Max"][y,t]*(dfGen[!,:Existing_Cap_MW][y] + EP[:vCAP][y] - EP[:vRETCAP][y]))
					   end

					   # AC power output from inverter = (DC Power flow from  PV + DC power from storage) * Inverter efficiency
					   @constraint(EP, dfGen[!,:EtaInverter][y]*(EP[:vPPVI][y,t] + EP[:vP][x,t]) == EP[:vP][y,t] )

					   # AC power capacity of the PV+ storage system is limited by the inverter loading ratio times DC power capacity of PV system
					   @constraint(EP, EP[:vP][y,t] <= (dfGen[!,:Existing_Cap_MW][y] + EP[:vCAP][y] - EP[:vRETCAP][y])/ dfGen[!,:InveterLoadRatio][y])

					   # DC coupled storage constraints to be modeled
					   # State of charge constraint (equals previous state + charge - discharge)
					   if h <= 1
						   # Links last time step with first time step, ensuring position in hour 1 is within eligible change from final hour position
						   #  Wrap-up Tolerance
						   # @constraint(EP, EP[:vS][x,t] >= (1 - setup["pTolerance"])*EP[:vS][x,tw_max]-(1/dfGen[!,:Eff_down][x]*EP[:vP][x,t])+(dfGen[!,:Eff_up][x]*EP[:vCHARGE][x,t])-(dfGen[!,:Self_disch][x]*EP[:vS][x,tw_max])) #energy stored for the next hour
						   # @constraint(EP, EP[:vS][x,t] <= (1 + setup["pTolerance"])*EP[:vS][x,tw_max]-(1/dfGen[!,:Eff_down][x]*EP[:vP][x,t])+(dfGen[!,:Eff_up][x]*EP[:vCHARGE][x,t])-(dfGen[!,:Self_disch][x]*EP[:vS][x,tw_max])) #energy stored for the next hour
						   @constraint(EP, EP[:vS][x,t] == (1 - setup["pTolerance"])*EP[:vS][x,tw_max]-(1/dfGen[!,:Eff_down][x]*EP[:vP][x,t])+(dfGen[!,:Eff_up][x]*EP[:vCHARGE][x,t])-(dfGen[!,:Self_disch][x]*EP[:vS][x,tw_max])) #energy stored for the next hour
					   else
						   @constraint(EP, EP[:vS][x,t] == EP[:vS][x,t-1]-(1/dfGen[!,:Eff_down][x]*EP[:vP][x,t])+(dfGen[!,:Eff_up][x]*EP[:vCHARGE][x,t])-(dfGen[!,:Self_disch][x]*EP[:vS][x,t]))  #energy stored for the next hour
					   end #END if
					   # Maximum energy stored must be less than energy capacity
					   @constraint(EP, EP[:vS][x,t] <= (1/dfGen[!,:Ratio_power_to_energy][x])*(dfGen[!,:Existing_Cap_MW][x] + EP[:vCAP][x] - EP[:vRETCAP][x]))
					   # Maximum charging rate must be less than power rating OR available storage capacity, whichever is less
					   # @constraint(EP, EP[:vCHARGE][x,t] <= ((1/dfGen[!,:Eff_up][x])*(dfGen[!,:Existing_Cap_MW][x] + EP[:vCAP][x] - EP[:vRETCAP][x])))
					   # @constraint(EP, EP[:vCHARGE][x,t] <= ((1/dfGen[!,:Ratio_power_to_energy][x])*(dfGen[!,:Existing_Cap_MW][x] + EP[:vCAP][x] - EP[:vRETCAP][x])-EP[:vS][x,t]))
					   @constraint(EP, EP[:vCHARGE][x,t] <= (dfGen[!,:Existing_Cap_MW][x] + EP[:vCAP][x] - EP[:vRETCAP][x]))
					   # Maximum discharging rate must be less than power rating OR available stored energy, whichever is less
					   # @constraint(EP, EP[:vP][x,t] <= ((dfGen[!,:Eff_down][x])*(dfGen[!,:Existing_Cap_MW][x] + EP[:vCAP][x] - EP[:vRETCAP][x])))
					   # @constraint(EP, EP[:vP][x,t] <= (EP[:vS][x,t]))
					   @constraint(EP, EP[:vP][x,t] <= (dfGen[!,:Existing_Cap_MW][x] + EP[:vCAP][x] - EP[:vRETCAP][x]))

					   # DC coupled storage cannot contribute operating reserves to the system
					   # NOTE: this needs to be revisited in the future
					   setupperbound(EP[:vREG_UP][x,t], 0.0)
					   setlowerbound(EP[:vREG_UP][x,t], 0.0)
					   setupperbound(EP[:vREG_DN][x,t], 0.0)
					   setlowerbound(EP[:vREG_DN][x,t], 0.0)
					   setupperbound(EP[:vRSV_DN][x,t], 0.0)
					   setlowerbound(EP[:vRSV_DN][x,t], 0.0)
					   setupperbound(EP[:vRSV_UP][x,t], 0.0)
					   setlowerbound(EP[:vRSV_UP][x,t], 0.0)

				   end #END "x" dim for dc coupled storage
			   end #END "y" dim for coupled PV system
		   end #END "for z"
		end
	end

	### End Constraints ###

	dModuleArgs["dExpressions"] = dExpressions
	dModuleArgs["dObjective"] = dObjective
	dModuleArgs["dPowerBalance"] = dPowerBalance

	return EP, dModuleArgs

end

# End Technologies Functions

# Begin Policies Functions

function trading(EP::Model, dModuleArgs::Dict)

	println("Hybrid Policies Module")

	inputs = dModuleArgs["inputs"]
	setup = dModuleArgs["setup"]
	dExpressions = dModuleArgs["dExpressions"]
	dObjective = dModuleArgs["dObjective"]
	dPowerBalance = dModuleArgs["dPowerBalance"]

	dfGen = inputs["dfGen"]

	G = inputs["G"]     # Number of resources (generators, storage, DR, and DERs)
	T = inputs["T"]     # Number of time steps (hours)
	Z = inputs["Z"]     # Number of zones
	W = inputs["W"]     # Number of subperiods

	Tw = div.(T,W) #total number of hours per subperiod

	### Variables ###

	# Auxiliary variables to track inports and exports into/out of each zone -
	# used for tracking emissions  changes due to exports and imports
	# Need to figure out if there is a particular setup case when these variables are activated
	@variable(EP, vEXPORTS[z=1:Z,t=1:T] >= 0)
	@variable(EP, vIMPORTS[z=1:Z,t=1:T] >= 0)

	### End Variables ###

	### Expressions ###

	## Objective Function Expressions ##

	# special cases for variable cost expression
	# Set variable costs of plants not within the region of interest, but supplying power to be zero
	if setup["Trading"]==1
		for y in (dfGen[(dfGen[!,:TRADER].=="E") | (dfGen[!,:TRADER].=="I"),:][!,:R_ID])
			for t in 1:T
				dExpressions["eCVar_out"][y,t] = 0
			end
		end
	end

	@expression(EP, eExportRev[y=1:G], 0)
	@expression(EP, eImportCost[y=1:G], 0)

	if setup["Trading"]==1
		for y in (dfGen[(dfGen[!,:TRADER].=="E"),:][!,:R_ID])
			eExportRev[y] = sum(inputs["omega"][t]* dfGen[!,:Var_OM_cost_per_MWh][y] * EP[:vP][y,t] for t in 1:T)
		end
		for y in (dfGen[(dfGen[!,:TRADER].=="I"),:][!,:R_ID])
			eImportCost[y] = sum(inputs["omega"][t]* dfGen[!,:Var_OM_cost_per_MWh][y] * EP[:vP][y,t] for t in 1:T)
		end
	end
	@expression(EP, eTotalexportRev, sum(eExportRev[y] for  y in 1:G ))
	@expression(EP, eTotalImportCost, sum(eImportCost[y] for  y in 1:G ))

	dObjective["eTotalexportRev"] = eTotalexportRev
	dObjective["eTotalImportCost"] = eTotalImportCost

	## End Objective Function Expressions ##

	## Power Balance Expressions ##
	if (setup["Trading"]==1) && (setup["HeatMarket"]!=1)
		@expression(EP, ePowerBalanceTradeI[t=1:T, z=1:Z],
			-sum(vP[y,t] for y=dfGen[(dfGen[!,:TRADER].=="I") .&  (dfGen[!,:zone].==z),:][!,:R_ID]))
		@expression(EP, ePowerBalanceTradeE[t=1:T, z=1:Z],
			sum(vP[y,t] for y=dfGen[(dfGen[!,:TRADER].=="E") .&  (dfGen[!,:zone].==z),:][!,:R_ID]))

		dPowerBalance["ePowerBalance"] = dPowerBalance["ePowerBalance"] + ePowerBalanceTradeI
		dPowerBalance["ePowerBalance"] = dPowerBalance["ePowerBalance"] + ePowerBalanceTradeE
	end
	## End Power Balance Expressions ##

	### End Expressions ###

	### Constratints ###

	# Define auxiliary variables to track imports and exports of electricity by zone
	# Not clear about the logic behind this constraint - vImports appears to be a free variable
	@constraint(EP, cImportsByZone[z=1:Z, t=1:T], vIMPORTS[z,t] >= -dExpressions["eNet_Export_Flows"][z,t] )
	@constraint(EP, cExportsByZone[z=1:Z, t=1:T], vEXPORTS[z,t] >= dExpressions["eNet_Export_Flows"][z,t] )

	for w in 1:W
		for h in 1:Tw
			t = Tw*(w-1)+h

			## Generator and resource constraints...
			## Integer Unit Commitment configuration for operational constraints

			## Generators in zones outside the region of interest
			## Simplified representation - hourly capacity limits, no other costs considered
			if setup["Trading"]==1
				for y in dfGen[(dfGen[!,:TRADER].== "I") .|(dfGen[!,:TRADER].== "E") ,:][!,:R_ID]
					# Maximum power generated per hour by  generators must be less than
					# hourly availability.
					@constraint(EP, EP[:vP][y,t] <= inputs["pP_Max"][y,t]*(dfGen[!,:Existing_Cap_MW][y] + EP[:vCAP][y] - EP[:vRETCAP][y]))
					# Dispatchable technologies cannot store energy
					setupperbound(EP[:vCHARGE][y,t], 0.0)
					setlowerbound(EP[:vCHARGE][y,t], 0.0)
					# No constraint on hourly ramp rate
				end #END "for y"
			end

			## Constraints applied to all AC connected energy storage technologies

		end # END for h in 1:Tw for operational constraints
	end # END for t in 1:W for operational constraints

	### End Constraints ###

	dModuleArgs["dExpressions"] = dExpressions
	dModuleArgs["dObjective"] = dObjective
	dModuleArgs["dPowerBalance"] = dPowerBalance

	return EP, dModuleArgs
end

function co2(EP::Model, dModuleArgs::Dict)

	println("C02 Policies Module")

	inputs = dModuleArgs["inputs"]
	setup = dModuleArgs["setup"]
	dExpressions = dModuleArgs["dExpressions"]
	dObjective = dModuleArgs["dObjective"]
	dPowerBalance = dModuleArgs["dPowerBalance"]

	dfGen = inputs["dfGen"]

	G = inputs["G"]     # Number of resources (generators, storage, DR, and DERs)
	T = inputs["T"]     # Number of time steps (hours)
	Z = inputs["Z"]     # Number of zones

	### Variables ###
	### End Variables ###

	### Expressions ###

	## Objective Function Expressions ##
	## End Objective Function Expressions ##

	## Power Balance Expressions ##
	## End Power Balance Expressions ##

	# CO2 emissions for resources "y" during hour "t" [tons]
	if (setup["HeatMarket"]==0)
		if (setup["UCommit"]>=1)
			# CO2 emissions for resources "y" during hour "t" [tons]
			@expression(EP, eEmissionsByPlant[y=1:G,t=1:T],(dfGen[!,:CO2_per_MWh][y]*(EP[:vP][y,t]+EP[:vCHARGE][y,t]) + dfGen[!,:CO2_per_Start][y]*EP[:vSTART][y,t]))
		else
			@expression(EP, eEmissionsByPlant[y=1:G,t=1:T],(dfGen[!,:CO2_per_MWh][y]*(EP[:vP][y,t]+EP[:vCHARGE][y,t])))
		end

		dExpressions["eEmissionsByPlant"] = eEmissionsByPlant

	elseif (setup["HeatMarket"]==1)
		if (setup["UCommit"]>=1)
			# CO2 emissions for resources "y" during hour "t" [tons]
			@expression(EP, eEmissionsByPlant[y=1:G,t=1:T],(dfGen[!,:CO2_per_MWh][y]*(EP[:vP][y,t]+EP[:vCHARGE][y,t])+ dfGen[!,:CO2_per_Start][y]*EP[:vSTART][y,t]))
			@expression(EP, eEmissionsNG[y=(dfGen[(dfGen[!,:NACC].==1) ,:][!,:R_ID]),t=1:T],inputs["pHeat_CO2"][1]*EP[:vHEAT_NG][y,t])
			for y in (dfGen[(dfGen[!,:NACC].==1) ,:][!,:R_ID])
				eEmissionsByPlant[y,:]  = eEmissionsByPlant[y,:] + eEmissionsNG[y,:]
			end
		else
			@expression(EP, eEmissionsByPlant[y=1:G,t=1:T],(dfGen[!,:CO2_per_MWh][y]*(EP[:vP][y,t]+EP[:vCHARGE][y,t])))
			@expression(EP, eEmissionsNG[y=(dfGen[(dfGen[!,:NACC].==1) ,:][!,:R_ID]),t=1:T],inputs["pHeat_CO2"][1]*EP[:vHEAT_NG][y,t])
			for y in (dfGen[(dfGen[!,:NACC].==1) ,:][!,:R_ID])
				eEmissionsByPlant[y,:]  = eEmissionsByPlant[y,:] + eEmissionsNG[y,:]
			end
		end

		dExpressions["eEmissionsByPlant"] = eEmissionsByPlant
		dExpressions["eEmissionsNG"] = eEmissionsNG
	end

	if setup["CO2Cap"] >= 1
	# flag = 1
	# if flag == 1
		# Emissions per zone = sum of emissions from each generator plus emissions rate applied to net imports minus emissions rate applied to net exports
		@expression(EP, eEmissionsByZone[z=1:Z, t=1:T], sum(eEmissionsByPlant[y,t] for y in dfGen[(dfGen[!,:zone].==z),:][!,:R_ID])
															+ inputs["pCO2ImportsRate"][z]*EP[:vIMPORTS][z,t]
															- inputs["pCO2ExportsRate"][z]*EP[:vEXPORTS][z,t] )
	else
		# Emissions per zone = sum of emissions from each generator
		@expression(EP, eEmissionsByZone[z=1:Z, t=1:T], sum(eEmissionsByPlant[y,t] for y in dfGen[(dfGen[!,:zone].==z),:][!,:R_ID]))
	end

	dExpressions["eEmissionsByZone"] = eEmissionsByZone

	### End Expressions ###

	### Constraints ###

	if (setup["CO2Cap"]==1)
		## Emissions constrains by zone/bus
		if(inputs["pMaxCO2Rate"][1] >= 0) # Apply emissions limit (tons/MWh)
			@constraint(EP, cCO2Emissions_zonal[z=1:Z], sum(inputs["omega"][t]*eEmissionsByZone[z,t] for t=1:T) <=
				inputs["pMaxCO2Rate"][z]*(sum(inputs["omega"][t]*inputs["pD"][t,z] for t=1:T) +
				setup["StorageLosses"]*sum(dExpressions["eELOSS"][y] for y=dfGen[(dfGen[!,:zone].==z),:][!,:R_ID]) ) )
		else # Apply absolute emissions limit (tons)
			@constraint(EP, cCO2Emissions_zonal[z=1:Z], sum(inputs["omega"][t]*eEmissionsByZone[z,t] for t=1:T) <=
				inputs["pMaxCO2"][z])
		end
	elseif (setup["CO2Cap"]==2)
		## Emissions constrains for entire system
		if(inputs["pMaxCO2Rate"][1] >= 0) # Apply emissions limit (tons/MWh)
			@constraint(EP, cCO2Emissions_systemwide, sum(inputs["omega"][t]*eEmissionsByZone[z,t] for z=1:Z, t=1:T) <=
				sum(inputs["pMaxCO2Rate"][z] * sum(inputs["omega"][t]*inputs["pD"][t,z] for t=1:T) for z=1:Z) +
				sum(inputs["pMaxCO2Rate"][z] * setup["StorageLosses"]*sum(dExpressions["eELOSS"][y] for y=dfGen[(dfGen[!,:zone].==z),:][!,:R_ID]) for z=1:Z) )
		else # Apply absolute emissions limit (tons)
			@constraint(EP, cCO2Emissions_systemwide, sum(inputs["omega"][t]*eEmissionsByZone[z,t] for z=1:Z, t=1:T) <=
				sum(inputs["pMaxCO2"][z] for z=1:Z) )
		end
	end #END CO2 emissions contraints

	@expression(EP, eCO2Emissions_power,  sum(inputs["omega"][t]*eEmissionsByZone[z,t] for z=1:Z, t=1:T))
	dExpressions["eCO2Emissions_power"] = eCO2Emissions_power


	# NOTE: The constraint above using was change to allow for between zone emissions 'trading'

	### End Constraints ###
	dModuleArgs["dExpressions"] = dExpressions
	dModuleArgs["dObjective"] = dObjective
	dModuleArgs["dPowerBalance"] = dPowerBalance

	return EP, dModuleArgs

end

function rps(EP::Model, dModuleArgs::Dict)

	println("RPS Policies Module")

	inputs = dModuleArgs["inputs"]
	setup = dModuleArgs["setup"]
	dExpressions = dModuleArgs["dExpressions"]

	dfGen = inputs["dfGen"]

	T = inputs["T"]     # Number of time steps (hours)
	Z = inputs["Z"]     # Number of zonests

	### Variables ###
	### End Variables ###

	### Expressions ###

	## Objective Function Expressions ##
	## End Objective Function Expressions ##

	## Power Balance Expressions ##
	## End Power Balance Expressions ##

	### End Expressions ###

	### Constratints ###

	## Renewable Portfolio Standard (minimum energy share from qualifying renewable resources) constraint
	if setup["RPS"] == 1
		# Sum across RPS eligible technologies and across time must be greater or equal than the qualifying renewables share times the total zonal demand
		@constraint(EP, cRPSShare[z=1:Z], sum(inputs["omega"][t]*EP[:vP][y,t] for y=dfGen[(dfGen[!,:RPS].==1) .&  (dfGen[!,:zone].==z),:][!,:R_ID], t=1:T)	>=
			inputs["pRPS"][z]*(sum(inputs["omega"][t]*inputs["pD"][t,z] for t=1:T) +
			setup["StorageLosses"]*sum(dExpressions["eELOSS"][y] for y=dfGen[(dfGen[!,:zone].==z),:][!,:R_ID]) ) )
	elseif setup["RPS"] == 2
		# Sum across RPS eligible technologies and across time must be greater or equal than the qualifying renewables share times the total system-wide demand
		@constraint(EP, cRPSShare, sum(inputs["omega"][t]*EP[:vP][y,t] for y=dfGen[(dfGen[!,:RPS].==1),:][!,:R_ID], t=1:T) >=
			sum(inputs["pRPS"][z]*inputs["omega"][t]*inputs["pD"][t,z] for t=1:T, z=1:Z ) +
			sum(inputs["pRPS"][z]*setup["StorageLosses"]*sum(dExpressions["eELOSS"][y] for y=dfGen[(dfGen[!,:zone].==z),:][!,:R_ID]) for z=1:Z)  )
	end

	### End Constraints ###

	dModuleArgs["dExpressions"] = dExpressions

	return EP, dModuleArgs
end

function ces(EP::Model, dModuleArgs::Dict)

	println("CES Policies Module")

	inputs = dModuleArgs["inputs"]
	setup = dModuleArgs["setup"]
	dExpressions = dModuleArgs["dExpressions"]

	dfGen = inputs["dfGen"]

	T = inputs["T"]     # Number of time steps (hours)
	Z = inputs["Z"]     # Number of zonests

	### Variables ###
	### End Variables ###

	### Expressions ###

	## Objective Function Expressions ##
	## End Objective Function Expressions ##

	## Power Balance Expressions ##
	## End Power Balance Expressions ##

	### End Expressions ###

	### Constratints ###

	## Clean Energy Standard (minimum energy share from qualifying clean energy resources) constraint
	if setup["CES"] == 1
		# Sum across CES eligible technologies and across time must be greater or equal than the qualifying renewables share times the total zonal demand
		@constraint(EP, cCESShare[z=1:Z], sum(inputs["omega"][t]*EP[:vP][y,t] for y=dfGen[(dfGen[!,:CES].==1) .&  (dfGen[!,:zone].==z),:][!,:R_ID], t=1:T)	>=
			inputs["pCES"][z]*(sum(inputs["omega"][t]*inputs["pD"][t,z] for t=1:T) +
			setup["StorageLosses"]*sum(dExpressions["eELOSS"][y] for y=dfGen[(dfGen[!,:zone].==z),:][!,:R_ID]) ) )
	elseif setup["CES"] == 2
		# Sum across CES eligible technologies and across time must be greater or equal than the qualifying renewables share times the total system-wide demand
		@constraint(EP, cCESShare, sum(inputs["omega"][t]*EP[:vP][y,t] for y=dfGen[(dfGen[!,:CES].==1),:][!,:R_ID], t=1:T) >=
			sum(inputs["pCES"][z]*inputs["omega"][t]*inputs["pD"][t,z] for t=1:T, z=1:Z ) +
			sum(inputs["pCES"][z]*setup["StorageLosses"]*sum(dExpressions["eELOSS"][y] for y=dfGen[(dfGen[!,:zone].==z),:][!,:R_ID]) for z=1:Z)  )
	end

	### End Constraints ###

	dModuleArgs["dExpressions"] = dExpressions

	return EP, dModuleArgs
end

function specificshare(EP::Model, dModuleArgs::Dict)

	println("Specific Share Policies Module")

	inputs = dModuleArgs["inputs"]
	setup = dModuleArgs["setup"]

	dfGen = inputs["dfGen"]

	T = inputs["T"]     # Number of time steps (hours)
	Z = inputs["Z"]     # Number of zonests

	### Variables ###
	### End Variables ###

	### Expressions ###

	## Objective Function Expressions ##
	## End Objective Function Expressions ##

	## Power Balance Expressions ##
	## End Power Balance Expressions ##

	### End Expressions ###

	### Constratints ###

	## technology specific energy share constraints
	if setup["Specific_Share"] == 1
		# Sum across VRE technologies and across time must be greater or equal than the VRE share times the total demand
		@constraint(EP, cMinSpecifcShare[y=dfGen[(dfGen[!,:Min_Share].==1),:][!,:R_ID]], sum(inputs["omega"][t]*EP[:vP][y,t] for  t=1:T)
			>= dfGen[!,:Min_Share_percent][y]*sum(inputs["omega"][t]*inputs["pD"][t,z] for t=1:T,z=dfGen[(dfGen[!,:R_ID].==y),:][!,:zone] ))
		@constraint(EP, cMaxSpecifcShare[y=dfGen[(dfGen[!,:Max_Share].==1),:][!,:R_ID]], sum(inputs["omega"][t]*EP[:vP][y,t] for  t=1:T)
			<= dfGen[!,:Max_Share_percent][y]*sum(inputs["omega"][t]*inputs["pD"][t,z] for t=1:T,z=dfGen[(dfGen[!,:R_ID].==y),:][!,:zone] ))

		# Enforcing total capacity of a certain type of technology across all regions to be fixed to a pre-defined value
		# Force total capacity of certain grouping of technologies (e.g. storage different locations) to sum to some pre-determined value
		if setup["ParameterScale"]==1 # Scaling input parameter values
			@constraint(EP, sum(dfGen[!,:Cap_size][y]*EP[:vCAP][y] for y =dfGen[dfGen[!,:Fixed_CAP].>0,:][!,:R_ID]) == mean(dfGen[dfGen[!,:Fixed_CAP].>0,:][!,:Fixed_CAP_MW])/1e+3)
		else
			@constraint(EP, sum(dfGen[!,:Cap_size][y]*EP[:vCAP][y] for y =dfGen[dfGen[!,:Fixed_CAP].>0,:][!,:R_ID]) == mean(dfGen[dfGen[!,:Fixed_CAP].>0,:][!,:Fixed_CAP_MW]))
		end

	end

	### End Constraints ###

  return EP, dModuleArgs
end

# End Policies Functions

end # End Module
