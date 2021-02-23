# Store the path of the current working directory

cd(dirname(@__FILE__))

# ENV["CPLEX_STUDIO_BINARIES"] = "/home/gridsan/gnhe/opt/ibm/ILOG/CPLEX_Studio129/cplex/bin/x86-64_linux"

working_path = pwd()

# GenX path
module_path = pwd()

# Load GenX modules
push!(LOAD_PATH, module_path)
println("Loading packages")


using DataFrames
using CSV
using JuMP
using MathProgBase #for fix_integers
using StatsBase
using LinearAlgebra
using Random
using JLD2, FileIO
using Dates

using Gurobi
using CPLEX

# using Revise

using GenX_Modular
using Hydrogen_Modular

# Import scenarios of interest and store in dataframe
Scenario_Param = CSV.read("Case_list_test.csv")

# Slurm command - if running as a job array on supercloud
# idx = ENV["SLURM_ARRAY_TASK_ID"]
# If running locally
idx = string(size(Scenario_Param,1))
idx = "1"
job_id = "1"
# Loop over each scenario and generate model outputs
c = parse(Int,idx)
# c = 5
println("Case:",c)
# for c = 1


datetimenow = Dates.format(now(),"mmddyy_HHMMSS")

# Load inputs
push!(LOAD_PATH, working_path)

# Set GenX parameters
mysetup=Dict() # Development Note: replace this soon with inputs from a Setup.csv file
mysetup["Trans_Loss_Segments"]=1 # Number of segments to use in piecewise linear approximation of losses; 1 = linear, >2 = piecewise quadratic
mysetup["Distr_Loss_Segments"]=0 # Number of segments to use in piecewise linear approximation of losses; 1 = linear, >2 = piecewise quadratic
mysetup["Distr_Margin_LV_Segments"]=0 # Number of segments to use in piecewise linear approximation of quadratic term in LV distribution margin equation
mysetup["Specific_Share"] =0 # 1 = active otherwise no active
mysetup["VRE_Share"]=0 #minimum VRE penetration; 0 = no active; 1 = active per node; 2 = active as system
mysetup["HeatMarket"]=0 # 1 = active otherwise no active
mysetup["Externalities"]=0 # 1 = active otherwise no active
mysetup["pTolerance"] = 0.001 # Wrap-up tolerance
mysetup["ParameterScale"] =1 # Flag to turn on parameter scaling wherein load, capacity and power variables defined in GW rather than MW. 1- scaling done, 0- use input as it is.
# NOTE: OUTPUTS ARE NOT SCALED BACK
mysetup["MIPGap"] = 0.01 # 1.0% MIP Gap
mysetup["MacOrWindows"] = "Mac"  # Set to either "Mac" (also works for Linux) or "Windows" to ensure use of proper file directory separator "\" or "/
mysetup["Dual_MIP"] = 0 # get dual of MIP
mysetup["Feasib_Tol"] = 1e-6 # All constraints must be satisfied to a tolerance of
mysetup["Pre_Solve"] = 2 # Controls the presolve level. A value of -1 corresponds to an automatic setting. Other options are off (0), conservative (1), or aggressive (2).

mysetup["solver_option"] = Scenario_Param[!,:solver_option][c]
mysetup["Crossover"] = Scenario_Param[!,:Crossover][c] # Determines the crossover strategy used to transform the interior solution produced by barrier into a basic solution.
mysetup["Method"] = Scenario_Param[!,:Method][c] # Algorithm used to solve continuous models or the root node of a MIP model. Options are: -1=automatic, 0=primal simplex, 1=dual simplex, 2=barrier, 3=concurrent, 4=deterministic concurrent.
mysetup["SolutionType"] = Scenario_Param[!,:SolutionType][c]
mysetup["Threads"] = 4

mysetup["NumericFocus"] = 1 #The NumericFocus parameter controls the degree to which the code attempts to detect and manage numerical issues. The default setting (0)
mysetup["Optimal_Tol"] = 1e-6 #Reduced costs must all be smaller than OptimalityTol in the improving direction in order for a model to be declared optimal.
mysetup["TimeLimit"] =2000000 # Time limit to terminate the solution algorithm, model could also terminate if it reaches MIPGap before this time
#NOTE: new setup option for time wrapping
mysetup["OperationWrapping"]=2 # 0 = Full Year First-Last Hour Year; 1 = Full Year First-Last Hour Week; 2 = Sample Weeks First-Last Hour Week;
mysetup["Trading"]=0
mysetup["StorageLosses"]=1 # 0 VRE and CO2 constraint DO NOT account for energy lost; 1 constraint DO  account for energy lost
mysetup["PrintModel"]=0 # 1 for including the model equation as an output; 0 for the model equation not being included as an output
mysetup["RPS"]=0 #minimum qualifying renewables penetration; 0 = no active; 1 = active per node; 2 = active as system
mysetup["CES"]=0 #minimum qualifying clean energy resource penetration; 0 = no active; 1 = active per node; 2 = active as system
mysetup["Fuel_Price_Var"]=0 # 1 = active otherwise no active. When active, scales fuel price over time based on Fuels_variability.csv input

# These are the parameters which will be tested...
# ...the values that they are initially set at will be overwritten
mysetup["CO2Cap"] = 0# 0 = no cap; 1 = zonal caps; 2 = overall cap;
mysetup["CO2_price_option"] = 1 # 0 = no CO2 price; 1 = CO2 price
mysetup["UCommit"]=2 # 0 = no UC, 1 integer clestering, 2 linearized clustering
mysetup["NetworkExpansion"]=1 # 1 = active otherwise no active
mysetup["Reserves"] = 0 # 0 = no reserves, 1 regulation (primary) and operating (secondary) reserves

mysetup["single_power_sector"] = 0
mysetup["single_hydrogen_sector"] = 0
mysetup["power_load_factor"] = 1

mysetup["battery_truck_option"] = 0
mysetup["H2_to_power_cost_factor"] = 1
mysetup["reliability_test"] = 0

mysetup["LDS"] = 1

# Initialize Optimizer
solver_option = mysetup["solver_option"]
if solver_option == 1
	Gurobi.GurobiSolver
	Opt = GurobiSolver(
	MIPGap = mysetup["MIPGap"],
	FeasibilityTol = mysetup["Feasib_Tol"],
	Presolve = mysetup["Pre_Solve"],
	Method = mysetup["Method"],
	Crossover = mysetup["Crossover"],
	NumericFocus = mysetup["NumericFocus"],
	# TimeLimit = mysetup["TimeLimit"],
	# 	Threads = mysetup["Threads"],
	OptimalityTol = mysetup["Optimal_Tol"]
	)
    println("Using Gurobi")
    println("Crossover:", mysetup["Crossover"])
    println("Method:", mysetup["Method"])
    # println("Threads:", mysetup["Threads"])
else
	CPLEX.Optimizer()
	Opt = CplexSolver(
	CPX_PARAM_EPGAP = mysetup["MIPGap"],
	CPX_PARAM_EPRHS = mysetup["Feasib_Tol"],
	CPXPARAM_LPMethod = mysetup["Method"],
	# CPX_PARAM_BARCROSSALG = mysetup["Crossover"],
	CPX_PARAM_SOLUTIONTYPE = mysetup["SolutionType"],
	# NumericFocus = mysetup["NumericFocus"],
	# CPX_PARAM_TILIM = mysetup["TimeLimit"],
	# 	CPX_PARAM_THREADS = mysetup["Threads"],
	# CPXPARAM_Benders_Strategy = 3,
	CPX_PARAM_EPOPT = mysetup["Optimal_Tol"]
	)
    println("Using CPLEX")
    println("SolutionType:", mysetup["SolutionType"])
    println("Method:", mysetup["Method"])
    # println("Threads:", mysetup["Threads"])
end

# Number of weeks to operate
FolderName = Scenario_Param[!,:FolderName][c]

# Define input folder to evaluate
inpath_power = string(working_path,"/Inputs_Wks_",FolderName)
inpath_H2 = string(working_path,"/Inputs_H2")
# Load inputs
println("Loading inputs")
inputs = load_power_inputs(mysetup,inpath_power)
# Unifying time horizon
mysetup["W"] = inputs["W"]
mysetup["Tw"] = inputs["H"]
mysetup["omega"] = inputs["omega"]
inputs_H2 = load_H2_inputs(mysetup,inpath_H2)

# Adjust parameters according to scenario setting
inputs, inputs_H2 = adjust_parameters(mysetup,inputs,inputs_H2,Scenario_Param, c)

# Print Settings
print_setting(mysetup,inputs,inputs_H2)


# Generate model
# creation of optimization instance and solution
# println("Generating Model")
EP, dModuleArgs = generate_Power_H2_model(mysetup,inputs,inputs_H2,Opt)


warm_start_option = inputs_H2["warm_start"]
if warm_start_option == 1
	EP, dModuleArgs = warm_start(mysetup,inputs,inputs_H2,EP,dModuleArgs)
end

fix_investment_option = inputs_H2["fix_investment_option"]
if fix_investment_option == 1
	EP, dModuleArgs = fix_planning_var(mysetup,inputs,inputs_H2,EP,dModuleArgs)
end

# Solve model
println("Solving Model")
# myresults = solve_power_system_model(EP,mysetup,myinputs,dModuleArgs)

## Solve Model
EP, dModuleArgs = solve_model(EP,dModuleArgs)

results_power = generate_Power_results(EP, mysetup, inputs, dModuleArgs)

results_H2 = generate_H2_results(EP, dModuleArgs)

#
# println("Results:")
# for key in keys(results_H2["H2PLAN"])
# 	println(key,":",round.(results_H2["H2PLAN"][key],digits = 1))
# end
# println(results_H2["H2COSTS"])

datetimenow = Dates.format(now(),"mmddyy_HHMMSS")

outpath = "$working_path/Outputs"
sep = "/"
output_folder = string(outpath,sep,"results_",job_id,"_",datetimenow)

if (isdir(output_folder)==false)
   mkdir(output_folder)
end

case_id= string("results_",job_id,"_case_",c)

case_path = string(output_folder,sep,case_id)
if (isdir(case_path)==false)
   mkdir(case_path)
end

save_H2_option = 1
if save_H2_option == 1
	println("Saving results of hydrogen sector")
	outpath_H2 = string(case_path,sep,case_id,"_H2")
	if (isdir(outpath_H2)==false)
	   mkdir(outpath_H2)
	end
	output_filename = string(outpath_H2,sep,case_id,"_H2.jld2")
	@save output_filename inputs_H2 results_H2
	save_H2_results(inputs_H2,results_H2,outpath_H2)
end

save_power_option = 1
if save_power_option == 1
	println("Saving results of power sector")
	# Creating outpath location to save
	outpath_power = string(case_path,sep,case_id,"_power")
# 	outpath_power = string(outpath,sep,"results_",job_id,"_case_",c)
	# Creating folder to write results, if it does not exist
	if (isdir(outpath_power)==false)
	   mkdir(outpath_power)
	end
	for result_key in collect(keys(results_power))
		CSV.write(string(outpath_power,sep,result_key,".csv"),results_power[result_key])
	end
	inputs_power = inputs
	output_filename = string(outpath_power,sep,case_id,"power.jld2")
	@save output_filename inputs_power results_power
end

println("------------------------------------------------------------------")
println("------------------------------------------------------------------")
println("Cost of Power Sector:")
println("------------------------------------------------------------------")
println(results_power["COSTS"])
println("------------------------------------------------------------------")
println("------------------------------------------------------------------")
println("Cost of Hydrogen Sector:")
println("------------------------------------------------------------------")
println(results_H2["H2COSTS"])
println("------------------------------------------------------------------")
println("Hydrogen Results:")
println("------------------------------------------------------------------")
for key in keys(results_H2["H2PLAN"])
	println(key,":",round.(results_H2["H2PLAN"][key],digits = 1))
end
# println(results_H2["H2COSTS"])
