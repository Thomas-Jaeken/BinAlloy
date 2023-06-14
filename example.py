import BinAlloy

##############
# Set simulation parameters
##############
T_start = 200
T_step_size = 10
N_steps = 10
Metro_step = 10000
Dimension = 10

##############
# Calculate the order parameter as function of temperature
# This creates a dictionary holding the results of the monte carlo simulation
##############
results = BinAlloy.run_monte_carlo(
    Dimension, Metro_step, N_steps, T_start, T_step_size)

##############
# Unpack the dict to illustrate what is contained.
##############

# List of temperatures for which the Monte Carlo Algorithm was run.
TempRange = results["temperature range"]

# Average of the Observables across the datasets generated at each temperature
OrderParam = results["order parameterr"]
Energy = results["energy"]
ShortOrder = results["short-range order"]

# Variance of the Observables across the datasets generated at each temperature
OrderVariance = results["order variance"]
EnergyVariance = results["energy variance"]
ShortOrdervariance = results["short order variancne"]

#  statistical inefficiencies as estimate of the error
OrderInefficiency = results["order inefficiency"]
EnergyInefficiency = results["energy inefficiency"]
ShortOrderInefficiency = results["short order inefficiency"]
