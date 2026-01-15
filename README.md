# TailCausalityTests
Code for the paper "General tests for pairwise causality in extremes" by Lisa Leimenstoll and Melanie Schienle

## External code and license

This project includes a modified version of code from:

Gnecco, Nicola, Nicolai Meinshausen, Jonas Peters, and Sebastian Engelke, 2019.
Original repository: https://github.com/nicolagnecco/causalXtreme
License: GPL-3.0

The adapted code is located in `Simulation_Study/k_simulation.R` and is distributed under
the same GPL-3.0 license, with minor modifications:
  - add second tail index and pareto distribution in simulate_data function
  - add pareto distribution in simulate_noise
  - in simulation_0 change arguments, add second tail index
  - add second tail index and pareto distribution in my_args
  - change of plot labels and adjust for two tail indices
  - add additional simualtion study for percentage of wrong causal inference 
    between two variables
