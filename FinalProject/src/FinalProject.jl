module FinalProject

using TrajectoryGamesBase: rollout, ProductDynamics, TrajectoryGame, state_bounds,
                           control_bounds, LinearDynamics, num_players, state_dim,
                           control_dim, stack_trajectories, unstack_trajectory,
                           TimeSeparableTrajectoryGameCost, GeneralSumCostStructure,
                           OpenLoopStrategy, JointStrategy, RecedingHorizonStrategy
using Symbolics: Symbolics, @variables
using ParametricMCPs: ParametricMCPs, ParametricMCP
using LinearAlgebra: norm, norm_sqr, Diagonal, I
using BlockArrays: BlockArray, Block, blocks
using PATHSolver: PATHSolver

include("parameters.jl")
export init_conds, get_P

include("dynamics.jl")
export orbital_double_integrator

include("parametric_game.jl")
export ParametricGame

include("parametric_optimization_problem.jl")
export ParametricOptimizationProblem, solve, total_dim

include("trajectory_game.jl")
export setup_trajectory_game, build_parametric_game

include("collision_avoid_game.jl")
export run, tests

end
