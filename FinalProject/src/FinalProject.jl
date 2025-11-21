module FinalProject

using SqpSolver
using Ipopt
using JuMP
using LinearAlgebra: norm, norm_sqr, Diagonal, I
using BlockArrays: BlockArray, Block
using TrajectoryGamesBase: rollout, ProductDynamics, TrajectoryGame, state_bounds,
                           control_bounds, LinearDynamics, num_players

include("parameters.jl")
export init_conds, get_P

include("dynamics.jl")
export orbital_double_integrator

include("trajectory_game.jl")
export setup_trajectory_game, get_objectives_and_constraints

include("collision_avoid_game.jl")
export run, tests

end
