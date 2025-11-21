module FinalProject

using SqpSolver
using Ipopt
using JuMP
using LinearAlgebra: norm, norm_sqr, Diagonal
using BlockArrays: BlockArray, Block
using TrajectoryGamesBase: rollout

include("parameters.jl")
export init_conds, get_P

include("dynamics.jl")
export orbital_double_integrator

include("collision_avoid_game.jl")
export run, tests

end
