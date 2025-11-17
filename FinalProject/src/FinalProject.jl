module FinalProject

using SqpSolver
using Ipopt
using JuMP
using LinearAlgebra: norm, norm_sqr, Diagonal
using BlockArrays: BlockArray, Block

include("dynamics.jl")
export custom_double_integrator

include("collision_avoid_game.jl")
export run, tests

end
