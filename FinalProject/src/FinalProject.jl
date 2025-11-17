module FinalProject

using SqpSolver
using Ipopt
using JuMP
using BlockArrays: BlockArray, Block

include("collision_avoid_game.jl")
export run

end
