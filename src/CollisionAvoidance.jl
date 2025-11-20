module CollisionAvoidance

using JuMP
using LinearAlgebra

export CollisionGameParams, build_model

struct CollisionGameParams
    v::Float64          # Vehicle speeds
    dt::Float64         # Time step
    N::Int              # Horizon length in time steps
    R_safe::Float64     # Safety distance goal
    x0_1::Vector{Float64} # Initial state agent 1
    x0_2::Vector{Float64} # Initial state agent 2
    w_path::Float64     # Weight for path deviation
    w_heading::Float64  # Weight for heading deviation
    w_control::Float64  # Weight for control effort
end

function build_model(params::CollisionGameParams)
    N = params.N
    dt = params.dt
    v = params.v
    R = params.R_safe
    
    model = Model()
    
    # State variables
    @variable(model, x1[1:3, 1:N+1])
    @variable(model, x2[1:3, 1:N+1])
    
    # Limit turning rate
    @variable(model, -1.0 <= u1[1:N] <= 1.0)
    @variable(model, -1.0 <= u2[1:N] <= 1.0)
    
    # Initial conditions
    # Fix initial state
    for i in 1:3
        @constraint(model, x1[i, 1] == params.x0_1[i])
        @constraint(model, x2[i, 1] == params.x0_2[i])
    end
    
    # Dynamics modeled as constraints
    for k in 1:N
        # Agent 1
        @NLconstraint(model, x1[1, k+1] == x1[1, k] + v * cos(x1[3, k]) * dt)
        @NLconstraint(model, x1[2, k+1] == x1[2, k] + v * sin(x1[3, k]) * dt)
        @constraint(model, x1[3, k+1] == x1[3, k] + u1[k] * dt)
        
        # Agent 2 
        @NLconstraint(model, x2[1, k+1] == x2[1, k] + v * cos(x2[3, k]) * dt)
        @NLconstraint(model, x2[2, k+1] == x2[2, k] + v * sin(x2[3, k]) * dt)
        @constraint(model, x2[3, k+1] == x2[3, k] + u2[k] * dt)
    end
    
    # Collision Avoidance 
    for k in 1:N+1
        @NLconstraint(model, (x1[1, k] - x2[1, k])^2 + (x1[2, k] - x2[2, k])^2 >= R^2)
    end
    
    # Objective
    # Minimize deviation from y=0, deviation from initial heading, and control effort
    target_th1 = params.x0_1[3]
    target_th2 = params.x0_2[3]

    @objective(model, Min, 
        sum(params.w_path * (x1[2, k]^2 + x2[2, k]^2) for k in 1:N+1) +
        sum(params.w_heading * ((x1[3, k] - target_th1)^2 + (x2[3, k] - target_th2)^2) for k in 1:N+1) +
        sum(params.w_control * (u1[k]^2 + u2[k]^2) for k in 1:N)
    )
    
    return model, x1, x2, u1, u2
end

end
