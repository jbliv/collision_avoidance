
mutable struct Init
    
    # S1 pos, S1 vel, S2 pos, S2 vel
    x::Vector{Float64}
    xnom::Vector{Float64}

    # number of time steps
    horizon::Float64

    # number of players, states, and control inputs
    num_players::Float64
    num_states::Float64
    num_control::Float64

end

# initial conditions
function init_conds()

    # initial true and nominal states
    x    = zeros(12)
    xnom = zeros(12)

    init = Init(x, xnom, 2, 2, length(x), 3)

    return init

end

# relevant parameters including probability of collision covariance matrix
function get_P(t)

    # TODO: ADD IN P_c LOGIC AS A FUNCTION OF TIME

end


