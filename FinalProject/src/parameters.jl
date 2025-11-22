
mutable struct Init
    
    # S1 pos, S1 vel, S2 pos, S2 vel
    x::AbstractVector{Float64}
    xnoms::Vector{AbstractVector{Float64}}

    # number of time steps
    horizon::Int

    # number of players, states, and control inputs
    num_players::Int
    num_states::Int
    num_control::Int

    # Q, R weights
    Q1::Matrix{Float64}
    Q2::Matrix{Float64}
    R1::Matrix{Float64}
    R2::Matrix{Float64}

end

# initial conditions
function init_conds()

    # initial true and nominal states
    horizon     = 1
    num_players = 2
    x           = BlockArray(zeros(12), [6, 6])
    num_control = 6

    # get nominal states
    xnom0 = zeros(12)
    xnoms = get_nominal_states(xnom0, horizon)

    # Q, R weights
    Q1 = I(6)
    Q2 = I(6)
    R1 = I(3)
    R2 = I(3)

    init = Init(x,
                xnoms,
                horizon,
                num_players,
                length(x),
                num_control,
                Q1,
                Q2,
                R1,
                R2)

    return init

end

# pre-compute nominal states
function get_nominal_states(xnom0, horizon)

    xnoms = Vector{AbstractVector{Float64}}(undef, horizon)

    for t = 1:horizon

        # TODO: INTEGRATE NOMINAL STATES HERE!!!

        # TEMPORARY
        xnoms[t] = BlockArray(zeros(12), [6, 6])

    end

    # NOTE: at each time, xnom needs to be a block array with x1nom, x2nom

    return xnoms

end

# relevant parameters including probability of collision covariance matrix
function get_P(t)

    # TODO: ADD IN P_c LOGIC AS A FUNCTION OF TIME

end


