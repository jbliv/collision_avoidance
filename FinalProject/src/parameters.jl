
mutable struct Init
    
    # S1 pos, S1 vel, S2 pos, S2 vel
    x::AbstractVector{Float64}
    xnoms::Vector{AbstractVector{Float64}}

    # how often to replan using next horizon (must be <= horizon)
    turn_length::Int

    # number of future time steps to use for current horizon (must be <= n_sim_steps)
    horizon::Int

    # total number of time steps for entire simulation
    n_sim_steps::Int

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
    turn_length = 2
    horizon     = 2
    n_sim_steps = 10
    num_players = 2
    x           = BlockArray(zeros(12), [6, 6])
    num_control = 6

    # get nominal states
    xnom0 = BlockArray(ones(12), [6, 6])
    xnoms = get_nominal_states(xnom0, n_sim_steps)

    # Q, R weights
    Q1 = I(6)
    Q2 = I(6)
    R1 = I(3)
    R2 = I(3)

    init = Init(x,
                xnoms,
                turn_length,
                horizon,
                n_sim_steps,
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

    # initialize nominal state array
    xnoms    = Vector{AbstractVector{Float64}}(undef, horizon)
    xnoms[1] = xnom0

    for t = 2:horizon

        # TODO: INTEGRATE NOMINAL STATES HERE!!!

        # TEMPORARY
        xnoms[t] = BlockArray(ones(12), [6, 6])

    end

    # NOTE: at each time, xnom needs to be a block array with x1nom, x2nom

    return xnoms

end

# relevant parameters including probability of collision covariance matrix
function get_P(t)

    # TODO: ADD IN P_c LOGIC AS A FUNCTION OF TIME

end


