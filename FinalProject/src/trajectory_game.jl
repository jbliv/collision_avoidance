
# set up two player/satellite trajectory game
# NOTE: currently set to have no environment and allow satellites to move anywhere in 3D space
function setup_trajectory_game(; init = init_conds(), environment = nothing)

    # locally define cost for each player
    cost = let
        function stage_cost(x, u, t, θ, xnom)
            x1, x2       = blocks(x)
            u1, u2       = blocks(u)
            x1nom, x2nom = blocks(xnom)

            # cost for player 1 and player 2
            [
                (x1-x1nom)'*init.Q1*(x1-x1nom) + u1'*init.R1*u1,
                (x2-x2nom)'*init.Q2*(x2-x2nom) + u2'*init.R2*u2,
            ]

            # condense into a single cost over time
            function reducer(stage_costs)
                reduce(.+, stage_costs) ./ length(stage_costs)
            end

            TimeSeparableTrajectoryGameCost(stage_cost, reducer, GeneralSumCostStructure(), 1.0)

        end

        function collision_avoid_constraint(xs, us, θ)

            # get collision avoidance constraint for each time (as a vector)
            mapreduce(vcat, xs) do x
                x1, x2 = blocks(x)

                # TODO: UPDATE THIS WITH CORRECT COLLISION AVOIDANCE CONSTRAINT!!!
                # TODO: INCORPORATE TIME HERE --> use get_P(t) to get P_c covariance matrix
                # NOTE: if necessary, can bring this outside of this setup_trajectory_game function
                # and just define it when shared inequality constraints are defined later

                # nonlinear collision avoidance constraint goes here
                # NOTE: for now, they must be 10 m away from each other
                (x2[1:3] - x1[1:3])' * (x2[1:3] - x1[1:3]) - 100

            end
        end

        # double integrator for three-dimensional Keplerian orbit
        dynamics = ProductDynamics(
                [orbital_double_integrator(init, i;
                state_bounds = (; lb = [-Inf, -Inf, -Inf, -Inf, -Inf, -Inf], ub = [Inf, Inf, Inf, Inf, Inf, Inf]),
                control_bounds = (; lb = [-Inf, -Inf, -Inf], ub = [Inf, Inf, Inf]),
            ) for i = 1:2]
        )

        TrajectoryGame(dynamics, cost, environment, collision_avoid_constraint)

    end

end


function build_parametric_game(; game = setup_trajectory_game(), horizon = 10)

    # check that the game only has two players
    N = 2
    N == num_players(game) || error("Should have only two players")

    # construct costs
    function player_cost(τ, θ, player_index)
        (; xs, us) = unpack_trajectory(τ; game.dynamics)
        ts = Iterators.eachindex(xs)
        Iterators.map(xs, us, ts) do x, u, t
            game.cost.stage_cost(x, u, t, θ)[player_index]
        end |> game.cost.reducer
    end

    # cost function for each player
    fs = [(τ, θ) -> player_cost(τ, θ, ii) for ii in 1:N]

    # individual equality constraints
    gs = Vector{Function}(undef, N)
    for i in 1:N
        gs[i] = (τ, θ) -> let
            (; xs, us) = unpack_trajectory(τ; init)
            
            # get all of the ith player's states
            xis    = xs[:][Block(i)]
            θi     = θ[Block(i)]

            # enforce the given initial condition
            g1 = xis[1] - θi

            # dynamics constraints
            ts = Iterators.eachindex(xs)
            g2 = mapreduce(vcat, ts[2:end]) do t
                xis[t] - game.dynamics(xis[t-1], us[t-1], t-1)
            end

            vcat(g1, g2)

        end
    end

    # dummy individual inequality constraint
    # TODO: POSSIBLY ADD INDIVIDUAL FUEL CONSTRAINTS HERE
    hs = [(τ, θ) -> [0] for _ in 1:N]

    # dummy shared equality constraints
    g̃ = (τ, θ) -> [0]

    # TODO: IMPLEMENT SHARED INEQUALITY CONSTRAINTS!!!
    # NOTE: the probability of collision/miss-distance constraint should be shared between the two
    # NOTE: use collision_avoid_constraint(xs, us, θ) here
    h̃ = (τ, θ) -> [0]

    # return parametric game
    ParametricGame(; 
        objectives = fs,
        equality_constraints = gs,
        inequality_constraints = hs,
        shared_equality_constraint = g̃,
        shared_inequality_constraint = h̃,
        parameter_dimension = state_dim(game.dynamics),
        primal_dimensions = [
            horizon * (state_dim(game.dynamics, ii) + control_dim(game.dynamics, ii)) for ii in 1:N
        ],
        equality_dimensions = [1 for _ in 1:N],
        inequality_dimensions = [1 for _ in 1:N],
        shared_equality_dimension = 0,
        shared_inequality_dimension = 0, # TODO: UPDATE THIS LATER
    )

end

# NOTE: need to pass in init.num_states, init.num_control, init.num_players into unpack_trajectories
# instead of passing in game.dynamics


# can use rollout which takes in a dynamics function which takes in x, u, t, so just internally in the dynamics
# function we need to determine what the xnom is to use in the linearization

