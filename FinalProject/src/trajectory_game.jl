
# set up two player/satellite trajectory game
# NOTE: currently set to have no environment and allow satellites to move anywhere in 3D space
function setup_trajectory_game(; init = init_conds(), environment = nothing)

    # locally define cost for each player
    cost = let
        function stage_cost(x, u, t)
            x1, x2       = blocks(x)
            u1, u2       = blocks(u)
            x1nom, x2nom = blocks(init.xnoms[t])

            # cost for player 1 and player 2
            [
                (x1-x1nom)'*init.Q1*(x1-x1nom) + u1'*init.R1*u1,
                (x2-x2nom)'*init.Q2*(x2-x2nom) + u2'*init.R2*u2,
            ]
        end

        # condense into a single cost over time
        function reducer(stage_costs)
            reduce(.+, stage_costs) ./ length(stage_costs)
        end

        TimeSeparableTrajectoryGameCost(stage_cost, reducer, GeneralSumCostStructure(), 1.0)

    end

    function collision_avoid_constraint(xs, us)

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

function get_objectives_and_constraints(; init = init_conds(), game = setup_trajectory_game(), horizon = 10)

    # check that the game only has two players
    N = 2
    N == num_players(game) || error("Should have only two players")

    # construct costs
    function player_cost(τ, player_index)
        (; xs, us) = unpack_trajectory(τ; game.dynamics)
        ts = Iterators.eachindex(xs)
        Iterators.map(xs, us, ts) do x, u, t
            game.cost.stage_cost(x, u, t)[player_index]
        end |> game.cost.reducer
    end

    # cost function for each player
    fs = [(τ) -> player_cost(τ, ii) for ii in 1:N]

    # individual equality constraints
    gs = Vector{Function}(undef, N)
    for i in 1:N
        gs[i] = (τ) -> let
            (; xs, us) = unpack_trajectory(τ; game.dynamics)

            # enforce the given initial condition
            g1 = xs[1][Block(i)] - init.x[Block(i)]
            println(g1)
            println(xs[1][Block(2)] - init.x[Block(2)])
            println(xs)

            # dynamics constraints
            ts = Iterators.eachindex(xs)
            g2 = mapreduce(vcat, ts[2:end]) do t
                xs[t][Block(i)] - game.dynamics(xs[t-1], us[t-1], t-1)
            end

            vcat(g1, g2)

        end
    end

    # dummy individual inequality constraint
    # TODO: POSSIBLY ADD INDIVIDUAL FUEL CONSTRAINTS HERE
    hs = [(τ) -> [0] for _ in 1:N]

    # dummy shared equality constraints
    g̃ = (τ) -> [0]

    # TODO: IMPLEMENT SHARED INEQUALITY CONSTRAINTS!!!
    # NOTE: the probability of collision/miss-distance constraint should be shared between the two
    # NOTE: use collision_avoid_constraint(xs, us, θ) here
    h̃ = (τ) -> [0]

    return fs, gs, hs, g̃, h̃

end

# NOTE: need to pass in init.num_states, init.num_control, init.num_players into unpack_trajectories
# instead of passing in game.dynamics


# unpack trajectory
function unpack_trajectory(flat_trajectory; dynamics::ProductDynamics)
    trajs = Iterators.map(1:num_players(dynamics), blocks(flat_trajectory)) do ii, τ
        horizon = Int(length(τ) / (state_dim(dynamics, ii) + control_dim(dynamics, ii)))
        num_states = state_dim(dynamics, ii) * horizon
        X = reshape(τ[1:num_states], (state_dim(dynamics, ii), horizon))
        U = reshape(τ[(num_states + 1):end], (control_dim(dynamics, ii), horizon))

        (; xs = eachcol(X) |> collect, us = eachcol(U) |> collect)
    end

    stack_trajectories(trajs)
end

# pack trajectory
function pack_trajectory(traj)
    trajs = unstack_trajectory(traj)
    mapreduce(vcat, trajs) do τ
        vcat(reduce(vcat, τ.xs), reduce(vcat, τ.us))
    end
end

