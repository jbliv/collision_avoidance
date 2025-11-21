
# set up two player/satellite trajectory game
# NOTE: currently set to have no environment and allow satellites to move anywhere in 3D space
function setup_trajectory_game(; environment = nothing)

    # locally define cost for each player
    cost = let
        function stage_cost(x, u, t, θ, xnom)
            x1, x2       = blocks(x)
            u1, u2       = blocks(u)
            x1nom, x2nom = blocks(xnom)

            # cost for player 1 and player 2
            # NOTE: currently set to only penalize fuel consumption
            [
                (x1-x1nom)'*Q1*(x1-x1nom) + u1'*R1*u1,
                (x2-x2nom)'*Q2*(x2-x2nom) + u2'*R2*u2,
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

                # nonlinear collision avoidance constraint goes here
                # NOTE: for now, they must be 10 m away from each other
                (x2[1:3] - x1[1:3])' * (x2[1:3] - x1[1:3]) - 100

            end
        end

        # double integrator for three-dimensional Keplerian orbit
        sat_dynamics = orbital_double_integrator(;
            state_bounds = (; lb = [-Inf, -Inf, -Inf, -Inf, -Inf, -Inf], ub = [Inf, Inf, Inf, Inf, Inf, Inf]),
            control_bounds = (; lb = [-Inf, -Inf, -Inf], ub = [Inf, Inf, Inf]),
        )
        # THIS NEEDS TO BE CHANGED ABOVE AND BELOW BECAUSE DYNAMICS ARE NOW DEPENDENT ON NOMINAL STATE
        dynamics = ProductDynamics([sat_dynamics for _ in 1:2])

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

    # dummy individual constraints
    # POSSIBLY UPDATE THIS!!!!!
    # NOTE: the dynamics should no longer be shared constraints but individual to each player
    # because they are linearized about the nominal trajectory

    gs = Vector{Function}(undef, N)
    for i in 1:N
        gs[i] = (τ, θ) -> let
            (; xs, us) = unpack_trajectory(τ; init)
            
            # get all of the ith player's states
            xis    = xs[:][Block(i)]
            xinoms = xnoms[:][Block(i)]
            θi     = θ[Block(i)]

            # enforce the given initial condition
            g1 = xis[1] - θi

            # dynamics constraints
            ts = Iterators.eachindex(xs)
            g2 = mapreduce(vcat, ts[2:end]) do t
                xis[t] - orbit_dynamics(xinoms[t-1], xis[t-1], us[t-1])
            end

            vcat(g1, g2)

        end

    end

    # dummy individual inequality constraint
    hs = [(τ, θ) -> [0] for _ in 1:N]

    # shared equality constraints
    # NOTE: the probability of collision/miss-distance constraint should be shared between the two
    g̃ = (τ, θ) -> let
        (; xs, us) = unpack_trajectory(τ; game.dynamics)

        # Force all players to start at the given initial condition.
        g̃1 = xs[1] - θ

        # Dynamics constraints.
        ts = Iterators.eachindex(xs)
        g̃2 = mapreduce(vcat, ts[2:end]) do t
            xs[t] - game.dynamics(xs[t - 1], us[t - 1])
        end

        vcat(g̃1, g̃2)
    end

end

# NOTE: need to pass in init.num_states, init.num_control, init.num_players into unpack_trajectories
# instead of passing in game.dynamics

