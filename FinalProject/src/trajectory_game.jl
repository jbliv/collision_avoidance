
using TrajectoryGamesBase
using ProgressMeter: ProgressMeter

# set up two player/satellite trajectory game
# NOTE: currently set to have no environment and allow satellites to move anywhere in 3D space
function setup_trajectory_game(; init = init_conds(), environment = nothing)

    # locally define cost for each player
    cost = let
        function stage_cost(x, u, t, θ)
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

# convert a TrajectoryGame to a ParametricGame 
function build_parametric_game(; game = setup_trajectory_game(), init = init_conds())
    N = 2
    N == num_players(game) || error("Should have only two players.")

    # Construct costs.
    function player_cost(τ, θ, player_index)
        (; xs, us) = unpack_trajectory(τ; game.dynamics)
        ts = Iterators.eachindex(xs)
        Iterators.map(xs, us, ts) do x, u, t
            game.cost.discount_factor^(t - 1) * game.cost.stage_cost(x, u, t, θ)[player_index]
        end |> game.cost.reducer
    end

    # cost function for each player
    fs = [(τ, θ) -> player_cost(τ, θ, ii) for ii in 1:N]

    # individual equality constraints
    gs = Vector{Function}(undef, N)
    for i in 1:N
        gs[i] = (τ, θ) -> let
            (; xs, us) = unpack_trajectory(τ; game.dynamics)

            # enforce the given initial condition
            g1 = xs[1][Block(i)] - init.x[Block(i)]

            # dynamics constraints
            ts = Iterators.eachindex(xs)
            g2 = mapreduce(vcat, ts[2:end]) do t
                xs[t][Block(i)] - game.dynamics(xs[t-1], us[t-1], t-1)[Block(i)]
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

    ParametricGame(;
        objectives = fs,
        equality_constraints = gs,
        inequality_constraints = hs,
        shared_equality_constraint = g̃,
        shared_inequality_constraint = h̃,
        parameter_dimension = state_dim(game.dynamics),
        primal_dimensions = [
            init.horizon * (state_dim(game.dynamics, ii) + control_dim(game.dynamics, ii)) for ii in 1:N
        ],
        equality_dimensions = [init.horizon * state_dim(game.dynamics, ii) for ii in 1:N],
        inequality_dimensions = [1 for _ in 1:N],
        shared_equality_dimension = 1,
        shared_inequality_dimension = 1,
    )
end

# generate an initial guess for primal variables following a zero input sequence
function generate_initial_guess(;
    game::TrajectoryGame{<:ProductDynamics},
    parametric_game::ParametricGame,
    horizon,
    initial_state,
)
    rollout_strategy =
        map(1:num_players(game)) do ii
            (x, t) -> zeros(control_dim(game.dynamics, ii))
        end |> JointStrategy

    zero_input_trajectory =
        rollout(game.dynamics, rollout_strategy, initial_state, horizon)

    vcat(
        pack_trajectory(zero_input_trajectory),
        zeros(total_dim(parametric_game) - sum(parametric_game.primal_dimensions)),
    )
end

# solve a parametric trajectory game, where the parameter is just the initial state
function TrajectoryGamesBase.solve_trajectory_game!(
    game::TrajectoryGame{<:ProductDynamics},
    horizon,
    initial_state,
    strategy;
    parametric_game = build_parametric_game(; game, horizon),
    verbose = false,
    solving_info = nothing,
)

    # Solve, maybe with warm starting.
    if !isnothing(strategy.last_solution) && strategy.last_solution.status == PATHSolver.MCP_Solved
        solution = solve(
            parametric_game,
            initial_state;
            initial_guess = strategy.last_solution.variables,
            verbose,
        )
    else
        solution = solve(
            parametric_game,
            initial_state;
            initial_guess = generate_initial_guess(; game, parametric_game, horizon, initial_state),
            verbose,
        )
    end

    if !isnothing(solving_info)
        push!(solving_info, solution.info)
    end

    # Update warm starting info.
    if solution.status == PATHSolver.MCP_Solved
        strategy.last_solution = solution
    end
    strategy.solution_status = solution.status

    # Pack solution into OpenLoopStrategy.
    trajs = unstack_trajectory(unpack_trajectory(mortar(solution.primals); game.dynamics))
    JointStrategy(map(traj -> OpenLoopStrategy(traj.xs, traj.us), trajs))
end

# receding horizon strategy that supports warm starting
Base.@kwdef mutable struct WarmStartRecedingHorizonStrategy
    game::TrajectoryGame
    parametric_game::ParametricGame
    receding_horizon_strategy::Any = nothing
    time_last_updated::Int = 0
    turn_length::Int
    horizon::Int
    last_solution::Any = nothing
    context_state::Any = nothing
    solution_status::Any = nothing
end

function (strategy::WarmStartRecedingHorizonStrategy)(state, time)
    plan_exists = !isnothing(strategy.receding_horizon_strategy)
    time_along_plan = time - strategy.time_last_updated + 1
    plan_is_still_valid = 1 <= time_along_plan <= strategy.turn_length

    update_plan = !plan_exists || !plan_is_still_valid
    if update_plan
        strategy.receding_horizon_strategy = TrajectoryGamesBase.solve_trajectory_game!(
            strategy.game,
            strategy.horizon,
            state,
            strategy;
            strategy.parametric_game,
        )
        strategy.time_last_updated = time
        time_along_plan = 1
    end

    strategy.receding_horizon_strategy(state, time_along_plan)
end

