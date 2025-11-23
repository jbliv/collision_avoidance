
using FinalProject

function run()

    # get initial conditions
    init = init_conds()

    # set up trajectory and parametric game
    game = setup_trajectory_game(; environment = nothing)
    parametric_game = build_parametric_game(; game, init)

    sim_steps = let
        progress = ProgressMeter.Progress(init.n_sim_steps)
        receding_horizon_strategy =
            WarmStartRecedingHorizonStrategy(; game, parametric_game, init.turn_length, init.horizon)

        rollout(
            game.dynamics,
            receding_horizon_strategy,
            init.x,
            init.n_sim_steps;
            get_info = (γ, x, t) ->
                (ProgressMeter.next!(progress); γ.receding_horizon_strategy),
        )

    end

    

end
 



