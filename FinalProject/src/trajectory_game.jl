
# set up two player/satellite trajectory game
# NOTE: currently set to have no environment and allow satellites to move anywhere in 3D space
function setup_trajectory_game(; environment = nothing)

    # locally define cost for each player
    cost = let
        function stage_cost(x, u, t, θ)
            x1, x2 = blocks(x)
            u1, u2 = blocks(u)

            # cost for player 1 and player 2
            # NOTE: currently set to only penalize fuel consumption
            [
                norm_sqr(u1),
                norm_sqr(u2),
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

        # three-dimensional double integrator
        sat_dynamics = custom_double_integrator(;
            state_bounds = (; lb = [-Inf, -Inf, -Inf, -Inf, -Inf, -Inf], ub = [Inf, Inf, Inf, Inf, Inf, Inf]),
            control_bounds = (; lb = [-Inf, -Inf, -Inf], ub = [Inf, Inf, Inf]),
        )
        # THIS NEEDS TO BE CHANGED ABOVE AND BELOW BECAUSE DYNAMICS ARE NOW DEPENDENT ON NOMINAL STATE
        dynamics = ProductDynamics([sat_dynamics for _ in 1:2])

        TrajectoryGame(dynamics, cost, environment, collision_avoid_constraint)

    end

end



