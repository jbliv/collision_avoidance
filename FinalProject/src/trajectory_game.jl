
using TrajectoryGamesBase
using ProgressMeter: ProgressMeter


# unpack trajectory
function unpack_trajectory(flat_trajectory; init)
    trajs = Iterators.map(1:init.num_players, blocks(flat_trajectory)) do ii, τ
        horizon = Int(length(τ) / ((init.num_states + init.num_control)/2))
        num_states = Int(init.num_states/2) * horizon
        X = reshape(τ[1:num_states], (Int(init.num_states/2), horizon))
        U = reshape(τ[(num_states + 1):end], (Int(init.num_control/2), horizon))

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
