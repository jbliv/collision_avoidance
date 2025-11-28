

function iterative_best_response(init; max_iters::Int=20, tol::Float64=1e-2)

    # use the nominal trajectory as the initial guess
    traj = get_initial_trajectories(init)

    # loop until converge to solution or reach maximum # of iterations
    for iter = 1:max_iters
        println("here")

        # store old state
        traj_old = deepcopy(traj)

        # loop through each player
        for i = 1:init.num_players

            # solve using SQP solver
            println("here")
            traj = solve_sqp(traj, i, init)

        end

        println("here")

        # check if converged
        # dx = norm(pack_trajectory(traj) - pack_trajectory(traj_old))

        # if dx < tol
        #     break
        # end

    end

    return traj

end

function solve_sqp(traj, i, init)

    # define model
    model = Model(optimizer_with_attributes(
        SqpSolver.Optimizer, 
        "external_optimizer" => Ipopt.Optimizer,
    ))

    num_vars = Int((init.num_states+init.num_control)/2)*init.horizon

    # TEMPORARY
    @variable(model, τi[1:Int((init.num_states+init.num_control)/2)*init.horizon])
    # @variable(model, τi[1:Int((init.num_states+init.num_control))*init.horizon])
    # τ       = Vector{Any}(undef, init.num_players * num_vars)
    # τ       = Vector{VariableRef}(undef, init.num_players * num_vars)

    τ_fixed = pack_trajectory(traj)

    # for j = 1:init.num_players
    #     inds = (j-1)*num_vars+1:j*num_vars
    #     if j == i
    #         τ[inds] = τi
    #     else
    #         # τ[inds] = τi[inds]
    #         τ[inds] = τ_fixed[inds]
    #     end
    #     # UPDATE IT LIKE BELOW
    # end

    τ = Vector{BlockArray}(undef, init.num_players*num_vars)
    if i == 1
        inds     = num_vars+1:2*num_vars
        τ = BlockArray(vcat(τi, τ_fixed[inds]), [num_vars, num_vars])
    else
        inds     = 1:num_vars
        τ = BlockArray(vcat(τ_fixed[inds], τi), [num_vars, num_vars])
    end

    # unpack trajectory
    (; xs, us) = unpack_trajectory(τ; init)

    # nonlinear objective
    f = get_objective(xs, us, i, init)
    @NLobjective(model, Min, f)

    # initial dynamic constraint
    for ix = 1:Int(init.num_states/2)
        @NLconstraint(model, xs[1][Block(i)][ix] - init.x[Block(i)][ix] == 0)
    end

    # nonlinear dynamic constraints
    for k = 2:init.horizon
        g = get_dynamics_constraints(xs, us, k, i, init)
        for ix = 1:Int(init.num_states/2)
            @NLconstraint(model, g[ix] == 0)
        end
    end

    optimize!(model)

    # Get optimal τi vector
    τi_solved = value.(τi)

    τ_solved = Vector{BlockArray}(undef, 2*num_vars)
    if i == 1
        inds     = num_vars+1:2*num_vars
        τ_solved = BlockArray(vcat(τi_solved, τ[inds]), [num_vars, num_vars])
    else
        inds     = 1:num_vars
        τ_solved = BlockArray(vcat(τ[inds], τi_solved), [num_vars, num_vars])
    end
    (; xs, us) = unpack_trajectory(τ_solved; init)

    return (xs = xs, us = us)

end


# dynamics constraint function
function get_dynamics_constraints(xs, us, k, i, init)

    xi   = xs[k][Block(i)]
    xim1 = xs[k-1][Block(i)]
    uim1 = us[k-1][Block(i)]

    g = xi - integrate(xim1, uim1, k, i, init)

end

# objective function
function get_objective(xs, us, i, init)

    if i == 1
        Q  = init.Q1
        R  = init.R1
    else
        Q  = init.Q2
        R  = init.R2
    end

    sum((xs[k][Block(i)] - init.xnoms[k][Block(i)])' * Q * (xs[k][Block(i)] - init.xnoms[k][Block(i)]) +
        (us[k][Block(i)]' * R * us[k][Block(i)] ) for k = 1:init.horizon)

end

function get_initial_trajectories(init)

    xs = Vector{BlockArray}(undef, init.horizon)
    us = Vector{BlockArray}(undef, init.horizon)

    for k = 1:init.horizon

        x_blocks = [init.xnoms[k][Block(1)], init.xnoms[k][Block(2)]]
        xs[k] = BlockArray(vcat(x_blocks...), map(length, x_blocks))

        u_blocks = [zeros(Int(init.num_control/2)), zeros(Int(init.num_control/2))]
        us[k] = BlockArray(vcat(u_blocks...), map(length, u_blocks))

    end

    return (xs = xs, us = us)

end