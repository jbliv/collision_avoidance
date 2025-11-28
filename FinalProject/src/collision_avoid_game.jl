
using FinalProject

function run()

    # get initial conditions
    init = init_conds()

    # set up trajectory and parametric game
    traj = iterative_best_response(init)

    return traj

end
 



