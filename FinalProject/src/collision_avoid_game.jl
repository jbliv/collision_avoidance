
using FinalProject

function run()

    # get initial conditions
    init = init_conds()

    # set up trajectory and parametric game
    game_setup       = setup_trajectory_game(; environment = nothing)
    fs, gs, hs, g̃, h̃ = get_objectives_and_constraints(; game = game_setup, horizon = init.horizon)

    # TESTING
    model = Model(optimizer_with_attributes(
        SqpSolver.Optimizer, 
        "external_optimizer" => Ipopt.Optimizer,
    ))

    # extract objectives and constraints
    @NLobjective(model, Min, fs)
    @NLconstraint(model, gs)
    @NLconstraint(model, g̃)
    @NLconstraint(model, hs)
    @NLconstraint(model, h̃)

    optimize!(model)
    value.(x)


end

function test()

    # set up model to use SQP solver
    model = Model(optimizer_with_attributes(
        SqpSolver.Optimizer, 
        "external_optimizer" => Ipopt.Optimizer,
    ))

    # define state vector (S1 pos, S1 vel, S2 pos, S2 vel)
    x̃  = @variable(model, x̃[1:12])
    x  = BlockArray(x̃, [6, 6])
    x₁ = x[Block(1)]
    x₂ = x[Block(2)]

    # define objective (cost function)
    @NLobjective(model, Min, x₁[1]^2 + x₂[1]^2)

    # define constraints
    @constraint(model, x₁[1] + x₂[1] ≥ 1)
    # @NLconstraint(model, x₁[1]^2 + x₂[1]^2 ≥ 1) # this breaks the solver without a good initial guess

    # solve
    optimize!(model)
    value.(x)

end



