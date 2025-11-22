
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

    @variable(model, τ[1:(init.num_states+init.num_control)*init.horizon])
    τ = BlockArray(τ, [Int(init.horizon*(init.num_states/2+init.num_control/2)), 
                       Int(init.horizon*(init.num_states/2+init.num_control/2))])
    

    # # extract objectives and constraints
    # f_new(τ) = fs[1](τ)
    # println(fs[1](τ))
    # JuMP.register(model, :f_new, 1, f_new; autodiff=false)

    function f_new_flat(x::T...) where T<:Real
        # x is now a tuple of numbers; convert to vector
        τ_vec = collect(x)
        
        # call your original fs[1], but you need to convert τ_vec to a BlockArray
        # with the same shape as fs[1] expects:
        τ_block = BlockArray(τ_vec, [Int(length(τ_vec)/2), Int(length(τ_vec)/2)]) 
    
        # compute the cost
        fs[1](τ_block)
    end

    register(model, :f_new_flat, length(τ), f_new_flat; autodiff = true)
    @NLobjective(model, Min, f_new_flat(τ...))

    # @NLobjective(model, Min, f_new(τ))
    println(gs[1](τ))
    # @constraint(model, gs[1](τ))
    # @constraint(model, g̃)
    # @constraint(model, hs)
    # @constraint(model, h̃)

    optimize!(model)
    # value.(x)

    # FOOD FOR THOUGHT: if I do MPC is that a feedback SQP approach
    # but without it and only doing one time horizon would be an open loop approach?

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



