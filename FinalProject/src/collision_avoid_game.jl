

# Steps:
# cd FinalProject
# julia
# ] activate .
# import Pkg
# Pkg.instantiate()
# using FinalProject
# FinalProject.run()

using FinalProject

function run()

    # set up model
    model = Model(Ipopt.Optimizer)

    # define state vector (S1 pos, S1 vel, S2 pos, S2 vel)
    x̃ = @variable(model, x[1:12])
    x = BlockArray(x̃, [3, 3, 3, 3])

    # define objective (cost function)
    @NLobjective(model, Min, x[1]^2 + x[2]^2)

    # define constraints
    @NLconstraint(model, x[1] + x[2] ≥ 1)

    # solve
    optimize!(model)
    value.(x)

end



