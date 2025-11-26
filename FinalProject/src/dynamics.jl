
# NOTE: this is based on planar_double_integrator from TrajectoryGamesExamples.jl but adapted for three-dimensional motion
function orbital_double_integrator(init, i; m = 1, mu = 3.986e5, kwargs...)

    # time-varying linear dynamics (nonlinear dynamics linearized about nominal)
    # NOTE: If m = 1, then input u := (ax, ay, az). If m != 1, then input u := (fx, fy, fz)
    A = Vector{Matrix{Float64}}(undef, init.n_sim_steps)
    for t = 1:init.n_sim_steps

        # get nominal position for ith player
        xnom = init.xnoms[t][Block(i)]
        r    = xnom[1:3]
        rmag = norm(r)

        # Jacobian of orbital dynamics w.r.t. state
        A[t] = I(6) + init.dt * [
            zeros(3,3)                                                       I(3)
            ((-mu / rmag^3) * I(3) + 3 * mu * (r * r') / rmag^5)  zeros(3,3)
        ]
    
    end

    # Jacobian of dynamics w.r.t. control input
    B = [
        init.dt * [
            zeros(3,3)
            I(3)
        ] / m
        for t in 1:init.n_sim_steps
    ]

    return LinearDynamics(; A, B, kwargs...)

end


