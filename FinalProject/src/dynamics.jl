
# NOTE: this is based on planar_double_integrator from TrajectoryGamesExamples.jl but adapted for three-dimensional motion
function orbital_double_integrator(init, i; dt = 0.1, m = 1, mu = 3.986e14, kwargs...)

    # get dt parameter
    dt2  = 0.5 * dt * dt

    # time-varying linear dynamics (nonlinear dynamics linearized about nominal)
    # NOTE: If m = 1, then input u := (ax, ay, az). If m != 1, then input u := (fx, fy, fz)
    A = Vector{Matrix{Float64}}(undef, init.horizon)
    for t = 1:init.horizon

        # get nominal position for ith player
        xnom = init.xnoms[t][Block(i)]
        r    = xnom[1:3]
        rmag = norm(r)

        # Jacobian of orbital dynamics w.r.t. state
        A[t] = [
            I(3)                               Diagonal([dt, dt, dt])
            (-mu / rmag^3) * I(3) + 3 * mu * (r * r') / rmag^5   I(3)
        ]
    
    end

    # Jacobian of dynamics w.r.t. control input
    B = [
        [
            dt2 0.0 0.0
            0.0 dt2 0.0
            0.0 0.0 dt2
            dt 0.0 0.0
            0.0 dt 0.0
            0.0 0.0 dt
        ] / m
        for t in 1:init.horizon
    ]

    return LinearDynamics(; A, B, kwargs...)

end


