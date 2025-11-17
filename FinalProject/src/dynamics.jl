
# NOTE: this is based on planar_double_integrator from TrajectoryGamesExamples.jl but adapted for three-dimensional motion
function custom_double_integrator(xnom; dt = 0.1, m = 1, mu = 3.986e14, kwargs...)

    # get position vector and magnitude 
    dt2  = 0.5 * dt * dt
    r    = xnom[1:3]
    rmag = norm(r)

    # Jacobian of dynamics with respect to the state
    A_lin           = I(6)
    A_lin[1:3, 4:6] = Diagonal([dt, dt, dt])
    A_lin[4:6, 1:3] = (-mu / rmag^3) * I(3) + 3 * mu * (r * r') / rmag^5

    # NOTE: If m = 1, then input u := (ax, ay, az). If m != 1, then input u := (fx, fy, fz)
    time_invariant_linear_dynamics(;
        A = A_lin,
        B = [
            dt2 0.0 0.0
            0.0 dt2 0.0
            0.0 0.0 dt2
            dt 0.0 0.0
            0.0 dt 0.0
            0.0 0.0 dt
        ] / m,
        kwargs...,
    )

end


