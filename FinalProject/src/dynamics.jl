
# # NOTE: this is based on planar_double_integrator from TrajectoryGamesExamples.jl but adapted for three-dimensional motion
# function orbital_double_integrator(init, i; m = 1, mu = 3.986e5, kwargs...)

#     # time-varying linear dynamics (nonlinear dynamics linearized about nominal)
#     # NOTE: If m = 1, then input u := (ax, ay, az). If m != 1, then input u := (fx, fy, fz)
#     A = Vector{Matrix{Float64}}(undef, init.n_sim_steps)
#     for t = 1:init.n_sim_steps

#         # get nominal position for ith player
#         xnom = init.xnoms[t][Block(i)]
#         r    = xnom[1:3]
#         rmag = norm(r)

#         # Jacobian of orbital dynamics w.r.t. state
#         A[t] = I(6) + init.dt * [
#             zeros(3,3)                                                       I(3)
#             ((-mu / rmag^3) * I(3) + 3 * mu * (r * r') / rmag^5)  zeros(3,3)
#         ]
    
#     end

#     # Jacobian of dynamics w.r.t. control input
#     B = [
#         init.dt * [
#             zeros(3,3)
#             I(3)
#         ] / m
#         for t in 1:init.n_sim_steps
#     ]

#     return LinearDynamics(; A, B, kwargs...)

# end

function orbital_double_integrator(init, i; xnom_start=nothing, m=1, mu=3.986e5, kwargs...)
    N = init.n_sim_steps
    A = Vector{Matrix{Float64}}(undef, N)

    for t = 1:N
        # If xnom_start is given, propagate it forward for linearization
        xnom = isnothing(xnom_start) ? init.xnoms[t][Block(i)] : xnom_start
        r = xnom[1:3]
        rmag = norm(r)

        A[t] = I(6) + init.dt * [
            zeros(3,3)                    I(3)
            (-mu/rmag^3*I(3) + 3*mu*(r*r')/rmag^5)  zeros(3,3)
        ]
    end

    B = [ init.dt * [zeros(3,3); I(3)] / m for t = 1:N ]

    LinearDynamics(; A, B, kwargs...)
end


