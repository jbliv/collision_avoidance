
# orbital dynamics
function dynamics(x, u; mu = 3.986e5)

    r    = x[1:3]
    # rmag = norm(r)
    rmag2 = sum(r[i]^2 for i=1:3)
    rmag = sqrt(rmag2 + 1e-6)

    xdot      = Vector{Any}(undef, 6)
    xdot[1:3] = x[4:6]
    xdot[4:6] = (-mu / rmag^3) * r + u

    return xdot

end

function integrate(x, u, k, i, init; mu = 3.986e5)

    # nominal position
    r    = init.xnoms[k][Block(i)][1:3]
    rmag = norm(r)

    # linear dynamics
    Ad = I(6) + init.dt * [
            zeros(3,3)                    I(3)
            (-mu/rmag^3*I(3) + 3*mu*(r*r')/rmag^5)  zeros(3,3)
    ]

    Bd = init.dt * [zeros(3,3); I(3)]

    return Ad * x + Bd * u

    # forward Euler
    # return x + init.dt * dynamics(x, u)

    # runge-kutta
    # k1 = dynamics(x, u)
    # k2 = dynamics(x + 0.5*dt*k1, u)
    # k3 = dynamics(x + 0.5*dt*k2, u)
    # k4 = dynamics(x + dt*k3, u)

    # return x + (dt/6) * (k1 + 2*k2 + 2*k3 + k4)

end

# function orbital_double_integrator(init, i; xnom_start=nothing, m=1, mu=3.986e5, kwargs...)
#     N = init.n_sim_steps
#     A = Vector{Matrix{Float64}}(undef, N)

#     for t = 1:N
#         # If xnom_start is given, propagate it forward for linearization
#         xnom = isnothing(xnom_start) ? init.xnoms[t][Block(i)] : xnom_start
#         r = xnom[1:3]
#         rmag = norm(r)

#         A[t] = I(6) + init.dt * [
#             zeros(3,3)                    I(3)
#             (-mu/rmag^3*I(3) + 3*mu*(r*r')/rmag^5)  zeros(3,3)
#         ]
#     end

#     B = [ init.dt * [zeros(3,3); I(3)] / m for t = 1:N ]

#     LinearDynamics(; A, B, kwargs...)
# end


