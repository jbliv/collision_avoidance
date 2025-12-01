
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

    # Note: when I expanded the simulation length I noticed a massive control spike at the beginning. I figured the dynamics could be causing it
    # the changes here fixed it

    xnom_prev = init.xnoms[k-1][Block(i)]
    xnom_curr = init.xnoms[k][Block(i)]

    r    = xnom_prev[1:3]
    rmag = norm(r)

    # linear dynamics (Jacobian at k-1)
    Ad = I(6) + init.dt * [
            zeros(3,3)                    I(3)
            (-mu/rmag^3*I(3) + 3*mu*(r*r')/rmag^5)  zeros(3,3)
    ]

    Bd = init.dt * [zeros(3,3); I(3)]

    return xnom_curr + Ad * (x - xnom_prev) + Bd * u

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


