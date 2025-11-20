using Pkg
Pkg.activate(joinpath(@__DIR__, ".."))

using CollisionAvoidance
using SqpSolver
using JuMP
using Ipopt
using Plots
gr() # Explicitly use GR backend

# Define parameters
params = CollisionGameParams(
    1.0,    # v
    0.2,    # dt
    100,    # Timesteps
    2.0,    # R_safe
    [-10.0, 0.0, 0.0], # x0_1
    [10.0, 0.0, pi],   # x0_2
    10.0,   # w_path
    10.0,   # w_heading
    0.1     # w_control
)

println("Build problem")
model, x1, x2, u1, u2 = build_model(params)

println("Set up optimizer")
set_optimizer(model, optimizer_with_attributes(
    SqpSolver.Optimizer, 
    "external_optimizer" => optimizer_with_attributes(Ipopt.Optimizer, "print_level" => 0),
    "max_iter" => 1000,
    "tol_residual" => 1e-4,
    "OutputFlag" => 0
))

println("Start solver")
optimize!(model)

println("Solve")
println("Termination status ", termination_status(model))
println("Primal status ", primal_status(model))

# Unpack solution
X1_val = value.(x1)
X2_val = value.(x2)
U1_val = value.(u1)
U2_val = value.(u2)

# Plot
p1 = plot(X1_val[1, :], X1_val[2, :], label="Agent 1", xlabel="x", ylabel="y", title="Trajectories", lw=2)
plot!(p1, X2_val[1, :], X2_val[2, :], label="Agent 2", lw=2)

# Plot safety distance
dists = sqrt.((X1_val[1, :] .- X2_val[1, :]).^2 .+ (X1_val[2, :] .- X2_val[2, :]).^2)
p2 = plot(dists, label="Distance", xlabel="Step", ylabel="Distance", title="Separation", lw=2)
plot!(p2, [1, params.N+1], [params.R_safe, params.R_safe], label="Safety Limit", linestyle=:dash, color=:red)

p3 = plot(U1_val, label="Control 1", xlabel="Step", ylabel="Turn Rate", title="Controls")
plot!(p3, U2_val, label="Control 2")

final_plot = plot(p1, p2, p3, layout=(3, 1), size=(800, 1000))
savefig(final_plot, "collision_avoidance.png")
println("Plot saved to collision_avoidance.png")
