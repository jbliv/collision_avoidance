
using Plots
using BlockArrays
using JuMP
using LinearAlgebra
using FinalProject

# Set backend
gr()

# run collision avoidance scenario
traj = FinalProject.run()
println("Starting plotting")

# get states and control inputs
xs = value.(traj.xs)
us = value.(traj.us)

# get x-, y-, and z-coordinates for both satellites over time
xA_x = [x[Block(1)][1] for x in xs]
xA_y = [x[Block(1)][2] for x in xs]
xA_z = [x[Block(1)][3] for x in xs]

xB_x = [x[Block(2)][1] for x in xs]
xB_y = [x[Block(2)][2] for x in xs]
xB_z = [x[Block(2)][3] for x in xs]

p1 = plot(xA_x, xA_y, label="Sat A", lw=2, xlabel="X-Position [km]", ylabel="Y-Position[km]")
plot!(p1, xB_x, xB_y, label="Sat B", lw=2)
savefig(p1, "trajectoryXY.png")

p2 = plot3d(xA_x, xA_y, xA_z, label="Sat A", lw=2)
plot!(p2, xB_x, xB_y, xB_z, label="Sat B", lw=2)
savefig(p2, "trajectory3D.png")

# get nominal states
init  = FinalProject.init_conds()
xnoms = init.xnoms

# get x-, y-, and z-coordinates for both nominal states over time
xnomA_x = [xnom[Block(1)][1] for xnom in xnoms]
xnomA_y = [xnom[Block(1)][2] for xnom in xnoms]
xnomA_z = [xnom[Block(1)][3] for xnom in xnoms]

xnomB_x = [xnom[Block(2)][1] for xnom in xnoms]
xnomB_y = [xnom[Block(2)][2] for xnom in xnoms]
xnomB_z = [xnom[Block(2)][3] for xnom in xnoms]

p3 = plot(xnomA_x, xnomA_y, label="Sat A", lw=2, xlabel="X-Position [km]", ylabel="Y-Position[km]")
plot!(p3, xnomB_x, xnomB_y, label="Sat B", lw=2)
savefig(p3, "nomtrajectoryXY.png")

p4 = plot3d(xnomA_x, xnomA_y, xnomA_z, label="Sat A", lw=2)
plot!(p4, xnomB_x, xnomB_y, xnomB_z, label="Sat B", lw=2)
savefig(p4, "nomtrajectory3D.png")

# get control input
# Convert acceleration (km/s^2) to Delta V (m/s) by multiplying by dt and 1000
dt = init.dt
to_mps = 1000.0 * dt

u_norm = [sqrt(u[Block(1)][1]^2 + u[Block(1)][2]^2 + u[Block(1)][3]^2) * to_mps for u in us]
println("Total Delta V (Sat A): $(sum(u_norm)) m/s")
# u_norm = cumsum(u_norm)
ts     = 1:length(us)
p5 = plot(ts, u_norm, lw=2, label="Sat A", xlabel="Time Step", ylabel="Delta V [m/s]")
savefig(p5, "control.png")

# plot separation distance
dists = [norm(x[Block(1)][1:3] - x[Block(2)][1:3]) for x in xs]
min_dist = minimum(dists)
p6 = plot(ts, dists, label="Separation Distance", lw=2, xlabel="Time Step", ylabel="Distance [km]", yaxis=:log)
hline!(p6, [min_dist], label="Min Dist: $(round(min_dist, digits=4)) km", color=:red, linestyle=:dash)
savefig(p6, "separation_distance.png")

# plot collision probability
pcs = Float64[]
HBR = 0.002 

for (i, x) in enumerate(xs)
    xA = x[Block(1)]
    xB = x[Block(2)]
    
    rA = xA[1:3]
    vA = xA[4:6]
    rB = xB[1:3]
    vB = xB[4:6]
    
    t_sec = i * init.dt

    PA = cov_ECI(rA, vA, t_sec)
    PB = cov_ECI(rB, vB, t_sec)
    P_combined = PA + PB
    
    rho = rA - rB

    try
        val = Pc(rho, P_combined, HBR)
        push!(pcs, val)
    catch e
        push!(pcs, 0.0)
    end
end

p7 = plot(ts, pcs, label="Probability of Collision", lw=2, xlabel="Time Step", ylabel="Pc")
savefig(p7, "collision_probability.png")

# directional control plots

uA_x = [u[Block(1)][1] * to_mps for u in us]
uA_y = [u[Block(1)][2] * to_mps for u in us]
uA_z = [u[Block(1)][3] * to_mps for u in us]

uB_x = [u[Block(2)][1] * to_mps for u in us]
uB_y = [u[Block(2)][2] * to_mps for u in us]
uB_z = [u[Block(2)][3] * to_mps for u in us]

p_ctrl_1 = plot(ts, [uA_x uB_x], label=["Sat A" "Sat B"], ylabel="dV_x [m/s]", lw=2)
p_ctrl_2 = plot(ts, [uA_y uB_y], label=["Sat A" "Sat B"], ylabel="dV_y [m/s]", lw=2)
p_ctrl_3 = plot(ts, [uA_z uB_z], label=["Sat A" "Sat B"], ylabel="dV_z [m/s]", lw=2)

p8 = plot(p_ctrl_1, p_ctrl_2, p_ctrl_3, layout=(3, 1), size=(800, 900), xlabel="Time Step")
savefig(p8, "detailed_control.png")

# 3D plot with Earth wireframe, (from ChatGPT scale is a bit off and doesn't show much of the orbit TODO: fix)
function sphere(r, n=20)
    u = range(0, 2π, length=n)
    v = range(0, π, length=n)
    x = r .* cos.(u) * sin.(v)'
    y = r .* sin.(u) * sin.(v)'
    z = r .* ones(n) * cos.(v)'
    return x, y, z
end

# Earth Radius ~ 6378 km
R_earth = 6378.0
sx, sy, sz = sphere(R_earth)

# Create 3D plot
# Use a better camera angle to see the polar orbits (isometric view)
p9 = plot3d(xA_x, xA_y, xA_z, label="Sat A", lw=4, color=:blue, legend=:topright)
plot!(p9, xB_x, xB_y, xB_z, label="Sat B", lw=4, color=:red)

# Add markers for start and end positions to make them visible
scatter!(p9, [xA_x[1]], [xA_y[1]], [xA_z[1]], color=:blue, markersize=4, label="Start A")
scatter!(p9, [xB_x[1]], [xB_y[1]], [xB_z[1]], color=:red, markersize=4, label="Start B")
scatter!(p9, [xA_x[end]], [xA_y[end]], [xA_z[end]], color=:blue, markershape=:square, markersize=4, label="End A")
scatter!(p9, [xB_x[end]], [xB_y[end]], [xB_z[end]], color=:red, markershape=:square, markersize=4, label="End B")

# Add Earth wireframe
# Use plot3d for wireframe lines manually to avoid GR surface issues
for i in 1:size(sx, 1)
    plot3d!(p9, sx[i, :], sy[i, :], sz[i, :], color=:lightgrey, alpha=0.3, label="", lw=0.5)
end
for j in 1:size(sx, 2)
    plot3d!(p9, sx[:, j], sy[:, j], sz[:, j], color=:lightgrey, alpha=0.3, label="", lw=0.5)
end

# Set camera to isometric view (Azimuth=45, Elevation=30)
plot!(p9, camera=(45, 30))

savefig(p9, "trajectory3D_earth.png")

