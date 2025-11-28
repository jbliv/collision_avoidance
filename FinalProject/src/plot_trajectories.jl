
using Plots
using BlockArrays
using JuMP

# run collision avoidance scenario
traj = FinalProject.run()

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

plot(xA_x, xA_y, label="Sat A", lw=2, xlabel="X-Position [km]", ylabel="Y-Position[km]")
plot!(xB_x, xB_y, label="Sat B", lw=2, xlabel="X-Position [km]", ylabel="Y-Position[km]")
savefig("trajectoryXY.png")

plot3d(xA_x, xA_y, xA_z, label="Sat A", lw=2)
plot!(xB_x, xB_y, xB_z, label="Sat B", lw=2)
savefig("trajectory3D.png")

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

plot(xnomA_x, xnomA_y, label="Sat A", lw=2, xlabel="X-Position [km]", ylabel="Y-Position[km]")
plot!(xnomB_x, xnomB_y, label="Sat B", lw=2, xlabel="X-Position [km]", ylabel="Y-Position[km]")
savefig("nomtrajectoryXY.png")

plot3d(xnomA_x, xnomA_y, xnomA_z, label="Sat A", lw=2)
plot!(xnomB_x, xnomB_y, xnomB_z, label="Sat B", lw=2)
savefig("nomtrajectory3D.png")

# get control input
u_norm = [sqrt(u[Block(1)][1]^2 + u[Block(1)][2]^2 + u[Block(1)][3]^2) for u in us]
println(sum(u_norm) / length(u_norm))
# u_norm = cumsum(u_norm)
ts     = 1:length(us)
plot(ts, u_norm, lw=2)
savefig("control.png")

