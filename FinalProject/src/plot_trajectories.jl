
using Plots
using BlockArrays

# run collision avoidance scenario
traj = FinalProject.run()

# get states and control inputs
xs = traj.xs
us = traj.us

# get x-, y-, and z-coordinates for both satellites over time
xA_x = [x[1] for x in xs]
xA_y = [x[2] for x in xs]
xA_z = [x[3] for x in xs]

xB_x = [x[7] for x in xs]
xB_y = [x[8] for x in xs]
xB_z = [x[9] for x in xs]

println(xA_x - xB_x)

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

