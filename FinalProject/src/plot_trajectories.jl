
using Plots

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
