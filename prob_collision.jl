include("orbit_utils.jl")
include("cov_utils.jl")
using DifferentialEquations
using LinearAlgebra

# orbit_utils.jl and cov_utils.jl hold supporting functions for this script
# Adjust TCA_days below to change number of days until time of closest approach
# of satellites 
# Adjust 

# Sat A = satellite in circular orbit
# Sat B = satellite in elliptical transfer orbit

# Initialize Parameters:
mu = 3.986e5        # km^3/s^2
re = 6378.1     # km
HBR = 0.002         # combined hard body radius (2 meters)
TCA_days = 7        # TCA in days
TCA = TCA_days*24*3600     # TCA in sec
tspan = (0.0, TCA)
arc_length_dist = 0.5   # distance along Sat A orbit between Sat A and Sat B at TCA

##------------------------- COMPUTE INIT STATE (ECI) -------------------------
# Smaller Circular orbit for Hohmann
a1 = re + 685           # km
e1 = 0.0
inc1 = deg2rad(98.2)    # sun synchronous orbit
argp1 = 0.0
RAAN1 = 0.0

# Larger circ orbit for Hohmann (Aura Spacecraft (Duncan & Long paper))
a2 = re + 705  # km
e2 = 0.0
inc2 = inc1    # sun synchronous orbit
argp2 = 0.0
RAAN2 = 0.0

# Hohmann transfer ellipse
a_t = (a1 + a2)/2
e_t = (a2 - a1) / (a1 + a2)
inc_t = inc1
argp_t = 0.0
RAAN_t = 0.0

# Relative phasing for collision in TCA_days
# Want Sat A (circular orbit) to be an arc length of 0.5km behind sat B at TCA
nu_meetB = pi                   # true anomaly at TCA to ensure sat B close approach
offset = arc_length_dist/a2     # central angle in radius required for assigned arc length
nu_meetA = pi - offset          # true anomaly at TCA to ensure sat A close approach

# Calculate initial true anom for circular orbit
n2  = sqrt(mu / a2^3)                 # outer circular mean motion
nu2_0 = mod(nu_meetA - n2*TCA, 2*pi)   # required initial true anomaly for sc in larger circ orbit

# Calculate initial true anom for elliptical transfer orbit
# First, must compute mean anomaly at nu_meetB on the transfer orbit
cosE_meet = (e_t + cos(nu_meetB)) / (1 + e_t*cos(nu_meetB))
E_meet = acos(cosE_meet)

# Account for fact that true anom, mean anom, and eccentric anom all in same half plane
if nu_meetB > pi
    E_meet = 2*pi - E_meet
end
M_meet = E_meet - e_t*sin(E_meet)
nt  = sqrt(mu / a_t^3)          # elliptical transfer mean motion

# Back out the initial mean anomaly for transfer at t=0
M_t0 = M_meet - nt*TCA
# Convert that back to true anomaly nut_0
nut_0 = true_anomaly_from_mean(M_t0, e_t)

# Get ECI states (X, Y, Z, Vx, Vy, Vz), km and km/s
satA_state0 = kepler2cart(a2, e2, inc2, argp2, RAAN2, nu2_0; mu=mu)
satB_state0 = kepler2cart(a_t, e_t, inc_t, argp_t, RAAN_t, nut_0; mu=mu)

println("\nSat A (circular) Initial circular state (r,v):")
println(satA_state0)

println("\nSat B (elliptical transfer) Initial State (r,v):")
println(satB_state0)

##------------------------- COMPUTE NOMINAL TRAJECTORY -------------------------
# Sat A (circular) propagation
prob_A = ODEProblem(two_body!, satA_state0, tspan, mu)
sol_A = solve(prob_A, reltol=1e-12, abstol=1e-14)

# Sat B (elliptical transfer) propagation
prob_B = ODEProblem(two_body!, satB_state0, tspan, mu)
sol_B = solve(prob_B, reltol=1e-12, abstol=1e-14)

# States at TCA
satA_stateTCA = sol_A(TCA)  #[r;v]
satB_stateTCA = sol_B(TCA)  #[r;v]

satA_posTCA = satA_stateTCA[1:3]      #[r]
satB_posTCA = satB_stateTCA[1:3]      #[r]
satA_velTCA = satA_stateTCA[4:6]      #[v]
satB_velTCA = satB_stateTCA[4:6]      #[v]

# Calculate miss distance (MD) vector (rho) in ECI km
rho_TCA_ECI = satA_posTCA .- satB_posTCA
MD = norm(rho_TCA_ECI)

# Using Sat A (circular) as reference for RTN frame
r_ref = satA_posTCA
v_ref = satA_velTCA

# Rotation matrix from ECI --> RTN
R_eci2rtn = RTN2ECI(r_ref, v_ref)

# Miss distance vector (rho) in RTN frame
rho_TCA_RTN = R_eci2rtn * rho_TCA_ECI
MD_RTN = norm(rho_TCA_RTN)

println("\nMiss Distance at TCA (km): ", MD_RTN)
println("Relative Position at TCA (km) in RTN: ", rho_TCA_RTN)

##--------------------------- COMPUTE COVARIANCE --------------------------------
P_satA = cov_ECI(satA_posTCA, satA_velTCA, TCA)
P_satB = cov_ECI(satB_posTCA, satB_velTCA, TCA)
P = P_satA + P_satB

## -------------------- COMPUTE PROBABILITY OF COLLISION ------------------------
# First determine encounter frame (following Dr. Jones Orbital Debris notes)
# T performs rotation from ECI to conjunction frame
r_rel = satB_posTCA .- satA_posTCA
v_rel = satB_velTCA .- satA_velTCA

u_hatx = r_rel / norm(r_rel)
u_haty = cross(r_rel, v_rel)/norm(cross(r_rel, v_rel))
T = vcat(u_hatx', u_haty')

rho_2D = T * rho_TCA_ECI
P_2D = T * P * T'

Pc_Alfriend = Pc(rho_2D, P_2D, HBR)
println("\nProbability of Collision Alfriend/Akella Method: ", Pc_Alfriend)

# Checking Alfriend/Akella 2D-Pc method above with Chan's implementation:
sig_x = sqrt(P_2D[1,1])
sig_y = sqrt(P_2D[2,2])

x = rho_2D[1]
y = rho_2D[2]

u = HBR^2 / (sig_x * sig_y)
v = (x^2 / sig_x^2) + (y^2 / sig_y^2)

Pc_Chan = exp(-v/2) * (1 - exp(-u/2))
println("\nProbability of Collision Chan Method: ", Pc_Chan)