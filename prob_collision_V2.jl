include("orbit_utils.jl")
include("cov_utils.jl")
using DifferentialEquations
using LinearAlgebra

# orbit_utils.jl and cov_utils.jl hold supporting functions for this script
# Adjust TCA_days below to change number of days until time of closest approach
# of satellites 
# Since set up is for a collision, TCA set by user must be < orbital period
# frac is multiplier to orbital period to set TCA, set frac < 1 to ensure 
# satellites collide only once during the time spacn of simulation

# Sat A = satellite in inclination 98.2 deg
# Sat B = satellite in inclination 110 deg

# Initialize Parameters:
mu = 3.986e5        # km^3/s^2
re = 6378.1         # km
HBR = 0.003         # combined hard body radius (3 meters)
frac = 0.5          # fraction of orbital period

##------------------------- COMPUTE INIT STATE (ECI) -------------------------
# Choose common circular orbit radius
alt = 700.0
a   = re + alt
e   = 0.0

# Different inclinations
incA = deg2rad(98.2)    # Sat A
incB = deg2rad(110)    # Sat B

# Same RAAN and argument of periapsis so they intersect at line of nodes
RAANA = 0.0
RAANB = 0.0
argp_A = 0.0
argp_B = 0.0

# Choose Time of Closest Approach (collision time)
period = 2*pi*sqrt(a^3/mu)      # T in sec
period_days = period/24/3600    # T in days

TCA_days = frac * period_days   # TCA in days
TCA = frac*period               # TCA in sec
tspan = (0.0, TCA)

# Mean motion (same for both)
n = sqrt(mu / a^3)

# We want both at ascending node (ν = 0) at t = TCA
# So M(TCA) = 0 mod 2π → M0 = -n*TCA mod 2π
M0 = mod(-n * TCA, 2π)

# Convert M0 → ν0 (for e=0, ν0 ≈ M0, but we’ll use your function)
true_anom0 = true_anomaly_from_mean(M0, e)

println("Initial true anomaly ν0 (rad) = ", true_anom0)

# Get states in ECI
satA_state0 = kepler2cart(a, e, incA, argp_A, RAANA, true_anom0; mu=mu)
satB_state0 = kepler2cart(a, e, incB, argp_B, RAANB, true_anom0; mu=mu)

println("Sat A initial state (r,v): ", satA_state0)
println("Sat B initial state (r,v): ", satB_state0)

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

# Assume baseline level of uncertainty at TCA (user defined)
# Nominal uncertainty is 10-50m (1-sigma) in radial, 5-20m (1-sigma) in Normal,
# and 50-500m (1-sigma) in-track for well tracked objects
baseline_unc = [(10/1000)^2, (50/1000)^2, (10/1000)^2]
P_add = diagm(baseline_unc)
P_add_ECI = R_eci2rtn' * P_add * R_eci2rtn      #Rotate RTN --> ECI

# Add baseline uncertainty to both sat A and sat B covariance
P = P_satA + P_satB + P_add_ECI + P_add_ECI

println("\nCombined Position Uncertainty: ", P)

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

# Typical Threshold used by NASA CARA team to perform CAM is 1e-4 Pc value
if Pc_Alfriend > 1e-4
    println("Pc Above Threshold: CAM Required")
else
    println("Pc Below Threshold: CAM Not Required")
end