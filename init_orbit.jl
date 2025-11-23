using LinearAlgebra
"""
Computes the initial states for 2 satellites to ensure a collision after a user-defined
time before Time of Closest Approach (TCA). TCA variable currently defined as 7 days,
so satellites assigned these initial states will "collide" after 7 days.
"""
##------------------------------ IMPORTANT FUNCTIONS ------------------------------
# Keplerian to Cartesian in ECI for circular & elliptical orbits
function kepler2cart(a, e, inc, w, RAAN, nu; mu = 3.986e5)

    p = a * (1 - e^2)
    r0 = p / (1 + e*cos(nu))

    # Position in perifocal coordinates (PQW)
    x = r0 * cos(nu)
    y = r0 * sin(nu)

    # Velocity in PQW
    vx_pqw = -sqrt(mu/p) * sin(nu)
    vy_pqw =  sqrt(mu/p) * (e + cos(nu))
    vz_pqw = 0.0

    # Rotation matrix PQW --> ECI
    cw, sw = cos(w), sin(w)
    cO, sO = cos(RAAN), sin(RAAN)
    ci, si = cos(inc),  sin(inc)

    Rz_O = [ cO -sO 0.0;
             sO  cO 0.0;
             0.0 0.0 1.0 ]

    Rx_i = [ 1.0 0.0 0.0;
             0.0 ci -si;
             0.0 si  ci ]

    Rz_w = [ cw -sw 0.0;
             sw  cw 0.0;
             0.0 0.0 1.0 ]

    Q = Rz_O * Rx_i * Rz_w
    r = Q * [x, y, 0.0]
    v = Q * [vx_pqw, vy_pqw, vz_pqw]
    return vcat(r, v)
end

function true_anomaly_from_mean(M, e; tol=1e-12, maxiter=20)
    # wrap M to [0, 2π)
    M = mod(M, 2*pi)
    E = 0.5 #Initial Guess

    # Solve Kepler's equation: M = E - e*sin(E)
    # Newton-Rhapson implementation from notes
    for _ in 1:maxiter
        m_mi = M - (E - e*sin(E))
        dm_dE = 1-e*cos(E)
        E_next = E + m_mi/dm_dE

        if abs(E_next - E) < tol
            E = E_next
            break
        end
        E = E_next
    end

    # Convert eccentric → true anomaly
    cos_nu = (cos(E) - e) / (1 - e*cos(E))
    nu = acos(cos_nu)

    if M > pi
        nu= 2*pi - nu
    end

    return nu
end

##------------------------------ COMPUTING INIT STATE ------------------------------
# Constants
mu = 3.986e5     # km^3/s^2
re = 6378.1     # km

# Inner orbit
a1 = re + 685           # 0 km altitude
e1 = 0.0
inc1 = deg2rad(98.2)    # sun synchronous orbit
argp1 = 0.0
RAAN1 = 0.0

# Outer orbit (Aura Spacecraft (Duncan & Long paper))
a2 = re + 705           # km
e2 = 0.0
inc2 = deg2rad(98.2)    # sun synchronous orbit
argp2 = 0.0
RAAN2 = 0.0

# Hohmann transfer ellipse
a_t = (a1 + a2)/2
e_t = (a2 - a1) / (a1 + a2)
inc_t = inc1
argp_t = 0.0
RAAN_t = 0.0

# Relative phasing for collision in 7 days
TCA = 7*24*3600      # seconds
nu_meet = pi         # meet at apogee of transfer

# Calculate initial true anom for circular orbit
n2  = sqrt(mu / a2^3)                 # outer circular mean motion
nu2_0 = mod(nu_meet - n2*TCA, 2*pi)   # required initial true anomaly for sc in larger circular orbit

# Calculate initial true anom for elliptical transfer orbit
# First, must compute mean anomaly at nu_meet on the transfer orbit
cosE_meet = (e_t + cos(nu_meet)) / (1 + e_t*cos(nu_meet))
E_meet = acos(cosE_meet)

# Account for fact that true anom, mean anom, and eccentric anom all in same half plane
if nu_meet > pi
    E_meet = 2*pi - E_meet
end
M_meet = E_meet - e_t*sin(E_meet)
nt  = sqrt(mu / a_t^3)          # elliptical transfer mean motion

# Back out the initial mean anomaly for transfer at t=0
M_t0 = M_meet - nt*TCA
# Convert that back to true anomaly nut_0
nut_0 = true_anomaly_from_mean(M_t0, e_t)

# Get ECI states (X, Y, Z, Vx, Vy, Vz), km and km/s
state_transfer0 = kepler2cart(a_t, e_t, inc_t, argp_t, RAAN_t, nut_0; mu=mu)
state_circ0     = kepler2cart(a2, e2, inc2, argp2, RAAN2, nu2_0; mu=mu)

println("Initial transfer state (r,v):")
println(state_transfer0)

println("\nInitial circular state (r,v):")
println(state_circ0)