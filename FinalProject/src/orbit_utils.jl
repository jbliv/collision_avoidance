
using LinearAlgebra
#
# Computes the initial states for 2 satellites to ensure a close approach after a 
# user-defined time before Time of Closest Approach (TCA). TCA variable currently 
# defined as 7 days, so satellites assigned these initial states will have a near
# miss after 7 days.
# Currently designed for a miss distance of .5 km at TCA (500m)

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

# Two-body dynamics (state = [r;v] in km, km/s)
function two_body!(dstate, state, mu, t)
    r = state[1:3]
    v = state[4:6]
    rnorm = norm(r)

    a = -mu*r / rnorm^3  # km/s^2

    dstate[1:3] = v
    dstate[4:6] = a
end