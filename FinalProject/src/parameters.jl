
mutable struct Init
    
    # S1 pos, S1 vel, S2 pos, S2 vel
    x::AbstractVector{Float64}
    xnoms::Vector{AbstractVector{Float64}}

    # time step
    dt::Float64

    # how often to replan using next horizon (must be <= horizon)
    turn_length::Int

    # number of future time steps to use for current horizon (must be <= n_sim_steps)
    horizon::Int

    # total number of time steps for entire simulation
    n_sim_steps::Int

    # number of time steps until reach closest approach
    TCA_sec::Int

    # number of players, states, and control inputs
    num_players::Int
    num_states::Int
    num_control::Int

    # Q, R weights
    Q1::Matrix{Float64}
    Q2::Matrix{Float64}
    R1::Matrix{Float64}
    R2::Matrix{Float64}

end

# initial conditions
function init_conds()

    # initial true and nominal states
    dt          = 0.1 # dt = 1 worked
    turn_length = 3
    horizon     = 200 # horizon = 20 worked
    n_sim_steps = horizon
    TCA_sec     = 10
    num_players = 2
    x           = get_init_states(TCA_sec)
    num_control = 6

    # NOTE: for now, the initial nominal and true states of the satellites
    # are set to be the same

    # get nominal states
    xnom0 = get_init_states(TCA_sec) # BlockArray(ones(12), [6, 6])
    xnoms = get_nominal_states(xnom0, n_sim_steps, dt)

    # Q, R weights
    Q1 = 100*I(6)
    Q2 = 100*I(6)
    R1 = 1*I(3)
    R2 = 1*I(3)

    init = Init(x,
                xnoms,
                dt,
                turn_length,
                horizon,
                n_sim_steps,
                TCA_sec,
                num_players,
                length(x),
                num_control,
                Q1,
                Q2,
                R1,
                R2)

    return init

end

# pre-compute nominal states
function get_nominal_states(xnom0, n_sim_steps, dt; mu = 3.986e5)

    # initialize nominal state array
    xnoms    = Vector{AbstractVector{Float64}}(undef, n_sim_steps)
    xnoms[1] = xnom0

    for t = 2:n_sim_steps

        # time span
        tspan = dt * [t-1, t]

        # Sat A (circular) propagation
        satA_state_prev = xnoms[t-1][Block(1)]
        prob_A = ODEProblem(two_body!, satA_state_prev, tspan, mu)
        sol_A  = DifferentialEquations.solve(prob_A, reltol=1e-12, abstol=1e-14)

        # Sat B (elliptical transfer) propagation
        satB_state_prev = xnoms[t-1][Block(2)]
        prob_B = ODEProblem(two_body!, satB_state_prev, tspan, mu)
        sol_B  = DifferentialEquations.solve(prob_B, reltol=1e-12, abstol=1e-14)

        # states at time t
        satA_statet = sol_A(tspan[end])  #[r;v]
        satB_statet = sol_B(tspan[end])  #[r;v]

        # add satellite state to nominal state array
        xnoms[t] = BlockArray(vcat(satA_statet, satB_statet), [6, 6])

    end

    return xnoms

end

# relevant parameters including probability of collision covariance matrix
function get_P(k, x; init = init_conds())

    # combined hard body radius (2 meters)
    HBR = 0.002 # km

    # get nominal positions and velocities
    satA_pos = init.xnoms[k][Block(1)][1:3]
    satA_vel = init.xnoms[k][Block(1)][4:6]
    satB_pos = init.xnoms[k][Block(2)][1:3]
    satB_vel = init.xnoms[k][Block(2)][4:6]
    satA_pos_true = x[Block(1)][1:3]
    satB_pos_true = x[Block(2)][1:3]

    # calculate miss distance (MD) vector (rho) in ECI km
    rho_TCA_ECI = satA_pos_true .- satB_pos_true

    ##--------------------------- COMPUTE COVARIANCE --------------------------------
    P_satA = cov_ECI(satA_pos, satA_vel, k*init.dt)
    P_satB = cov_ECI(satB_pos, satB_vel, k*init.dt)
    P      = P_satA + P_satB

    ## -------------------- COMPUTE PROBABILITY OF COLLISION ------------------------
    # First determine encounter frame (following Dr. Jones Orbital Debris notes)
    # T performs rotation from ECI to conjunction frame
    r_rel = satB_pos .- satA_pos
    v_rel = satB_vel .- satA_vel

    u_hatx = r_rel / norm(r_rel)
    u_haty = cross(r_rel, v_rel)/norm(cross(r_rel, v_rel))
    T      = vcat(u_hatx', u_haty')

    rho_2D = T * rho_TCA_ECI
    P_2D   = T * P * T'

    # Pc_Alfriend = Pc(rho_2D, P_2D, HBR)

    # # # # # # # # # # # # # # # # #
    # Logic for h constraint below  #
    # # # # # # # # # # # # # # # # #

    # Pc < 1e-4
    # exponent < ln(1e-4 * 2 * sqrt(detP) / (HBR^2))
    # rho' * P_inv * rho > -2 * ln(1e-4 * 2 * sqrt(detP) / (HBR^2))

    P_2D_inv = inv(P_2D)
    h = 2 * log(1e-4 * 2 * sqrt(det(P_2D)) / (HBR^2)) + rho_2D' * P_2D_inv * rho_2D

    return h

end

# get initial states
function get_init_states(TCA_sec; mu = 3.986e5, re = 6378.1)

    # TODO: MAKE IT SO THAT dt IS CHANGED IN DYNAMICS AND n_sim_steps, etc... in init_conds() 
    # ARE CHANGED TO REFLECT TCA (WANT n_sim_steps > TCA)

    # TCA_sec = number of seconds to get to TCA
    TCA_hours = TCA_sec/60/60

    # Initialize Parameters:
    arc_length_dist = 0.5   # distance along Sat A orbit between Sat A and Sat B at TCA
    TCA_days        = TCA_hours/24        # TCA in days
    TCA             = TCA_days*24*3600     # TCA in sec

    # # # # # # # # # # #
    # HOHMANN TRANSFER  #
    # # # # # # # # # # #

    # # Smaller Circular orbit for Hohmann
    # a1    = re + 685           # km
    # e1    = 0.0
    # inc1  = deg2rad(98.2)    # sun synchronous orbit
    # argp1 = 0.0
    # RAAN1 = 0.0

    # # Larger circ orbit for Hohmann (Aura Spacecraft (Duncan & Long paper))
    # a2    = re + 705  # km
    # e2    = 0.0
    # inc2  = inc1    # sun synchronous orbit
    # argp2 = 0.0
    # RAAN2 = 0.0

    # # Hohmann transfer ellipse
    # a_t    = (a1 + a2)/2
    # e_t    = (a2 - a1) / (a1 + a2)
    # inc_t  = inc1
    # argp_t = 0.0
    # RAAN_t = 0.0

    # # Relative phasing for collision in TCA_days
    # # Want Sat A (circular orbit) to be an arc length of 0.5km behind sat B at TCA
    # nu_meetB = pi                   # true anomaly at TCA to ensure sat B close approach
    # offset   = arc_length_dist/a2     # central angle in radius required for assigned arc length
    # nu_meetA = pi - offset          # true anomaly at TCA to ensure sat A close approach

    # # Calculate initial true anom for circular orbit
    # n2    = sqrt(mu / a2^3)                 # outer circular mean motion
    # nu2_0 = mod(nu_meetA - n2*TCA, 2*pi)   # required initial true anomaly for sc in larger circ orbit

    # # Calculate initial true anom for elliptical transfer orbit
    # # First, must compute mean anomaly at nu_meetB on the transfer orbit
    # cosE_meet = (e_t + cos(nu_meetB)) / (1 + e_t*cos(nu_meetB))
    # E_meet    = acos(cosE_meet)

    # # Account for fact that true anom, mean anom, and eccentric anom all in same half plane
    # if nu_meetB > pi
    #     E_meet = 2*pi - E_meet
    # end
    # M_meet = E_meet - e_t*sin(E_meet)
    # nt     = sqrt(mu / a_t^3)          # elliptical transfer mean motion

    # # Back out the initial mean anomaly for transfer at t=0
    # M_t0 = M_meet - nt*TCA
    # # Convert that back to true anomaly nut_0
    # nut_0 = true_anomaly_from_mean(M_t0, e_t)

    # # Get ECI states (X, Y, Z, Vx, Vy, Vz), km and km/s
    # satA_state0 = kepler2cart(a2, e2, inc2, argp2, RAAN2, nu2_0; mu=mu)
    # satB_state0 = kepler2cart(a_t, e_t, inc_t, argp_t, RAAN_t, nut_0; mu=mu)

    # # # # # # # # # # # # # # # # # # # # #
    # SAME ORBIT BUT DIFFERENT INCLINATIONS #
    # # # # # # # # # # # # # # # # # # # # #

    # Choose common circular orbit radius
    alt = 700.0
    a   = re + alt
    e   = 0.0

    # Different inclinations
    incA = deg2rad(98.2)    # Sat A
    # incA = deg2rad(109)
    incB = deg2rad(110)    # Sat B

    # Same RAAN and argument of periapsis so they intersect at line of nodes
    RAANA = 0.0
    RAANB = 0.0
    argp_A = 0.0
    argp_B = 0.0

    # Choose Time of Closest Approach (collision time)
    period = 2*pi*sqrt(a^3/mu)      # T in sec
    period_days = period/24/3600    # T in days

    # frac = 0.01667
    # TCA_days = frac * period_days   # TCA in days
    # TCA = frac*period               # TCA in sec
    tspan = (0.0, TCA)

    # Mean motion (same for both)
    n = sqrt(mu / a^3)

    # We want both at ascending node (ν = 0) at t = TCA
    # So M(TCA) = 0 mod 2π → M0 = -n*TCA mod 2π
    M0 = mod(-n * TCA, 2π)

    # Convert M0 → ν0 (for e=0, ν0 ≈ M0, but we’ll use your function)
    true_anom0 = true_anomaly_from_mean(M0, e)

    # Get states in ECI
    satA_state0 = kepler2cart(a, e, incA, argp_A, RAANA, true_anom0; mu=mu)
    satB_state0 = kepler2cart(a, e, incB, argp_B, RAANB, true_anom0; mu=mu)

    # concatenate states
    x = BlockArray(vcat(satA_state0, satB_state0), [6, 6])

    return x

end


