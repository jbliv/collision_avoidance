
using LinearAlgebra
#
# Provides position covariance as function of time to model uncertainty 
# evolution as satellites approach TCA.

# cov_ECI() - Given position, velocity (r,v) in ECI and time in seconds, 
# computes an estimated position covariance matrix in the ECI frame at a 
# given timein seconds prior to the Time of Closest Approach (TCA)
#

##------------------------------ COMPUTING COVARIANCE ------------------------------
#Prediction Error Polynomials from Duncan/Long Paper (specific to satellites w/ altitude
# of ~700km)
#Provide fitted estimates of mean positional prediction error as function of time 
#in days along each RTN axis
rad_sigma_m(t)   = -0.222*t^2 + 2.444*t + 1.668     #Radial
tang_sigma_m(t)  =  6.217*t^2 + 29.851*t + 5.747    #In-Track
cross_sigma_m(t) = -0.203*t^2 + 4.589*t + 0.246     #Cross-Track

# Converting to kilometers
rad_sigma_km(t)   = rad_sigma_m(t)   / 1000.0
tang_sigma_km(t)  = tang_sigma_m(t)  / 1000.0
cross_sigma_km(t) = cross_sigma_m(t) / 1000.0

##------------------------------ IMPORTANT FUNCTIONS ------------------------------
# Build 3x3 Position Covariance Matrix in km^2
function cov_RTN(t_days)
    sigma_R = rad_sigma_km(t_days)
    sigma_T = tang_sigma_km(t_days)
    sigma_N = cross_sigma_km(t_days)

    return Diagonal([sigma_R^2, sigma_T^2, sigma_N^2])
end

# Rotatinging from RTN-->ECI Frame
function RTN2ECI(r, v)
    r_hat = r / norm(r)
    h_hat = cross(r, v) / norm(cross(r, v)) #aka n hat
    t_hat = cross(h_hat, r_hat)

    return [r_hat'; t_hat'; h_hat']   # rows = R,T,N axes in ECI
end

# Transform position covariance RTN-->ECI
function cov_ECI(r, v, t_seconds)
    t_days = t_seconds / (24*3600)

    pos_cov_RTN = cov_RTN(t_days)
    C = RTN2ECI(r, v)

    return C' * pos_cov_RTN * C
end

# Approximate value for Pc assuming constant probability density
# over collision sphere
# Pc = R^2/(2*det(P)^.5) * exp(-(rho^T * P^-1 * rho)/2)
# 2D-Pc method uses a 2x2 covariance
function Pc(rho::AbstractVector, P::AbstractMatrix, HBR::Real)
    detP = det(P)
    @assert detP > 0 "Covariance matrix Σ must be positive definite"

    P_inv = inv(P)
    exponent = -0.5 * dot(rho, P_inv * rho)

    return (HBR^2) / (2 * sqrt(detP)) * exp(exponent)
end