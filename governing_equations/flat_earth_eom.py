import math
import numpy as np

def flat_earth_eom(t, x, aircraft_model):
    '''
    Arguments:
        t: time [s]
        x: state vector at time t (numpy array)
            x[0] = u_b_ms
            x[1] = v_b_ms
            x[2] = w_b_ms
            x[3] = p_b_rads
            x[4] = q_b_rads
            x[5] = r_b_rads
            x[6] = phi_rad
            x[7] = theta_rad
            x[8] = psi_rad
            x[9] = x_n_m
            x[10] = y_n_m
            x[11] = z_n_m
        aircraft_model: model data stored as a dictionary containing parameters

    Returns:
        dx: time derivative of state vector
    '''

    dx = np.zeros((12,), dtype=float)

    u_b_ms = x[0]
    v_b_ms = x[1]
    w_b_ms = x[2]
    p_b_rads = x[3]
    q_b_rads = x[4]
    r_b_rads = x[5]
    phi_rad = x[6]
    theta_rad = x[7]
    psi_rad = x[8]
    x_n_m = x[9]
    y_n_m = x[10]
    z_n_m = x[11]

    # get mass and moments of inertia 
    m_kg = aircraft_model["m_kg"]
    Jxx_b_kgm2 = aircraft_model["Jxx_b_kgm2"]
    Jxz_b_kgm2 = aircraft_model["Jxz_b_kgm2"]
    Jyy_b_kgm2 = aircraft_model["Jyy_b_kgm2"]
    Jzz_b_kgm2 = aircraft_model["Jzz_b_kgm2"]

    # TODO: air data calculation(Mach, altitude, AoA, AoS)

    # TODO: atmosphere model

    # gravity
    gz_n_mps2 = 9.81

    # gravity (body-fixed)
    gx_b_mps2 = -math.sin(theta_rad) * gz_n_mps2
    gy_b_mps2 = math.sin(phi_rad) * math.cos(theta_rad) * gz_n_mps2
    gz_b_mps2 = math.cos(phi_rad) * math.cos(theta_rad) * gz_n_mps2

    # TODO: external forces
    Fx_b_kgmps2 = 0.0
    Fy_b_kgmps2 = 0.0
    Fz_b_kgmps2 = 0.0

    # TODO: external moments
    l_b_kgm2ps2 = 0.0
    m_b_kgm2ps2 = 0.0
    n_b_kgm2ps2 = 0.0

    # Denominator for rpy equations
    denom = Jxx_b_kgm2 * Jzz_b_kgm2 - Jxz_b_kgm2**2


    # compute dx

    #udot
    dx[0] = (1 / m_kg * Fx_b_kgmps2 
             + gx_b_mps2
             - w_b_ms * q_b_rads
             + v_b_ms * r_b_rads)
    #vdot 
    dx[1] = (1 / m_kg * Fy_b_kgmps2
             + gy_b_mps2
             - u_b_ms * r_b_rads
             + w_b_ms * p_b_rads)
    #wdot
    dx[2] = (1 / m_kg * Fz_b_kgmps2
             + gz_b_mps2
             - v_b_ms * p_b_rads
             + u_b_ms * q_b_rads)
    #pdot
    dx[3] = ((Jxz_b_kgm2 * (Jxx_b_kgm2 - Jyy_b_kgm2 + Jzz_b_kgm2)
             * p_b_rads * q_b_rads)
             - (Jzz_b_kgm2 * (Jzz_b_kgm2 - Jyy_b_kgm2) + Jxz_b_kgm2**2)
             * q_b_rads * r_b_rads
             + (Jzz_b_kgm2 * l_b_kgm2ps2)
             + (Jxz_b_kgm2 * n_b_kgm2ps2)) / denom
    #qdot
    dx[4] = ((Jzz_b_kgm2 - Jxx_b_kgm2) * r_b_rads * p_b_rads
             - Jxz_b_kgm2 * (p_b_rads**2 - r_b_rads**2) + m_b_kgm2ps2) / Jyy_b_kgm2

    #rdot
    dx[5] = (-Jxz_b_kgm2 * (Jxx_b_kgm2 - Jyy_b_kgm2 + Jzz_b_kgm2) * q_b_rads * r_b_rads
             + (Jxx_b_kgm2 * (Jxx_b_kgm2 - Jyy_b_kgm2) + Jxz_b_kgm2**2) * p_b_rads * q_b_rads
             + Jxz_b_kgm2 * l_b_kgm2ps2
             + Jxx_b_kgm2 * n_b_kgm2ps2) / denom

    # euler angle kinematics
    #phidot
    dx[6] = (p_b_rads 
             + q_b_rads * math.sin(phi_rad) * math.tan(theta_rad) 
             + r_b_rads * math.cos(phi_rad) * math.tan(theta_rad))
    #thetadot
    dx[7] = (q_b_rads * math.cos(phi_rad)
             - r_b_rads * math.sin(phi_rad))
    #psidot
    dx[8] = (q_b_rads * math.sin(phi_rad) / math.cos(theta_rad)
             + r_b_rads * math.cos(phi_rad) / math.cos(theta_rad))

    # TODO: Position (navigation equations)
    # Just use velocity measurements, rotate from body -> NED (123)
    c_theta = math.cos(theta_rad)
    s_theta = math.sin(theta_rad)
    c_phi = math.cos(phi_rad)
    s_phi = math.sin(phi_rad)
    c_psi = math.cos(psi_rad)
    s_psy = math.sin(psi_rad)

    dx[9] = (c_theta * c_phi * u_b_ms
             + (-c_phi * s_psi + s_phi * s_theta * c_psi) * v_b_ms
             + (s_phi * s_psi + c_phi * s_theta * c_psi) * w_b_ms)

    dx[10] = (c_phi * s_psi * u_b_ms
              + (c_phi * c_psi + s_phi * s_theta * s_psi) * v_b_ms
              + (-s_phi * c_psi + c_phi * s_theta * s+psi) * w_b_ms)

    dx[11] = (-s_theta * u_b_ms 
              + s_phi * c_theta * v_b_ms
              + c_phi * c_theta * w_b_ms)

    return dx
