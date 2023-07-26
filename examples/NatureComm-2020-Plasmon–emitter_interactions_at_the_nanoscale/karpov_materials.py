import numpy as np


def evalDrudeParams(R, WLs):
    GammaB = 3.5 * 10**13
    Gamma = GammaB + ((1.4 * 10**6) / (R * 10**-9))
    w = (2 * np.pi * 3 * 1E8) / (WLs * 1E-9)
    wp = 1.4
    return w, wp, Gamma


def eps_r(w, wp, Gamma):
    return 5 - (w**2 * (wp * 1E16)**2) / (w**4 - w**2 * Gamma**2)


def eps_i(w, wp, Gamma):
    return (w * (wp * 1E16)**2 * Gamma) / (w**4 - w**2 * Gamma**2)


def eps_drude(w, wp, Gamma):
    return eps_r(w, wp, Gamma) + 1j * eps_i(w, wp, Gamma)


def get_d_params(WL_min, WL_max, d_param_filename):
    from_disk = np.loadtxt(d_param_filename[0])
    omega_star_ratio = from_disk[0, :]
    d_perp = from_disk[1, :] + 1j*from_disk[2, :]
    from_disk = np.loadtxt(d_param_filename[1])
    d_parl = from_disk[1, :] + 1j*from_disk[2, :]

    c = 299792458  # m/s
    h_reduced = 6.5821e-16  # eV s
    omega_p_star = 3.81  # eV 

    # min_lim_omega_star_ratio = 0.87
    # max_lim_omega_star_ratio = 0.99

    d_perp_in, d_parl_in = [],[]
    WLs_nm = []
    om_rat_plot = []
    WL_d = 2*np.pi/((omega_star_ratio * omega_p_star/c)/h_reduced)*1e9  # nm
    for i in range(len(omega_star_ratio)):
        if WL_d[i] > WL_max: 
            continue
        if WL_d[i] < WL_min:
            continue
        # om_star_rat = omega_star_ratio[i]
        # if (om_star_rat < min_lim_omega_star_ratio
        #         or om_star_rat > max_lim_omega_star_ratio):
        #     continue
        # om_rat_plot.append(om_star_rat)
        WLs_nm.append(WL_d[i])
        d_perp_in.append(d_perp[i])
        d_parl_in.append(d_parl[i])

    WLs_nm = np.array(WLs_nm)
    d_perp_in = np.array(d_perp_in)
    d_parl_in = np.array(d_parl_in)
    return WLs_nm, d_perp_in, d_parl_in
