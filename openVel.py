import h5py as h5
import numpy as np
filename = '01122022_avg_vels.h5'


def pos_finite_vels_p5p7():
    """ Returns a truth table corresponding to shots that are postive definite at
        P5P7
        Inputs:
            None
        Outputs:
            positive_finite_vels: Boolean array indicating if the vels is positive and finite.
    """
    vels = get_p5p7()[0]
    vel_plus_bool = np.where(vels > 0, True, False)
    vel_finite_bool = np.isfinite(vels)

    return [vel_plus_bool*vel_finite_bool]


def pos_finite_vels_p19p21():
    """ Returns a truth table corresponding to shots that are postive definite at
        P19P21
        Inputs:
            None
        Outputs:
            positive_finite_vels: Boolean array indicating if the vels is positive and finite.
    """
    vels = get_p19p21()[0]
    vel_plus_bool = np.where(vels >= 0, True, False)
    vel_finite_bool = np.isfinite(vels)

    return [vel_plus_bool*vel_finite_bool]


def pos_finite_vels_p33p35():
    """ Returns a truth table corresponding to shots that are postive definite at
        P33P35
        Inputs:
            None
        Outputs:
            positive_finite_vels: Boolean array indicating if the vels is positive and finite.
    """
    vels = get_p33p35()[0]
    vel_plus_bool = np.where(vels >= 0, True, False)
    vel_finite_bool = np.isfinite(vels)

    return [vel_plus_bool*vel_finite_bool]


def get_p5p7():
    """ Returns the time-averaged velocities at P5P7
    Inputs:
        None
    Outputs:
        vels (N,) numpy array (units km/s)
    """
    datFil = h5.File(filename, 'r')
    vels = datFil['P5P7/vels'][()]
    datFil.close()

    return [vels]


def get_p19p21():
    """ Returns the time-averaged velocities at P19P21
    Inputs:
        None
    Outputs:
        vels (N,) numpy array (units km/s)
    """
    datFil = h5.File(filename, 'r')
    vels = datFil['P19P21/vels'][()]
    datFil.close()

    return [vels]


def get_p33p35():
    """ Returns the time-averaged velocities at P33P35
    Inputs:
        None
    Outputs:
        vels (N,) numpy array (units km/s)
    """
    datFil = h5.File(filename, 'r')
    vels = datFil['P33P35/vels'][()]
    datFil.close()

    return [vels]
