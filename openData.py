import h5py as h5


filename = "Dataset_01122022.h5"


def get_b_time():
    """Function that outputs the time (microsec) data.
    Input:
        None
    Output:
        time - 1D numpy array
    """
    dfile = h5.File(filename, "r")
    time = dfile["time/time_us"][()]
    dfile.close()

    return [time[1:]]


def get_bdot_time():
    """Function that outputs the time (microsec) data.
    Input:
        None
    Output:
        time - 1D numpy array
    """
    dfile = h5.File(filename, "r")
    time = dfile["time/time_us"][()]
    dfile.close()

    return [time]


def get_pos5_b():
    """Function that outputs the magnetic fields of P5.
    Input:
        None
    Output:
        Br, Bt, Bz - (X, X, X) all 2D numpy arrays
    """
    dfile = h5.File(filename, "r")
    Br = dfile["mag_probe/positions/probe5/r/b"][()]
    Bt = dfile["mag_probe/positions/probe5/t/b"][()]
    Bz = dfile["mag_probe/positions/probe5/z/b"][()]
    dfile.close()

    return [Br, Bt, Bz]


def get_pos5_bdot():
    """Function that outputs the bdot fields of P5.
    Input:
        None
    Output:
        Brdot, Btdot, Bzdot - (X, X, X) all 2D numpy arrays
    """
    dfile = h5.File(filename, "r")
    Brdot = dfile["mag_probe/positions/probe5/r/bdot"][()]
    Btdot = dfile["mag_probe/positions/probe5/t/bdot"][()]
    Bzdot = dfile["mag_probe/positions/probe5/z/bdot"][()]
    dfile.close()

    return [Brdot, Btdot, Bzdot]


def get_pos7_b():
    """Function that outputs the magnetic fields of P7.
    Input:
        None
    Output:
        Br, Bt, Bz - (X, X, X) all 2D numpy arrays
    """
    dfile = h5.File(filename, "r")
    Br = dfile["mag_probe/positions/probe7/r/b"][()]
    Bt = dfile["mag_probe/positions/probe7/t/b"][()]
    Bz = dfile["mag_probe/positions/probe7/z/b"][()]
    dfile.close()

    return [Br, Bt, Bz]


def get_pos7_bdot():
    """Function that outputs the bdot fields of P7.
    Input:
        None
    Output:
        Brdot, Btdot, Bzdot - (X, X, X) all 2D numpy arrays
    """
    dfile = h5.File(filename, "r")
    Brdot = dfile["mag_probe/positions/probe7/r/bdot"][()]
    Btdot = dfile["mag_probe/positions/probe7/t/bdot"][()]
    Bzdot = dfile["mag_probe/positions/probe7/z/bdot"][()]
    dfile.close()

    return [Brdot, Btdot, Bzdot]


def get_pos19_b():
    """Function that outputs the magnetic fields of P19.
    Input:
        None
    Output:
        Br, Bt, Bz - (X, X, X) all 2D numpy arrays
    """
    dfile = h5.File(filename, "r")
    Br = dfile["mag_probe/positions/probe19/r/b"][()]
    Bt = dfile["mag_probe/positions/probe19/t/b"][()]
    Bz = dfile["mag_probe/positions/probe19/z/b"][()]
    dfile.close()

    return [Br, Bt, Bz]


def get_pos19_bdot():
    """Function that outputs the bdot fields of P19.
    Input:
        None
    Output:
        Brdot, Btdot, Bzdot - (X, X, X) all 2D numpy arrays
    """
    dfile = h5.File(filename, "r")
    Brdot = dfile["mag_probe/positions/probe19/r/bdot"][()]
    Btdot = dfile["mag_probe/positions/probe19/t/bdot"][()]
    Bzdot = dfile["mag_probe/positions/probe19/z/bdot"][()]
    dfile.close()

    return [Brdot, Btdot, Bzdot]


def get_pos21_b():
    """Function that outputs the magnetic fields of P21.
    Input:
        None
    Output:
        Br, Bt, Bz - (X, X, X) all 2D numpy arrays
    """
    dfile = h5.File(filename, "r")
    Br = dfile["mag_probe/positions/probe21/r/b"][()]
    Bt = dfile["mag_probe/positions/probe21/t/b"][()]
    Bz = dfile["mag_probe/positions/probe21/z/b"][()]
    dfile.close()

    return [Br, Bt, Bz]


def get_pos21_bdot():
    """Function that outputs the bdot fields of P21.
    Input:
        None
    Output:
        Brdot, Btdot, Bzdot - (X, X, X) all 2D numpy arrays
    """
    dfile = h5.File(filename, "r")
    Brdot = dfile["mag_probe/positions/probe21/r/bdot"][()]
    Btdot = dfile["mag_probe/positions/probe21/t/bdot"][()]
    Bzdot = dfile["mag_probe/positions/probe21/z/bdot"][()]
    dfile.close()

    return [Brdot, Btdot, Bzdot]


def get_pos33_b():
    """Function that outputs the magnetic fields of P33.
    Input:
        None
    Output:
        Br, Bt, Bz - (X, X, X) all 2D numpy arrays
    """
    dfile = h5.File(filename, "r")
    Br = dfile["mag_probe/positions/probe33/r/b"][()]
    Bt = dfile["mag_probe/positions/probe33/t/b"][()]
    Bz = dfile["mag_probe/positions/probe33/z/b"][()]
    dfile.close()

    return [Br, Bt, Bz]


def get_pos33_bdot():
    """Function that outputs the bdot fields of P33.
    Input:
        None
    Output:
        Brdot, Btdot, Bzdot - (X, X, X) all 2D numpy arrays
    """
    dfile = h5.File(filename, "r")
    Brdot = dfile["mag_probe/positions/probe33/r/bdot"][()]
    Btdot = dfile["mag_probe/positions/probe33/t/bdot"][()]
    Bzdot = dfile["mag_probe/positions/probe33/z/bdot"][()]
    dfile.close()

    return [Brdot, Btdot, Bzdot]


def get_pos35_b():
    """Function that outputs the magnetic fields of P35.
    Input:
        None
    Output:
        Br, Bt, Bz - (X, X, X) all 2D numpy arrays
    """
    dfile = h5.File(filename, "r")
    Br = dfile["mag_probe/positions/probe35/r/b"][()]
    Bt = dfile["mag_probe/positions/probe35/t/b"][()]
    Bz = dfile["mag_probe/positions/probe35/z/b"][()]
    dfile.close()

    return [Br, Bt, Bz]


def get_pos35_bdot():
    """Function that outputs the bdot fields of P35.
    Input:
        None
    Output:
        Brdot, Btdot, Bzdot - (X, X, X) all 2D numpy arrays
    """
    dfile = h5.File(filename, "r")
    Brdot = dfile["mag_probe/positions/probe35/r/bdot"][()]
    Btdot = dfile["mag_probe/positions/probe35/t/bdot"][()]
    Bzdot = dfile["mag_probe/positions/probe35/z/bdot"][()]
    dfile.close()

    return [Brdot, Btdot, Bzdot]


def get_discharge_data():
    """Function that outputs the discharge data.
    Input:
        None
    Output:
        Voltage, Current - (X, X) all 2D numpy arrays
    """
    dfile = h5.File(filename, "r")
    volt = -dfile["discharge/dis_V"][()] * 1e-3
    curr = dfile["discharge/dis_I"][()] * 1e-2
    dfile.close()

    return [volt, curr]
