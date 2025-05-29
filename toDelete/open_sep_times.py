import numpy as np


filename = 'separation_times_by_shot.csv'


def get_vel_arrays():
    my_data = np.transpose(np.genfromtxt(
        filename, delimiter=',', usecols=(10, 11, 12), skip_header=2))
    return [my_data[0], my_data[1], my_data[2]]
