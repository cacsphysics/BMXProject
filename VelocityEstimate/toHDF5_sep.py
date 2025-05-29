import h5py as h5
import open_sep_times as oz

""" The difference between the two toHDF5.py scripts is the source data.
"""


def get_data():

    P5P7_vel, P19P21_vel, P33P35_vel = oz.get_vel_arrays()

    return [P5P7_vel, P19P21_vel, P33P35_vel]


h5f = h5.File('01122022_avg_vels.h5', 'w')
h5f.attrs["date"] = 1122022
h5f.attrs["time window"] = '60 to 160 microseconds'
h5f.attrs["bandpass"] = '50 to 500kHz'

# Creating data groups
P5P7 = h5f.create_group('P5P7')
#P5P7.attrs['units'] = 'km/s'
P19P21 = h5f.create_group('P19P21')
#P19P21.attrs['units'] = 'km/s'
P33P35 = h5f.create_group('P33P35')
#P33P35.attrs['units'] = 'km/s'

P5P7_vel, P19P21_vel, P33P35_vel = get_data()

vels5 = P5P7.create_dataset('vels', data=P5P7_vel)
vels5.attrs['units'] = 'km/s'

vels19 = P19P21.create_dataset('vels', data=P19P21_vel)
vels19.attrs['units'] = 'km/s'

vels33 = P33P35.create_dataset('vels', data=P33P35_vel)
vels33.attrs['units'] = 'km/s'

h5f.close()
