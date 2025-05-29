import numpy as np
import matplotlib.pyplot as plt

import openData as od
import basic_tools as bt
import compute_pdf as cp
import index_finders as infi
import cacs_library as cl
import openVel as ov


norm_struct = cl.structure_function.normalized_structure_function
plt.style.use("square_presentation.mplstyle")


def determine_limits(vels, min, max):
    """Returns a truth table corresponding the (min, max)-range
    P5P7
    Inputs:
        vels - 1D numpy array
        min - float
        max -float
    Outputs:
        band_vels: Boolean array indicating if the vels are within (min, max)
    """
    vel_min_bool = np.where(vels > min, True, False)
    vel_max_bool = np.where(vels < max, True, False)

    return [vel_min_bool * vel_max_bool]


def get_B5():
    vels = ov.get_p5p7()[0]

    # removing shots with bad velocities

    vels_correct = ov.pos_finite_vels_p5p7()[0]
    Br, Bt, Bz = od.get_pos5_b()

    vels = vels[vels_correct]
    Br = Br[vels_correct]
    Bt = Bt[vels_correct]
    Bz = Bz[vels_correct]

    vels_correct_2 = determine_limits(vels, 0, 120)[0]
    Br = Br[vels_correct_2]
    Bt = Bt[vels_correct_2]
    Bz = Bz[vels_correct_2]

    # Bmag = bt.get_mag(Br, Bt, Bz)

    return Br, Bt, Bz


def get_B19():
    vels = ov.get_p19p21()[0]

    # removing shots with bad velocities

    vels_correct = ov.pos_finite_vels_p19p21()[0]
    Br, Bt, Bz = od.get_pos19_b()

    vels = vels[vels_correct]
    Br = Br[vels_correct]
    Bt = Bt[vels_correct]
    Bz = Bz[vels_correct]

    vels_correct_2 = determine_limits(vels, 0, 120)[0]
    Br = Br[vels_correct_2]
    Bt = Bt[vels_correct_2]
    Bz = Bz[vels_correct_2]

    Bmag = bt.get_mag(Br, Bt, Bz)

    return Br, Bt, Bz


def get_B33():
    vels = ov.get_p33p35()[0]

    # removing shots with bad velocities

    vels_correct = ov.pos_finite_vels_p33p35()[0]
    Br, Bt, Bz = od.get_pos33_b()

    vels = vels[vels_correct]
    Br = Br[vels_correct]
    Bt = Bt[vels_correct]
    Bz = Bz[vels_correct]

    vels_correct_2 = determine_limits(vels, 0, 120)[0]
    Br = Br[vels_correct_2]
    Bt = Bt[vels_correct_2]
    Bz = Bz[vels_correct_2]

    Bmag = bt.get_mag(Br, Bt, Bz)

    return Br, Bt, Bz


def avg_comp_incrs(x1, x2, x3, delay_index):

    inc_r = cp.get_increment(x1, x1, delay_index)[0]
    inc_t = cp.get_increment(x2, x2, delay_index)[0]
    inc_z = cp.get_increment(x3, x3, delay_index)[0]
    inc = np.asarray([inc_r, inc_t, inc_z])
    return [np.mean(inc, axis=0)]


def get_correct_bins(bins):
    """Output a list of values corresponding to the mean of bin edges."""
    bin_list = []
    for i in range(bins.shape[0] - 1):
        bin_list.append(np.mean(bins[i : i + 1]))
    bin_list = np.asarray(bin_list)
    return [bin_list]


find_index = infi.find_Index
HPF = cl.tools.HPF


Br5, Bt5, Bz5 = get_B5()
Br19, Bt19, Bz19 = get_B19()
Br33, Bt33, Bz33 = get_B33()
# B5mag = bt.get_mag(Br5, Bt5, Bz5)


time = od.get_b_time()[0]

index1 = find_index(time, 60)
index2 = find_index(time, 85)
index3 = find_index(time, 110)
index4 = find_index(time, 135)
index5 = find_index(time, 160)

delay = 0.2 * 1e-6  # must be in units of seconds
delay_array = np.linspace(0, 50, 101) * delay + 0.2 * 1e-6
# bins = 300

# Here I split the time domain into four different regions using the indices {index1, index2, index3, index4, index5}
# incr_regionI = []
# kurtosis_regionI = []
# incr_regionII = []
# incr_regionIII = []
# incr_regionIV = []
kurt_list_regionI = []
kurt_std_list_regionI = []
kurt_list_regionII = []
kurt_std_list_regionII = []
kurt_list_regionIII = []
kurt_std_list_regionIII = []
kurt_list_regionIV = []
kurt_std_list_regionIV = []
# kurt_list_rand_regionI = []
# kurt_std_list_rand_regionI = []


for dt in delay_array:

    kurt_temp_list_I = []
    kurt_temp_list_II = []
    kurt_temp_list_III = []
    kurt_temp_list_IV = []
    # kurt_temp_list_rand = []
    delay_index = cp.get_delay_index(dt)[0]

    shot_lim_list = [
        [0, 10],
        [10, 20],
        [20, 30],
        [30, 40],
        [40, 50],
        [50, 60],
        [60, 70],
        [70, 80],
        [80, 93],
    ]
    for shot_range in shot_lim_list:

        incrs_I = np.asarray([])
        incrs_II = np.asarray([])
        incrs_III = np.asarray([])
        incrs_IV = np.asarray([])

        for num in np.arange(shot_range[0], shot_range[1]):

            # Region I [60us, 85us]
            incrs_r = cp.get_increment(
                Br5[num][index1:index2], Br5[num][index1:index2], delay_index
            )[0]
            incrs_t = cp.get_increment(
                Bt5[num][index1:index2], Bt5[num][index1:index2], delay_index
            )[0]
            incrs_z = cp.get_increment(
                Bz5[num][index1:index2], Bz5[num][index1:index2], delay_index
            )[0]
            base_incrs = 1 / 3 * (incrs_r + incrs_t + incrs_z)
            standard_I = (base_incrs - np.mean(base_incrs)) / np.std(base_incrs)
            incrs_I = np.append(incrs_I, standard_I)

            # Region II [85us, 110us]
            incrs_r = cp.get_increment(
                Br5[num][index2:index3], Br5[num][index2:index3], delay_index
            )[0]
            incrs_t = cp.get_increment(
                Bt5[num][index2:index3], Bt5[num][index2:index3], delay_index
            )[0]
            incrs_z = cp.get_increment(
                Bz5[num][index2:index3], Bz5[num][index2:index3], delay_index
            )[0]
            base_incrs = 1 / 3 * (incrs_r + incrs_t + incrs_z)
            standard_II = (base_incrs - np.mean(base_incrs)) / np.std(base_incrs)
            incrs_II = np.append(incrs_II, standard_II)

            # Region III [110us, 135us]
            incrs_r = cp.get_increment(
                Br5[num][index3:index4], Br5[num][index3:index4], delay_index
            )[0]
            incrs_t = cp.get_increment(
                Bt5[num][index3:index4], Bt5[num][index3:index4], delay_index
            )[0]
            incrs_z = cp.get_increment(
                Bz5[num][index3:index4], Bz5[num][index3:index4], delay_index
            )[0]
            base_incrs = 1 / 3 * (incrs_r + incrs_t + incrs_z)
            standard_III = (base_incrs - np.mean(base_incrs)) / np.std(base_incrs)
            incrs_III = np.append(incrs_III, standard_III)

            # Region IV [135us, 160us]
            incrs_r = cp.get_increment(
                Br5[num][index4:index5], Br5[num][index4:index5], delay_index
            )[0]
            incrs_t = cp.get_increment(
                Bt5[num][index4:index5], Bt5[num][index4:index5], delay_index
            )[0]
            incrs_z = cp.get_increment(
                Bz5[num][index4:index5], Bz5[num][index4:index5], delay_index
            )[0]
            base_incrs = 1 / 3 * (incrs_r + incrs_t + incrs_z)
            standard_IV = (base_incrs - np.mean(base_incrs)) / np.std(base_incrs)
            incrs_IV = np.append(incrs_IV, standard_IV)

        kurtosis_I = norm_struct(incrs_I, power=4)
        kurt_temp_list_I.append(kurtosis_I)
        kurtosis_II = norm_struct(incrs_II, power=4)
        kurt_temp_list_II.append(kurtosis_II)
        kurtosis_III = norm_struct(incrs_III, power=4)
        kurt_temp_list_III.append(kurtosis_III)
        kurtosis_IV = norm_struct(incrs_IV, power=4)
        kurt_temp_list_IV.append(kurtosis_IV)

    kurt_temp_list_I = np.asarray(kurt_temp_list_I)
    kurt_list_regionI.append(np.mean(kurt_temp_list_I))
    kurt_std_list_regionI.append(np.std(kurt_temp_list_I))

    kurt_temp_list_II = np.asarray(kurt_temp_list_II)
    kurt_list_regionII.append(np.mean(kurt_temp_list_II))
    kurt_std_list_regionII.append(np.std(kurt_temp_list_II))

    kurt_temp_list_III = np.asarray(kurt_temp_list_III)
    kurt_list_regionIII.append(np.mean(kurt_temp_list_III))
    kurt_std_list_regionIII.append(np.std(kurt_temp_list_III))

    kurt_temp_list_IV = np.asarray(kurt_temp_list_IV)
    kurt_list_regionIV.append(np.mean(kurt_temp_list_IV))
    kurt_std_list_regionIV.append(np.std(kurt_temp_list_IV))


fig = plt.figure()
gs = fig.add_gridspec(2, 3)
ax0 = fig.add_subplot(gs[0:2, 0:2])
ax1 = fig.add_subplot(gs[0, 2])
ax2 = fig.add_subplot(gs[1, 2])

ax0.axhline(y=3, color="black", zorder=1)
ax0.errorbar(
    delay_array * 1e6,
    kurt_list_regionI,
    yerr=kurt_std_list_regionI,
    fmt="o",
    ecolor="orange",
    label=r"[$60\mu s$, $85\mu s$]",
)

ax0.errorbar(
    delay_array * 1e6,
    kurt_list_regionII,
    yerr=kurt_std_list_regionII,
    fmt="1",
    ecolor="green",
    label=r"[$85\mu s$, $110\mu s$]",
)

ax0.errorbar(
    delay_array * 1e6,
    kurt_list_regionIII,
    yerr=kurt_std_list_regionIII,
    fmt="s",
    ecolor="red",
    label=r"[$110\mu s$, $135\mu s$]",
)
ax0.errorbar(
    delay_array * 1e6,
    kurt_list_regionIV,
    yerr=kurt_std_list_regionIV,
    fmt="p",
    ecolor="purple",
    label=r"[$135\mu s$, $160\mu s$]",
)
ax0.legend(loc="center right", frameon=False)
ax0.set_ylabel("Kurtosis")
ax0.set_xlabel(r"Time Delay ($\mu s$)")

# Second plot

ax1.axhline(y=3, color="black", zorder=1)
ax1.errorbar(
    delay_array * 1e6,
    kurt_list_regionI,
    yerr=kurt_std_list_regionI,
    fmt="o",
    ecolor="orange",
    label=r"[$60\mu s$, $85\mu s$]",
)

ax1.errorbar(
    delay_array * 1e6,
    kurt_list_regionII,
    yerr=kurt_std_list_regionII,
    fmt="1",
    ecolor="green",
    label=r"[$85\mu s$, $110\mu s$]",
)

ax1.errorbar(
    delay_array * 1e6,
    kurt_list_regionIII,
    yerr=kurt_std_list_regionIII,
    fmt="s",
    ecolor="red",
    label=r"[$110\mu s$, $135\mu s$]",
)
ax1.errorbar(
    delay_array * 1e6,
    kurt_list_regionIV,
    yerr=kurt_std_list_regionIV,
    fmt="p",
    ecolor="purple",
    label=r"[$135\mu s$, $160\mu s$]",
)

ax1.set_xlim([0, 1.5])

# Third plot

ax2.axhline(y=3, color="black", zorder=1)
ax2.errorbar(
    delay_array * 1e6,
    kurt_list_regionI,
    yerr=kurt_std_list_regionI,
    fmt="o",
    ecolor="orange",
    label=r"[$60\mu s$, $85\mu s$]",
)

ax2.errorbar(
    delay_array * 1e6,
    kurt_list_regionII,
    yerr=kurt_std_list_regionII,
    fmt="1",
    ecolor="green",
    label=r"[$85\mu s$, $110\mu s$]",
)

ax2.errorbar(
    delay_array * 1e6,
    kurt_list_regionIII,
    yerr=kurt_std_list_regionIII,
    fmt="s",
    ecolor="red",
    label=r"[$110\mu s$, $135\mu s$]",
)
ax2.errorbar(
    delay_array * 1e6,
    kurt_list_regionIV,
    yerr=kurt_std_list_regionIV,
    fmt="p",
    ecolor="purple",
    label=r"[$135\mu s$, $160\mu s$]",
)
ax2.set_xlim([0, 4.5])
ax2.set_ylim([2, 3.5])


for ax in [ax1, ax2]:
    ax.yaxis.tick_right()
    ax.yaxis.set_label_position("right")
    ax.set_ylabel("Kurtosis")
    ax.set_xlabel(r"Time Delay ($\mu s$)")


ax0.text(
    0.95,
    0.9605,
    "(a)",
    horizontalalignment="center",
    verticalalignment="center",
    transform=ax0.transAxes,
)
ax1.text(
    0.9,
    0.9,
    "(b)",
    horizontalalignment="center",
    verticalalignment="center",
    transform=ax1.transAxes,
)
ax2.text(
    0.9,
    0.9,
    "(c)",
    horizontalalignment="center",
    verticalalignment="center",
    transform=ax2.transAxes,
)
plt.savefig(
    "Figures/Manuscript/Kurtosis/ensem_componentavg_kurt_10shot_ranges_three.pdf"
)
plt.savefig(
    "Figures/Manuscript/Kurtosis/ensem_componentavg_kurt_10shot_ranges_three.png"
)

# plt.show()
plt.close()
