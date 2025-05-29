import numpy as np
import matplotlib.pyplot as plt

import openData as od
import basic_tools as bt
import compute_pdf as cp
import index_finders as infi
import cacs_library as cl
import openVel as ov

plt.style.use("tall_presentation.mplstyle")


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
# B5mag = bt.get_mag(Br5, Bt5, Bz5)


time = od.get_b_time()[0]

index1 = find_index(time, 60)
index2 = find_index(time, 85)
index3 = find_index(time, 110)
index4 = find_index(time, 135)
index5 = find_index(time, 160)

delay = 0.5 * 1e-6  # must be in units of seconds
delay_index = cp.get_delay_index(delay)[0]
bins = 300

# Here I split the time domain into four different regions using the indices {index1, index2, index3, index4, index5}
incr_regionI = []
incr_regionII = []
incr_regionIII = []
incr_regionIV = []


for num in np.arange(0, Br5.shape[0]):
    # B5mag = HPF(Br5, 5e4)

    incrs_r = cp.get_increment(
        Br5[num][index1:index2], Br5[num][index1:index2], delay_index
    )[0]

    incrs_t = cp.get_increment(
        Bt5[num][index1:index2], Bt5[num][index1:index2], delay_index
    )[0]

    incrs_z = cp.get_increment(
        Bz5[num][index1:index2], Bz5[num][index1:index2], delay_index
    )[0]

    incrs = (incrs_r + incrs_t + incrs_z) / 3

    standard = (incrs - np.mean(incrs)) / np.std(incrs)
    incr_regionI.append(standard)

    # Why am I computing the histogram and the bins?
    # hist, bin_edges = np.histogram(standard, bins=bins, density=True)

    incrs_r = cp.get_increment(
        Br5[num][index2:index3], Br5[num][index2:index3], delay_index
    )[0]

    incrs_t = cp.get_increment(
        Bt5[num][index2:index3], Bt5[num][index2:index3], delay_index
    )[0]

    incrs_z = cp.get_increment(
        Bz5[num][index2:index3], Bz5[num][index2:index3], delay_index
    )[0]

    incrs = (incrs_r + incrs_t + incrs_z) / 3

    standard = (incrs - np.mean(incrs)) / np.std(incrs)
    incr_regionII.append(standard)

    incrs_r = cp.get_increment(
        Br5[num][index3:index4], Br5[num][index3:index4], delay_index
    )[0]

    incrs_t = cp.get_increment(
        Bt5[num][index3:index4], Bt5[num][index3:index4], delay_index
    )[0]

    incrs_z = cp.get_increment(
        Bz5[num][index3:index4], Bz5[num][index3:index4], delay_index
    )[0]

    incrs = (incrs_r + incrs_t + incrs_z) / 3

    standard = (incrs - np.mean(incrs)) / np.std(incrs)
    incr_regionIII.append(standard)

    incrs_r = cp.get_increment(
        Br5[num][index4:index5], Br5[num][index4:index5], delay_index
    )[0]

    incrs_t = cp.get_increment(
        Bt5[num][index4:index5], Bt5[num][index4:index5], delay_index
    )[0]

    incrs_z = cp.get_increment(
        Bz5[num][index4:index5], Bz5[num][index4:index5], delay_index
    )[0]

    incrs = (incrs_r + incrs_t + incrs_z) / 3

    standard = (incrs - np.mean(incrs)) / np.std(incrs)
    incr_regionIV.append(standard)

plt.figure()
# plt.title(r"$\Delta \tau$ = $0.1\mu s$")

incr_regionI = np.asarray(incr_regionI).flatten()
incr_regionII = np.asarray(incr_regionII).flatten()
incr_regionIII = np.asarray(incr_regionIII).flatten()
incr_regionIV = np.asarray(incr_regionIV).flatten()

hist, bin_edges = np.histogram(incr_regionI, bins=bins, density=True)
bin_centers = 0.5 * (bin_edges[1:] + bin_edges[:-1])
bin_edges = get_correct_bins(bins=bin_edges)[0]
plt.semilogy(
    bin_edges, hist / 0.000001, "o", markersize=4, label=r"[60$\mu s$, 85$\mu s$]"
)

hist, bin_edges = np.histogram(incr_regionII, bins=bins, density=True)
bin_centers = 0.5 * (bin_edges[1:] + bin_edges[:-1])
bin_edges = get_correct_bins(bins=bin_edges)[0]
plt.semilogy(
    bin_edges, hist / 0.0001, "o", markersize=4, label=r"[85$\mu s$, 110$\mu s$]"
)

hist, bin_edges = np.histogram(incr_regionIII, bins=bins, density=True)
bin_centers = 0.5 * (bin_edges[1:] + bin_edges[:-1])
bin_edges = get_correct_bins(bins=bin_edges)[0]
plt.semilogy(
    bin_edges, hist / 0.01, "o", markersize=4, label=r"[110$\mu s$, 135$\mu s$]"
)

hist, bin_edges = np.histogram(incr_regionIV, bins=bins, density=True)
bin_centers = 0.5 * (bin_edges[1:] + bin_edges[:-1])
bin_edges = get_correct_bins(bins=bin_edges)[0]
plt.semilogy(bin_edges, hist, "o", markersize=4, label=r"[135$\mu s$, 160$\mu s$]")

plt.xlim([-6, 6])
plt.text(
    -4,
    1e5,
    r"$\tau = 0.5\mu s$",
    horizontalalignment="center",
    verticalalignment="center",
)
plt.legend(frameon=False)
plt.ylabel("PDF")
plt.xlabel(r"$\frac{\Delta(B_{avg}) - <\Delta(B_{avg})>}{\sigma_{\Delta B_{avg}}}$")
plt.savefig(
    f"Figures/Manuscript/incr_standardized_density_ensemble_componentavg_0.5us.pdf"
)
plt.savefig(
    f"Figures/Manuscript/incr_standardized_density_ensemble_componentavg_0.5us.png"
)
# plt.show()
plt.close()
