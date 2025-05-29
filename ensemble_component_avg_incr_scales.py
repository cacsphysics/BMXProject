import numpy as np
import matplotlib.pyplot as plt

import openData as od
import basic_tools as bt
import compute_pdf as cp
import index_finders as infi
import cacs_library as cl
import openVel as ov

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

# index1 = find_index(time, 60)
index2 = find_index(time, 85)
index3 = find_index(time, 110)
# index4 = find_index(time, 135)
# index5 = find_index(time, 160)

delay = np.asarray([0.5, 1.0, 2.0, 5.0]) * 1e-6  # must be in units of seconds
bins = 300

# Here I split the time domain into four different regions using the indices {index1, index2, index3, index4, index5}


fig = plt.figure()
gs = fig.add_gridspec(2, 2)


axI = fig.add_subplot(gs[0, 0])
axII = fig.add_subplot(gs[0, 1], sharex=axI, sharey=axI)
axIII = fig.add_subplot(gs[1, 0], sharex=axI, sharey=axI)
axIV = fig.add_subplot(gs[1, 1], sharex=axI, sharey=axI)
ax = [axI, axII, axIII, axIV]
for count, tau in enumerate(delay):
    incr_regionII = []
    delay_index = cp.get_delay_index(tau)[0]
    for num in np.arange(0, Br5.shape[0]):

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

    incr_regionII = np.asarray(incr_regionII).flatten()
    hist, bin_edges = np.histogram(incr_regionII, bins=bins, density=True)
    bin_centers = 0.5 * (bin_edges[1:] + bin_edges[:-1])
    bin_edges = get_correct_bins(bins=bin_edges)[0]
    if count == 1:
        ax[count].semilogy(
            bin_edges, hist, "o", markersize=6, label=r"$\Delta B_s(\tau)$"
        )
    else:
        ax[count].semilogy(bin_edges, hist, "o", markersize=6)


ax[1].legend(frameon=False)

ax[1].yaxis.tick_right()
ax[3].yaxis.tick_right()
ax[0].set_ylabel("PDF")
ax[2].set_ylabel("PDF")
ax[2].set_xlabel(
    r"$\frac{\Delta(B_{avg}) - <\Delta(B_{avg})>}{\sigma_{\Delta B_{avg}}}$"
)
ax[3].set_xlabel(
    r"$\frac{\Delta(B_{avg}) - <\Delta(B_{avg})>}{\sigma_{\Delta B_{avg}}}$"
)

ax_text_list = [
    r"$\tau$ = $0.5\mu s$",
    r"$\tau$ = $1.0\mu s$",
    r"$\tau$ = $2.0\mu s$",
    r"$\tau$ = $5.0\mu s$",
]
for num, ax in enumerate(ax):
    ax.text(
        0.145,
        0.875,
        ax_text_list[num],
        horizontalalignment="center",
        verticalalignment="center",
        transform=ax.transAxes,
    )

axI.set_xlim([-6, 6])
plt.savefig(
    f"Figures/Manuscript/incr_standardized_density_ensemble_componentavg_scales.pdf"
)
plt.savefig(
    f"Figures/Manuscript/incr_standardized_density_ensemble_componentavg_scales.png"
)
# plt.show()
plt.close()
