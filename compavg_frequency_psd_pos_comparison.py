import numpy as np
import matplotlib.pylab as plt
import tools_basic as tb
import openVel as ov
import openData as od
import cacs_library as cl

sw = cl.spectrum_wwind.spectrum_wwind
find_Index = cl.tools.finding_Index_Time
plt.style.use("presentation.mplstyle")


def determine_limits(vels, min, max):
    """Returns a truth table corresponding the (min, max)-range
    Inputs:
        vels - 3D numpy array
        min - float
        max -float
    Outputs:
        band_vels: Boolean array indicating if the vels are within (min, max)
    """
    vel_min_bool = np.where(vels > min, True, False)
    vel_max_bool = np.where(vels < max, True, False)

    return [vel_min_bool * vel_max_bool]


time = od.get_bdot_time()[0]

### Position 5
Brdot, Btdot, Bzdot = od.get_pos5_bdot()

index1 = find_Index(time * 1e-6, 60)
index2 = find_Index(time * 1e-6, 160)

vels = ov.get_p5p7()[0]

# removing shots with bad velocities

vels_correct = ov.pos_finite_vels_p5p7()[0]

vels = vels[vels_correct]
Brdot = Brdot[vels_correct]
Btdot = Btdot[vels_correct]
Bzdot = Bzdot[vels_correct]

vels_correct_2 = determine_limits(vels, 0, 120)[0]
vels = vels[vels_correct_2]
Brdot = Brdot[vels_correct_2]
Btdot = Btdot[vels_correct_2]
Bzdot = Bzdot[vels_correct_2]


for num in range(0, Brdot.shape[0]):
    # Computing Br-fft
    freq, _, _, pwr_rden_dot, _, _, _, _, _ = sw(
        Brdot[num][index1:index2], time[index1:index2] * 1e-6
    )
    # Computing Bt-fft
    freq, _, _, pwr_tden_dot, _, _, _, _, _ = sw(
        Btdot[num][index1:index2], time[index1:index2] * 1e-6
    )
    # Computing Bz-fft
    freq, _, _, pwr_zden_dot, _, _, _, _, _ = sw(
        Bzdot[num][index1:index2], time[index1:index2] * 1e-6
    )

    pwr_rden_dot = pwr_rden_dot[1:-2] / freq[1:-2] ** 2
    pwr_tden_dot = pwr_tden_dot[1:-2] / freq[1:-2] ** 2
    pwr_zden_dot = pwr_zden_dot[1:-2] / freq[1:-2] ** 2

    pwr_den_dot = (pwr_rden_dot + pwr_tden_dot + pwr_zden_dot) / 3
    if num == 0:
        avg_P5_psd_dot = pwr_den_dot / Brdot.shape[0]
    else:
        avg_P5_psd_dot += pwr_den_dot / Brdot.shape[0]

# Pos 19
Brdot, Btdot, Bzdot = od.get_pos19_bdot()


vels = ov.get_p19p21()[0]

# removing shots with bad velocities

vels_correct = ov.pos_finite_vels_p19p21()[0]

vels = vels[vels_correct]
Brdot = Brdot[vels_correct]
Btdot = Btdot[vels_correct]
Bzdot = Bzdot[vels_correct]

vels_correct_2 = determine_limits(vels, 0, 120)[0]
vels = vels[vels_correct_2]
Brdot = Brdot[vels_correct_2]
Btdot = Btdot[vels_correct_2]
Bzdot = Bzdot[vels_correct_2]


for num in range(0, Brdot.shape[0]):
    # Computing Br-fft
    freq, _, _, pwr_rden_dot, _, _, _, _, _ = sw(
        Brdot[num][index1:index2], time[index1:index2] * 1e-6
    )
    # Computing Bt-fft
    freq, _, _, pwr_tden_dot, _, _, _, _, _ = sw(
        Btdot[num][index1:index2], time[index1:index2] * 1e-6
    )
    # Computing Bz-fft
    freq, _, _, pwr_zden_dot, _, _, _, _, _ = sw(
        Bzdot[num][index1:index2], time[index1:index2] * 1e-6
    )
    pwr_rden_dot = pwr_rden_dot[1:-2] / freq[1:-2] ** 2
    pwr_tden_dot = pwr_tden_dot[1:-2] / freq[1:-2] ** 2
    pwr_zden_dot = pwr_zden_dot[1:-2] / freq[1:-2] ** 2

    pwr_den_dot = (pwr_rden_dot + pwr_tden_dot + pwr_zden_dot) / 3
    if num == 0:
        avg_P19_psd_dot = pwr_den_dot / Brdot.shape[0]
    else:
        avg_P19_psd_dot += pwr_den_dot / Brdot.shape[0]

# Position 33

Brdot, Btdot, Bzdot = od.get_pos33_bdot()


vels = ov.get_p33p35()[0]

# removing shots with bad velocities

vels_correct = ov.pos_finite_vels_p33p35()[0]

vels = vels[vels_correct]
Brdot = Brdot[vels_correct]
Btdot = Btdot[vels_correct]
Bzdot = Bzdot[vels_correct]

vels_correct_2 = determine_limits(vels, 0, 120)[0]
vels = vels[vels_correct_2]
Brdot = Brdot[vels_correct_2]
Btdot = Btdot[vels_correct_2]
Bzdot = Bzdot[vels_correct_2]


for num in range(0, Brdot.shape[0]):
    # Computing Br-fft
    freq, _, _, pwr_rden_dot, _, _, _, _, _ = sw(
        Brdot[num][index1:index2], time[index1:index2] * 1e-6
    )
    # Computing Bt-fft
    freq, _, _, pwr_tden_dot, _, _, _, _, _ = sw(
        Btdot[num][index1:index2], time[index1:index2] * 1e-6
    )
    # Computing Bz-fft
    freq, _, _, pwr_zden_dot, _, _, _, _, _ = sw(
        Bzdot[num][index1:index2], time[index1:index2] * 1e-6
    )

    pwr_rden_dot = pwr_rden_dot[1:-2] / freq[1:-2] ** 2
    pwr_tden_dot = pwr_tden_dot[1:-2] / freq[1:-2] ** 2
    pwr_zden_dot = pwr_zden_dot[1:-2] / freq[1:-2] ** 2

    pwr_den_dot = (pwr_rden_dot + pwr_tden_dot + pwr_zden_dot) / 3
    if num == 0:
        avg_P33_psd_dot = pwr_den_dot / Brdot.shape[0]
    else:
        avg_P33_psd_dot += pwr_den_dot / Brdot.shape[0]

fig = plt.figure()
# fig = plt.figure(dpi=600)
gs = fig.add_gridspec(1, 1)

ax = fig.add_subplot(gs[0, 0])


ax.loglog(freq[1:-2], avg_P5_psd_dot, label="z = 41.6cm")
ax.loglog(freq[1:-2], avg_P19_psd_dot, label="z = 59.8cm")
ax.loglog(freq[1:-2], avg_P33_psd_dot, label="z = 78.0cm")

K41 = freq[1:-2] ** (-5 / 3)
ax.loglog(
    freq[1:-2],
    1e1 * K41 / K41[0] * np.max(avg_P5_psd_dot),
    color="black",
    label=r"K41 ($f^{-5/3}$)",
    linestyle="--",
    alpha=0.75,
    zorder=1,
)


# ax.axvline(x=5e4, color="grey", linestyle=":") # 20us windows
ax.axvline(x=6.5e4, color="red", zorder=1)  # feature frequency
ax.axvline(x=4.8e6, color="grey", linestyle=":")  # cyclotron frequency
# ax.axvline(x=1 / (25e-6), color="grey", linestyle=":")  # 25us windows
# ax.set_ylim([10**-8, 10**8])

ax.set_ylabel(r"PSD (arb.)")
ax.set_xlabel(r"Frequency ($Hz$)")
ax.legend(frameon=False)
# ax.legend(frameon=False)
plt.show()
# plt.savefig(f"Figures/Manuscript/psd_frequency_compavg_comparison.png")
# plt.savefig(f"Figures/Manuscript/psd_frequency_compavg_comparison.pdf")
plt.close()
