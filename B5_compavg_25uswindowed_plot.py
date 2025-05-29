import numpy as np
import matplotlib.pylab as plt
import tools_basic as tb
import openData as od
import cacs_library as cl
import openVel as ov


sw = cl.spectrum_wwind.spectrum_wwind
find_Index = cl.tools.finding_Index_Time
plt.style.use("presentation.mplstyle")


def determine_limits(vels, min, max):
    """Returns a truth table corresponding the (min, max)-range
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


time = od.get_bdot_time()[0]
Brdot, Btdot, Bzdot = od.get_pos5_bdot()

Bmag = tb.compute_mag(Brdot, Btdot, Bzdot)

# I do not use the velocities to convert to wavenumber but I want to be consist with the analysis for the manuscript; on what shot numbers I use.

vels = ov.get_p5p7()[0]

# removing shots with bad velocities
vels_correct = ov.pos_finite_vels_p5p7()[0]  # All shots for the p5p7 are okay.
# print(vels_correct)
vels = vels[vels_correct]
Bmag = Bmag[vels_correct]
Brdot = Brdot[vels_correct]
Btdot = Btdot[vels_correct]
Bzdot = Bzdot[vels_correct]
# print(Bmag.shape)
# removing shots with bad velocities
vels_correct_2 = determine_limits(vels, 0, 120)[0]
vels = vels[vels_correct_2]
Bmag = Bmag[vels_correct_2]

Brdot = Brdot[vels_correct_2]
Btdot = Btdot[vels_correct_2]
Bzdot = Bzdot[vels_correct_2]

# The removal of the bad velocities does not preserve the shot number.


quarter_1_psd = []
quarter_2_psd = []
quarter_3_psd = []
quarter_4_psd = []
for num in range(0, Brdot.shape[0]):

    # First 25us
    index1 = find_Index(time * 1e-6, 60)
    index2 = find_Index(time * 1e-6, 85)
    quarter_freq, _, _, pwr_r, _, _, _, _, _ = sw(
        Brdot[num][index1:index2], time * 1e-6
    )
    quarter_freq, _, _, pwr_t, _, _, _, _, _ = sw(
        Btdot[num][index1:index2], time * 1e-6
    )
    quarter_freq, _, _, pwr_z, _, _, _, _, _ = sw(
        Bzdot[num][index1:index2], time * 1e-6
    )
    pwr_den = 1 / 3 * (pwr_r + pwr_t + pwr_z)[1:-1] / (quarter_freq[1:-1] ** 2)

    quarter_1_psd.append(pwr_den)

    # Second 25us
    index1 = find_Index(time * 1e-6, 85)
    index2 = find_Index(time * 1e-6, 110)
    quarter_freq, _, _, pwr_r, _, _, _, _, _ = sw(
        Brdot[num][index1:index2], time * 1e-6
    )
    quarter_freq, _, _, pwr_t, _, _, _, _, _ = sw(
        Btdot[num][index1:index2], time * 1e-6
    )
    quarter_freq, _, _, pwr_z, _, _, _, _, _ = sw(
        Bzdot[num][index1:index2], time * 1e-6
    )
    pwr_den = 1 / 3 * (pwr_r + pwr_t + pwr_z)[1:-1] / (quarter_freq[1:-1] ** 2)

    quarter_2_psd.append(pwr_den)
    # Third 25us
    index1 = find_Index(time * 1e-6, 110)
    index2 = find_Index(time * 1e-6, 135)
    quarter_freq, _, _, pwr_r, _, _, _, _, _ = sw(
        Brdot[num][index1:index2], time * 1e-6
    )
    quarter_freq, _, _, pwr_t, _, _, _, _, _ = sw(
        Btdot[num][index1:index2], time * 1e-6
    )
    quarter_freq, _, _, pwr_z, _, _, _, _, _ = sw(
        Bzdot[num][index1:index2], time * 1e-6
    )
    pwr_den = 1 / 3 * (pwr_r + pwr_t + pwr_z)[1:-1] / (quarter_freq[1:-1] ** 2)

    quarter_3_psd.append(pwr_den)

    # Fourth 25us
    index1 = find_Index(time * 1e-6, 135)
    index2 = find_Index(time * 1e-6, 160)
    quarter_freq, _, _, pwr_r, _, _, _, _, _ = sw(
        Brdot[num][index1:index2], time * 1e-6
    )
    quarter_freq, _, _, pwr_t, _, _, _, _, _ = sw(
        Btdot[num][index1:index2], time * 1e-6
    )
    quarter_freq, _, _, pwr_z, _, _, _, _, _ = sw(
        Bzdot[num][index1:index2], time * 1e-6
    )
    pwr_den = 1 / 3 * (pwr_r + pwr_t + pwr_z)[1:-1] / (quarter_freq[1:-1] ** 2)

    quarter_4_psd.append(pwr_den)

quarter_1_psd = np.asarray(quarter_1_psd)
quarter_1_psd = np.average(quarter_1_psd, axis=0)
quarter_2_psd = np.asarray(quarter_2_psd)
quarter_2_psd = np.average(quarter_2_psd, axis=0)
quarter_3_psd = np.asarray(quarter_3_psd)
quarter_3_psd = np.average(quarter_3_psd, axis=0)
quarter_4_psd = np.asarray(quarter_4_psd)
quarter_4_psd = np.average(quarter_4_psd, axis=0)

# If the magentic field is in Gaussian units then the outputs is ergs/cm^3


f, ax = plt.subplots()
plt.loglog(quarter_freq[1:-1], quarter_1_psd, label=r"[$60\mu s,~85\mu s$]")
plt.loglog(quarter_freq[1:-1], quarter_2_psd, label=r"[$85\mu s,~110\mu s$]")
plt.loglog(quarter_freq[1:-1], quarter_3_psd, label=r"[$110\mu s,~135\mu s$]")
plt.loglog(quarter_freq[1:-1], quarter_4_psd, label=r"[$135\mu s,~160\mu s$]")

K41 = quarter_freq[1:-1] ** (-5 / 3)  # K41 region

Ktrans = quarter_freq[1:-1] ** (-2.72)  # Transition region

Kdis = quarter_freq[1:-1] ** (-4.5)  # Dissipation region

# index1 = find_Index(quarter_freq[1:-1] * 1e-6, 2 * 1e5)  # 1.0 MHz


plt.loglog(
    quarter_freq[1:-1][:6],
    10 ** (0.75) * K41[:6] / K41[0] * quarter_2_psd[0],
    color="black",
    linestyle="--",
    alpha=0.75,
)

plt.text(
    0.1,
    0.86,
    r"$\alpha_{K}\sim - 5/3$",
    transform=ax.transAxes,
    rotation=-20,
    color="black",
)

index1 = 6
index2 = find_Index(quarter_freq[1:-1] * 1e-6, 1e6)  # 10 MHz

plt.loglog(
    quarter_freq[1:-1][index1 : index2 - 1],
    10 ** (1.5) * Ktrans[index1 : index2 - 1] / Ktrans[0] * quarter_2_psd[0],
    color="blue",
    linestyle="--",
    alpha=0.65,
)

plt.text(
    0.3,
    0.7,
    r"$\alpha_{T}\sim - 2.7$",
    transform=ax.transAxes,
    rotation=-30,
    color="blue",
)

index1 = find_Index(quarter_freq[1:-1] * 1e-6, 1e6)  # 1.0 MHz
index2 = find_Index(quarter_freq[1:-1] * 1e-6, 1e7)  # 10 MHz

plt.loglog(
    quarter_freq[1:-1][index1:index2],
    1e4 * Kdis[index1:index2] / Kdis[0] * quarter_2_psd[0],
    color="red",
    linestyle="--",
    alpha=0.65,
)

plt.text(
    0.525,
    0.425,
    r"$\alpha_{dis}\sim - 4.5$",
    transform=ax.transAxes,
    rotation=-43,
    color="red",
)

plt.axvline(
    x=1 / (2.5e-6), alpha=0.75, linestyle="-."
)  # tranistion point from high kurtosis to low kurtosis
plt.axvline(x=4.9e6, alpha=0.75, linestyle=":")  # ion cyclotron frequency
plt.ylim([0.5e-21, 0.5e-10])
plt.ylabel("PWR (arb)")
plt.xlabel("Frequency (Hz)")
plt.legend(frameon=False)
# plt.show()
plt.savefig("Figures/Manuscript/freq_psd_sum_p5_window_quaters.pdf")
plt.close()


"""axI = fig.add_subplot(gs[0, 0])
axII = fig.add_subplot(gs[1, 0])

axI.set_ylabel(r"$f^{2.98}$ PWR ($ergs/cm^3$)")
# axI.set_xlabel(f"frequency ($Hz$)")
axI.loglog(
    quarter_freq[1:-1] * 10 ** (-6),
    quarter_freq[1:-1] ** (factI) * full_psd[1:-1],
    label=r"[$60\mu s$, $160\mu s$]",
    zorder=2,
)
axI.loglog(
    quart_freq[1:-1] * 10 ** (-6),
    quart_freq[1:-1] ** factI * quarter_1_psd[1:-1],
    label=r"[$60\mu s$, $85\mu s$]",
    zorder=2,
)
axI.loglog(
    quart_freq[1:-1] * 10 ** (-6),
    quart_freq[1:-1] ** factI * quarter_2_psd[1:-1],
    label=r"[$85\mu s$, $110\mu s$]",
    zorder=2,
)
axI.loglog(
    quart_freq[1:-1] * 10 ** (-6),
    quart_freq[1:-1] ** factI * quarter_3_psd[1:-1],
    label=r"[$110\mu s$, $135\mu s$]",
    zorder=2,
)
axI.loglog(
    quart_freq[1:-1] * 10 ** (-6),
    quart_freq[1:-1] ** factI * quarter_4_psd[1:-1],
    label=r"[$135\mu s$, $160\mu s$]",
    zorder=2,
)


axII.set_ylabel(r"$f^{4.83}$ PWR ($ergs/cm^3$)")
axII.set_xlabel(f"frequency ($MHz$)")
axII.loglog(
    quarter_freq[1:-1] * 10 ** (-6),
    quarter_freq[1:-1] ** (factII) * full_psd[1:-1],
    label=r"[$60\mu s$, $160\mu s$]",
    zorder=2,
)
axII.loglog(
    quart_freq[1:-1] * 10 ** (-6),
    quart_freq[1:-1] ** factII * quarter_1_psd[1:-1],
    label=r"[$60\mu s$, $85\mu s$]",
    zorder=2,
)
axII.loglog(
    quart_freq[1:-1] * 10 ** (-6),
    quart_freq[1:-1] ** factII * quarter_2_psd[1:-1],
    label=r"[$85\mu s$, $110\mu s$]",
    zorder=2,
)
axII.loglog(
    quart_freq[1:-1] * 10 ** (-6),
    quart_freq[1:-1] ** factII * quarter_3_psd[1:-1],
    label=r"[$110\mu s$, $135\mu s$]",
    zorder=2,
)
axII.loglog(
    quart_freq[1:-1] * 10 ** (-6),
    quart_freq[1:-1] ** factII * quarter_4_psd[1:-1],
    label=r"[$135\mu s$, $160\mu s$]",
    zorder=2,
)

# index1 = find_Index(quart_freq * 1e-6, 0.2e6)
# index2 = find_Index(quart_freq * 1e-6, 20e6)

# f5 = 10 ** (2) * quart_freq[index1:index2] ** (-5 / 3)

# plt.loglog(
#    quart_freq[index1:index2],
#    quart_freq[index1:index2] ** fact * f5 * 10**8,
#    color="tab:blue",
#    linestyle="--",
# )

# plt.text(
#    1.25e6,
#    1.25,
#    rf"$\alpha\sim 5/3$",
#    fontsize=16,
#    horizontalalignment="center",
#    verticalalignment="center",
#    rotation=-30.5,
#    color="tab:blue",
# )


axI.legend(frameon=False)
for ax in [axI, axII]:
    ax.minorticks_on()
    ax.tick_params(axis="both", which="both", direction="in")
    ax.axvline(x=5, zorder=1, color="black")  # Cyclotron Frequency (cyc/sec)
    ax.axvline(
        x=(1 / 0.62), zorder=1, color="black", linestyle=":"
    )  # 1/Tau_c in manuscript; zero-time delay curvature.
    ax.grid(which="both")
axI.set_xticklabels("")
axI.axvspan(0.3, 0.85, alpha=0.5, color="orange")
axII.axvspan(2, 5, alpha=0.5, color="orange")

# plt.show()
plt.savefig(f"Figures/QuickPlots/P5/pwr_25uswindows_compensatedboth_shotavg.pdf")
plt.close()
"""
