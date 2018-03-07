import matplotlib.pyplot as plt
from astropy.table import Table
from matplotlib.ticker import MaxNLocator
from matplotlib.cm import Dark2 as Set1
import numpy as np

ngc2419 = Table.read("ngc2419-mucciarelli.csv")
ngc2808 = Table.read("ngc2808-mucciarelli.csv")

candidates = Table.read("/Users/arc/research/projects/inactive/Mg-K-Hollyweirdos/article/Total112listdr2kna.csv")

translate_label = lambda label: "[{}]".format("/".join([ea.strip() for ea in label.strip("[]").split("/")]))
all_labels = [
  ("[Mg I/Fe]", "[K I/Fe]"),
#  ("[Mg I/Fe]", "[Ca I/Fe I]"),
  ("[Mg I/Fe]", "[Sc II/Fe I]"),
#  ("[Mg I/Fe]", "[Ti I/Fe I]"),
#  ("[Mg I/Fe]", "[Ce II/Fe I]"),
  ("[Mg I/Fe]", "[Ba II/Y II]"),
  ("[Mg I/Fe]", "[Eu II/K I]")
]


fig, axes = plt.subplots(len(all_labels), 1, figsize=(3, 8))
ax = axes[0]

ax.scatter(ngc2419["mg_fe"], ngc2419["k_fe"], facecolor="k",
          label="NGC 2419")
ax.errorbar(ngc2419["mg_fe"], ngc2419["k_fe"],
            xerr=ngc2419["u_mg_fe"], yerr=ngc2419["u_k_fe"],
            fmt=None, ecolor="k", zorder=-1,
            label="_")

#ax.scatter(ngc2808["mg_fe"], ngc2808["k_fe"],
#           facecolor="#FFFFFF", marker="s", edgecolor="k", linewidth=1.5, zorder=100, label="NGC 2808")
#ax.errorbar(ngc2808["mg_fe"], ngc2808["k_fe"],
#            xerr=0.2, yerr=0.2, fmt=None, ecolor="k", zorder=50,
#            label="_")

apogee = Table.read("/Users/arc/research/projects/active/the-battery-stars/catalogs/allStar-l31c.2.fits")

qc = (apogee["SNR"] > 200) \
   * (apogee["TEFF"] > 0) \
   * (apogee["MG_FE"] > -5000) \
   * (apogee["K_FE"] > -5000) \
   * (apogee["FE_H"] > -5000) \
   * (apogee["MG_FE"] > -0.6) \
   * (apogee["K_FE"] < 1.0) \
   * ((apogee["K_FE"] + apogee["FE_H"]) > -2) \
   * (apogee["K_FE"] > -1)

ax.scatter(apogee["MG_FE"][qc], apogee["K_FE"][qc],
    s=1, c="#666666", zorder=-1, alpha=0.5, label=r"Milky Way")
ax.errorbar(apogee["MG_FE"][qc], apogee["K_FE"][qc],
    xerr=apogee["MG_FE_ERR"][qc], yerr=apogee["K_FE_ERR"][qc],
    fmt=None, ecolor="#666666", label="_", elinewidth=0.25, zorder=-100)

ax.scatter(
  candidates["cannon_alpha_m_1"], candidates["K_FE"], 
  facecolor=Set1(0), edgecolor="k", marker="^", s=50)


ax.set_xticks([-2, -1, 0, 1])
ax.set_yticks([0, 1, 2])

ax.set_xlim(-2, 1.25)
ax.set_ylim(-0.5, 2.5)



pop1 = Table.read("cohen-kirby-pop1.csv")
pop2 = Table.read("cohen-kirby-pop2.csv")


def get_values(data, x_y, fill_error=0.10):

    xlabel, ylabel = x_y.strip("[]").split("/")

    xrow = (data["element"] == xlabel)
    yrow = (data["element"] == ylabel)

    xvalues = []
    yvalues = []
    xerrors = []
    yerrors = []
    for star in data.dtype.names[1::2]:
        xvalues.append(data[star][xrow])
        xerrors.append(data["{}_error".format(star)][xrow])

        if ylabel.upper() != "H" and not ylabel.upper().startswith("FE"):
          yvalues.append(data[star][yrow])
          yerrors.append(data["{}_error".format(star)][yrow])

    xvalues = np.array(xvalues)
    xerrors = np.array(xerrors)

    yvalues = np.array(yvalues)
    yerrors = np.array(yerrors)

    xerrors[~np.isfinite(xerrors)] = fill_error
    yerrors[~np.isfinite(yerrors)] = fill_error

    if yvalues.size == 0:
        yvalues = np.zeros_like(xvalues)
        yerrors = np.zeros_like(xerrors)

    values = xvalues - yvalues
    errors = np.sqrt(xerrors**2 + yerrors**2)

    return (values.flatten(), errors.flatten())


for ax, labels in zip(axes, all_labels):

    xlabel, ylabel = labels

    for pop, color in zip((pop1, pop2), (Set1(1), Set1(2))):

        x, xerr = get_values(pop, xlabel)
        y, yerr = get_values(pop, ylabel)

        ax.scatter(x, y, c=color, edgecolor="k", s=50)
        ax.errorbar(x, y, xerr=xerr, yerr=yerr, fmt=None, zorder=-1, ecolor=color)

    ax.set_xlim(axes[0].get_xlim())
    ax.set_xticks(axes[0].get_xticks())

    if ax.is_last_row():
        ax.set_xlabel(translate_label(xlabel))
    else:
        ax.set_xticklabels([])

    if not ax.is_first_row():
        ax.tick_params(direction="inout", axis="both", bottom=True, top=True)


    ax.set_ylabel(translate_label(ylabel))
    ax.yaxis.set_major_locator(MaxNLocator(3))

    
fig.tight_layout()
fig.subplots_adjust(hspace=0)

fig.axes[1].set_yticks([0, 1])
fig.axes[2].set_yticks([0, 1])
fig.axes[3].set_yticks([-2, -1, 0])

fig.savefig("figure2.pdf", dpi=300)

