import matplotlib.pyplot as plt
import numpy as np
from pymatgen.core.structure import Structure
from pymatgen.electronic_structure.plotter import DosPlotter
from pymatgen.io.vasp.outputs import Vasprun
from triboflow.phys.shaper import Shaper
plt.rcParams.update({'axes.titlesize': 'x-large'})

struct = Structure.from_file(f"Si_recon_5.vasp")
layers = Shaper.get_layers(struct)
c_values = sorted(layers.keys())
layer_size = 2
top_c = c_values[-layer_size:]
bot_c = c_values[:layer_size]

layer_bot = [item for sublist in [layers[i] for i in bot_c] for item in sublist]
layer_top = [item for sublist in [layers[i] for i in top_c] for item in sublist]

area = Shaper.get_surface_area(struct)
num_sites = struct.num_sites

bulk_en_ref = -5.41378666
evtoj = 16.02176
energies = []
slab_ens = []
for angle in np.arange(0, 50, 5):
    vr = Vasprun(f"{angle}.xml")

    slab = Structure.from_file(f"Si_recon_{angle}.vasp")
    slab_en = vr.final_energy
    slab_ens.append((angle, slab_en))
    surfen = 0.5 * evtoj * (slab_en - num_sites * bulk_en_ref) / area
    energies.append((angle, surfen))

    tdos = vr.complete_dos
    plotter = DosPlotter()
    plotter.add_dos(label='total', dos=tdos)
    dos_dict = tdos.get_spd_dos()
    for orbital, dos in dos_dict.items():
        if not orbital.name == 'd':
            plotter.add_dos(label=f"{orbital}", dos=dos)

    plt = plotter.get_plot()
    plt.title(f"DOS of Si(100) 2x1 reconstruction for u={angle} degrees")
    leg = plt.gca().get_legend()
    ltext = leg.get_texts()  # all the text.Text instance in the legend
    plt.setp(ltext, fontsize=30)
    plt.savefig(f"DOS_{angle}.png", dpi=300, bbox_inches='tight')

slab_ens = np.asarray(slab_ens)

# energies = np.asarray(energies)
# xarr, yarr = map(np.array, zip(*energies))
# xarr_norm = (xarr - xarr.mean()) / xarr.std()
# fit = np.polyfit(xarr_norm, yarr, 6)
# f = np.poly1d(fit)
# fig, ax = plt.subplots()
# crit = f.deriv().r
# r_crit = crit[crit.imag == 0].real
# r_crit = np.round(r_crit*xarr.std() + xarr.mean(), 3)
#
# ax.scatter(xarr, energies[:, 1])
# ax.plot(xarr, f(xarr_norm), label='fit')
# ax.set_title("Surface energy vs u")
# ax.set_xlabel("asdsa")
# plt.grid(True)
# plt.xlabel("u (degrees)")
# plt.ylabel("$\\gamma$ ($J/m^2$)")
# plt.text(0.4, 0.9, f'Optimized u: {r_crit[-1]} degrees',
#          transform=plt.gca().transAxes)
# plt.legend(loc="best", shadow=True)
#
# plt.savefig('u_vs_surfen.png', dpi=300, bbox_inches='tight')
