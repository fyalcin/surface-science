import matplotlib.pyplot as plt
import numpy as np
from pymatgen.core.structure import Structure
from pymatgen.electronic_structure.plotter import DosPlotter
from pymatgen.io.vasp.outputs import Vasprun

plt.rcParams.update({'axes.titlesize': 'x-large'})
bulk_en_ref = -5.41378666

# for angle in np.arange(0, 50, 5):
#     vr = Vasprun(f"{angle}.xml")
#     tdos = vr.complete_dos
#     plotter = DosPlotter()
#     slab = tdos.structure
#
#     # dos_dict = tdos.get_spd_dos()
#     index = 0
#     dos_dict = tdos.get_site_spd_dos(slab[index])
#     for orbital, dos in dos_dict.items():
#         plotter.add_dos(label=f"{orbital}", dos=dos)
#
#     plt = plotter.get_plot()
#     plt.title(f"DOS of Si(100) 2x1 reconstruction for u={angle} degrees")
#     leg = plt.gca().get_legend()
#     ltext = leg.get_texts()  # all the text.Text instance in the legend
#     plt.setp(ltext, fontsize=30)
#     plt.savefig(f"DOS_{angle}_site_{index}.png", dpi=300, bbox_inches='tight')

plotter = DosPlotter()
index = 0
for angle in [0, 35]:
    vr = Vasprun(f"{angle}.xml")
    tdos = vr.complete_dos
    slab = tdos.structure

    plotter.add_dos(label=f'tot_{angle}', dos=tdos)
    # dos_dict = tdos.get_spd_dos()
    # #
    # # dos_dict = tdos.get_site_spd_dos(slab[index])
    # for orbital, dos in dos_dict.items():
    #     if not orbital.name == 'd':
    #         plotter.add_dos(label=f"{orbital}_{angle}", dos=dos)

plt = plotter.get_plot()
plt.title(f"DOS of Si(100) 2x1 reconstruction for u={angle} degrees from site {index}")
leg = plt.gca().get_legend()
ltext = leg.get_texts()  # all the text.Text instance in the legend
plt.setp(ltext, fontsize=30)
plt.savefig(f"DOS_compare_site_{index}.png", dpi=300, bbox_inches='tight')