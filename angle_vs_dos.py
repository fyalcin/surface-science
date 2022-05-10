import matplotlib.pyplot as plt
import numpy as np
from pymatgen.electronic_structure.plotter import DosPlotter
from pymatgen.io.vasp.outputs import Vasprun

plt.rcParams.update({'axes.titlesize': 'x-large'})
bulk_en_ref = -5.41378666

for angle in np.arange(0, 50, 5):
    vr = Vasprun(f"{angle}.xml")
    tdos = vr.complete_dos
    plotter = DosPlotter()

    dos_dict = tdos.get_spd_dos()
    for orbital, dos in dos_dict.items():
        plotter.add_dos(label=f"{orbital}", dos=dos)

    plt = plotter.get_plot()
    plt.title(f"DOS of Si(100) 2x1 reconstruction for u={angle} degrees")
    leg = plt.gca().get_legend()
    ltext = leg.get_texts()  # all the text.Text instance in the legend
    plt.setp(ltext, fontsize=30)
    plt.savefig(f"DOS_{angle}.png", dpi=300, bbox_inches='tight')
