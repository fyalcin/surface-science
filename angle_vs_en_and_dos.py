import matplotlib.pyplot as plt
import numpy as np
from pymatgen.electronic_structure.plotter import DosPlotter
from pymatgen.io.vasp.outputs import Vasprun


energies = []
for angle in np.arange(0, 50, 5):
    plotter = DosPlotter()
    vr = Vasprun(f"{angle}.xml")
    energies.append((angle, vr.final_energy))
    cdos = vr.complete_dos

    plotter.add_dos(f'{angle}_tdos', cdos)

    plotter.show()

energies = np.asarray(energies)
xarr, yarr = map(np.array, zip(*energies))
xarr = (xarr - xarr.mean()) / xarr.std()
fit = np.polyfit(xarr, yarr, 6)
f = np.poly1d(fit)
fig, ax = plt.subplots()
crit = f.deriv().r
r_crit = crit[crit.imag == 0].real

ax.scatter(xarr, energies[:, 1])
ax.plot(xarr, f(xarr), label='fit')

plt.show()
