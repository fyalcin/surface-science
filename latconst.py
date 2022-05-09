import numpy as np
from pymatgen.analysis.eos import EOS

latvecs = np.array([[0.5, 0.5, 0], [0, 0.5, 0.5], [0.5, 0, 0.5]])
cp = np.dot(latvecs[0], np.cross(latvecs[1], latvecs[2]))

avsE = np.asarray([(3.6, -.46750841E+01),
                   (3.64285714, -.47421398E+01),
                   (3.68571429, -.47943187E+01),
                   (3.72857143, -.48328264E+01),
                   (3.77142857, -.48592685E+01),
                   (3.81428571, -.48751213E+01),
                   (3.85714286, -.48815703E+01),
                   (3.9, -.48795946E+01),
                   (3.94285714, -.48701487E+01),
                   (3.98571429, -.48538949E+01),
                   (4.02857143, -.48316337E+01),
                   (4.07142857, -.48041168E+01),
                   (4.11428571, -.47718587E+01),
                   (4.15714286, -.47349530E+01),
                   (4.2, -.46936060E+01)])

fit = np.polyfit(x=avsE[:, 0], y=avsE[:, 1], deg=3)

eos = EOS()

volumes = cp * avsE[:, 0] ** 3
energies = avsE[:, 1]

fit = eos.fit(volumes, energies)
plt = fit.plot(bbox_inches='tight')
latconst = np.round((fit.eos_params[-1] * (1 / cp)) ** (1 / 3), 4)

plt.text(0.1, 0.9, f'Lattice constant at minimum volume: {latconst}',
         transform=plt.gca().transAxes)
plt.savefig('latconst_converge.png', dpi=300, bbox_inches='tight')
