import numpy as np
from pymatgen.analysis.eos import EOS
from pymatgen.core.structure import Structure
from pymatgen.transformations.standard_transformations import DeformStructureTransformation

ouc = Structure.from_file('Si_ouc.vasp')
orig_volume = ouc.volume

latvecs = np.array([[0.5, 0.5, 0], [0, 0.5, 0.5], [0.5, 0, 0.5]])
cp = np.dot(latvecs[0], np.cross(latvecs[1], latvecs[2]))

# avsE = np.asarray([(3.6, -.46750841E+01),
#                    (3.62857143, -.47214613E+01),
#                    (3.65714286, -.47611578E+01),
#                    (3.6857142, -.47943186E+01),
#                    (3.71428571, -.48214271E+01),
#                    (3.74285714, -.48429245E+01),
#                    (3.77142857, -.48592685E+01),
#                    (3.8, -.48710310E+01),
#                    (3.82857143, -.48783123E+01),
#                    (3.85714286, -.48815703E+01),
#                    (3.88571429, -.48811563E+01),
#                    (3.91428571, -.48772424E+01),
#                    (3.94285714, -.48701487E+01),
#                    (3.97142857, -.48600275E+01),
#                    (4., -.48470712E+01)])
avsE = np.asarray([(0.85, -.14098878E+02),
                   (0.86578947, -.15922672E+02),
                   (0.88157895, -.17429093E+02),
                   (0.89736842, -.18653624E+02),
                   (0.91315789, -.19628393E+02),
                   (0.92894737, -.20381982E+02),
                   (0.94473684, -.20940706E+02),
                   (0.96052632, -.21328818E+02),
                   (0.97631579, -.21568465E+02),
                   (0.99210526, -.21679578E+02),
                   (1.00789474, -.21680243E+02),
                   (1.02368421, -.21586827E+02),
                   (1.03947368, -.21413893E+02),
                   (1.05526316, -.21174435E+02),
                   (1.07105263, -.20879995E+02),
                   (1.08684211, -.20540724E+02),
                   (1.10263158, -.20165767E+02),
                   (1.11842105, -.19762988E+02),
                   (1.13421053, -.19339070E+02),
                   (1.15, -.18899604E+02)])
volumes = []
latconst_vol = {}
for deform in avsE:
    scale_factor = deform[0]
    scale_array = [(scale_factor, 0, 0), (0, scale_factor, 0), (0, 0, scale_factor)]
    deformer = DeformStructureTransformation(deformation=scale_array)
    deformed_ouc = deformer.apply_transformation(ouc.copy())
    volume = deformed_ouc.volume
    volumes.append(volume)
    latconst_vol[deformed_ouc.lattice.a] = volume

# avsE = np.asarray([(3.6, -.46750841E+01),
#                    (3.64285714, -.47421398E+01),
#                    (3.68571429, -.47943187E+01),
#                    (3.72857143, -.48328264E+01),
#                    (3.77142857, -.48592685E+01),
#                    (3.81428571, -.48751213E+01),
#                    (3.85714286, -.48815703E+01),
#                    (3.9, -.48795946E+01),
#                    (3.94285714, -.48701487E+01),
#                    (3.98571429, -.48538949E+01),
#                    (4.02857143, -.48316337E+01),
#                    (4.07142857, -.48041168E+01),
#                    (4.11428571, -.47718587E+01),
#                    (4.15714286, -.47349530E+01),
#                    (4.2, -.46936060E+01)])

eos = EOS()

# volumes = cp * avsE[:, 0] ** 3
energies = avsE[:, 1]

fit = eos.fit(volumes, energies)
plt = fit.plot(bbox_inches='tight')
latconst = np.round(ouc.lattice.a*(fit.eos_params[-1]/orig_volume)**(1/3), 4)

plt.text(0.1, 0.9, f'Optimized lattice constant: {latconst}A',
         transform=plt.gca().transAxes)
plt.text(0.1, 0.8, f'Optimized energy per atom: {fit.eos_params[0]/ouc.num_sites}eV',
         transform=plt.gca().transAxes)
plt.savefig('latconst_converge.png', dpi=300, bbox_inches='tight')
