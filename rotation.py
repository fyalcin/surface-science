import numpy as np
from scipy.spatial.transform import Rotation as R
from triboflow.utils.structure_manipulation import slab_from_file
from pymatgen.core.structure import Structure


def rotate_and_shift(vec, shift, angle):
    # angle = angle * np.pi / 180
    vec = np.asarray(vec)
    shift = np.asarray(shift)

    vec_shifted = vec - shift
    rotmat = R.from_euler('y', angle, degrees=True)
    vec_rotated = rotmat.apply(vec_shifted)
    # rotmat = np.array(((np.cos(angle), np.sin(angle)), (-np.sin(angle), np.cos(angle))))
    # vec_rotated = np.dot(rotmat, vec_shifted)
    vec_final = vec_rotated + shift

    return vec_final


shift = np.array((23.93000, 3.83959, 0.00000))

vec = np.array((25.20563, 1.91979, 0.46429))

positions = []
for angle in [0, -20, -40, -60, -80]:
    final_pos = rotate_and_shift(vec, shift, angle) - vec
    positions.append(final_pos)

struct = Structure.from_file('Si_100_new_2x1.vasp')

slab = Slab(lattice=struct.lattice, species=struct.species, coords=struct.frac_coords, miller_index=miller,
            oriented_unit_cell=struct, shift=0, scale_factor=np.array((1, 0, 0), (0, 1, 0), (0, 0, 1)))
indices = [1, 2, 21, 22]
