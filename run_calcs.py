import numpy as np
from pymatgen.core.structure import Structure
from scipy.spatial.transform import Rotation as R
from atomate.vasp.workflows.presets.core import wf_bandstructure_no_opt
from fireworks.core.rocket_launcher import rapidfire
from fireworks import LaunchPad, Workflow, Firework

from triboflow.firetasks.start_swfs import FT_StartSurfaceEnergySWF
from triboflow.utils.database import Navigator
from datetime import datetime

lpad = LaunchPad.auto_load()



rapidfire(lpad)

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


indices_to_move = [[1, 7], [2, 8], [28, 30], [27, 29]]

for angle in np.arange(0, 50, 5):
    slab = Structure.from_file('Si100_2x1_unrecon.vasp')
    for pair in indices_to_move:
        site, ref = pair[0] - 1, pair[1] - 1
        shift = slab.cart_coords[ref]
        vec = slab.cart_coords[site]
        coeff = 1 if site % 2 else -1
        diff = rotate_and_shift(vec, shift, coeff * angle) - vec
        slab.translate_sites(indices=site, vector=diff, frac_coords=False, to_unit_cell=True)
    # slab.to('poscar', f"Si_recon_{angle}.vasp")
    wf = wf_bandstructure_no_opt(slab)
    lpad.add_wf(wf)



