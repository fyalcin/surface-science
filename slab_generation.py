from pymatgen.core.surface import SlabGenerator
from pymatgen.ext.matproj import MPRester


def generate_slabs(mpid,
                   miller_index,
                   slab_thickness,
                   vacuum_thickness,
                   symmetrize=False,
                   to_file=False,
                   filter_slabs=False):
    with MPRester() as m:
        conv_bulk = m.get_structure_by_material_id(mpid, conventional_unit_cell=True)

    sg = SlabGenerator(initial_structure=conv_bulk,  # conventional bulk structure
                       miller_index=miller_index,  # miller index for slabs
                       min_slab_size=slab_thickness,  # minimum slab thickness
                       min_vacuum_size=vacuum_thickness,  # minimum vacuum region thickness
                       center_slab=True,  # whether to center the slabs in the c-direction
                       in_unit_planes=True,  # number of layers or Angstroms for thicknesses
                       primitive=True,  # whether to apply cell reduction on generated slabs
                       max_normal_search=max([abs(m) for m in miller_index]))

    slabs = sg.get_slabs(symmetrize=symmetrize)

    # if filter_slabs:
    #     filtered_slabs = []
    #     for slab in slabs:
    #         slab.add_oxidation_state_by_guess()
    #         if not slab.is_polar():
    #             filtered_slabs.append(slab)
    #     slabs = filtered_slabs
    #     if len(slabs) == 0:
    #         return None

    if to_file:
        millerstr = ''.join([str(m) for m in miller_index])
        for index, slab in enumerate(slabs):
            formula = slab.composition.reduced_formula
            slab.to('poscar', f'{formula}_{millerstr}_term_{index}.vasp')

    return slabs


au_mpid = 'mp-81'
lamno3_mpid = 'mp-19025'

au_slabs = generate_slabs(mpid=au_mpid,
                          miller_index=(1, 1, 1),
                          slab_thickness=10,
                          vacuum_thickness=20,
                          symmetrize=False,
                          to_file=True,
                          filter_slabs=False)
