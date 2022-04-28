from pymatgen.core.surface import SlabGenerator
from pymatgen.ext.matproj import MPRester


def generate_slabs(mpid,
                   miller_index,
                   slab_thickness,
                   vacuum_thickness,
                   conv=True,
                   symmetrize=False,
                   to_file=False,
                   filter_slabs=False):
    with MPRester() as m:
        conv_bulk = m.get_structure_by_material_id(mpid, conventional_unit_cell=conv)

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
        bulk_formula = conv_bulk.composition.reduced_formula
        conv_bulk.to('poscar', f'{bulk_formula}_conv_bulk.vasp')
        for index, slab in enumerate(slabs):
            formula = slab.composition.reduced_formula
            slab.to('poscar', f'{formula}_{millerstr}_conv_{conv}_term_{index}.vasp')

    return conv_bulk, slabs


au_mpid = 'mp-81'
lamno3_mpid = 'mp-19025'

au_slabs = generate_slabs(mpid=au_mpid,
                          miller_index=(1, 0, 0),
                          slab_thickness=10,
                          vacuum_thickness=20,
                          conv=False,
                          symmetrize=True,
                          to_file=True,
                          filter_slabs=False)
