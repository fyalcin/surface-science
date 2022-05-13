from pymatgen.core.surface import SlabGenerator, ReconstructionGenerator
from pymatgen.ext.matproj import MPRester
from triboflow.phys.shaper import Shaper


def generate_slabs(mpid,
                   miller_index,
                   slab_thickness,
                   vacuum_thickness,
                   prim=True,
                   symmetrize=False,
                   to_file=False,
                   filter_slabs=False,
                   recon=None):
    with MPRester() as m:
        conv_bulk = m.get_structure_by_material_id(mpid, conventional_unit_cell=True)

    sg = SlabGenerator(initial_structure=conv_bulk,  # conventional bulk structure
                       miller_index=miller_index,  # miller index for slabs
                       min_slab_size=slab_thickness,  # minimum slab thickness
                       min_vacuum_size=vacuum_thickness,  # minimum vacuum region thickness
                       center_slab=True,  # whether to center the slabs in the c-direction
                       in_unit_planes=False,  # number of layers or Angstroms for thicknesses
                       primitive=prim,  # whether to apply cell reduction on generated slabs
                       max_normal_search=max([abs(m) for m in miller_index]))

    slabs = sg.get_slabs(symmetrize=symmetrize)
    for slab in slabs:
        slab.add_oxidation_state_by_guess()

    slabs = [slab for slab in slabs if not slab.is_polar()]
    millerstr = ''.join([str(m) for m in miller_index])

    if recon:
        rg = ReconstructionGenerator(initial_structure=conv_bulk,
                                     min_slab_size=slab_thickness,
                                     min_vacuum_size=vacuum_thickness,
                                     reconstruction_name=recon)

        recon_slabs = rg.build_slabs()
        unrecon_slabs = rg.get_unreconstructed_slabs()
        if to_file:
            for index, slab in enumerate(recon_slabs):
                formula = slab.composition.reduced_formula
                slab.to('poscar', f'{formula}_{millerstr}_prim_{prim}_recon_term_{index}.vasp')
            for index, slab in enumerate(unrecon_slabs):
                formula = slab.composition.reduced_formula
                slab.to('poscar', f'{formula}_{millerstr}_prim_{prim}_unrecon_term_{index}.vasp')

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
        bulk_formula = conv_bulk.composition.reduced_formula
        conv_bulk.to('poscar', f'{bulk_formula}_conv_bulk.vasp')
        slab = slabs[0]
        slab.oriented_unit_cell.to('poscar', f"{bulk_formula}_{slab.miller_index}_ouc.vasp")

        for index, slab in enumerate(slabs):
            formula = slab.composition.reduced_formula
            slab.to('poscar', f'{formula}_{millerstr}_prim_{prim}_term_{index}.vasp')

    if recon:
        return {'bulk': conv_bulk, 'slabs': slabs, 'recon_slabs': recon_slabs, 'unrecon_slabs': unrecon_slabs}
    else:
        return {'bulk': conv_bulk, 'slabs': slabs}


mpid = 'mp-2657'
lamno3_mpid = 'mp-19025'

struct_dict = generate_slabs(mpid=mpid,
                             miller_index=(1, 1, 0),
                             slab_thickness=16,
                             vacuum_thickness=15,
                             prim=True,
                             symmetrize=True,
                             to_file=True,
                             filter_slabs=False,
                             recon=None)

# print(f"Unrecon slabs have {len(Shaper.get_layers(unrecon_slabs[0]))} layers")
# print(f"1x1 slabs have {len(Shaper.get_layers(slabs[0]))} layers")
#
# sd_array = []
# for i in range(len(recon_slabs[0].sites)):
#     sd_array.append([False, False, True])
# recon_slabs[0].add_site_property('selective_dynamics', sd_array)
# recon_slabs[0].to('poscar', 'Si100_recon.vasp')
