import os
from os import path

import matplotlib.pyplot as plt
import matplotlib.image as mpimg
from matplotlib.patches import FancyArrow

from raichu.reactions.pks_elongation_reactions import PKS_SUBUNIT_TO_MONOMER, pks_elongation
from raichu.reactions.pks_tailoring_reactions import *
from raichu.drawing.RaichuDrawer import RaichuDrawer
from raichu.data.attributes import ATTRIBUTES, ALL_PKS_ELONGATION_UNITS
from raichu.reactions.nrps_condensation import condensation_nrps, label_nrp_central_chain
from raichu.reactions.nrps_tailoring_reactions import nrps_epimerization, nrps_methylation
from raichu.central_chain_detection.label_central_chain import label_pk_central_chain

import interactive.flatfiles

FLATFILES = os.path.dirname(interactive.flatfiles.__file__)
PARAS_SMILES = os.path.join(FLATFILES, "PARAS_smiles.txt")


KR_NO_KETOREDUCTASE_ACTIVITY = ['KR_C1', 'KR_C2', 'KR_inactive']

ELONGATION_UNIT_TO_TEXT = {'malonylcoa': 'Malonyl-CoA',
                           'methylmalonylcoa': 'Methylmalonyl-CoA',
                           'methoxymalonylacp': 'Methoxymalonyl-ACP',
                           'ethylmalonylcoa': 'Ethylmalonyl-CoA',
                           'pk': 'Unknown elongation unit'}

TAILOR_DOMAIN_SHORT_TO_LONG = {'E': 'Epimerization', 'nMT' : 'N-methylation'}

TAILOR_DOMAIN_TO_COLOUR = {'KR': 'red', 'KR_inactive': 'red', 'KR*': 'red',
                           'KR_A1': 'red', 'KR_A2': 'red', 'KR_B1': 'red',
                           'KR_B2': 'red', 'KR_C1': 'red', 'KR_C2': 'red',
                           'DH': 'blue', 'ER': 'green'}


def get_substrate_structure(substrate, aa_dict):
    if substrate.upper() in aa_dict:
        structure = read_smiles(aa_dict[substrate.upper()])
    else:
        structure = read_smiles(substrate)

    return structure


def cluster_to_structure(modules, visualization_mechanism=False,
                         draw_structures_per_module=False,
                         attach_to_acp=False,
                         draw_mechanism_per_module=False):
    """If visualization_mechanism == False: returns the final polyketide
    structure resulting from the input PKS cluster
    If visualization_mechanism == True: returns a visualization of the elong-
    ation and/or tailoring reactions performed by each PKS module
    If attach_to_acp == True: the polyketide will be attached to an ACP domain
    object of type Domain. Can be combined with previous optional flags.
    If draw_structures_per_module == True: saves the polyketide intermediates
    as 'module_name'.png, doesn't display any visualisations

    modules: [Starter module, elongation module, ..., terminator_module] with
    Starter module: ['module name','starter_module_pks','SMILES starter unit']
    or ['module name','starter_module_nrps','SMILES starter amino acid'],
    elongation modules: ['module name', 'elongation_module_pks',
    'elongation_monomer', ['KR', 'DH', 'ER']] or ['module name',
    'elongation_module_nrps', 'amino acid'],
    terminator module: ['module name', 'terminator_module_pks',
    'elongation_monomer', ['KR', 'DH', 'ER']] or ['module name',
    'terminator_module_nrps', 'amino acid'],
    """
    last_module_nrps = False
    # Check if last module is an NRPS module:
    if modules[-1][1] == 'elongation_module_nrps' or \
            modules[-1][1] == 'starter_module_nrps' or \
            modules[-1][1] == 'terminator_module_nrps':
        last_module_nrps = True

    # Construct dict to find the SMILES of all amino acids recognized by PARAS
    dict_aa_smiles = make_dict_aa_smiles()

    modules = modules[::]
    list_drawings_per_module = []

    # Retrieve and remove starter module
    starter_module = modules[0]
    if starter_module[1] == 'starter_module_pks':
        modules = modules[1:]

    assert len(starter_module) >= 3

    starter_module_type = starter_module[1]
    starter_module_smiles = starter_module[2]

    assert starter_module_type == 'starter_module_pks' or starter_module_type == 'starter_module_nrps'
    if starter_module_type == 'starter_module_pks':
        assert len(starter_module) == 4
    elif starter_module_type == 'starter_module_nrps':
        assert len(starter_module) == 4 or len(starter_module) == 5

    if len(starter_module) > 3:
        tailoring_domains = starter_module[3]

    # If starter module = PKS module: construct starter unit from
    # supplied SMILES string

    if starter_module_type == 'starter_module_pks':
        starter_unit = read_smiles(starter_module_smiles)
        starter_unit.add_attributes(ATTRIBUTES, boolean=True)
        starter_unit = label_pk_central_chain (starter_unit)

        if attach_to_acp:
            chain_intermediate = attach_to_domain_pk(starter_unit)
        else:
            chain_intermediate = starter_unit
        if draw_structures_per_module:
            drawing = RaichuDrawer(chain_intermediate, dont_show=True)
            list_drawings_per_module.append([drawing])

    # Iterate over remaining modules

    for module in modules:
        list_filenames = []

        # If module is a PKS module:

        if module[1] == 'elongation_module_pks' or module[1] == 'terminator_module_pks':
            # If last reaction was an NRPS reaction
            locations = chain_intermediate.find_substructures(read_smiles('C(=O)(O)CN'))
            b_locations = chain_intermediate.find_substructures(read_smiles('C(=O)(O)CCN'))
            asp_locations = chain_intermediate.find_substructures(read_smiles('C(=O)(O)CC(C=O)N'))
            mal_amino_c_locations = chain_intermediate.find_substructures(read_smiles('NC(=O)CC(O)=O'))

            location_nr = len(locations) + len(b_locations) - len(asp_locations) - len(mal_amino_c_locations)
            if location_nr > 0:
                chain_intermediate = attach_to_domain_nrp(chain_intermediate)

            elif any(atom.annotations.c2_acid for atom in chain_intermediate.graph):
                chain_intermediate = attach_to_domain_nrp(chain_intermediate)

            module_name, module_type, elongation_unit, list_domains, gene_name = module
            assert elongation_unit in ALL_PKS_ELONGATION_UNITS

            # Reset atom colours in the first structure depicted to black
            for atom in chain_intermediate.graph:
                atom.draw.colour = 'black'

            elongation_monomer = PKS_SUBUNIT_TO_MONOMER[elongation_unit]

            if visualization_mechanism:
                # If a previous process was stopped prematurely, make sure the
                # image files are removed before continuing
                if path.exists('1.png'):
                    os.remove('1.png')
                RaichuDrawer(chain_intermediate, save_png='1.png', dpi=500)
                list_filenames.append('1.png')

            new_chain_intermediate = pks_elongation(chain_intermediate, elongation_monomer)
            chain_intermediate = new_chain_intermediate

            if visualization_mechanism:
                if path.exists('2.png'):
                    os.remove('2.png')
                RaichuDrawer(chain_intermediate, save_png='2.png', dpi=500)
                list_filenames.append('2.png')

            # If the module does not contain any tailoring domains:
            if len(list_domains) == 0:
                if visualization_mechanism:
                    display_reactions(['1.png', '2.png'],
                                      list_domains, elongation_unit,
                                      module_name,
                                      draw_mechanism_per_module)

            # If module contains tailoring domains, perform tailoring reactions
            for domain in list_domains:
                inactive_kr_domain = False

                if domain.startswith('KR'):
                    if domain == 'KR':
                        new_chain_intermediate = ketoreduction(chain_intermediate)
                        chain_intermediate = new_chain_intermediate
                    elif domain == 'KR_inactive':
                        inactive_kr_domain = True
                    else:
                        kr_type = domain[-2:]
                        new_chain_intermediate = \
                            ketoreduction(chain_intermediate, kr_type = kr_type)
                        chain_intermediate = new_chain_intermediate
                    if visualization_mechanism and not inactive_kr_domain:
                        if path.exists('3.png'):
                            os.remove('3.png')
                        RaichuDrawer(chain_intermediate, save_png='3.png', dpi=500)
                        list_filenames.append('3.png')
                        # Check if the module also contains a DH domain
                        if 'DH' not in list_domains:
                            display_reactions(['1.png', '2.png', '3.png'],
                                              list_domains, elongation_unit,
                                              module_name,
                                              draw_mechanism_per_module)
                    elif inactive_kr_domain and visualization_mechanism:
                        list_domains = []

                        display_reactions(['1.png', '2.png'], list_domains,
                                          elongation_unit, module_name,
                                          draw_mechanism_per_module)

                if domain == 'DH':
                    assert any(domain.startswith('KR') for domain in list_domains)
                    executable = True
                    for domain_type in list_domains:
                        if domain_type in KR_NO_KETOREDUCTASE_ACTIVITY:
                            executable = False
                    if executable:
                        new_chain_intermediate = dehydration(chain_intermediate)
                        chain_intermediate = new_chain_intermediate
                        if visualization_mechanism:
                            if path.exists('4.png'):
                                os.remove('4.png')
                            RaichuDrawer(chain_intermediate, save_png='4.png', dpi=500)
                            list_filenames.append('4.png')
                            if 'ER' not in list_domains:
                                display_reactions(['1.png', '2.png', '3.png','4.png'],
                                                  list_domains, elongation_unit,
                                                  module_name, draw_mechanism_per_module)

                if domain == 'ER':
                    assert any(domain.startswith('KR') for domain in list_domains)
                    assert 'DH' in list_domains
                    executable = True
                    for domain_type in list_domains:
                        if domain_type in KR_NO_KETOREDUCTASE_ACTIVITY:
                            executable = False
                    if executable:
                        new_chain_intermediate = enoylreduction(chain_intermediate)
                        chain_intermediate = new_chain_intermediate
                        if visualization_mechanism:
                            if path.exists('5.png'):
                                os.remove('5.png')
                            RaichuDrawer(chain_intermediate, save_png='5.png', dpi=500)
                            list_filenames.append('5.png')
                            display_reactions(['1.png', '2.png', '3.png','4.png','5.png'],
                                              list_domains, elongation_unit,
                                              module_name,
                                              draw_mechanism_per_module)

            if draw_structures_per_module and attach_to_acp:
                drawing = RaichuDrawer(chain_intermediate, dont_show=True)
                list_drawings_per_module.append([drawing])

        # If the module is an NRPS module:
        elif module[1] == 'elongation_module_nrps' or \
                module[1] == 'terminator_module_nrps' or \
                module[1] == 'starter_module_nrps':
            module_name, module_type, aa_specifity, list_tailoring_domains, gene_name = module
            if module[1] == 'starter_module_nrps':
                starter_unit = get_substrate_structure(starter_module_smiles, dict_aa_smiles)
                starter_unit.add_attributes(ATTRIBUTES, boolean=True)
                label_nrp_central_chain(starter_unit, module_type='starter')
                chain_intermediate = starter_unit

            if visualization_mechanism and attach_to_acp:
                for atom in chain_intermediate.graph:
                    atom.draw.colour = 'black'
                # If a previous process was stopped prematurely, make sure the
                # image files are removed before continuing
                if path.exists('1.png'):
                    os.remove('1.png')
                copy_chain_intermediate = chain_intermediate.deepcopy()
                copy_chain_intermediate.find_cycles()
                if not any(atom.annotations.domain_type for
                           atom in copy_chain_intermediate.graph):
                    copy_attached = attach_to_domain_nrp(copy_chain_intermediate)
                else:
                    copy_attached = copy_chain_intermediate
                RaichuDrawer(copy_attached, save_png='1.png', dpi=500)
                list_filenames.append('1.png')

            # Reset atom colours in the first structure depicted to black
            for atom in chain_intermediate.graph:
                atom.draw.colour = 'black'

            # Build amino acid structure from dict containig PARAS SMILES
            aa_specifity = aa_specifity.upper()

            aa_structure = get_substrate_structure(aa_specifity, dict_aa_smiles)

            # Check that the structure is an amino acid
            if module_type != 'starter_module_nrps':
                if not (len(aa_structure.find_substructures(
                        read_smiles('CN'))) > 0
                        and len(aa_structure.find_substructures(
                            read_smiles('C(O)=O'))) > 0):
                    raise ValueError(
                        f'The structure: {aa_specifity}, is not an amino acid')

            # Perform condensation reaction
            if module_type != 'starter_module_nrps':
                chain_intermediate = condensation_nrps(aa_structure, chain_intermediate)

            if len(list_tailoring_domains) == 0:
            # Save drawings if necessary
                if draw_structures_per_module and attach_to_acp:
                    copy_chain_intermediate = chain_intermediate.deepcopy()
                    # copy_chain_intermediate.find_cycles()
                    copy_attached = attach_to_domain_nrp(copy_chain_intermediate)
                    copy_attached.refresh_structure()
                    copy_attached.set_connectivities()
                    copy_attached.find_cycles()
                    drawing = RaichuDrawer(copy_attached, dont_show=True)
                    list_drawings_per_module.append([drawing])

                elif visualization_mechanism and attach_to_acp:
                        copy_chain_intermediate = chain_intermediate.deepcopy()
                        copy_chain_intermediate.find_cycles()
                        copy_attached = attach_to_domain_nrp(copy_chain_intermediate)
                        if path.exists('2.png'):
                            os.remove('2.png')
                        RaichuDrawer(copy_attached, save_png='2.png', dpi=500)
                        list_filenames.append('2.png')
                        elongation_unit = aa_specifity.lower()
                        display_reactions(['1.png', '2.png'], list_tailoring_domains, elongation_unit, module_name, draw_mechanism_per_module)
                        # Necessary for thioesterase reactions when last module is an
                        # NRPS module
                        if module == modules[-1]:
                            chain_intermediate = copy_attached
            else:
                for domain in list_tailoring_domains:
                    if domain == 'E':
                        if draw_structures_per_module and attach_to_acp:
                            chain_intermediate = nrps_epimerization(chain_intermediate)
                            if domain == list_tailoring_domains[-1]:
                                copy_chain_intermediate = chain_intermediate.deepcopy()
                                # copy_chain_intermediate.find_cycles()
                                copy_attached = attach_to_domain_nrp(
                                    copy_chain_intermediate, 'PCP')
                                copy_attached.refresh_structure()
                                copy_attached.set_connectivities()
                                copy_attached.find_cycles()
                                drawing = RaichuDrawer(copy_attached,
                                                       dont_show=True)
                                list_drawings_per_module.append([drawing])
                        elif visualization_mechanism and attach_to_acp:
                            copy_chain_intermediate = chain_intermediate.deepcopy()
                            copy_chain_intermediate.find_cycles()
                            copy_attached = attach_to_domain_nrp(
                                copy_chain_intermediate)
                            after_reaction_filename = None

                            if domain == list_tailoring_domains[0]:
                                before_reaction_filename = '2.png'
                                after_reaction_filename = '3.png'
                                if path.exists(before_reaction_filename):
                                    os.remove(before_reaction_filename)
                                RaichuDrawer(copy_attached,
                                             save_png=before_reaction_filename,
                                             dpi=500)
                                list_filenames.append(before_reaction_filename)
                            elif len(list_tailoring_domains) == 2 and domain == list_tailoring_domains[1]:
                                after_reaction_filename = '4.png'

                            chain_intermediate = nrps_epimerization(
                                chain_intermediate)
                            copy_chain_intermediate = chain_intermediate.deepcopy()
                            copy_chain_intermediate.find_cycles()
                            copy_attached = attach_to_domain_nrp(
                                copy_chain_intermediate)

                            if path.exists(after_reaction_filename):
                                os.remove(after_reaction_filename)
                            RaichuDrawer(copy_attached, save_png=after_reaction_filename, dpi=500)
                            list_filenames.append(after_reaction_filename)
                            if domain == list_tailoring_domains[-1]:
                                elongation_unit = aa_specifity.lower()
                                display_reactions(list_filenames,
                                                  list_tailoring_domains,
                                                  elongation_unit, module_name,
                                                  draw_mechanism_per_module)
                                # Necessary for thioesterase reactions when last module is an
                                # NRPS module
                                if module == modules[-1]:
                                    chain_intermediate = copy_attached
                        elif attach_to_acp:
                            chain_intermediate = nrps_epimerization(chain_intermediate)

                    elif domain == 'nMT':
                        if draw_structures_per_module and attach_to_acp:
                            chain_intermediate = nrps_methylation(
                                chain_intermediate)
                            if domain == list_tailoring_domains[-1]:
                                copy_chain_intermediate = chain_intermediate.deepcopy()
                                # copy_chain_intermediate.find_cycles()
                                copy_attached = attach_to_domain_nrp(
                                    copy_chain_intermediate)
                                copy_attached.refresh_structure()
                                copy_attached.set_connectivities()
                                copy_attached.find_cycles()
                                drawing = RaichuDrawer(copy_attached,
                                                       dont_show=True)
                                list_drawings_per_module.append([drawing])
                        elif visualization_mechanism and attach_to_acp:
                            copy_chain_intermediate = chain_intermediate.deepcopy()
                            copy_chain_intermediate.find_cycles()
                            copy_attached = attach_to_domain_nrp(
                                copy_chain_intermediate)
                            if domain == list_tailoring_domains[0]:
                                before_reaction_filename = '2.png'
                                after_reaction_filename = '3.png'
                                if path.exists(before_reaction_filename):
                                    os.remove(before_reaction_filename)
                                RaichuDrawer(copy_attached,
                                             save_png=before_reaction_filename,
                                             dpi=500)
                                list_filenames.append(before_reaction_filename)
                            elif len(list_tailoring_domains) == 2 and domain == \
                                    list_tailoring_domains[1]:
                                after_reaction_filename = '4.png'

                            chain_intermediate = nrps_methylation(
                                chain_intermediate)
                            copy_chain_intermediate = chain_intermediate.deepcopy()
                            copy_chain_intermediate.find_cycles()
                            copy_attached = attach_to_domain_nrp(
                                copy_chain_intermediate)

                            if path.exists(after_reaction_filename):
                                os.remove(after_reaction_filename)
                            RaichuDrawer(copy_attached,
                                         save_png=after_reaction_filename,
                                         dpi=500)
                            list_filenames.append(after_reaction_filename)
                            if domain == list_tailoring_domains[-1]:
                                elongation_unit = aa_specifity.lower()
                                display_reactions(list_filenames,
                                                  list_tailoring_domains,
                                                  elongation_unit, module_name,
                                                  draw_mechanism_per_module)
                                # Necessary for thioesterase reactions when last module is an
                                # NRPS module
                                if module == modules[-1]:
                                    chain_intermediate = copy_attached

                        elif attach_to_acp:
                            chain_intermediate = nrps_methylation(chain_intermediate)

    # Reset the atom color in the final structure to black
    for atom in chain_intermediate.graph:
        atom.draw.colour = 'black'

    if not visualization_mechanism and not draw_structures_per_module:
        if last_module_nrps and attach_to_acp:
            # chain_intermediate.find_cycles()
            chain_intermediate = attach_to_domain_nrp(chain_intermediate)
            chain_intermediate.find_cycles()
            RaichuDrawer(chain_intermediate)
        else:
            # chain_intermediate.find_cycles()
            RaichuDrawer(chain_intermediate)

    # Remove all intermediate structure files
    for i in range(1, 6):
        if path.exists(f"{i}.png"):
            os.remove(f"{i}.png")

    if draw_structures_per_module:
        return(list_drawings_per_module)

    return chain_intermediate


def display_reactions(structures, tailoring_domains, elongation_unit,
                      module_name, draw_mechanism_per_module):
    """

    """
    # Create plot
    fig = plt.figure(constrained_layout=True, figsize=[15, 8],
                     frameon=False,
                     num=f"Quick reaction mechanism {module_name}")
    nr_subplots = 5 + (len(tailoring_domains) * 3)
    ax = fig.subplots(1, nr_subplots)
    ax = fig.add_gridspec(1, nr_subplots)

    # Read in all images
    images = []
    for im_filename in structures:
        image = mpimg.imread(im_filename)
        images.append(image) # Prevent disposal by garbage collector

    # position_1 and position_2 describe the subgraphs used to place the arrows
    # and structures in each reaction mechanism
    position_1 = 0
    position_2 = 2

    # Iterate over all image filenames where structure drawings are saved
    nr_of_structures = len(structures)
    for i in range(nr_of_structures):
        # Add chemical structures
        image = images[i]
        ax1 = fig.add_subplot(ax[0, position_1:position_2])
        ax1.imshow(image)

        # Add arrows
        if i != nr_of_structures - 1:
            ax1 = fig.add_subplot(ax[0, position_2])
            arrow = FancyArrow(0, 0.5, 0.89, 0, overhang=0.3, head_width=0.025,
                               head_length=0.1, color='black')
            images.append(arrow)
            ax1.add_patch(arrow)

            # Add text above arrows, default is black text, font size 12
            arrow_text_size = 12
            arrow_text_colour = 'black'

            # First arrow describes the PKS elongation unit or amino acid
            if i == 0:
                if elongation_unit in ALL_PKS_ELONGATION_UNITS:
                    arrow_text = ELONGATION_UNIT_TO_TEXT[elongation_unit]
                else:
                    arrow_text = elongation_unit

            # Other arrows describe tailoring reactions
            else:
                if elongation_unit in ALL_PKS_ELONGATION_UNITS:
                    arrow_text = tailoring_domains[i-1]
                    arrow_text_size = 17
                    arrow_text_colour = TAILOR_DOMAIN_TO_COLOUR[tailoring_domains[i-1]]
                else:
                    arrow_text = TAILOR_DOMAIN_SHORT_TO_LONG[tailoring_domains[i-1]]

            ax1.text(0.5, 0.55, arrow_text, ha='center',
             fontdict={'size': arrow_text_size, 'color': arrow_text_colour})

            # Move to the next subplots in the graph
            position_1 += 3
            position_2 += 3

    # Remove all axes from the figure
    for ax in fig.get_axes():
        ax.axis('off')

    # Save reaction mechanism drawing
    if draw_mechanism_per_module:
        plt.savefig(f'{module_name}_quick_mechanism.png')
        plt.clf()
        plt.close()
    else:
        plt.show()


def make_dict_aa_smiles():
    """
    Returns a dict as {name_aa : SMILES_aa} from the amino acids in the
    PARAS_smiles.txt file
    """
    # Parse list SMILES amino acids attached to PCP to dict name -> Structure
    lines_aa = open(PARAS_SMILES, 'r', encoding='utf8').readlines()
    dict_aa_smiles = {}

    for line in lines_aa:
        line = line.strip()
        name, smiles = line.split()
        name = name.upper()
        dict_aa_smiles[name] = smiles

    return dict_aa_smiles







