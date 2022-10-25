import os

from pikachu.chem.structure import Structure
from pikachu.reactions.functional_groups import find_bonds, find_atoms
from pikachu.reactions.basic_reactions import hydrolysis, internal_condensation
from pikachu.general import structure_to_smiles

from raichu.data.molecular_moieties import SC_BOND, O_OH, O_BETAPROPRIOLACTONE_O,\
    O_BETAPROPRIOLACTONE_TERMINAL_O
from raichu.reactions.general import initialise_atom_attributes
from raichu.drawing.drawer import RaichuDrawer


def release_linear_reduction(chain_intermediate: Structure) -> Structure:
    sc_bonds = find_bonds(SC_BOND, chain_intermediate)
    if len(sc_bonds) != 1:
        raise ValueError("Cannot release product from carrier domain as no SC bond is present")

    sc_bond = sc_bonds[0]

    carbon = sc_bond.get_neighbour('C')
    sulphur = sc_bond.get_neighbour('S')

    chain_intermediate.break_bond(sc_bond)

    chain_intermediate.add_atom('H', [carbon])
    chain_intermediate.add_atom('H', [sulphur])

    initialise_atom_attributes(chain_intermediate)

    structures = chain_intermediate.split_disconnected_structures()
    for structure in structures:
        if carbon in structure.graph:
            structure.refresh_structure()
            return structure

    raise RuntimeError("Could not release chain through reduction.")


def release_linear_thioesterase(chain_intermediate):
    """Carries out the thioesterase reaction on the polyketide chain
    intermediate, returning the free acid linear polyketide product as a
    PIKAChU Structure object

    chain_intermediate: PIKAChU Structure object of the polyketide chain
    intermediate, either attached to a PKS domain or not
    """
    # Find S-H bond in chain intermediate and break the bond and define atoms
    sc_bonds = find_bonds(SC_BOND, chain_intermediate)
    assert len(sc_bonds) == 1
    carbon = sc_bonds[0].get_neighbour('C')

    linear_product = None

    structures = hydrolysis(chain_intermediate, sc_bonds[0])
    for structure in structures:
        if carbon in structure.graph:
            linear_product = structure
            break

    assert linear_product

    initialise_atom_attributes(linear_product)

    carbon = linear_product.get_atom(carbon)
    carbon.annotations.set_annotation('terminal_c', True)

    oxygens = carbon.get_neighbours('O')
    terminal_oxygen = None

    for oxygen in oxygens:
        if oxygen.has_neighbour('H'):
            terminal_oxygen = oxygen
            break

    assert terminal_oxygen

    terminal_oxygen.annotations.set_annotation('terminal_o', True)

    return linear_product


def cyclic_release(linear_product, o_oh_n_amino):
    """Performs the thioesterase reactions on the input chain_intermediate
     using the -OH group defined by the input O-atom participating in that
     internal -OH group, returns the circular product as PIKAChU Structure
     object.

     chain_intermediate: PIKAChU Structure object of a polyketide
     o_oh_n_amino: PIKAChU Atom object of the O-atom in the -OH group or N-atom
     in the amino group that the function should use to perform the
     thioesterase reaction.
    """

    # Remove H from internal -OH group and refresh structure
    cyclisation_site = linear_product.get_atom(o_oh_n_amino)
    h_atom = cyclisation_site.get_neighbour('H')
    assert h_atom
    h_bond = cyclisation_site.get_bond(h_atom)
    assert h_bond
    terminal_oxygen = None
    terminal_carbon = None

    for atom in linear_product.graph:
        if atom.annotations.terminal_o:
            terminal_oxygen = atom
        elif atom.annotations.terminal_c:
            terminal_carbon = atom

    assert terminal_oxygen and terminal_carbon

    oh_bond = terminal_oxygen.get_bond(terminal_carbon)

    cyclic_product, water = internal_condensation(linear_product, oh_bond, h_bond)

    return cyclic_product


def find_o_betapropriolactone(polyketide):
    """
    Finds and returns the oxygen atom (PIKAChU atom object) in the -OH group
    that shouldn't be used by the thioesterase_circular_product function, as
    this will create a beta-propriolactone compound, which does not occur in
    polyketide synthesis, if present. Otherwise the function returns None

    polyketide: PIKAChU structure object of a polyketide
    """
    o_propriolactone = find_atoms(O_BETAPROPRIOLACTONE_O, polyketide)
    o_propriolactone_terminal = find_atoms(O_BETAPROPRIOLACTONE_TERMINAL_O, polyketide)

    assert len(o_propriolactone) == len(o_propriolactone_terminal)

    o_beta_propriolactones = []

    for i, atom in enumerate(o_propriolactone):
        terminal_o = o_propriolactone_terminal[i]
        if terminal_o.annotations.terminal_o:
            o_beta_propriolactones.append(atom)

    if len(o_beta_propriolactones) == 0:
        return None
    elif len(o_beta_propriolactones) == 1:
        return o_beta_propriolactones[0]
    else:
        raise ValueError('Error: this molecule is not a polyketide, as the \
        carbon in the beta ketone/hydroxyl group is bound to an additional \
         oxygen atom')


def thioesterase_all_products(chain_intermediate, out_folder=None):
    """Performs all thioesterase reactions on the input chain_intermediate
     using all internal amino and -OH groups except for the -OH group that
     leads to the formation of a beta-propriolactone compound, which does not
     occur in polyketide synthesis. Returns a list of PIKAChU Structure objects
     of all possible thioesterase products.

     chain_intermediate: PIKAChU Structure object of a polyketide/NRP
    """
    # Perform first thioesterase reaction, generating linear polyketide/NRP

    chain_intermediate_copy = chain_intermediate.deepcopy()
    linear_product = release_linear_thioesterase(chain_intermediate_copy)
    if not out_folder:
        RaichuDrawer(linear_product)

    # Find OH groups in polyketide/NRP, perform cyclization for each -OH group
    o_oh_atoms = find_atoms(O_OH, linear_product)

    o_oh_atoms_filtered = []
    for atom in o_oh_atoms:
        if not atom.annotations.terminal_o and atom not in o_oh_atoms_filtered and atom.has_neighbour('H'):
            o_oh_atoms_filtered.append(atom)

    # Find amino groups in polyketide/NRP, perform cyclization for each group

    amino_n_atoms_filtered = []

    for atom in linear_product.graph:
        if atom.type == 'N':
            if len(atom.get_neighbours('H')) == 2 and atom.has_neighbour('C') and not atom.aromatic and atom \
                    not in amino_n_atoms_filtered:
                amino_n_atoms_filtered.append(atom)

    # Define -OH group that should not be used to carry out the thioesterase
    # reaction (distance -S and internal -OH group)
    o_not_to_use = find_o_betapropriolactone(linear_product)

    list_product_drawings = []
    circular_smiles = []

    # Perform all possible thioesterase reactions leading to the formation of
    # circular products using the internal amino groups, save Structure objects
    # to list

    for n_amino in amino_n_atoms_filtered:
        linear_product_copy = linear_product.deepcopy()
        n_atom = linear_product_copy.get_atom(n_amino)
        product = cyclic_release(linear_product_copy, n_atom)
        circular_smiles.append(structure_to_smiles(product))
        list_product_drawings.append(product)

    # Perform all possible thioesterase reactions leading to the formation of
    # circular products using the internal -OH groups, save Structure objects
    # to list
    for o_oh in o_oh_atoms_filtered:
        linear_product_copy = linear_product.deepcopy()
        oh_atom = linear_product_copy.get_atom(o_oh)

        if oh_atom != o_not_to_use:
            product = cyclic_release(linear_product_copy, oh_atom)
            circular_smiles.append(structure_to_smiles(product))
            list_product_drawings.append(product)

    if out_folder:
        smiles_path = os.path.join(out_folder, "product_smiles.txt")
        smiles_file = open(smiles_path, 'w')

        file_path = os.path.join(out_folder, f"product_0.png")

        if os.path.exists(file_path):
            os.remove(file_path)

        linear_smiles = structure_to_smiles(linear_product)
        smiles_file.write(f"product_0\t{linear_smiles}\n")

        # png_from_smiles(linear_smiles, file_path)

        drawing = RaichuDrawer(linear_product, save_png=file_path)
        drawing.draw_structure()

        for i, product in enumerate(list_product_drawings):
            smiles = circular_smiles[i]
            smiles_file.write(f"product_{i}\t{smiles}\n")

            file_path = os.path.join(out_folder, f"product_{i + 1}.png")
            if os.path.exists(file_path):
                os.remove(file_path)
            drawing = RaichuDrawer(product, save_png=file_path)
            drawing.draw_structure()

            # png_from_smiles(smiles, file_path)

        smiles_file.close()

    else:

        # Draw all products
        for product in list_product_drawings:
            RaichuDrawer(product)

    return list_product_drawings
