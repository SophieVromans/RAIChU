from pikachu.general import read_smiles, draw_structure
from pikachu.reactions.functional_groups import combine_structures

from raichu.data.attributes import ATTRIBUTES


def epimerization(chiral_centre):
    new_chirality = None
    if chiral_centre.chiral == 'clockwise':
        new_chirality = 'counterclockwise'
    elif chiral_centre.chiral == 'counterclockwise':
        new_chirality = 'clockwise'

    chiral_centre.chiral = new_chirality


def epimerize(nrp):
    """

    """
    did_reaction = False
    # Define epimerization reaction target
    chiral_c_ep_atoms = []
    for atom in nrp.graph:
        if atom.annotations.chiral_c_ep:
            chiral_c_ep_atoms.append(atom)

    # If the substrate is a (fatty) acid, epimerization is not possible
    if len(chiral_c_ep_atoms) == 0:
       print('Warning: Cannot perform epimerization reaction on non-amino acid substrate!')

    assert len(chiral_c_ep_atoms) < 2

    # Carry out epimerization
    if len(chiral_c_ep_atoms) == 1:
        chiral_c = chiral_c_ep_atoms[0]
        epimerization(chiral_c)
        if chiral_c.chiral:
            did_reaction = True

    return nrp, did_reaction


def methylation(target_atom, structure):

    methyl_group = read_smiles('C')
    methyl_group.add_attributes(ATTRIBUTES, boolean=True)
    carbon = methyl_group.atoms[0]
    hydrogen_1 = methyl_group.atoms[1]
    bond_1 = carbon.get_bond(hydrogen_1)

    hydrogen_2 = target_atom.get_neighbour('H')
    if not hydrogen_2:
        raise Exception("Can't methylate this atom!")

    bond_2 = target_atom.get_bond(hydrogen_2)

    methylated_structure = combine_structures([structure, methyl_group])

    methylated_structure.break_bond(bond_1)
    methylated_structure.break_bond(bond_2)

    methylated_structure.make_bond(carbon, target_atom, methylated_structure.find_next_bond_nr())
    methylated_structure.make_bond(hydrogen_1, hydrogen_2, methylated_structure.find_next_bond_nr())

    structures = methylated_structure.split_disconnected_structures()

    for s in structures:
        if carbon.nr in s.atoms:
            return s


def n_methylate(nrp):
    """

    """
    # Define N-methylation reaction target
    n_meth_locations = []
    for atom in nrp.graph:
        if atom.annotations.n_atom_nmeth:
            n_meth_locations.append(atom)

    # If the substrate is not an amino acid, N-methylation is not possible
    if len(n_meth_locations) == 0:
        print('Warning: Cannot perform N-methylation on a non-amino acid substrate!')

    assert len(n_meth_locations) < 2

    product = None

    if len(n_meth_locations) == 1:

        n_meth = n_meth_locations[0]

        # Check if the N atom has a hydrogen group necessary for the reaction, and
        # not a cyclic amiono acid such as proline
        executable = True
        if not n_meth.has_neighbour('H'):
            executable = False
            product = nrp

        # Carry out N-methylation
        if executable:
            product = methylation(n_meth, nrp)

    else:
        product = nrp

    assert product

    return product, True


if __name__ == "__main__":
    s = read_smiles(r"NCC(=O)NCC(=O)O")
    target_atom = s.atoms[4]
    structure = methylation(target_atom, s)
    draw_structure(structure)





