from pikachu.reactions.functional_groups import combine_structures
from pikachu.general import draw_structure
from raichu.data.molecular_moieties import THIOESTERBOND
from raichu.reactions.general import label_rest_groups


def pks_elongation(chain_intermediate, elongation_monomer):
    """
    Returns a Structure object of the PK chain after a single elongation step
    elongation step.

    structure: PIKAChU Structure object of the RC(=O)S PK intermediate before
    the elongation step
    elongation_monomer: ['SMILES_elongation_unit', index_c_to_c, index_c_to_s']
    """

    # If this is the first elongation reaction on the starter unit, define
    # central chain atoms in the starter unit
    assert any(atom.annotations.in_central_chain for atom in chain_intermediate.graph)

    h_to_remove_1 = None
    h_to_remove_2 = None

    for atom in elongation_monomer.c_to_pk_intermediate.neighbours:
        if atom.type == 'H':
            h_to_remove_1 = atom
            break

    for atom in elongation_monomer.c_to_s.neighbours:
        if atom.type == 'H':
            h_to_remove_2 = atom
            break

    assert h_to_remove_1 and h_to_remove_2

    # If the elongation unit contains a wildcard, add a number to the
    # '*' atom to display in the structure drawing
    label_rest_groups(elongation_monomer.structure, chain_intermediate)

    thioester_bonds = find_bonds(chain_intermediate, THIOESTERBOND)

    assert len(thioester_bonds) == 1
    for bond in thioester_bonds:
        if bond.atom_1.type == 'S':
            s_pkchain = bond.atom_1
            c_pkchain = bond.atom_2
        elif bond.atom_2.type == 'S':
            s_pkchain = bond.atom_2
            c_pkchain = bond.atom_1

    new_structure = combine_structures([elongation_monomer.structure, chain_intermediate])

    # Remove the Hs in the malonyl derivative in order to add it to the pk chain
    elongation_monomer.structure.remove_atom(new_structure.get_atom(h_to_remove_1))
    elongation_monomer.structure.remove_atom(new_structure.get_atom(h_to_remove_2))

    # Find thioester bond in growing chain, define Cs to attach new unit

    # Breaking the thioester bond in the PK chain
    for atom in chain_intermediate.graph:
        atom.hybridise()
    for bond in thioester_bonds[:]:
        chain_intermediate.break_bond(bond)
    chain_intermediate.get_connectivities()
    chain_intermediate.set_connectivities()
    chain_intermediate.set_atom_neighbours()
    chain_intermediate.make_bond_lookup()

    # Refresh malonylunit

    elongation_monomer.structure.refresh_structure()

    # Combining structures PK chain and elongation unit into one Structure object
    pk_chain_and_malonyl = (elongation_monomer.structure, chain_intermediate)

    combined = combine_structures(pk_chain_and_malonyl)

    combined.get_connectivities()
    combined.set_connectivities()
    combined.set_atom_neighbours()
    combined.make_bond_lookup()

    # Adding the bonds to form the new molecule after the single elongation step
    new_bond_nr = combined.find_next_bond_nr()
    combined.make_bond(elongation_monomer.c_to_pk_intermediate, c_pkchain, new_bond_nr)
    new_bond_nr = combined.find_next_bond_nr()
    combined.make_bond(elongation_monomer.c_to_s, s_pkchain, new_bond_nr)
    combined.get_connectivities()
    combined.set_connectivities()
    combined.set_atom_neighbours()
    combined.make_bond_lookup()
    for bond_nr, bond in combined.bonds.items():
        bond.set_bond_summary()

    for atom in combined.graph:
        atom.annotations.c2_acid = False

    return combined


def find_bonds(structure, bond_type):
    """
    Returns a list of Pikachu.bond objects of the bonds of the indicated type
    in the structure of interest

    structure: Pikachu.structure object of the molecule of interest
    bond_type: BondDefiner object indicating the bond type that is searched
    in the indicated structure
    """
    for atom in structure.graph:
        atom.get_connectivity()
    locations = structure.find_substructures(bond_type.structure)
    bonds = []
    for match in locations:
        atom_1 = match.atoms[bond_type.atom_1]
        atom_2 = match.atoms[bond_type.atom_2]
        bond = structure.bond_lookup[atom_1][atom_2]
        bonds.append(bond)

    return bonds
