from raichu.class_domain import *
from pikachu.reactions.functional_groups import find_atoms, GroupDefiner, combine_structures
from pikachu.reactions.basic_reactions import condensation
from raichu.reactions.general import initialise_atom_attributes
from pikachu.general import read_smiles

POLYKETIDE_S = GroupDefiner('Sulphur atom polyketide', 'SC(C)=O', 0)
POLYKETIDE_S_INSERTED_O = GroupDefiner(
    'Sulphur atom polyketide inserted O', 'SC(O)=O', 0)
POLYKETIDE_S_REDUCED_STARTER = GroupDefiner('Sulphur atom polyketide reduced O', 'SC(C)O', 0)
NRP_C = GroupDefiner('C atom to attach to PCP domain', 'NCC(O)=O', 2)
RIPP_N = GroupDefiner('N atom to attach to leader', 'NCC=O', 0)
AMINO_ACID_BACKBONE = read_smiles('NCC(=O)O')
AMINO_ACID_BACKBONE_N_TERMINUS = read_smiles('NCC(=O)N')
B_AMINO_ACID_BACKBONE = read_smiles('NCCC(=O)O')
ASP_ACID_C = GroupDefiner('Asp', 'NC(C=O)CC(O)=O', 5)
B_NRP_C = GroupDefiner('Beta C atom to attach to PCP domain', 'NCCC(O)=O', 3)
MAL_AMINO = GroupDefiner('malonyl_starter_amino', 'NC(=O)CC(O)=O', 4)


def condensation_with_reverse_numbering(structure_1, structure_2, oh_bond, h_bond):
    """
    Returns condensed product of structure 1 and structure 2, and a water molecule, where the numbering of structure 1 is kept

    Input
    ----------
    structure_1: Structure object, structure containing -OH leaving group
    structure_2: Structure object, structure containing -H leaving group
    oh_bond: Bond object, bond attached to -OH leaving group
    h_bond: Bond object, bond attached to -H leaving group

    Returns
    -------
    list of [product, water], both Structure objects, with product the condensed product of structure 1 and structure 2
        and water a water molecule.

    """

    # Ensure that the defined OH-bond actually has an -OH leaving group
    o_atom = None
    o_neighbour = None

    for atom in oh_bond.neighbours:
        if not o_atom and atom.type == 'O' and atom.has_neighbour('H'):
            o_atom = atom
        else:
            o_neighbour = atom

    if not o_atom:
        raise Exception(
            f"The selected bond {oh_bond} does not attach to an -OH leaving group.")

    # Ensure that the defined H-bond actually has an -H leaving group

    h_atom = None
    h_neighbour = None

    for atom in h_bond.neighbours:
        if not h_atom and atom.type == 'H':
            h_atom = atom
        else:
            h_neighbour = atom

    if not h_atom:
        raise Exception(
            f"The selected bond {h_bond} does not attach to an -H leaving group.")

    # Break the bonds between the leaving groups and the rest of the product

    structure_1.break_bond(oh_bond)
    structure_2.break_bond(h_bond)

    # Put the two structures into the same object
    structure = combine_structures([structure_2, structure_1])
    # Retrieve the atoms that have to form bonds: h_atom with o_atom to form water, and h_neighbour with o_neighbour
    # to form the product

    h_atom = structure.get_atom(h_atom)
    o_atom = structure.get_atom(o_atom)
    h_neighbour = structure.get_atom(h_neighbour)
    o_neighbour = structure.get_atom(o_neighbour)

    # Create the bonds

    structure.make_bond(h_atom, o_atom, structure.find_next_bond_nr())
    structure.make_bond(h_neighbour, o_neighbour,
                        structure.find_next_bond_nr())

    # Put the water and the product into different Structure instances

    structures = structure.split_disconnected_structures()

    water = None
    product = None

    # Find out which of the structures is your product and which is your water

    for structure in structures:
        if h_atom in structure.graph:
            water = structure
        elif h_neighbour in structure.graph:
            structure.refresh_structure(find_cycles=True)
            product = structure

    return [product, water]


def attach_to_domain_pk(polyketide):
    """
    Attaches the sulphur atom in the input polyketide to a PKS domain and
    returns the attached structure as a PIKAChU Structure object

    polyketide: PIKAChU Structure object, to-be attached structure
    """
    # Add attributes to input molecule
    initialise_atom_attributes(polyketide)

    # Create domain
    domain = make_scaffold_domain('ACP')
    sh_bond = domain.bond_lookup[domain.atoms[1]][domain.atoms[2]]
    hydrogen = domain.atoms[2]
    sulphur_1 = domain.atoms[1]

    locations_sulphur = find_atoms(POLYKETIDE_S, polyketide)
    if len(locations_sulphur) == 0:
        locations_sulphur = find_atoms(
            POLYKETIDE_S_INSERTED_O, polyketide)
        if len(locations_sulphur) == 0:
            locations_sulphur = find_atoms(POLYKETIDE_S_REDUCED_STARTER, polyketide)
    assert len(locations_sulphur) == 1
    sulphur_2 = locations_sulphur[0]
    carbon = sulphur_2.get_neighbour('C')
    sc_bond = sulphur_2.get_bond(carbon)

    structure = combine_structures([domain, polyketide])
    structure.break_bond(sh_bond)
    structure.break_bond(sc_bond)

    structure.make_bond(hydrogen, sulphur_2, structure.find_next_bond_nr())
    structure.make_bond(carbon, sulphur_1, structure.find_next_bond_nr())

    split = structure.split_disconnected_structures()

    tethered_polyketide = None

    for structure in split:
        if carbon in structure.graph:
            tethered_polyketide = structure
            break

    assert tethered_polyketide
    initialise_atom_attributes(tethered_polyketide)

    return tethered_polyketide


def attach_to_domain_nrp(nrp):
    """
    Attaches the input NRP to a PCP domain and returns the product as a
    PIKAChU Structure object

    nrp: PIKAChU Structure object, to-be attached NRP
    """

    initialise_atom_attributes(nrp)

    # Create domain
    domain = make_scaffold_domain('PCP')

    # Remove OH group from carboxylic acid in NRP, to allow attachment
    # to domain

    if nrp.find_substructures(AMINO_ACID_BACKBONE):
        locations_c_to_domain = find_atoms(NRP_C, nrp)
        assert len(locations_c_to_domain) == 1
        c_atom_to_domain = locations_c_to_domain[0]
    elif nrp.find_substructures(B_AMINO_ACID_BACKBONE):
        locations_c_to_domain = find_atoms(B_NRP_C, nrp)
        asp_acid_cs = find_atoms(ASP_ACID_C, nrp)
        mal_amino_cs = find_atoms(MAL_AMINO, nrp)

        for asp_acid_c in asp_acid_cs:
            if asp_acid_c in locations_c_to_domain:
                locations_c_to_domain.remove(asp_acid_c)
        for mal_amino_c in mal_amino_cs:
            if mal_amino_c in locations_c_to_domain:
                locations_c_to_domain.remove(mal_amino_c)
        assert len(locations_c_to_domain) == 1
        c_atom_to_domain = locations_c_to_domain[0]
    else:
        c_atom_to_domain = None
        for atom in nrp.graph:
            if atom.annotations.c2_acid:
                c_atom_to_domain = atom

    assert c_atom_to_domain

    oxygens = c_atom_to_domain.get_neighbours('O')

    hydroxyl_oxygen = None
    hydroxyl_bond = None

    hydrogen_bond = domain.bond_lookup[domain.atoms[1]][domain.atoms[2]]

    assert hydrogen_bond.has_neighbour('S')
    assert hydrogen_bond.has_neighbour('H')

    for oxygen in oxygens:
        bond = c_atom_to_domain.get_bond(oxygen)
        if bond.type == 'single' and oxygen.has_neighbour('H'):
            hydroxyl_oxygen = oxygen
            hydroxyl_bond = bond
            break

    assert hydroxyl_oxygen and hydroxyl_bond

    structure = condensation(nrp, domain, hydroxyl_bond, hydrogen_bond)[0]

    initialise_atom_attributes(structure)

    return structure


def attach_to_follower_ripp(ripp):
    initialise_atom_attributes(ripp)

    # Create domain
    domain = make_scaffold_peptide('Follower')
    # Remove OH group from carboxylic acid in NRP, to allow attachment
    # to domain

    if ripp.find_substructures(AMINO_ACID_BACKBONE):
        locations_c_to_domain = find_atoms(NRP_C, ripp)
        assert len(locations_c_to_domain) == 1
        c_atom_to_domain = locations_c_to_domain[0]
    elif ripp.find_substructures(B_AMINO_ACID_BACKBONE):
        locations_c_to_domain = find_atoms(B_NRP_C, ripp)
        asp_acid_cs = find_atoms(ASP_ACID_C, ripp)
        mal_amino_cs = find_atoms(MAL_AMINO, ripp)

        for asp_acid_c in asp_acid_cs:
            if asp_acid_c in locations_c_to_domain:
                locations_c_to_domain.remove(asp_acid_c)
        for mal_amino_c in mal_amino_cs:
            if mal_amino_c in locations_c_to_domain:
                locations_c_to_domain.remove(mal_amino_c)
        assert len(locations_c_to_domain) == 1
        c_atom_to_domain = locations_c_to_domain[0]
    else:
        c_atom_to_domain = None
        for atom in ripp.graph:
            if atom.annotations.c2_acid:
                c_atom_to_domain = atom

    assert c_atom_to_domain

    oxygens = c_atom_to_domain.get_neighbours('O')

    hydroxyl_oxygen = None
    hydroxyl_bond = None

    hydrogen_bond = domain.bond_lookup[domain.atoms[1]][domain.atoms[2]]

    assert hydrogen_bond.has_neighbour('N')
    assert hydrogen_bond.has_neighbour('H')

    for oxygen in oxygens:
        bond = c_atom_to_domain.get_bond(oxygen)
        if bond.type == 'single' and oxygen.has_neighbour('H'):
            hydroxyl_oxygen = oxygen
            hydroxyl_bond = bond
            break

    assert hydroxyl_oxygen and hydroxyl_bond
    structure = condensation_with_reverse_numbering(ripp, domain, hydroxyl_bond, hydrogen_bond)[0]
    initialise_atom_attributes(structure)
    return structure


def attach_to_leader_ripp(ripp):
    initialise_atom_attributes(ripp)

    # Create domain
    domain = make_scaffold_peptide('Leader')
    n_atom_to_domain = None
    if ripp.find_substructures(AMINO_ACID_BACKBONE_N_TERMINUS):
        locations_n_to_domain = find_atoms(RIPP_N, ripp)
        for nitrogen in locations_n_to_domain:
            if len(nitrogen.get_neighbours('H')) == 2:
                n_atom_to_domain = nitrogen
                break
    assert n_atom_to_domain
    oxygen_bond = domain.bond_lookup[domain.atoms[0]][domain.atoms[1]]

    assert oxygen_bond.has_neighbour('O')
    hydrogen = n_atom_to_domain.get_neighbours('H')[0]
    nitrogen_bond = n_atom_to_domain.get_bond(hydrogen)
    assert nitrogen_bond
    structure = condensation(
        domain, ripp, oxygen_bond, nitrogen_bond)[0]
    initialise_atom_attributes(structure)

    return structure
