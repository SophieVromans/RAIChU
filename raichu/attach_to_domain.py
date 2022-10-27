from raichu.class_domain import *
from pikachu.reactions.functional_groups import find_atoms, GroupDefiner, combine_structures
from pikachu.reactions.basic_reactions import condensation
from raichu.reactions.general import initialise_atom_attributes
from pikachu.general import read_smiles

POLYKETIDE_S = GroupDefiner('Sulphur atom polyketide', 'SC(C)=O', 0)
NRP_C = GroupDefiner('C atom to attach to PCP domain', 'NCC(O)=O', 2)
AMINO_ACID_BACKBONE = read_smiles('NCC(=O)O')
B_AMINO_ACID_BACKBONE = read_smiles('NCCC(=O)O')
ASP_ACID_C = GroupDefiner('Asp', 'NC(C=O)CC(O)=O', 5)
B_NRP_C = GroupDefiner('Beta C atom to attach to PCP domain', 'NCCC(O)=O', 3)
MAL_AMINO = GroupDefiner('malonyl_starter_amino', 'NC(=O)CC(O)=O', 4)


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
