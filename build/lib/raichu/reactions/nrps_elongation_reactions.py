from pikachu.reactions.basic_reactions import condensation, hydrolysis
from pikachu.reactions.functional_groups import find_bonds, find_atoms
from raichu.class_domain import ATTRIBUTES
from raichu.central_chain_detection.label_central_chain import label_nrp_central_chain
from raichu.data.molecular_moieties import THIOESTERBOND, THIOESTER_CARBON, N_AMINO_ACID, B_N_AMINO_ACID
from raichu.reactions.general import label_rest_groups, initialise_atom_attributes, reset_nrp_annotations
from raichu.attach_to_domain import attach_to_domain_nrp


def nrps_elongation(amino_acid, chain_intermediate):
    """
    Returns the NRPS condensation product as a PIKAChU Structure object

    amino_acid: PIKAChU Structure object of the amino acid
    nrp_intermediate: PIKAChU Structure object of the NRP intermediate
    """

    # Initialize annotations in all atoms in amino acid
    initialise_atom_attributes(amino_acid)

    chain_intermediate, oh_bond = sulphur_to_hydroxyl(chain_intermediate)

    # Label rest groups
    label_rest_groups(amino_acid, chain_intermediate)

    # Ensure a central chain has been defined
    assert any(atom.annotations.in_central_chain for atom in chain_intermediate.graph)

    # Find amino group of amino acid or beta-amino acid
    n_atoms_aa = find_atoms(N_AMINO_ACID, amino_acid)
    if not n_atoms_aa:
        n_atoms_aa = find_atoms(B_N_AMINO_ACID, amino_acid)

    assert len(n_atoms_aa) == 1

    n_atom = n_atoms_aa[0]
    h_bond = None

    for bond in n_atom.bonds:
        for neighbour in bond.neighbours:
            if neighbour.type == 'H':
                h_bond = bond
                break

    assert h_bond

    # Label atoms in amino acid that end up in central peptide chain
    label_nrp_central_chain(amino_acid)

    # Reset chiral_c_ep and n_atom_nmeth AtomAnnotations attribute for all atoms in the intermediate;
    # This way, only the newly incorporated amino acid is labelled for tailoring
    reset_nrp_annotations(chain_intermediate)

    # Carry out condensation reaction using build-in PIKAChU function
    condensation_product = condensation(chain_intermediate, amino_acid, oh_bond, h_bond)[0]

    # Refresh condensation product
    condensation_product.refresh_structure(find_cycles=True)

    # Initialize annotations for atoms that don't have any yet
    initialise_atom_attributes(condensation_product)
    chain_intermediate = attach_to_domain_nrp(condensation_product)

    return chain_intermediate


def sulphur_to_hydroxyl(thioester_structure):
    """
    Identifies and removes the sulphur (and attached atoms/domains) in the
    thioester structure and replaces it with a hydroxyl group, creating a
    carboxylic acid group. Returns the product of this reaction, also as a
    PIKAChU Structure object.

    thioester_structure: PIKAChU Structure object of a thioester (R-C(S)=O)
    """
    # Find thioester bond in the input structure
    found_bonds_thioester = find_bonds(THIOESTERBOND, thioester_structure)
    found_carbon_thioester = find_atoms(THIOESTER_CARBON, thioester_structure)

    assert len(found_carbon_thioester) == 1
    assert len(found_bonds_thioester) == 1

    carbon_thioester = found_carbon_thioester[0]
    sh_bond = found_bonds_thioester[0]

    structures = hydrolysis(thioester_structure, sh_bond)
    combined = None
    for structure in structures:
        if carbon_thioester in structure.graph:
            combined = structure
            break

    assert combined

    oxygens = carbon_thioester.get_neighbours('O')

    oxygen_hydroxyl = None

    for oxygen in oxygens:
        if oxygen.get_bond(carbon_thioester).type == 'single' and oxygen.has_neighbour('H'):
            oxygen_hydroxyl = oxygen
            break

    assert oxygen_hydroxyl

    oh_bond = combined.bond_lookup[oxygen_hydroxyl][carbon_thioester]

    for atom in combined.graph:
        if not hasattr(atom.annotations, 'in_central_chain'):
            for attribute in ATTRIBUTES:
                atom.annotations.add_annotation(attribute, False)

    return combined, oh_bond
