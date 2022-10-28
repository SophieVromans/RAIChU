from pikachu.reactions.functional_groups import BondDefiner, GroupDefiner, find_atoms, find_bonds, combine_structures
from pikachu.reactions.basic_reactions import condensation, hydrolysis, internal_condensation
from pikachu.chem.chirality import same_chirality
from pikachu.chem.structure import Structure
from raichu.domain.domain_types import KRDomainSubtype, ERDomainSubtype
from raichu.reactions.general import initialise_atom_attributes
from pikachu.general import read_smiles
from raichu.data.attributes import ATTRIBUTES

FIRST_C = GroupDefiner('first_c', 'CC(S)=O', 1)
RECENT_ALPHA_C = GroupDefiner('recent_alpha_c', 'CCC(S)=O', 1)
RECENT_BETA_C = GroupDefiner('recent_beta_c', 'CCC(S)=O', 0)
RECENT_BETA_C_OH = GroupDefiner('recent_beta_c_oh', 'OCCC(S)=O', 0)
RECENT_ALPHA_METHYL_C = GroupDefiner('recent_alpha_methyl_c', 'CC(C)C(S)=O', 3)
RECENT_ELONGATION = BondDefiner('recent_elongation', 'O=C(C)CC(=O)S', 0, 1)
RECENT_ELONGATION_CC = BondDefiner('recent_elongation_CC', 'O=C(C)CC(=O)S', 1, 2)
RECENT_REDUCTION_COH = BondDefiner('recent_reduction_C-OH', 'OC(C)CC(=O)S', 0, 1)
RECENT_REDUCTION_MMAL_CHIRAL_C = GroupDefiner('recent_reduction_mmal_chiral_c', 'CCC(C(=O)S)C', 2)
RECENT_REDUCTION_C = GroupDefiner('recent_reduction_mal', 'OCCC(=O)S', 2)
RECENT_REDUCTION_CC = BondDefiner('recent_reduction_C-C', 'OCCC(=O)S', 1, 2)
RECENT_REDUCTION_CC_SHIFTED = BondDefiner('recent_reduction_C-C_shifted', 'CC(O)CC(S)=O', 0, 1)
RECENT_DEHYDRATION = BondDefiner('recent_dehydration', 'SC(C=CC)=O', 2, 3)
RECENT_EONYL_REDUCTION=BondDefiner('recent_eonyl_reduction',"CCCC(S)=O",2,3)
RECENT_EONYL_REDUCTION_CC=BondDefiner('recent_eonyl_reduction',"CCCC(S)=O",2,1)
RECENT_REDUCTION_CH=BondDefiner('recent_reduction_C-H', '[H]C(C)(O)CC(S)=O', 0, 1)
S_KR = GroupDefiner('C1 atom before KR reaction', 'SC(C)=O', 0)
ER_MMAL_CARBON = GroupDefiner('Chiral carbon atom after enoylreduction of mmal', 'SC(=O)C(C)CC', 3)
ER_S_CARBON = GroupDefiner('S-carbon atom after enoylreduction of mmal', 'SC(=O)C(C)CC', 1)
HYDROXYL_GROUP_TWO_MODULES_UPSTREAM_ALPHA_WITH_DOUBLE_BOND=BondDefiner("hydroxyl_group_two_module_upstream_alpha_with_double_bond","OC\C=C\CCC(S)=O",0,1)
HYDROXYL_GROUP_TWO_MODULES_UPSTREAM_BETA_WITH_DOUBLE_BOND=BondDefiner("hydroxyl_group_two_module_upstream_beta_with_double_bond","OCC\C=C\CCC(S)=O",0,1)
HYDROXYL_GROUP_TWO_MODULES_UPSTREAM_ALPHA_WITH_DOUBLE_BOND_SHIFTED=BondDefiner("hydroxyl_group_two_module_upstream_alpha_with_double_bond_shifted","O\C=C\CCCC(S)=O",0,1)
HYDROXYL_GROUP_TWO_MODULES_UPSTREAM_BETA_WITH_DOUBLE_BOND_SHIFTED=BondDefiner("hydroxyl_group_two_module_upstream_beta_with_double_bond_shifted","OC\C=C\CCCC(S)=O",0,1)
HYDROXYL_GROUP_TWO_MODULES_UPSTREAM_ALPHA=BondDefiner("hydroxyl_group_two_module_upstream_alpha","OCCCCCC(S)=O",0,1)
HYDROXYL_GROUP_TWO_MODULES_UPSTREAM_BETA=BondDefiner("hydroxyl_group_two_module_upstream_beta","OCCCCCCC(S)=O",0,1)


#STILL MISSING: E/Z-configured double bonds, E/Z-Gamma-beta-dehydrogenase

def ketoreduction(chain_intermediate: Structure, kr_type: KRDomainSubtype) -> Tuple[Structure, bool]:
    """
    Performs the ketoreductase reaction on the PKS chain intermediate, returns
    the reaction product as a PIKAChU Structure object

    chain_intermediate: PIKAChU Structure object, PKS chain intermediate just
    after a PKS elongation step
    """

    if kr_type.name == 'C1':
        return chain_intermediate, False

    chiral_c = None

    # Find chiral carbon atom
    chiral_c_locations = find_atoms(RECENT_REDUCTION_MMAL_CHIRAL_C, chain_intermediate)

    # Identify beta-ketone bond, identify O- and C-atom participating in bond

    beta_ketone_bonds = find_bonds(RECENT_ELONGATION, chain_intermediate)

    if not len(beta_ketone_bonds) == 1:
        return chain_intermediate, False

    beta_ketone_bond = beta_ketone_bonds[0]

    carbonyl_oxygen = beta_ketone_bond.get_neighbour('O')
    carbonyl_carbon = beta_ketone_bond.get_neighbour('C')

    if not carbonyl_carbon or not carbonyl_oxygen:
        raise ValueError("Beta-ketone bond does not have carbon or oxygen neighbour.")

    # Identify beta-ketone bond, identify O- and C-atom participating in bond
    if kr_type.name != "C2":

        # Change carbonyl bond to single bond
        beta_ketone_bond.make_single()

        # Add H atom to form hydroxyl group and another H to the C
        for neighbour in beta_ketone_bond.neighbours:
            chain_intermediate.add_atom('H', [neighbour])

        # Add atom annotations to new hydrogen atoms
        initialise_atom_attributes(chain_intermediate)

        chain_intermediate.refresh_structure()
        carbonyl_carbon = chain_intermediate.get_atom(carbonyl_carbon)

        h_atom = carbonyl_carbon.get_neighbour('H')
        o_atom = carbonyl_carbon.get_neighbour('O')
        c_atoms = carbonyl_carbon.get_neighbours('C')

        # Define the two carbons adjacent to the carbonyl atom as 'top_c' and 'bottom_c'

        top_cs = find_atoms(RECENT_REDUCTION_C, chain_intermediate)
        assert len(top_cs) == 1

        top_c = top_cs[0]
        assert top_c in c_atoms

        bottom_c = None
        for atom in c_atoms:
            if atom != top_c:
                bottom_c = atom

        assert bottom_c

        # Set the chirality of the carbonyl carbon

        counterclockwise_order = [h_atom, o_atom, top_c, bottom_c]

        if kr_type.name.startswith('A'):

            if same_chirality(counterclockwise_order, carbonyl_carbon.neighbours):
                carbonyl_carbon.chiral = 'clockwise'
            else:
                carbonyl_carbon.chiral = 'counterclockwise'

        elif kr_type.name.startswith('B'):
            if same_chirality(counterclockwise_order, carbonyl_carbon.neighbours):
                carbonyl_carbon.chiral = 'counterclockwise'
            else:
                carbonyl_carbon.chiral = 'clockwise'

        # If the type of KR domain is not known, the chirality of the carbonyl carbon is set to None
        elif kr_type.name == "UNKNOWN":
            carbonyl_carbon.chiral = None

        else:
            raise ValueError(f'KR domain of type {kr_type.name} is not supported by RAIChU or does not exist')

    # See if the previous elongation step was performed using methylmalonyl-CoA,
    # perform epimerization if required

    if chiral_c_locations:

        chiral_c = chiral_c_locations[0]

    if chiral_c:
        sulphur_locations = find_atoms(S_KR, chain_intermediate)
        assert len(sulphur_locations) == 1
        sulphur = sulphur_locations[0]
        assert sulphur

        c_1 = None
        for neighbour in sulphur.neighbours:
            if neighbour.annotations.in_central_chain and neighbour.type == 'C':
                c_1 = neighbour
                break
        assert c_1

        c_3 = None
        for neighbour in chiral_c.neighbours:
            if neighbour.annotations.in_central_chain and neighbour != c_1:
                c_3 = neighbour
        assert c_3

        hydrogen = None
        first_sidechain_atom = None

        for atom in chiral_c.neighbours:
            if atom != c_1 and atom != c_3 and atom.type != 'H' and \
                    not atom.annotations.in_central_chain:
                first_sidechain_atom = atom
            elif atom.type == 'H':
                hydrogen = atom

        assert hydrogen
        assert first_sidechain_atom

        clockwise_order = [first_sidechain_atom, c_1, c_3, hydrogen]

        if kr_type.name.endswith('1'):
            if same_chirality(clockwise_order, chiral_c.neighbours):
                chiral_c.chiral = 'clockwise'
            else:
                chiral_c.chiral = 'counterclockwise'
        elif kr_type.name.endswith('2'):
            if same_chirality(clockwise_order, chiral_c.neighbours):
                chiral_c.chiral = 'counterclockwise'
            else:
                chiral_c.chiral = 'clockwise'
        elif kr_type.name == "UNKNOWN":
            chiral_c.chiral = None
        else:
            raise ValueError(f'KR domain of type {kr_type.name} is not supported by RAIChU or does not exist')

    chain_intermediate.refresh_structure()
    return chain_intermediate, True


def dehydration(chain_intermediate: Structure) -> Tuple[Structure, bool]:
    """
    Performs the dehydratase reaction on the PKS chain intermediate, returns
    the reaction product as a PIKAChU Structure object

    chain_intermediate: PIKAChU Structure object, PKS chain intermediate where
    the beta ketone group has been recently reduced by the KR domain
    """

    # Find and define the atoms that participate in the bond changes
    co_bonds = find_bonds(RECENT_REDUCTION_COH, chain_intermediate)
    cc_bonds = find_bonds(RECENT_REDUCTION_CC, chain_intermediate)

    if not len(co_bonds) == 1 or not len(cc_bonds) == 1:
        return chain_intermediate, False

    c1 = None
    c2 = None
    o_oh = None

    co_bond = co_bonds[0]

    for neighbour in co_bond.neighbours:
        if neighbour.type == 'C':
            c2 = neighbour
        elif neighbour.type == 'O':
            o_oh = neighbour

    cc_bond = cc_bonds[0]

    for neighbour in cc_bond.neighbours:
        if neighbour != c2:
            c1 = neighbour

    # Remove hydroxyl group from c2
    chain_intermediate.break_bond(co_bond)

    # Remove H-atom from c1
    bond_to_break = None
    hydrogen = None

    for atom in chain_intermediate.graph:
        if atom == c1:
            for neighbour in atom.neighbours:
                if neighbour.type == 'H':
                    hydrogen = neighbour
                    for bond in neighbour.bonds:
                        bond_to_break = bond
                    break

    assert bond_to_break and hydrogen

    chain_intermediate.break_bond(bond_to_break)

    # Patch. When the bond is broken, the electron is not removed from
    # the orbitals of the C atom.
    c1 = chain_intermediate.get_atom(c1)

    # Make double c1=c2 bond
    chain_intermediate.bond_lookup[c1][c2].make_double()

    # Make bond between removed hydroxyl group and hydrogen atom to form water
    next_bond_nr = chain_intermediate.find_next_bond_nr()
    chain_intermediate.make_bond(hydrogen, o_oh, next_bond_nr)

    # Remove separate water structure from the Structure object
    structure_1, structure_2 = chain_intermediate.split_disconnected_structures()
    if c1 in structure_1.graph:
        chain_intermediate = structure_1
    else:
        chain_intermediate = structure_2

    chain_intermediate.refresh_structure()

    return chain_intermediate, True


def enoylreduction(chain_intermediate: Structure, er_subtype: ERDomainSubtype) -> Tuple[Structure, bool]:
    """
    Performs the enoylreductase reaction on the PKS chain intermediate, returns
    the reaction product as a PIKAChU Structure object

    chain_intermediate: PIKAChU Structure object, PKS chain intermediate where
    the beta ketone group has been recently reduced and dehydrated by the KR
    and ER domains, respectively
    er_subtype: ERDomainSubtype
    """

    for atom in chain_intermediate.graph:
        atom.set_connectivity()

    # Find double bond, change to single bond
    double_cc_bonds = find_bonds(RECENT_DEHYDRATION, chain_intermediate)
    if not len(double_cc_bonds) == 1:
        return chain_intermediate, False

    double_cc_bond = double_cc_bonds[0]
    double_cc_bond.make_single()

    # Add H-atom to the C atoms participating in the new single bond
    for neighbour in double_cc_bond.neighbours:
        chain_intermediate.add_atom('H', [neighbour])

    # Give annotation to added H-atom
    initialise_atom_attributes(chain_intermediate)
    chain_intermediate.refresh_structure()

    er_carbons = set(find_atoms(ER_MMAL_CARBON, chain_intermediate))
    if len(er_carbons) == 1 and er_subtype.name != "UNKNOWN":
        er_carbon = list(er_carbons)[0]
        if len(er_carbon.get_neighbours('H')) == 1 and len(er_carbon.get_neighbours('C')) == 3:
            s_carbons = set(find_atoms(ER_S_CARBON, chain_intermediate))
            assert len(s_carbons) == 1

            s_carbon = list(s_carbons)[0]
            bottom_carbon = None
            sidechain_carbon = None
            for carbon_neighbour in er_carbon.get_neighbours('C'):
                if carbon_neighbour.annotations.in_central_chain and carbon_neighbour != s_carbon:
                    bottom_carbon = carbon_neighbour
                elif not carbon_neighbour.annotations.in_central_chain:
                    sidechain_carbon = carbon_neighbour

            hydrogen_neighbour = er_carbon.get_neighbour('H')
            assert bottom_carbon and sidechain_carbon and hydrogen_neighbour

            clockwise_s_order = [s_carbon, hydrogen_neighbour, sidechain_carbon, bottom_carbon]
            has_s_chirality = same_chirality(clockwise_s_order, er_carbon.neighbours)
            if er_subtype.name == 'S' and has_s_chirality:
                er_carbon.chiral = 'clockwise'
            elif er_subtype.name == 'S' and not has_s_chirality:
                er_carbon.chiral = 'counterclockwise'
            elif er_subtype.name == 'R' and not has_s_chirality:
                er_carbon.chiral = 'clockwise'
            elif er_subtype.name == 'R' and has_s_chirality:
                er_carbon.chiral = 'counterclockwise'
            else:
                raise ValueError(f"RAIChU does not support ER domain subtype {er_subtype.name}")

    return chain_intermediate
#trans-AT-PKS tailoringReactions
def hydroxylation(target_atom, structure: Structure) -> Structure:
    """
    Returns the hydroxylated structure thats hydroxylated at the target atom.

    structure: PIKAChU Structure object
    target_atom:  PIKAChU atom object
    """

    hydroxyl_group = read_smiles('o')
    hydroxyl_group.add_attributes(ATTRIBUTES, boolean=True)
    oxygen = hydroxyl_group.atoms[0]
    hydrogen_1 = hydroxyl_group.atoms[1]
    bond_1 = oxygen.get_bond(hydrogen_1)
    hydrogen_2 = target_atom.get_neighbour('H')
    if not hydrogen_2:
        raise Exception("Can't oxidate this atom!")

    bond_2 = target_atom.get_bond(hydrogen_2)

    hydroxylated_structure = combine_structures([structure, hydroxyl_group])

    hydroxylated_structure.break_bond(bond_1)
    hydroxylated_structure.break_bond(bond_2)

    hydroxylated_structure.make_bond(oxygen, target_atom, hydroxylated_structure.find_next_bond_nr())
    hydroxylated_structure.make_bond(hydrogen_1, hydrogen_2,hydroxylated_structure.find_next_bond_nr())

    structures = hydroxylated_structure.split_disconnected_structures()

    for s in structures:
        if oxygen.nr in s.atoms:
            return s
def methylation(target_atom, structure: Structure) -> Structure:
    """
    Returns the structure thats methylated at the target atom.

    structure: PIKAChU Structure object
    target_atom:  PIKAChU atom object
    """


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

def find_OH_two_modules_upstream(structure: Structure):
    """
    Returns the O atom that is connected to the first hydroxygroup two modules upstream (also if chain contains double bonds)
    structure: PIKAChU Structure object
    """
    bonds= find_bonds(HYDROXYL_GROUP_TWO_MODULES_UPSTREAM_ALPHA,structure)
    if not bonds:
        bonds= find_bonds(HYDROXYL_GROUP_TWO_MODULES_UPSTREAM_BETA,structure)
    if not bonds:
        bonds= find_bonds(HYDROXYL_GROUP_TWO_MODULES_UPSTREAM_BETA_WITH_DOUBLE_BOND,structure)
    if not bonds:
        bonds= find_bonds(HYDROXYL_GROUP_TWO_MODULES_UPSTREAM_BETA_WITH_DOUBLE_BOND_SHIFTED,structure)
    return bonds

def smallest_cyclisation(structure: Structure) -> Structure:
    """
    Cyclises with the OH group two modules upstream to the smalles possible ring

    structure: PIKAChU Structure object
    """
    #reduce Ketogroup first
    structure=ketoreduction(structure, KRDomainSubtype(1))
    oh_bond=find_OH_two_modules_upstream(structure)
    assert len(oh_bond)>0, "No hydroxy group two modules upstream availiable"
    oh_bond=oh_bond[0]
    h_bond=find_bonds(RECENT_REDUCTION_CH, structure)[0]
    structure=internal_condensation(structure, oh_bond, h_bond)[0]
    return structure


def alpha_methyl_transferase(structure: Structure) -> Structure:
    """
    Returns the structure thats methylated at the alpha-c.

    structure: PIKAChU Structure object
    target_atom:  PIKAChU atom object
    """
    #find atom to add methylgroup
    alpha_c=find_atoms(RECENT_ALPHA_C, structure)[0]
    structure=methylation(alpha_c,structure)
    return structure

def beta_methyl_transferase(structure: Structure) -> Structure:
    """
    Returns the structure thats methylated at the beta-c.

    structure: PIKAChU Structure object
    target_atom:  PIKAChU atom object
    """
    #find atom to add methylgroup
    beta_c=find_atoms(RECENT_BETA_C, structure)[0]
    structure=methylation(beta_c,structure)
    return structure

def beta_hydroxy_methyl_transferase(structure: Structure) -> Structure:
    """
    Returns the structure thats methylated at the beta-oh-goup.

    structure: PIKAChU Structure object
    target_atom:  PIKAChU atom object
    """
    #find atom to add methylgroup

    beta_c_oh=find_atoms(RECENT_BETA_C_OH, structure)[0]
    structure=methylation(beta_c_oh,structure)
    return structure

def alpha_hydroxylase(structure: Structure) -> Structure:
    """
    Returns the structure thats hydroxylated at the alpha-c.

    structure: PIKAChU Structure object
    target_atom:  PIKAChU atom object
    """
    #find atom to add methylgroup
    alpha_c=find_atoms(RECENT_ALPHA_C, structure)[0]
    structure=hydroxylation(alpha_c,structure)
    return structure

def alpha_L_methyl_transferase(structure: Structure) -> Structure:
    """
    Returns the structure thats methylated at the alpha- c in L configuration.

    structure: PIKAChU Structure object
    target_atom:  PIKAChU atom object
    """
    #find atom to add methylgroup
    alpha_c=find_atoms(RECENT_ALPHA_C, structure)[0]
    structure=methylation(alpha_c,structure)
    structure.refresh_structure()
    h_atom = alpha_c.get_neighbour('H')
    c_atoms = alpha_c.get_neighbours('C')

    # Define the two carbons adjacent to the carbonyl atom as 'top_c' and 'bottom_c'
    methyl_c = find_atoms(RECENT_ALPHA_METHYL_C, structure)[0]
    assert methyl_c
    top_cs = find_atoms(FIRST_C, structure)
    assert len(top_cs) == 1

    top_c = top_cs[0]
    assert top_c in c_atoms

    bottom_c = None
    for atom in c_atoms:
        if atom != top_c and atom != methyl_c:
            bottom_c = atom

    assert bottom_c

    # Set the chirality of the carbonyl carbon

    counterclockwise_order = [h_atom, methyl_c, top_c, bottom_c]

    if same_chirality(counterclockwise_order, alpha_c.neighbours):
        alpha_c.chiral = 'clockwise'
    else:
        alpha_c.chiral = 'counterclockwise'
    alpha_c.chirality="counterclockwise"
    initialise_atom_attributes(structure)
    structure.refresh_structure()

    return structure

def gamma_beta_dehydratase(chain_intermediate: Structure) -> Structure:
    """
    Performs the dehydratase reaction on the PKS chain intermediate wuth shifted double bonds, returns
    the reaction product as a PIKAChU Structure object

    chain_intermediate: PIKAChU Structure object, PKS chain intermediate where
    the beta ketone group has been recently reduced by the KR domain
    """
    # Find and define the atoms that participate in the bond changes
    co_bonds = find_bonds(RECENT_REDUCTION_COH, chain_intermediate)
    cc_bonds = find_bonds(RECENT_REDUCTION_CC_SHIFTED, chain_intermediate)

    if not len(co_bonds) == 1 or not len(cc_bonds) == 1:
        return chain_intermediate

    c1 = None
    c2 = None
    o_oh = None

    co_bond = co_bonds[0]

    for neighbour in co_bond.neighbours:
        if neighbour.type == 'C':
            c2 = neighbour
        elif neighbour.type == 'O':
            o_oh = neighbour

    cc_bond = cc_bonds[0]

    for neighbour in cc_bond.neighbours:
        if neighbour != c2:
            c1 = neighbour

    # Remove hydroxyl group from c2
    chain_intermediate.break_bond(co_bond)

    # Remove H-atom from c1
    bond_to_break = None
    hydrogen = None

    for atom in chain_intermediate.graph:
        if atom == c1:
            for neighbour in atom.neighbours:
                if neighbour.type == 'H':
                    hydrogen = neighbour
                    for bond in neighbour.bonds:
                        bond_to_break = bond
                    break

    assert bond_to_break and hydrogen

    chain_intermediate.break_bond(bond_to_break)

    # Patch. When the bond is broken, the electron is not removed from
    # the orbitals of the C atom.
    c1 = chain_intermediate.get_atom(c1)

    # Make double c1=c2 bond
    chain_intermediate.bond_lookup[c1][c2].make_double()

    # Make bond between removed hydroxyl group and hydrogen atom to form water
    next_bond_nr = chain_intermediate.find_next_bond_nr()
    chain_intermediate.make_bond(hydrogen, o_oh, next_bond_nr)

    # Remove separate water structure from the Structure object
    structure_1, structure_2 = chain_intermediate.split_disconnected_structures()
    if c1 in structure_1.graph:
        chain_intermediate = structure_1
    else:
        chain_intermediate = structure_2

    chain_intermediate.refresh_structure()

    return chain_intermediate
