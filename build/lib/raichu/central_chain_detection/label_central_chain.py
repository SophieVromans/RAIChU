from pikachu.reactions.functional_groups import GroupDefiner, find_atoms, find_bonds
from pikachu.general import read_smiles
from raichu.data.attributes import ATTRIBUTES
from raichu.reactions.general import initialise_atom_attributes
from raichu.attach_to_domain import attach_to_domain_pk
from raichu.data.molecular_moieties import AMINO_FATTY_ACID, AMINO_ACID_BACKBONE, N_AMINO_ACID, C1_AMINO_ACID, \
    C2_AMINO_ACID, BETA_AMINO_ACID_BACKBONE, B_N_AMINO_ACID, B_C1_AMINO_ACID, B_C2_AMINO_ACID, B_C3_AMINO_ACID, \
    B_LEAVING_BOND, ACID_C1

POLYKETIDE_S = GroupDefiner('Sulphur atom polyketide', 'SC(C)=O', 0)


def label_pk_central_chain(pk_starter_unit):
    """Finds the the atoms in the central chain of the polyketide starter unit,
    sets the in_central_chain Atom object attribute to True/False accordingly,
    after which the pk_starter_unit Structure object is returned.

    pk_starter_unit: PIKAChU Structure object of the PK starter unit
    """

    initialise_atom_attributes(pk_starter_unit)
    central_chain = []
    visited = []

    # Check if structure is a polyketide, define sulphur attached to domain
    locations_sulphur = find_atoms(POLYKETIDE_S, pk_starter_unit)
    assert len(locations_sulphur) == 1
    sulphur_polyketide = locations_sulphur[0]
    central_chain.append(sulphur_polyketide)

    # Iterate over atom neighbours, check if the atom belongs to the central
    # carbon chain (i.e., not a methyl/ethyl/methoxy sidechain) and add to list
    for neighbour in sulphur_polyketide.neighbours:
        if neighbour.type == 'C':
            chain_carbon = neighbour
            end_carbon = False
            central_chain.append(chain_carbon)
            visited.append(chain_carbon)

            while not end_carbon:
                for next_atom in chain_carbon.neighbours:
                    ethyl_branch = False
                    inside_cycle = False
                    methyl_group = False

                    # Build lists of neighbouring atom types
                    if next_atom.type == 'C' and next_atom not in visited:

                        # Keep track of carbon neighbours

                        # Confirm carbon doesn't belong to ethyl sidebranch
                        types = []
                        c_neighbours = next_atom.get_neighbours('C')
                        if len(c_neighbours) == 2 and next_atom.get_neighbour('O') is None:

                            for atom in c_neighbours:
                                for neighbour in atom.neighbours:
                                    types.append(neighbour.type)

                            if types.count('H') == 4 and types.count('C') == 4:
                                ethyl_branch = True

                            for atom in c_neighbours:
                                if atom not in visited and atom != next_atom:
                                    ethyl_branch = False

                            visited.append(chain_carbon)

                        # Carbon in a terminal carboxylic acid group is the
                        # final carbon in the central chain
                        if len(next_atom.get_neighbours('O')) == 2:
                            central_chain.append(next_atom)
                            end_carbon = True
                            visited.append(chain_carbon)

                        # Confirm carbon is not part of cycle
                        if len(c_neighbours) == 2:
                            if all(atom.in_ring(pk_starter_unit) for atom in c_neighbours):
                                print(15)
                                visited.append(next_atom)
                                inside_cycle = True
                                end_carbon = True

                        # Case if carbon is part of a benzene ring
                        if len(c_neighbours) == 3 and len(next_atom.neighbours) == 3:
                            # central_chain.append(next_atom) #added 10/04/2021
                            visited.append(next_atom)
                            inside_cycle = True
                            end_carbon = True

                        # Case where the PK starter ends in two methyl branches
                        if len(next_atom.get_neighbours('C')) == 3 and len(next_atom.get_neighbours('H')) == 1:
                            nr_methyl_brances = 0
                            methyl_carbon = None
                            for c_neighbour in c_neighbours:
                                c_neighbour_neighbour_types = []
                                for c_neighbour_neighbour \
                                        in c_neighbour.neighbours:
                                    c_neighbour_neighbour_types.append(
                                        c_neighbour_neighbour.type)
                                if c_neighbour_neighbour_types.count('H') == 3:
                                    methyl_carbon = c_neighbour
                                    nr_methyl_brances += 1

                            if nr_methyl_brances == 2:
                                assert methyl_carbon is not None
                                methyl_group = True
                                visited.append(next_atom)
                                central_chain.append(next_atom)
                                visited.append(methyl_carbon)
                                central_chain.append(methyl_carbon)
                                end_carbon = True

                        # Confirm carbon doesn't belong to methyl sidebranch
                        if len(next_atom.get_neighbours('H')) == 3 or \
                                (len(next_atom.get_neighbours('H')) == 2 and len(next_atom.get_neighbours('*')) == 1):
                            methyl_group = True

                            # Terminal methyl group is not a side branch:
                            count_c_neighbours_not_visited = 0
                            for neighbouring_atom in chain_carbon.neighbours:
                                if neighbouring_atom.type == 'C' and neighbouring_atom not in visited:
                                    count_c_neighbours_not_visited += 1
                            if count_c_neighbours_not_visited < 2:
                                methyl_group = False
                                visited.append(chain_carbon)
                                for neighbour in chain_carbon.neighbours:
                                    if neighbour.type == 'C' and neighbour not in visited:
                                        next_neighbours = []
                                        for next_neighbour in neighbour.neighbours:
                                            next_neighbours.append(
                                                next_neighbour.type)
                                        if next_neighbours.count('H') == 3 or (
                                                next_neighbours.count(
                                                        'H') == 2 and (
                                                next_neighbours.count(
                                                        '*')) == 1):
                                            central_chain.append(neighbour)
                                end_carbon = True

                            visited.append(chain_carbon)

                        # Identification of central chain carbon atoms through
                        # exclusion
                        if not ethyl_branch and not methyl_group and not \
                                end_carbon and not inside_cycle:
                            central_chain += [next_atom]
                            chain_carbon = next_atom
                            visited.append(chain_carbon)

    # Set in_central_chain Atom attribute
    for atom in pk_starter_unit.graph:
        if atom in central_chain:
            atom.annotations.in_central_chain = True
        else:
            atom.annotations.in_central_chain = False


def label_nrp_central_chain(peptide, module_type='elongation'):
    initialise_atom_attributes(peptide)

    as_normal = False
    if module_type == 'starter':
        if peptide.find_substructures(AMINO_FATTY_ACID):
            acid = peptide
            label_acid_central_chain(acid)
        else:
            as_normal = True

    if as_normal or module_type == 'elongation':

        if peptide.find_substructures(AMINO_ACID_BACKBONE):
            n_atoms_aa = find_atoms(N_AMINO_ACID, peptide)
            c1_atoms_aa = find_atoms(C1_AMINO_ACID, peptide)
            c2_atoms_aa = find_atoms(C2_AMINO_ACID, peptide)

            assert len(n_atoms_aa) == 1
            assert len(c1_atoms_aa) == 1
            assert len(c2_atoms_aa) == 1

            n_atoms_aa[0].annotations.in_central_chain = True
            n_atoms_aa[0].annotations.n_atom_nmeth = True
            c1_atoms_aa[0].annotations.in_central_chain = True
            c1_atoms_aa[0].annotations.chiral_c_ep = True
            c2_atoms_aa[0].annotations.in_central_chain = True
        elif peptide.find_substructures(BETA_AMINO_ACID_BACKBONE):
            n_atoms_aa = find_atoms(B_N_AMINO_ACID, peptide)
            c1_atoms_aa = find_atoms(B_C1_AMINO_ACID, peptide)
            c2_atoms_aa = find_atoms(B_C2_AMINO_ACID, peptide)
            c3_atoms_aa = find_atoms(B_C3_AMINO_ACID, peptide)

            assert len(n_atoms_aa) == 1
            assert len(c1_atoms_aa) == 1
            assert len(c2_atoms_aa) == 1
            assert len(c3_atoms_aa) == 1

            if n_atoms_aa[0].in_ring(peptide) and c1_atoms_aa[0].in_ring(peptide) and c2_atoms_aa[0].in_ring(peptide):
                label_acid_central_chain(peptide)

            else:

                n_atoms_aa[0].annotations.in_central_chain = True
                n_atoms_aa[0].annotations.n_atom_nmeth = True
                c1_atoms_aa[0].annotations.in_central_chain = True
                c1_atoms_aa[0].annotations.chiral_c_ep = True
                c2_atoms_aa[0].annotations.in_central_chain = True
                c3_atoms_aa[0].annotations.in_central_chain = True

                bonds = find_bonds(B_LEAVING_BOND, peptide)
                assert len(bonds) == 1
                bond = bonds[0]
                if bond.atom_1.type == 'O':
                    bond.atom_1.annotations.leaving_oh_o = True

                    bond.atom_2.annotations.leaving_oh_h = True

                elif bond.atom_2.type == 'O':
                    bond.atom_2.annotations.leaving_oh_o = True

                    bond.atom_1.annotations.leaving_oh_h = True

        else:
            # If NRPS starter unit is a (fatty) acid instead of an amino acid
            acid = peptide
            label_acid_central_chain(acid)


def label_acid_central_chain(acid):
    c1_atom = None
    c2_atom = None
    oh_o_atom = None
    oh_h_atom = None
    visited = []

    c1_atoms_aa = find_atoms(ACID_C1, acid)
    if c1_atoms_aa:
        c1_atom = c1_atoms_aa[0]

    assert c1_atom

    for neighbour in c1_atom.neighbours:
        neighbours = []

        if neighbour.type == 'C':
            for next_atom in neighbour.neighbours:
                neighbours.append(next_atom.type)
            if neighbours.count('O') == 2 and len(neighbours) == 3:
                c2_atom = neighbour

                # Find the leaving -OH group and annotate it

                for bond in c2_atom.bonds:
                    if bond.type == 'single' and (bond.has_neighbour('O')):
                        if bond.atom_1.type == 'O':
                            bond.atom_1.annotations.leaving_oh_o = True
                            oh_o_atom = bond.atom_1

                            bond.atom_2.annotations.leaving_oh_h = True
                            oh_h_atom = bond.atom_2

                        elif bond.atom_2.type == 'O':
                            bond.atom_2.annotations.leaving_oh_o = True
                            oh_o_atom = bond.atom_2

                            bond.atom_1.annotations.leaving_oh_h = True
                            oh_h_atom = bond.atom_1

    assert c2_atom
    assert oh_o_atom
    assert oh_h_atom

    c1_atom.annotations.in_central_chain = True
    c2_atom.annotations.in_central_chain = True
    c2_atom.annotations.c2_acid = True

    visited.append(c1_atom)
    visited.append(c2_atom)

    starting_point = c1_atom
    endpoint = False
    ring_members = 0
    if starting_point.in_ring(acid):
        ring_members += 1

    # Add atoms in chain with exactly 2 non-H neighbours to the central chain
    while not endpoint:
        for neighbour in starting_point.get_non_hydrogen_neighbours():
            neighbour_in_ring = neighbour.in_ring(acid)
            if neighbour_in_ring:
                ring_members += 1
            neighbour_found = False
            if len(neighbour.get_non_hydrogen_bonds()) == 2 and neighbour.type == 'C' and neighbour not in visited \
                    and not (ring_members > 2 and neighbour_in_ring):
                neighbour.annotations.in_central_chain = True
                visited.append(neighbour)
                starting_point = neighbour
                neighbour_found = True

            if not neighbour_found and neighbour.type == 'C' and neighbour not in visited and \
                    len(neighbour.get_non_hydrogen_bonds()) >= 2 and not (ring_members > 2 and neighbour_in_ring):
                terminal_count = 0
                for atom in neighbour.neighbours:
                    if atom not in visited and len(atom.get_non_hydrogen_bonds()) == 1:
                        terminal_count += 1

                if len(neighbour.get_non_hydrogen_bonds()) - terminal_count == 2:
                    neighbour.annotations.in_central_chain = True
                    visited.append(neighbour)
                    starting_point = neighbour
                    neighbour_found = True

                if neighbour_found:
                    break

            if not neighbour_found and len(neighbour.get_non_hydrogen_bonds()) == 2 and neighbour not in visited \
                    and not (ring_members > 2 and neighbour_in_ring):
                neighbour.annotations.in_central_chain = True
                visited.append(neighbour)
                starting_point = neighbour
                break

            neighbours = []
            for nb in starting_point.get_non_hydrogen_neighbours():
                if nb not in visited:
                    neighbours.append(nb)

            if any(atom.in_ring(acid) for atom in neighbours) and neighbour_in_ring:
                endpoint = True
                for atom in neighbours:
                    if len(atom.get_non_hydrogen_bonds()) == 1 and atom.get_non_hydrogen_bonds()[0].type == 'single':
                        atom.annotations.in_central_chain = True
                break

            elif not any(len(atom.get_non_hydrogen_bonds()) == 2 for atom in neighbours):

                endpoint = True

                for atom in neighbours:
                    terminal_count = 0
                    for atom_2 in atom.neighbours:
                        if atom_2 not in visited and len(atom_2.get_non_hydrogen_bonds()) == 1:
                            terminal_count += 1

                    if len(atom.get_non_hydrogen_bonds()) - terminal_count == 2:
                        endpoint = False
                        break

                if endpoint:

                    # Add the last atom in the chain to the central chain
                    for atom in neighbours:
                        if len(atom.get_non_hydrogen_bonds()) == 1 and atom.get_non_hydrogen_bonds()[0].type == 'single':
                            atom.annotations.in_central_chain = True
                    break


if __name__ == "__main__":
    starter_units_antismash = ['SC(=O)CC', 'SC(CC(O)=O)=O', 'SC(CC(O)=O)=O',
                               'SC(C(C(O)=O)CC)=O', 'SC(C(C(O)=O)OC)=O',
                               'SC(C*)=O',
                               'SC(C(C)CC)=O', 'SC(C1C(CCC1)C(=O)O)=O',
                               'SC(C)=O',
                               'SC(C1=CC=CC=C1)=O', 'SC(CC(C)C)=O',
                               'SC(C(C(=O)O)CC[Cl])=O']
    starter_unit_error = 'SC(C1=CC=CC=C1)=O'
    struct = read_smiles(starter_unit_error)
    struct.add_attributes(ATTRIBUTES, boolean=True)
    struct = attach_to_domain_pk(struct)
    label_pk_central_chain(struct)
    struct.print_graph()
    for atom in struct.graph:
        print(atom, atom.annotations.in_central_chain)


