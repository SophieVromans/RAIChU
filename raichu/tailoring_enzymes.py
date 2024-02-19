from enum import Enum, unique
import itertools

from pikachu.drawing.drawing import Drawer
from raichu.reactions.general_tailoring_reactions import (
    proteolytic_cleavage,
    find_atoms_for_tailoring,
    remove_atom,
    single_bond_oxidation,
    addition,
    oxidative_bond_formation,
    epoxidation,
    double_bond_reduction,
    double_bond_shift,
    macrolactam_formation,
    cyclodehydration,
    change_chirality,
    excise_from_structure,
    reductive_bond_breakage,
)
from raichu.data.attributes import PRENYL_TRANSFERASE_SUBSTRATES_TO_SMILES
from raichu.data.molecular_moieties import (
    CO_BOND,
    CC_DOUBLE_BOND,
    PEPTIDE_BOND,
    CC_SINGLE_BOND,
    KETO_GROUP,
    C_CARBOXYL,
    ASPARTIC_ACID,
    GLUTAMIC_ACID,
    CYSTEINE,
    SERINE,
    THREONINE,
    REDUCED_SERINE,
    REDUCED_THREONINE,
    C1_AMINO_ACID_ATTACHED,
    ARGININE_SECONDARY_N1,
    ARGININE_SECONDARY_N2,
    ARGININE_SECONDARY_N3,
    ESTER_BOND,
)
from pikachu.reactions.functional_groups import find_atoms, find_bonds
from pikachu.reactions.basic_reactions import hydrolysis


@unique
class TailoringEnzymeType(Enum):
    # Group transfer reactions
    METHYLTRANSFERASE = 1
    C_METHYLTRANSFERASE = 2
    N_METHYLTRANSFERASE = 3
    O_METHYLTRANSFERASE = 4
    HYDROXYLASE = 5
    EPOXIDASE = 6
    PRENYLTRANSFERASE = 7
    ACETYLTRANSFERASE = 8
    ACYLTRANSFERASE = 9
    AMINOTRANSFERASE = 10
    HALOGENASE = 11
    METHYL_MUTASE = 12

    # Oxidoreduction
    DOUBLE_BOND_REDUCTASE = 13
    DOUBLE_BOND_ISOMERASE = 14
    DEHYDROGENASE = 15
    KETO_REDUCTION = 16
    ALCOHOL_DEHYDROGENASE = 17

    # Elimination
    PEPTIDASE = 18
    PROTEASE = 19
    MONOAMINE_OXIDASE = 20
    DEHYDRATASE = 21
    THREONINE_SERINE_DEHYDRATASE = 22
    DECARBOXYLASE = 23
    SPLICEASE = 24
    ARGINASE = 25

    # Cyclization
    OXIDATIVE_BOND_SYNTHASE = 26
    MACROLACTAM_SYNTHETASE = 27
    CYCLODEHYDRASE = 28
    LANTHIPEPTIDE_CYCLASE = 29
    LANTHIONINE_SYNTHETASE = 30
    THIOPEPTIDE_CYCLASE = 31

    # Epimerization
    AMINO_ACID_EPIMERASE = 32

    # Bond breakage
    HYDROLASE = 33
    REDUCTIVE_LYASE = 34

    @staticmethod
    def from_string(label: str) -> "TailoringEnzymeType":
        for value in TailoringEnzymeType:
            if str(value.name) == label:
                return value
        raise ValueError(f"Unknown tailoring enzyme: {label}")


class TailoringEnzyme:
    def __init__(
        self,
        gene_name,
        enzyme_type,
        modification_sites: list = None,
        substrate: str = None,
    ) -> None:
        self.gene_name = gene_name
        self.type = TailoringEnzymeType.from_string(enzyme_type)
        self.modification_sites = modification_sites
        self.substrate = substrate

    def do_tailoring(self, structure):
        """
        Performs tailoring reaction
        """
        if len(self.modification_sites) == 0:
            return structure
        if self.type.name == "HYDROXYLASE":
            for atom in self.modification_sites:
                if len(atom) == 0:
                    continue
                atom = atom[0]  # only one atom is hydroxylated at a time
                atom = structure.get_atom(atom)
                structure = addition(atom, "O", structure)
        elif self.type.name in [
            "METHYLTRANSFERASE",
            "C_METHYLTRANSFERASE",
            "N_METHYLTRANSFERASE",
            "O_METHYLTRANSFERASE",
        ]:
            for atom in self.modification_sites:
                if len(atom) == 0:
                    continue
                atom = atom[0]  # only one atom is methylated at a time
                atom = structure.get_atom(atom)
                structure = addition(atom, "C", structure)
        elif self.type.name == "PRENYLTRANSFERASE":
            for atom in self.modification_sites:
                if len(atom) == 0:
                    continue
                atom = atom[0]  # only one atom is methylated at a time
                atom = structure.get_atom(atom)
                if self.substrate not in PRENYL_TRANSFERASE_SUBSTRATES_TO_SMILES:
                    raise ValueError(
                        f"Not implemented prenyltransferase substrate: {self.substrate}"
                    )
                substrate = PRENYL_TRANSFERASE_SUBSTRATES_TO_SMILES[self.substrate]
                structure = addition(atom, substrate, structure)
        elif self.type.name == "ACETYLTRANSFERASE":
            for atom in self.modification_sites:
                if len(atom) == 0:
                    continue
                atom = atom[0]  # only one atom is methylated at a time
                atom = structure.get_atom(atom)
                structure = addition(atom, "[H]C(C)=O", structure)
        elif self.type.name == "ACYLTRANSFERASE":
            for atom in self.modification_sites:
                if len(atom) == 0:
                    continue
                atom = atom[0]  # only one atom is methylated at a time
                atom = structure.get_atom(atom)
                if self.substrate:
                    structure = addition(atom, self.substrate, structure)
        elif self.type.name == "OXIDATIVE_BOND_SYNTHASE":
            for atoms in self.modification_sites:
                if len(atoms) < 2:
                    continue
                atom1 = structure.get_atom(atoms[0])
                carbon_1 = structure.get_atom(atoms[1])
                structure = oxidative_bond_formation(atom1, carbon_1, structure)
        elif self.type.name == "EPOXIDASE":
            for atoms in self.modification_sites:
                if len(atoms) < 2:
                    continue
                atom1 = structure.get_atom(atoms[0])
                carbon_1 = structure.get_atom(atoms[1])
                structure = epoxidation(atom1, carbon_1, structure)
        elif self.type.name == "DOUBLE_BOND_REDUCTASE":
            for atoms in self.modification_sites:
                if len(atoms) < 2:
                    continue
                atom1 = structure.get_atom(atoms[0])
                carbon_1 = structure.get_atom(atoms[1])
                structure = double_bond_reduction(atom1, carbon_1, structure)
        elif self.type.name == "DEHYDROGENASE":
            for atoms in self.modification_sites:
                if len(atoms) < 2:
                    continue
                atom1 = structure.get_atom(atoms[0])
                carbon_1 = structure.get_atom(atoms[1])
                structure = single_bond_oxidation(atom1, carbon_1, structure)
        elif self.type.name == "DOUBLE_BOND_ISOMERASE":
            for atoms in self.modification_sites:
                if len(atoms) < 4:
                    continue
                old_double_bond_atom1 = structure.get_atom(atoms[0])
                old_double_bond_atom2 = structure.get_atom(atoms[1])
                new_double_bond_atom1 = structure.get_atom(atoms[2])
                new_double_bond_atom2 = structure.get_atom(atoms[3])
                if len(set(atoms)) == len(atoms):
                    raise ValueError(
                        "The bonds need to be adjacent to perform a dauble bond shift."
                    )
                structure = double_bond_shift(
                    structure,
                    old_double_bond_atom1,
                    old_double_bond_atom2,
                    new_double_bond_atom1,
                    new_double_bond_atom2,
                )
        elif self.type.name == "AMINOTRANSFERASE":
            for atom in self.modification_sites:
                if len(atom) == 0:
                    continue
                atom1 = atom[0]  # only one atom is modified at a time
                atom1 = structure.get_atom(atom1)
                oxygen = atom1.get_neighbour("O")
                structure = double_bond_reduction(atom1, oxygen, structure)
                oxygen = structure.get_atom(oxygen)
                structure = remove_atom(oxygen, structure)
                atom = structure.get_atom(atom)
                structure = addition(atom, "N", structure)
        elif self.type.name == "KETO_REDUCTION":
            for atom in self.modification_sites:
                if len(atom) == 0:
                    continue
                atom1 = atom[0]  # only one atom is modified at a time
                atom1 = structure.get_atom(atom1)
                if atom1.type != "O":
                    raise ValueError(
                        f"Can not perform KETO_REDUCTION on atom {atom1}, since there is no oxygen to be reduced."
                    )
                carbon_1 = atom1.get_neighbour("C")
                structure = double_bond_reduction(atom1, carbon_1, structure)
        elif self.type.name == "ALCOHOL_DEHYDROGENASE":
            for atom in self.modification_sites:
                if len(atom) == 0:
                    continue
                atom1 = atom[0]  # only one atom is modified at a time
                atom1 = structure.get_atom(atom1)
                if atom1.type != "O":
                    raise ValueError(
                        f"Can not perform ALCOHOL_DEHYDROGENASE on atom {atom1}, since there is no oxygen to be reduced."
                    )
                carbon_1 = atom1.get_neighbour("C")
                structure = single_bond_oxidation(atom1, carbon_1, structure)
        elif self.type.name == "DECARBOXYLASE":
            for atom in self.modification_sites:
                if len(atom) == 0:
                    continue
                atom1 = atom[0]  # only one atom is modified at a time
                atom1 = structure.get_atom(atom1)
                structure = remove_atom(atom1, structure)
        elif self.type.name == "DEHYDRATASE":
            for atoms in self.modification_sites:
                if len(atoms) < 2:
                    continue
                atom1 = structure.get_atom(atoms[0])
                carbon_1 = structure.get_atom(atoms[1])
                oxygen = atom1.get_neighbour("O")
                if not oxygen:
                    oxygen = carbon_1.get_neighbour("O")
                if not oxygen:
                    raise ValueError(
                        f"Can not perform DEHYDRATASE on atoms {atom1} and {carbon_1}, since there is no hydroxygroup to be removed."
                    )
                structure = remove_atom(oxygen, structure)
                if atom1.get_bond(carbon_1):
                    structure = single_bond_oxidation(atom1, carbon_1, structure)
                else:
                    structure = oxidative_bond_formation(atom1, carbon_1, structure)
        elif self.type.name == "MONOAMINE_OXIDASE":
            for atom in self.modification_sites:
                if len(atom) == 0:
                    continue
                atom1 = atom[0]  # only one atom is modified at a time
                atom1 = structure.get_atom(atom1)
                if atom1.type != "N":
                    raise ValueError(
                        f"Can not perform MONOAMINE_OXYDASE on atom {atom1}, since there is no nitrogen to be removed."
                    )
                structure = remove_atom(atom1, structure)
        elif self.type.name == "HALOGENASE":
            for atom in self.modification_sites:
                if len(atom) == 0:
                    continue
                atom1 = atom[0]  # only one atom is modified at a time
                atom1 = structure.get_atom(atom1)
                if self.substrate in ["F", "Cl", "Br", "I"]:
                    structure = addition(atom1, self.substrate, structure)
        elif self.type.name in ["PROTEASE", "PEPTIDASE"]:
            for atoms in self.modification_sites:
                if len(atoms) != 2:
                    continue
                atom1 = structure.get_atom(atoms[0])
                carbon_1 = structure.get_atom(atoms[1])
                structure = proteolytic_cleavage(atom1.get_bond(carbon_1), structure)
        elif self.type.name == "MACROLACTAM_SYNTHETASE":
            for atom in self.modification_sites:
                if len(atom) == 0:
                    continue
                atom1 = atom[0]  # only one atom is modified at a time
                atom1 = structure.get_atom(atom1)
                structure = macrolactam_formation(structure, atom1)
        elif self.type.name == "CYCLODEHYDRASE":
            for atom in self.modification_sites:
                if len(atom) == 0:
                    continue
                atom1 = atom[0]  # only one atom is modified at a time
                atom1 = structure.get_atom(atom1)
                break_all = False
                for neighbour_1 in atom1.neighbours:
                    if neighbour_1.type == "C":
                        for neighbour_2 in neighbour_1.neighbours:
                            if neighbour_2.type == "C":
                                nitrogen = neighbour_2.get_neighbour("N")
                                if nitrogen:
                                    carbon = nitrogen.get_neighbour("C")
                                    oxygen = carbon.get_neighbour("O")
                                    if carbon and nitrogen:
                                        structure = cyclodehydration(
                                            structure, atom1, oxygen
                                        )
                                        # breaking out of all loops
                                        break_all = True
                                        break
                        if break_all:
                            break
                else:
                    raise ValueError(
                        "No downstream amino acid for cyclodehydration availiable."
                    )
        elif self.type.name == "THIOPEPTIDE_CYCLASE":
            for atoms in self.modification_sites:
                if len(atoms) != 2:
                    continue
                carbon_1 = structure.get_atom(atoms[0])
                carbon_2_candidates = carbon_1.get_neighbours("C")
                carbon_2 = None
                for carbon_2_candidate in carbon_2_candidates:
                    if carbon_2_candidate.has_neighbour("N"):
                        carbon_2 = carbon_2_candidate
                assert carbon_2
                carbon_3 = structure.get_atom(atoms[1])
                carbon_4_candidates = carbon_3.get_neighbours("C")
                carbon_4 = None
                for carbon_4_candidate in carbon_4_candidates:
                    if carbon_4_candidate.has_neighbour("N"):
                        carbon_4 = carbon_4_candidate
                assert carbon_4
                nitrogen = carbon_4.get_neighbour("N")
                carbon_5 = [
                    carbon
                    for carbon in nitrogen.get_neighbours("C")
                    if carbon != carbon_4
                ][0]
                assert nitrogen, carbon_5
                oxygen = carbon_5.get_neighbour("O")
                structure = double_bond_reduction(carbon_1, carbon_2, structure)
                structure = double_bond_reduction(carbon_5, oxygen, structure)
                structure = double_bond_shift(
                    structure, carbon_3, carbon_4, carbon_4, nitrogen
                )
                structure = oxidative_bond_formation(carbon_1, carbon_3, structure)
                structure = oxidative_bond_formation(carbon_2, carbon_5, structure)
        elif self.type.name == "THREONINE_SERINE_DEHYDRATASE":
            for atom in self.modification_sites:
                if len(atom) == 0:
                    continue
                atom1 = atom[0]  # only one atom is modified at a time
                atom1 = structure.get_atom(atom1)
                carbon_1 = atom1.get_neighbour("C")
                carbon_2_candidates = carbon_1.get_neighbours("C")
                carbon_2 = None
                for carbon_2_candidate in carbon_2_candidates:
                    if carbon_2_candidate.has_neighbour("N"):
                        carbon_2 = carbon_2_candidate
                assert carbon_2
                structure = remove_atom(atom1, structure)
                structure = single_bond_oxidation(carbon_2, carbon_1, structure)
        elif self.type.name == "LANTHIPEPTIDE_CYCLASE":
            for atoms in self.modification_sites:
                if len(atoms) != 2:
                    continue
                if structure.get_atom(atoms[0]).type == "S":
                    sulfur = structure.get_atom(atoms[0])
                    carbon_1 = structure.get_atom(atoms[1])
                    carbon_2_candidates = carbon_1.get_neighbours("C")
                    carbon_2 = None
                    for carbon_2_candidate in carbon_2_candidates:
                        if carbon_2_candidate.has_neighbour("N"):
                            carbon_2 = carbon_2_candidate
                    assert carbon_2
                    structure = double_bond_reduction(carbon_1, carbon_2, structure)
                    structure = oxidative_bond_formation(sulfur, carbon_1, structure)
                if structure.get_atom(atoms[0]).type == "C":
                    carbon_1_1 = structure.get_atom(atoms[0])
                    carbon_1_2_candidates = carbon_1_1.get_neighbours("C")
                    carbon_1_2 = None
                    for carbon_1_2_candidate in carbon_1_2_candidates:
                        if carbon_1_2_candidate.has_neighbour("N"):
                            carbon_1_2 = carbon_1_2_candidate
                    assert carbon_1_2
                    carbon_2_1 = structure.get_atom(atoms[1])  # already cyclized carbon
                    carbon_2_2_candidates = carbon_2_1.get_neighbours("C")
                    carbon_2_2 = None
                    for carbon_2_2_candidate in carbon_2_2_candidates:
                        if carbon_2_2_candidate.has_neighbour("N"):
                            carbon_2_2 = carbon_2_2_candidate
                    assert carbon_2_2
                    if carbon_2_1.get_bond(carbon_2_2).type == "double":
                        structure = double_bond_reduction(
                            carbon_2_1, carbon_2_2, structure
                        )
                    structure = double_bond_reduction(carbon_1_1, carbon_1_2, structure)
                    structure = oxidative_bond_formation(
                        carbon_1_1, carbon_2_2, structure
                    )
        elif self.type.name == "LANTHIONINE_SYNTHETASE":
            for atoms in self.modification_sites:
                if len(atoms) != 2:
                    continue
                if structure.get_atom(atoms[0]).type == "S":
                    sulfur = structure.get_atom(atoms[0])
                    carbon_1 = structure.get_atom(atoms[1])
                    oxygen = carbon_1.get_neighbour("O")
                    if oxygen:
                        structure = remove_atom(oxygen, structure)
                    structure = oxidative_bond_formation(sulfur, carbon_1, structure)
                if structure.get_atom(atoms[0]).type == "C":
                    carbon = structure.get_atom(atoms[0])
                    oxygen_2 = carbon.get_neighbour("O")
                    carbon_1 = structure.get_atom(atoms[1])
                    oxygen = carbon_1.get_neighbour("O")
                    carbon_2_candidates = carbon_1.get_neighbours("C")
                    carbon_2 = None
                    for carbon_2_candidate in carbon_2_candidates:
                        if carbon_2_candidate.has_neighbour("N"):
                            carbon_2 = carbon_2_candidate
                    assert carbon_2
                    if oxygen:
                        structure = remove_atom(oxygen, structure)
                    if oxygen_2:
                        structure = remove_atom(oxygen_2, structure)
                    structure = oxidative_bond_formation(carbon, carbon_2, structure)
        elif self.type.name == "AMINO_ACID_EPIMERASE":
            # atom needs to be amino acid alpha atom
            for atom in self.modification_sites:
                if len(atom) == 0:
                    continue
                carbon = structure.get_atom(atom[0])
                assert carbon.type == "C"
                structure = change_chirality(carbon, structure)

        elif self.type.name == "SPLICEASE":
            for atoms in self.modification_sites:
                if len(atoms) != 2:
                    continue
                carbon_1 = structure.get_atom(atoms[0])
                carbon_2 = structure.get_atom(atoms[1])
                structure = excise_from_structure(carbon_1, carbon_2, structure)
        elif self.type.name == "HYDROLASE":
            for atoms in self.modification_sites:
                if len(atoms) != 2:
                    continue
                hetero_atom = structure.get_atom(atoms[0])
                carbon = structure.get_atom(atoms[1])  # gets oh attached
                bond = hetero_atom.get_bond(carbon)
                structure = hydrolysis(structure, bond)
        elif self.type.name == "REDUCTIVE_LYASE":
            for atoms in self.modification_sites:
                if len(atoms) != 2:
                    continue
                carbon_1 = structure.get_atom(atoms[0])
                carbon_2 = structure.get_atom(atoms[1])
                structure = reductive_bond_breakage(carbon_1, carbon_2, structure)
        elif self.type.name == "ARGINASE":
            # atom needs to be arginine secondary nitrogen
            for atom in self.modification_sites:
                if len(atom) == 0:
                    continue
                nitrogen = structure.get_atom(atom[0])
                assert nitrogen.type == "N"
                carbon = [
                    carbon
                    for carbon in nitrogen.get_neighbours("C")
                    if [atom.type for atom in carbon.neighbours].count("N") == 3
                ][0]
                bond = nitrogen.get_bond(carbon)
                assert bond
                if bond.type == "double":
                    bond.make_single()
                structure.break_bond(bond)
                structure_1, structure_2 = structure.split_disconnected_structures()
                if nitrogen in structure_1.graph:
                    structure = structure_1
                else:
                    structure = structure_2
                structure.add_atom("H", [nitrogen])
                structure.refresh_structure(find_cycles=True)
        elif self.type.name == "METHYL_MUTASE":
            for atoms in self.modification_sites:
                if len(atoms) != 2:
                    continue
                transferred_c = atoms[0]
                # Assert its actually a methyl group
                assert [atom.type for atom in transferred_c.neighbours].count("H") == 3
                source = None
                source = [
                    atom for atom in transferred_c.neighbours if atom.type != "H"
                ][0]
                assert source
                destination_c = atoms[1]
                assert destination_c.has_neighbour("H")
                assert transferred_c.type == "C" and destination_c.type == "C"

                structure = reductive_bond_breakage(source, transferred_c, structure)
                structure = addition(destination_c, "C", structure)
                structure.refresh_structure(find_cycles=True)

        return structure

    def get_possible_sites(self, structure, out_file=None):
        possible_sites = []
        if self.type.name in [
            "HYDROXYLASE",
        ]:
            possible_sites.extend(
                [
                    [atom]
                    for atom in find_atoms_for_tailoring(structure, "C")
                    if atom.has_neighbour("H")
                ]
            )
        elif self.type.name in [
            "C_METHYLTRANSFERASE",
            "N_METHYLTRANSFERASE",
            "O_METHYLTRANSFERASE",
        ]:
            atom = self.type.name.split("_")[0]
            possible_sites.extend(
                [
                    [atom]
                    for atom in find_atoms_for_tailoring(structure, atom)
                    if atom.has_neighbour("H")
                ]
            )
        elif self.type.name in [
            "METHYLTRANSFERASE",
            "PRENYLTRANSFERASE",
            "ACETYLTRANSFERASE",
            "ACYLTRANSFERASE",
            "OXIDATIVE_BOND_SYNTHASE",
            "HALOGENASE",
            "METHYL_MUTASE",
        ]:
            possible_sites.extend(
                [
                    [atom]
                    for atom in find_atoms_for_tailoring(structure, "C")
                    if atom.has_neighbour("H")
                ]
            )
            possible_sites.extend(
                [
                    [atom]
                    for atom in find_atoms_for_tailoring(structure, "N")
                    if atom.has_neighbour("H")
                ]
            )
            possible_sites.extend(
                [
                    [atom]
                    for atom in find_atoms_for_tailoring(structure, "O")
                    if atom.has_neighbour("H")
                ]
            )
            possible_sites.extend(
                [
                    [atom]
                    for atom in find_atoms_for_tailoring(structure, "S")
                    if atom.has_neighbour("H")
                ]
            )
        elif self.type.name in ["SPLICEASE"]:
            possible_sites.extend(
                [[atom] for atom in find_atoms_for_tailoring(structure, "C")]
            )
            possible_sites.extend(
                [[atom] for atom in find_atoms_for_tailoring(structure, "N")]
            )
            possible_sites.extend(
                [[atom] for atom in find_atoms_for_tailoring(structure, "O")]
            )
            possible_sites.extend(
                [[atom] for atom in find_atoms_for_tailoring(structure, "S")]
            )

        elif self.type.name in ["EPOXIDASE", "DOUBLE_BOND_REDUCTASE"]:
            peptide_bonds = find_bonds(CC_DOUBLE_BOND, structure)
            for bond in peptide_bonds:
                possible_sites.append(bond.neighbours)

        elif self.type.name == "DEHYDROGENASE":
            peptide_bonds = find_bonds(CC_SINGLE_BOND, structure)
            for bond in peptide_bonds:
                possible_sites.append(bond.neighbours)

        elif self.type.name == "REDUCTIVE_LYASE":
            possible_sites.extend(
                [
                    bond.neighbours
                    for bond in structure.bonds.values()
                    if "H" not in [atom.type for atom in bond.neighbours]
                ]
            )

        elif self.type.name == "DOUBLE_BOND_ISOMERASE":
            peptide_bonds = find_bonds(CC_DOUBLE_BOND, structure)
            for bond in peptide_bonds:
                neighbouring_bonds = bond.get_neighbouring_bonds()
                for neighbouring_bond in neighbouring_bonds:
                    if not "H" in [atom.type for atom in neighbouring_bond.neighbours]:
                        possible_sites.append(
                            bond.neighbours + neighbouring_bond.neighbours
                        )

        elif self.type.name == "AMINOTRANSFERASE":
            possible_sites.extend(
                [[atom] for atom in find_atoms(KETO_GROUP, structure)]
            )

        elif self.type.name == "KETO_REDUCTION":
            oxygens = find_atoms(KETO_GROUP, structure)
            possible_sites.extend([[atom] for atom in oxygens])

        elif self.type.name == "ALCOHOL_DEHYDROGENASE":
            possible_sites.extend(
                [[atom] for atom in find_atoms_for_tailoring(structure, "O")]
            )

        elif self.type.name == "DECARBOXYLASE":
            possible_sites.extend(
                [[atom] for atom in find_atoms(C_CARBOXYL, structure)]
            )

        elif self.type.name == "DEHYDRATASE":
            co_bonds = find_bonds(CO_BOND, structure)
            for co_bond in co_bonds:
                neighbouring_bonds = co_bond.get_neighbouring_bonds()
                for neighbouring_bond in neighbouring_bonds:
                    if (
                        not "H" in [atom.type for atom in neighbouring_bond.neighbours]
                        and neighbouring_bond.type == "single"
                    ):
                        for neighbouring_atom in neighbouring_bond.neighbours:
                            if neighbouring_atom != co_bond.get_neighbour(
                                "C"
                            ) and neighbouring_atom.has_neighbour("H"):
                                possible_sites.append(neighbouring_bond.neighbours)

        elif self.type.name == "MONOAMINE_OXIDASE":
            n_atoms_with_one_h = find_atoms_for_tailoring(structure, "N")
            for n_atom in n_atoms_with_one_h:
                if [atom.type for atom in n_atom.neighbours].count("H") == 2:
                    possible_sites.append([n_atom])
        elif self.type.name in ["PROTEASE", "PEPTIDASE"]:
            peptide_bonds = find_bonds(PEPTIDE_BOND, structure)
            for bond in peptide_bonds:
                possible_sites.append(bond.neighbours)
        elif self.type.name in ["HYDROLASE"]:
            ester_bonds = find_bonds(ESTER_BOND, structure)
            for bond in ester_bonds:
                possible_sites.append(bond.neighbours)
        elif self.type.name == "MACROLACTAM_SYNTHETASE":
            asp_glu_oxygen = find_atoms(ASPARTIC_ACID, structure) + find_atoms(
                GLUTAMIC_ACID, structure
            )
            possible_sites.extend(asp_glu_oxygen)
        elif self.type.name == "CYCLODEHYDRASE":
            cys_ser_thr_x = (
                find_atoms(CYSTEINE, structure)
                + find_atoms(SERINE, structure)
                + find_atoms(THREONINE, structure)
            )
            possible_sites.extend([[atom] for atom in cys_ser_thr_x])
        elif self.type.name == "THREONINE_SERINE_DEHYDRATASE":
            ser_thr_x = find_atoms(SERINE, structure) + find_atoms(THREONINE, structure)
            possible_sites.extend([[atom] for atom in ser_thr_x])
        elif self.type.name == "LANTHIPEPTIDE_CYCLASE":
            cys_x = find_atoms(CYSTEINE, structure)
            ser_thr_c = find_atoms(REDUCED_SERINE, structure) + find_atoms(
                REDUCED_THREONINE, structure
            )
            combinations = [list(t) for t in itertools.product(cys_x, ser_thr_c)]
            combinations.extend(
                [list(t) for t in itertools.product(ser_thr_c, ser_thr_c)]
            )
            possible_sites.extend(combinations)
        elif self.type.name == "LANTHIONINE_SYNTHETASE":
            cys_x = find_atoms(CYSTEINE, structure)
            ser_thr_x = find_atoms(SERINE, structure) + find_atoms(THREONINE, structure)
            ser_thr_c = [atom.get_neighbour("C") for atom in ser_thr_x]
            combinations = [list(t) for t in itertools.product(cys_x, ser_thr_c)]
            combinations.extend(
                [list(t) for t in itertools.product(ser_thr_c, ser_thr_c)]
            )
            possible_sites.extend(combinations)
        elif self.type.name == "AMINO_ACID_EPIMERASE":
            alpha_cs_amino_acid_backbone = find_atoms(C1_AMINO_ACID_ATTACHED, structure)
            possible_sites.extend([[atom] for atom in alpha_cs_amino_acid_backbone])
        elif self.type.name == "ARGINASE":
            arginine_n1 = find_atoms(ARGININE_SECONDARY_N_1, structure)
            arginine_n2 = find_atoms(ARGININE_SECONDARY_N_2, structure)
            arginine_n3 = find_atoms(ARGININE_SECONDARY_N_3, structure)
            arginine_n = arginine_n1 + arginine_n2 + arginine_n3
            possible_sites.extend([[atom] for atom in arginine_n])
        elif self.type.name == "THIOPEPTIDE_CYCLASE":
            ser_thr_c = find_atoms(REDUCED_SERINE, structure) + find_atoms(
                REDUCED_THREONINE, structure
            )
            possible_sites.extend(
                [list(t) for t in itertools.product(ser_thr_c, ser_thr_c)]
            )

        if out_file:
            drawing = Drawer(structure)
            site_list = []

            for possible_site in possible_sites:
                if type(possible_site) == list:
                    site_list += possible_site
                else:
                    site_list.append(possible_site)
            drawing.write_svg(out_file, numbered_atoms=site_list)
        return possible_sites
