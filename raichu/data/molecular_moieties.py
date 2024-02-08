from raichu.class_domain import make_scaffold_domain
from pikachu.general import read_smiles
from pikachu.reactions.functional_groups import BondDefiner, GroupDefiner, combine_structures
from raichu.reactions.general import initialise_atom_attributes


_STARTER_PKS_TO_SULPHUR = {"PROPIONYL_COA": 0,
                           "ACETYL_COA": 0,
                           "BENZOYL_COA": 0,
                           "METHYL_BUTYRYL_COA_3": 0,
                           "METHYL_BUTYRYL_COA_2": 0,
                           "TRANS_CYCLOPENTANE_DICARBOXYL_COA": 0,
                           "CYCLOHEXANE_CARBOXYL_COA": 0,
                           "HYDROXY_MALONYL_COA_2": 0,
                           "HYDROXY_MALONYL_COA_2R": 0,
                           "HYDROXY_MALONYL_COA_2S": 0,
                           "CHLOROETHYL_MALONYL_COA": 0,
                           "ISOBUTYRYL_COA": 0,
                           "GLYCINE": 0,
                           "HYDROXY_PROPENOYL_COA_3_23E": 5,
                           "HYDROXY_BUTENOYL_COA_3_23E": 5,
                           "DIHYDROXY_BUTANOLYL_COA_2R3": 8,
                           "TRIHYDROXY_PROPANOLYL_COA_233": 7,
                           "O_METHYLACETYL_COA": 3,
                           "HYDROXY_PROPENOYL_COA_3_23Z": 5,
                           "OXOMALONYL_COA_2": 6,
                           "METHYL_HYDROXY_PROPENOYL_COA_2_3_23Z": 6,
                           "DIHYDROXY_BUTANOLYL_COA_23": 7,
                           "DIHYDROXY_BUTANOLYL_COA_2S3S": 9,
                           "HEPTATRIENOYL_COA": 0,
                           "HYDROXYPROPIONYL_COA_2R": 5,
                           "DIHYDROXY_PROPANOLYL_COA_33": 6,
                           "LACTYL_COA": 5,
                           "PHENYLACETYLCOA": 0,
                           "METHOXYFORMYL_COA": 3
                           }


_STARTER_PKS_TO_BACKBONE = {"PROPIONYL_COA": [1, 3, 4],
                            "ACETYL_COA": [1, 2],
                            "BENZOYL_COA": [1, 2, 3],
                            "METHYL_BUTYRYL_COA_3": [1, 3, 4, 5],
                            "METHYL_BUTYRYL_COA_2": [1, 3, 5, 6],
                            "TRANS_CYCLOPENTANE_DICARBOXYL_COA": [1, 2, 4, 9, 11],
                            "CYCLOHEXANE_CARBOXYL_COA": [1, 3, 4, 5],
                            "HYDROXY_MALONYL_COA_2": [1, 3, 4, 6],
                            "HYDROXY_MALONYL_COA_2R": [1, 3, 5, 7],
                            "HYDROXY_MALONYL_COA_2S": [1, 3, 5, 7],
                            "CHLOROETHYL_MALONYL_COA": [1, 2, 6, 7, 8],
                            "ISOBUTYRYL_COA": [1, 3, 4],
                            "GLYCINE": [1, 2, 3],
                            "HYDROXY_PROPENOYL_COA_3_23E": [1, 2, 3, 4],
                            "HYDROXY_BUTENOYL_COA_3_23E": [1, 2, 3, 4],
                            "DIHYDROXY_BUTANOLYL_COA_2R3": [0, 1, 3, 6],
                            "TRIHYDROXY_PROPANOLYL_COA_233": [0, 1, 3, 5],
                            "O_METHYLACETYL_COA": [0, 1, 2],
                            "HYDROXY_PROPENOYL_COA_3_23Z": [1, 3, 4],
                            "OXOMALONYL_COA_2": [0, 1, 3, 5],
                            "METHYL_HYDROXY_PROPENOYL_COA_2_3_23Z": [1, 3, 5],
                            "DIHYDROXY_BUTANOLYL_COA_23": [0, 1, 3, 5],
                            "DIHYDROXY_BUTANOLYL_COA_2S3S": [0, 1, 4, 7],
                            "HEPTATRIENOYL_COA": [1, 3, 4, 5, 6, 7, 8],
                            "HYDROXYPROPIONYL_COA_2R": [0, 1, 4],
                            "DIHYDROXY_PROPANOLYL_COA_33": [0, 1, 3, 4],
                            "LACTYL_COA": [0, 1, 4],
                            "PHENYLACETYLCOA": [1, 3, 4, 5, 6],
                            "METHOXYFORMYL_COA": [0, 1, 2]
                            }


class PksStarterUnit:
    def __init__(self, name, smiles, central_chain_indices):
        self.name = name
        self.smiles = smiles
        self.structure = read_smiles(self.smiles)
        initialise_atom_attributes(self.structure)
        self.structure.refresh_structure()
        self.sulphur = self.structure.atoms[_STARTER_PKS_TO_SULPHUR[self.name]]

        for atom_index in central_chain_indices:
            self.structure.atoms[atom_index].annotations.in_central_chain = True

    def attach_to_acp(self):
        domain = make_scaffold_domain('ACP')
        sh_bond = domain.bond_lookup[domain.atoms[1]][domain.atoms[2]]
        hydrogen = domain.atoms[2]
        sulphur_1 = domain.atoms[1]
        carbon = self.sulphur.get_neighbour('C')
        sc_bond = self.sulphur.get_bond(carbon)

        structure = combine_structures([domain, self.structure])
        structure.break_bond(sh_bond)
        structure.break_bond(sc_bond)

        structure.make_bond(hydrogen, self.sulphur, structure.find_next_bond_nr())
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



class PksElongationUnit:
    def __init__(self, name, smiles, c_idx, carbonyl_idx):
        """
        name: str, elongation unit name
        smiles: str, elongation unit SMILES
        c_idx: int, index of carbon that will be attached to the thioester carbon of the chain intermediate
        carbonyl_idx: int, index of the carbon of the carbonyl group
        """
        self.name = name
        self.smiles = smiles
        self.structure = read_smiles(self.smiles)
        initialise_atom_attributes(self.structure)
        self.structure.refresh_structure()

        self.c_to_pk_intermediate = None
        self.c_to_s = None

        for atom in self.structure.graph:
            if atom.nr == c_idx:
                self.c_to_pk_intermediate = atom
            elif atom.nr == carbonyl_idx:
                self.c_to_s = atom

        assert self.c_to_pk_intermediate is not None and self.c_to_s is not None

        self.c_to_pk_intermediate.annotations.in_central_chain = True
        self.c_to_s.annotations.in_central_chain = True

        if self.smiles == 'O=CCC':
            # TODO: CHECK CHIRALITY ACTUALLY WORKS!!!!
            self.c_to_pk_intermediate.chiral = 'clockwise'


def make_starter_monomer(name, smiles):
    assert name in _STARTER_PKS_TO_BACKBONE

    monomer = PksStarterUnit(name, smiles, _STARTER_PKS_TO_BACKBONE[name])
    return monomer


def make_elongation_monomer(name):
    if name == 'MALONYL_COA':
        monomer = PksElongationUnit('Malonyl CoA', 'CC=O', 0, 1)
    elif name == 'METHYLMALONYL_COA':
        monomer = PksElongationUnit('Methylmalonyl CoA', 'O=CCC', 2, 1)
    elif name == 'METHOXYMALONYL_ACP':
        monomer = PksElongationUnit("Methoxymalonyl ACP", 'O=CCOC', 2, 1)
    elif name == 'ETHYLMALONYL_COA':
        monomer = PksElongationUnit("Ethylmalonyl CoA", 'O=CCCC', 2, 1)
    elif name == 'WILDCARD':
        monomer = PksElongationUnit("Wildcard", 'O=CC*', 2, 1)
    else:
        raise ValueError(f"{name} not recognised by RAIChU as PKS elongation unit.")

    return monomer


ESTER_BOND = BondDefiner('thioester_bond', 'COC=O', 0, 1)
THIOESTERBOND = BondDefiner('thioester_bond', 'SC(C)=O', 0, 1)
THIOESTERBOND_OXYGEN_INSERTED = BondDefiner(
    'thioester_bond_oxygen_inserted', 'SC(O)=O', 0, 1)
THIOESTERBOND_REDUCED_O = BondDefiner('thioester_bond_reduced_O', 'SC(C)O', 0, 1)
THIOESTER_CARBON = GroupDefiner('thioester carbon', 'SC(C)=O', 1)
THIOESTER_CARBON_REDUCED_O = GroupDefiner('thioester carbon reduced oxygen', 'SC(C)O', 1)
THIOESTER_CARBON_OXYGEN_INSERTED = GroupDefiner('thioester carbon oxygen inserted', 'SC(O)=O', 1)
LEAVING_OH_BOND = BondDefiner('Leaving -OH group bond', 'C(=O)(O)CN', 0, 2)
N_AMINO_ACID = GroupDefiner('Nitrogen atom amino acid', 'NCC(=O)O', 0)
C1_AMINO_ACID = GroupDefiner('C1 atom amino acid', 'NCC(=O)O', 1)
C2_AMINO_ACID = GroupDefiner('C2 atom amino acid', 'NCC(=O)O', 2)
N_AMINO_ACID_ATTACHED = GroupDefiner('Nitrogen atom amino acid', 'NCC=O', 0)
C1_AMINO_ACID_ATTACHED = GroupDefiner('C1 atom amino acid', 'NCC=O', 1)
C2_AMINO_ACID_ATTACHED = GroupDefiner('C2 atom amino acid', 'NCC=O', 2)
OH_AMINO_ACID = GroupDefiner('OH atom amino acid', 'NCC(O)=O', 3)
AMINO_ACID_BACKBONE = read_smiles('NCC(=O)O')
AMINO_ACID_BACKBONE_ATTACHED = read_smiles('NCC=O')
BETA_AMINO_ACID_BACKBONE = read_smiles('NCCC(=O)O')
B_N_AMINO_ACID = GroupDefiner('Nitrogen atom beta amino acid', 'NCCC(=O)O', 0)
B_C1_AMINO_ACID = GroupDefiner('C1 atom beta amino acid', 'NCCC(=O)O', 1)
B_C2_AMINO_ACID = GroupDefiner('C2 atom beta amino acid', 'NCCC(=O)O', 2)
B_C3_AMINO_ACID = GroupDefiner('C3 atom beta amino acid', 'NCCC(=O)O', 3)
B_LEAVING_BOND = BondDefiner('beta leaving bond', 'NCCC(=O)O', 3, 5)
AMINO_FATTY_ACID = read_smiles('CCCCCCCC(=O)O')
ACID_C1 = GroupDefiner('C1 atom (fatty) acid', 'CC(=O)O', 0)
ACID_C2 = GroupDefiner('C2 atom (fatty) acid', 'CC(=O)O', 1)
COABOND = BondDefiner('CoA_bond', 'CC(NCCC(NCCSC)=O)=O', 8, 9)
CYCLIC_CYS_C1_ATTACHED = GroupDefiner('oxidation_target', r"SC(=O)C1CSC=N1", 3)
CYCLIC_CYS_C2_ATTACHED = GroupDefiner('oxidation_target', r"SC(=O)C1CSC=N1", 4)
CYCLIC_SER_C1_ATTACHED = GroupDefiner('oxidation_target', r"SC(=O)C1COC=N1", 3)
CYCLIC_SER_C2_ATTACHED = GroupDefiner('oxidation_target', r"SC(=O)C1COC=N1", 4)

PEPTIDE_BOND = BondDefiner('peptide_bond', 'CC(=O)NC', 1, 3)
CC_DOUBLE_BOND = BondDefiner('cc_double_bond', 'C=C', 0, 1)
CO_DOUBLE_BOND = BondDefiner('co_double_bond', 'C=O', 0, 1)
CC_SINGLE_BOND = BondDefiner('cc_single_bond', 'CC', 0, 1)
C_TERMINAL_OH = BondDefiner('c_terminal_oh', 'NCC(O)=O', 2, 3)
SH_BOND = BondDefiner('recent_elongation', 'SC(C)=O', 0, 1)
SC_BOND = BondDefiner('recent_elongation', 'SC(C)=O', 0, 1)
CO_BOND = BondDefiner('recent_elongation', 'CO', 0, 1)
N_AMINO = GroupDefiner('N_amino', 'CN', 1) 
C_CARBOXYL = GroupDefiner('c_carboxyl', 'C(O)=O', 0)
O_OH = GroupDefiner('O_oh', 'CO', 1)
C_CH = GroupDefiner("C_Ch", "CC", 0)
N_NH = GroupDefiner('N_nh', 'CN', 1)
KETO_GROUP = BondDefiner('keto_group', 'C=O', 0, 1)
KETO_GROUP_C = GroupDefiner("keto_group_c", "C=O", 0)
O_BETAPROPRIOLACTONE = GroupDefiner('o_betapropriolactone', 'SC(CCO)=O', 4)
O_BETAPROPRIOLACTONE_O = GroupDefiner('o_betapropriolactone', 'OC(CCO)=O', 4)
O_BETAPROPRIOLACTONE_KETO_C = GroupDefiner('o_betapropriolactone_keto', 'C(CCO)=O', 1)
O_BETAPROPRIOLACTONE_KETO_OH = GroupDefiner('o_betapropriolactone_keto', 'C(CCO)=O', 4)
O_BETAPROPRIOLACTONE_TERMINAL_O = GroupDefiner('o_betapropriolactone', 'OC(CCO)=O', 0)
O_BETAPROPRIOLACTONE_TERMINAL_S = GroupDefiner('o_betapropriolactone', 'SC(CCO)=O', 0)
ASPARTIC_ACID = GroupDefiner('aspartic_acid', 'NC(CC(O)=O)C=O', 4)
GLUTAMIC_ACID = GroupDefiner('glutamic_acid', 'NC(CCC(O)=O)C=O', 5)
ARGININE_SECONDARY_N = GroupDefiner(
    "arginine_secondary_n", "NC(CCCNC(N)=N)C(O)=O", 5)
CYSTEINE = GroupDefiner('cysteine', 'NC(CS)C=O', 3)
SERINE = GroupDefiner('serine', 'NC(CO)C=O', 3)
THREONINE = GroupDefiner('threonine', 'CC(O)C(N)C=O', 2)
REDUCED_SERINE = GroupDefiner('reduced_serine', 'NC(=C)C=O', 2)
REDUCED_THREONINE = GroupDefiner('reduced_threonine', 'C=CC(N)C=O', 0)
PYROPHOSPHATE_BOND = BondDefiner(
    "pyrophosphate_bond", "COP", 0, 1)
ATTACHED_SERINE_OX = GroupDefiner('attached_serine', 'O=CNC(CO)C(S)=O', 0)
ATTACHED_CYSTEINE_OX = GroupDefiner('attached_cysteine', 'O=CNC(CS)C(S)=O', 0)
ATTACHED_SERINE_O = GroupDefiner('attached_serine', 'O=CNC(CO)C(S)=O', 5)
ATTACHED_CYSTEINE_S = GroupDefiner('attached_cysteine', 'O=CNC(CS)C(S)=O', 5)

