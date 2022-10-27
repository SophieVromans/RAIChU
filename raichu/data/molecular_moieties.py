import os

from pikachu.general import read_smiles
from pikachu.reactions.functional_groups import BondDefiner, GroupDefiner
from raichu.reactions.general import initialise_atom_attributes

import raichu.data

FLATFILES = os.path.dirname(raichu.data.__file__)
AA_SMILES = os.path.join(FLATFILES, "PARAS_smiles.txt")


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


def parse_smiles():
    """
    Returns a dict as {name_aa : SMILES_aa} from the amino acids in the
    PARAS_smiles.txt file
    """
    # Parse list SMILES amino acids attached to PCP to dict name -> Structure
    name_to_smiles = {}
    with open(AA_SMILES, 'r', encoding='utf8') as paras_smiles:
        for line in paras_smiles:
            line = line.strip()
            print(line)
            name, smiles = line.split('\t')
            name = name.upper()
            name_to_smiles[name] = smiles

    return name_to_smiles


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


THIOESTERBOND = BondDefiner('thioester_bond', 'SC(C)=O', 0, 1)
THIOESTER_CARBON = GroupDefiner('thioester carbon', 'SC(C)=O', 1)
LEAVING_OH_BOND = BondDefiner('Leaving -OH group bond', 'C(=O)(O)CN', 0, 2)
N_AMINO_ACID = GroupDefiner('Nitrogen atom amino acid', 'NCC(=O)O', 0)
C1_AMINO_ACID = GroupDefiner('C1 atom amino acid', 'NCC(=O)O', 1)
C2_AMINO_ACID = GroupDefiner('C2 atom amino acid', 'NCC(=O)O', 2)
AMINO_ACID_BACKBONE = read_smiles('NCC(=O)O')
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

SC_BOND = BondDefiner('recent_elongation', 'SC(C)=O', 0, 1)
CO_BOND = BondDefiner('recent_elongation', 'CO', 0, 1)
N_AMINO = GroupDefiner('N_amino', 'CN', 1)
O_OH = GroupDefiner('O_oh', 'CO', 1)

O_BETAPROPRIOLACTONE = GroupDefiner('o_betapropriolactone', 'SC(CCO)=O', 4)
O_BETAPROPRIOLACTONE_O = GroupDefiner('o_betapropriolactone', 'OC(CCO)=O', 4)
O_BETAPROPRIOLACTONE_TERMINAL_O = GroupDefiner('o_betapropriolactone', 'OC(CCO)=O', 0)
