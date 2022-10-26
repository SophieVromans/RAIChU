from paras.features import get_smiles
from raichu.data.molecular_moieties import NAME_TO_ELONGATION_MONOMER
from enum import Enum, unique

_PKS_TO_SMILES = {"WILDCARD": "SC(C([*])C(O)=O)=O",
                  "MALONYL_COA": "SC(CC(O)=O)=O",
                  "METHYLMALONYL_COA": "SC([C@H](C)C(O)=O)=O",
                  "METHOXYMALONYL_ACP": "SC(C(C(O)=O)OC)=O",
                  "METHYLBUTYRYL_COA_2S": "SC(=O)[C@@H](C)CC",
                  "METHYLBUTYRYL_COA_2R": "SC(=O)[C@H](C)CC",
                  "ETHYLMALONYL_COA": "SC(C(C(O)=O)CC)=O",
                  "PROPIONYL_COA": "SC(=O)CC",
                  "ACETYL_COA": "SC(C)=O",
                  "BENZOYL_COA": "SC(C1=CC=CC=C1)=O",
                  "METHYL_BUTYRYL_COA_3": "SC(=O)CC(C)C",
                  "METHYL_BUTYRYL_COA_2": "SC(=O)C(C)CC",
                  "TRANS_CYCLOPENTANE_DICARBOXYL_COA": "SC([C@H]1[C@@H](CCC1)C(=O)O)=O",
                  "CYCLOHEXANE_CARBOXYL_COA": "SC(=O)C1CCCCC1",
                  "HYDROXY_MALONYL_COA_2": "SC(=O)C(C(=O)O)O",
                  "HYDROXY_MALONYL_COA_2S": "SC(=O)[C@@H](C(=O)O)O",
                  "HYDROXY_MALONYL_COA_2R": "SC(=O)[C@H](C(=O)O)O",
                  "CHLOROETHYL_MALONYL_COA": "SC(C(C(=O)O)CC[Cl])=O",
                  "ISOBUTYRYL_COA": "SC(=O)C(C)C"}


@unique
class PksStarterSubstrate(Enum):
    PROPIONYL_COA = 1
    ACETYL_COA = 2
    BENZOYL_COA = 3
    METHYL_BUTYRYL_COA_3 = 4
    METHYL_BUTYRYL_COA_2 = 5
    TRANS_CYCLOPENTANE_DICARBOXYL_COA = 6
    CYCLOHEXANE_CARBOXYL_COA = 7
    HYDROXY_MALONYL_COA_2 = 8
    HYDROXY_MALONYL_COA_2S = 9
    HYDROXY_MALONYL_COA_2R = 10
    CHLOROETHYL_MALONYL_COA = 11
    ISOBUTYRYL_COA = 12

    @staticmethod
    def from_string(label: str):
        for value in PksElongationSubstrate:
            if str(value) == label:
                return value
        raise ValueError(f"Unknown PKS starter substrate: {label}")

    @staticmethod
    def get_smiles(label: str):
        for value in PksElongationSubstrate:
            if str(value) == label:
                return _PKS_TO_SMILES[value]

        raise ValueError(f"Unknown PKS starter substrate: {label}")

@unique
class PksElongationSubstrate(Enum):
    WILDCARD = 1
    MALONYL_COA = 2
    METHYLMALONYL_COA = 3
    METHOXYMALONYL_ACP = 4
    METHYLBUTYRYL_COA_2S = 5
    METHYLBUTYRYL_COA_2R = 6
    ETHYLMALONYL_COA = 7

    @staticmethod
    def from_string(label: str):
        for value in PksElongationSubstrate:
            if str(value) == label:
                return value
        raise ValueError(f"Unknown PKS elongation substrate: {label}")

    @staticmethod
    def get_smiles(label: str):
        for value in PksElongationSubstrate:
            if str(value) == label:
                return _PKS_TO_SMILES[value]

        raise ValueError(f"Unknown PKS elongation substrate: {label}")


class Substrate:
    def __init__(self, name: str, smiles: str) -> None:
        self.name = name
        self.smiles = smiles


class NRPSSubstrate(Substrate):
    def __init__(self, name: str) -> None:
        smiles = get_smiles(name)
        super().__init__(name, smiles)


class PKSSubstrate(Substrate):
    def __init__(self, name: str) -> None:
        smiles = _PKS_TO_SMILES.get(name, None)
        if smiles is None:
            raise ValueError(f"Cannot fetch SMILES string for PKS substrate {name}.")
        super().__init__(name, smiles)
        if name in [v.name for v in PksElongationSubstrate]:
            self.elongation_monomer = NAME_TO_ELONGATION_MONOMER[self.name]
        elif name in [v.name for v in PksStarterSubstrate]:
            self.elongation_monomer = None
        else:
            raise ValueError(f"PKS substrate {self.name} is not recognised by RAIChU.")
