from paras.features import get_smiles
from raichu.data.molecular_moieties import make_elongation_monomer
from enum import Enum, unique

_PKS_TO_SMILES = {"WILDCARD": r"SC(C([*])C(O)=O)=O",
                  "MALONYL_COA": r"SC(CC(O)=O)=O",
                  "METHYLMALONYL_COA": r"SC([C@H](C)C(O)=O)=O",
                  "METHOXYMALONYL_ACP": r"SC(C(C(O)=O)OC)=O",
                  "METHYLBUTYRYL_COA_2S": r"SC(=O)[C@@H](C)CC",
                  "METHYLBUTYRYL_COA_2R": r"SC(=O)[C@H](C)CC",
                  "ETHYLMALONYL_COA": r"SC(C(C(O)=O)CC)=O",
                  "PROPIONYL_COA": r"SC(=O)CC",
                  "ACETYL_COA": r"SC(C)=O",
                  "BENZOYL_COA": r"SC(C1=CC=CC=C1)=O",
                  "METHYL_BUTYRYL_COA_3": r"SC(=O)CC(C)C",
                  "METHYL_BUTYRYL_COA_2": r"SC(=O)C(C)CC",
                  "TRANS_CYCLOPENTANE_DICARBOXYL_COA": r"SC([C@H]1[C@@H](CCC1)C(=O)O)=O",
                  "CYCLOHEXANE_CARBOXYL_COA": r"SC(=O)C1CCCCC1",
                  "HYDROXY_MALONYL_COA_2": r"SC(=O)C(C(=O)O)O",
                  "HYDROXY_MALONYL_COA_2R": r"SC(=O)[C@@H](C(=O)O)O",
                  "HYDROXY_MALONYL_COA_2S": r"SC(=O)[C@H](C(=O)O)O",
                  "CHLOROETHYL_MALONYL_COA": r"SC(C(C(=O)O)CC[Cl])=O",
                  "ISOBUTYRYL_COA": r"SC(=O)C(C)C",
                  "GLYCINE": "C(C(=O)O)N",
                  "HYDROXY_PROPENOYL_COA_3_23E": r"[H]\C(O)=C/C(S)=O",
                  "HYDROXY_BUTENOYL_COA_3_23E": r"C\C(O)=C/C(S)=O",
                  "DIHYDROXY_BUTANOLYL_COA_2R3": r"CC(O)[C@@H](O)C(O)S",
                  "TRIHYDROXY_PROPANOLYL_COA_233": r"OC(O)C(O)C(O)S",
                  "O_METHYLACETYL_COA": r"COC(S)=O",
                  "HYDROXY_PROPENOYL_COA_3_23Z": r"[H]\C(O)=C\C(S)=O",
                  "OXOMALONYL_COA_2": r"OC(=O)C(=O)C(S)=O",
                  "METHYL_HYDROXY_PROPENOYL_COA_2_3_23Z": r"[H]\C(O)=C(/C)C(S)=O",
                  "DIHYDROXY_BUTANOLYL_COA_23": r"CC(O)C(O)C(O)S",
                  "DIHYDROXY_BUTANOLYL_COA_2S3S": r"C[C@H](O)[C@H](O)C(O)S",
                  "HEPTATRIENOYL_COA": r"SC(=O)C1=CC=CC=C3",
                  "HYDROXYPROPIONYL_COA_2R": r"C[C@@H](O)C(S)=O",
                  "DIHYDROXY_PROPANOLYL_COA_33": r"OC(O)CC(O)S",
                  "LACTYL_COA": r"C[C@@H](O)C(S)=O",
                  "PHENYLACETYLCOA": r"SC(=O)CC1=CC=CC=C1",
                  "METHOXYFORMYL_COA": r"COC(S)=O",
                  "GLYCERYL_COA": r"OCC(O)C(S)=O"
                  }


@unique
class PksStarterSubstrate(Enum):
    PROPIONYL_COA = 1
    ACETYL_COA = 2
    BENZOYL_COA = 3
    METHYL_BUTYRYL_COA_3 = 4  # 3-methylbutyryl-CoA
    METHYL_BUTYRYL_COA_2 = 5  # 2-methylbutyryl-CoA
    TRANS_CYCLOPENTANE_DICARBOXYL_COA = 6
    CYCLOHEXANE_CARBOXYL_COA = 7
    HYDROXY_MALONYL_COA_2 = 8  # 2-hydroxymalonyl-CoA
    HYDROXY_MALONYL_COA_2S = 9  # 2S-hydroxymalonyl-CoA
    HYDROXY_MALONYL_COA_2R = 10  # 2R-hydroxymalonyl-CoA
    CHLOROETHYL_MALONYL_COA = 11
    ISOBUTYRYL_COA = 12
    GLYCINE = 13
    HYDROXY_PROPENOYL_COA_3_23E = 14  # 3-hydroxy-2,3-E-propenoyl-CoA
    HYDROXY_BUTENOYL_COA_3_23E = 15  # 3-hydroxy-2,3-E-butenoyl-CoA
    DIHYDROXY_BUTANOLYL_COA_2R3 = 16  # 2R,3-dihydroxy-butanolyl-CoA
    TRIHYDROXY_PROPANOLYL_COA_233 = 17  # 2,3,3-trihydroxy-propanolyl-CoA
    O_METHYLACETYL_COA = 18  # O-methylacetyl-CoA
    HYDROXY_PROPENOYL_COA_3_23Z = 19  # 3-hydroxy-2,3-Z-propenoyl-CoA
    OXOMALONYL_COA_2 = 20  # 2-oxomalonyl-CoA
    METHYL_HYDROXY_PROPENOYL_COA_2_3_23Z = 21  # 2-methyl-3-hydroxy-2,3-Z-propenoyl-CoA
    DIHYDROXY_BUTANOLYL_COA_23 = 22  # 2,3-dihydroxy-butanolyl-CoA
    DIHYDROXY_BUTANOLYL_COA_2S3S = 23  # 2S,3S-dihydroxy-butanolyl-CoA
    HEPTATRIENOYL_COA = 24  # 2,4,6-heptatrienoyl-CoA
    HYDROXYPROPIONYL_COA_2R = 25  # 2R-hydroxypropionyl-CoA
    DIHYDROXY_PROPANOLYL_COA_33 = 26  # 3,3-dihydroxy-propanolyl-CoA
    LACTYL_COA=27
    PHENYLACETYLCOA=28
    METHOXYFORMYL_COA=29
    GLYCERYL_COA=30

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
            self.elongation_monomer = make_elongation_monomer(self.name)
        elif name in [v.name for v in PksStarterSubstrate]:
            self.elongation_monomer = None
        else:
            raise ValueError(f"PKS substrate {self.name} is not recognised by RAIChU.")
