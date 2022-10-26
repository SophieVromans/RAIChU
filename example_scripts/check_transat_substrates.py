from pikachu.general import read_smiles, svg_from_smiles
from pikachu.fingerprinting.similarity import get_jaccard_index

TRANSATOR_CLADE_TO_STARTER_SUBSTRATE = {'TRANS_AT_PKS_NON_ELONGATING_DB': 'OC(=O)CC(S)=O',
                                        # TODO: Figure out what is going on with this SMILES string
                                        'TRANS_AT_PKS_BETA_OH_KETO': 'OC(C(O)=O)C(S)=O OC(=O)C(=O)C(S)=O',
                                        'TRANS_AT_PKS_NON_ELONGATING_BETA_OH': 'OC(=O)CC(S)=O',
                                        'TRANS_AT_PKS_ALPHA_OH': 'OC(O)CC(O)S',
                                        'TRANS_AT_PKS_EDB': '[H]\\C(O)=C/C(S)=O',
                                        'TRANS_AT_PKS_ARST': 'SC(=O)C1=CC=CC=C3', 'TRANS_AT_PKS_PYR': 'OC(=O)CC(S)=O',
                                        'TRANS_AT_PKS_ALPHAME_EDB': 'C\\C(O)=C/C(S)=O',
                                        'TRANS_AT_PKS_OXI': 'OC(=O)CC(S)=O', 'TRANS_AT_PKS_AA': 'C(C(=O)O)N',
                                        'TRANS_AT_PKS_ALPHAME_BETA_L_OH': 'CC(O)[C@@H](O)C(O)S',
                                        'TRANS_AT_PKS_NON_ELONGATING_BETA_L_OH': 'O[C@H](C(O)=O)C(S)=O',
                                        'TRANS_AT_PKS_RED_SHDB': 'OC(=O)CC(S)=O',
                                        'TRANS_AT_PKS_DB': '[H]\\C(O)=C/C(S)=O',
                                        'TRANS_AT_PKS_KETO': 'OC(=O)C(=O)C(S)=O',
                                        'TRANS_AT_PKS_NON_ELONGATING': 'OC(=O)CC(S)=O',
                                        'TRANS_AT_PKS_BETA_D_OH': 'OC(C(O)=O)C(S)=O',
                                        'TRANS_AT_PKS_BETA_L_OH': 'O[C@H](C(O)=O)C(S)=O',
                                        'TRANS_AT_PKS_BR': 'OC(=O)CC(S)=O',
                                        'TRANS_AT_PKS_ALPHABETA_OH': 'OC(O)C(O)C(O)S',
                                        'TRANS_AT_PKS_UNST': 'OC(=O)CC(S)=O', 'TRANS_AT_PKS_ST': 'CC(S)=O',
                                        'TRANS_AT_PKS_MEOST': 'COC(S)=O', 'TRANS_AT_PKS_ACST': 'CCC(S)=O CC(S)=O',
                                        'TRANS_AT_PKS_ALPHAME': 'CC(O)CC(S)=O CC(O)C(O)C(S)=O CC(O)C(=O)C(S)=O',
                                        'TRANS_AT_PKS_BETA_D_OME': 'COC(C(O)=O)C(S)=O',
                                        'TRANS_AT_PKS_NON_ELONGATING_PYR': 'OC(=O)CC(S)=O',
                                        'TRANS_AT_PKS_BETA_ME': 'OC(=O)C(=C)C(S)=O CC(C(O)=O)C(S)=O',
                                        'TRANS_AT_PKS_BETA_OH_EDB': 'OC(C(O)=O)C(S)=O O\\C=C\\C(S)=O',
                                        'TRANS_AT_PKS_NON_ELONGATING_OXA': 'CC1=NC(C)=CO1 CC1=NC(C)=CS8',
                                        'TRANS_AT_PKS_LACST': 'C[C@@H](O)C(S)=O', 'TRANS_AT_PKS_SHDB': 'OC(=O)CC(S)=O',
                                        'TRANS_AT_PKS_OUT': 'OC(=O)CC(S)=O', 'TRANS_AT_PKS_BETA_OH': 'OC(C(O)=O)C(S)=O',
                                        'TRANS_AT_PKS_ZDB': '[H]\\C(O)=C\\C(S)=O',
                                        'TRANS_AT_PKS_BETA_MEDB': '[H]\\C(O)=C(/C)C(S)=O',
                                        'TRANS_AT_PKS': 'OC(=O)CC(S)=O', 'TRANS_AT_PKS_OXA': 'OC(=O)CC(S)=O',
                                        'TRANS_AT_PKS_ALPHAME_BETA_D_OH': 'C[C@H](O)[C@H](O)C(O)S',
                                        'TRANS_AT_PKS_ALPHAME_BETAOH': 'CC(O)C(O)C(O)S',
                                        'TRANS_AT_PKS_NON_ELONGATING_ALPHAME_EDB': 'C\\C(O)=C/C(S)=O',
                                        'TRANS_AT_PKS_RED': 'OC(=O)CC(S)=O'}

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

unique_trans_smiles = set()
trans_name_to_name = {}

for trans_name, trans_smiles in TRANSATOR_CLADE_TO_STARTER_SUBSTRATE.items():
    match_found = False

    try:
        s2 = read_smiles(trans_smiles)
    except Exception:
        print(f"Invalid SMILES {trans_smiles} for substrate {trans_name}.")
        continue

    for name, smiles in _PKS_TO_SMILES.items():
        try:
            s1 = read_smiles(smiles)
        except Exception:
            print(f"Invalid SMILES {smiles} for substrate {name}.")
            continue

        jaccard_index = get_jaccard_index(s1, s2)
        if jaccard_index == 1:
            if match_found:
                print(f"Warning: match already found for {trans_name}: {trans_name_to_name[trans_name]}")
            trans_name_to_name[trans_name] = name
            print(trans_smiles, smiles)
            match_found = True

    if not match_found:
        unique_trans_smiles.add(trans_smiles)

unique_trans_smiles = '.'.join(unique_trans_smiles)
svg_from_smiles(unique_trans_smiles, 'trans_at.svg')
print(unique_trans_smiles)

print(trans_name_to_name)
print(get_jaccard_index(read_smiles(r"OC(C(O)=O)C(S)=O"), read_smiles(r"SC(=O)C(C(=O)O)O")))




