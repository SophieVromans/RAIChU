def parse_smiles(paras_file):
    name_to_smiles = {}
    with open(paras_file, 'r') as smiles_file:
        for line in smiles_file:
            line = line.strip()
            if line:
                name, smiles = line.split()
                name_to_smiles[name] = smiles

    return name_to_smiles


