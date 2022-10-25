import os

from interactive.gene import Gene
from interactive.substrate import Substrate
from interactive.parsers import parse_smiles
import interactive.flatfiles

FLATFILES = os.path.dirname(interactive.flatfiles.__file__)
PARAS_SMILES = os.path.join(FLATFILES, "PARAS_smiles.txt")
PKS_SMILES = os.path.join(FLATFILES, "at_specificities.txt")

pks_name_to_smiles = {'malonylcoa': "CC(=O)S",
                      'methylmalonylcoa': "CCC(=O)S",
                      'methoxymalonylacp': "COCC(=O)S",
                      'ethylmalonylcoa': "CCCC(=O)S",
                      'wildcard': "[*]CC(=O)S"}

nrps_name_to_smiles = parse_smiles(PARAS_SMILES)

pks_starter = parse_smiles(PKS_SMILES)
pks_starter_reverse = {}

for name, smiles in pks_starter.items():
    pks_starter_reverse[smiles] = name


def load_tabular(tabular_file, screen):
    genes = []
    previous_gene_name = None

    with open(tabular_file, 'r') as tabular:
        tabular.readline()
        for line in tabular:
            line_info = line.split('\t')
            for i, item in enumerate(line_info[:]):
                line_info[i] = item.strip()

            if len(line_info) == 5:
                bgc_type, gene_name, module_id, module_type, specificity = line_info
                tailoring = ""
            else:
                bgc_type, gene_name, module_id, module_type, specificity, tailoring = line_info
            tailoring_domains = tailoring.split('|')
            if not specificity:
                specificity = None

            if previous_gene_name != gene_name:
                previous_gene_name = gene_name
                gene = Gene(screen, len(genes))
                gene.name = gene_name
                genes.append(gene)

            module_type = module_type.split('_')[-1].upper()

            genes[-1].add_module(len(genes[-1].modules), module_type)

            for domain in tailoring_domains:
                if domain:
                    subtype = None
                    if domain.startswith('KR'):
                        domain_and_subtype = domain.split('_')
                        if len(domain_and_subtype) == 2:
                            domain, subtype = domain_and_subtype

                    genes[-1].modules[-1].add_domain(domain)
                    if subtype:
                        for dom in genes[-1].modules[-1].domains:
                            if dom.type == 'KR':
                                dom.subtype = subtype

                if specificity:
                    for dom in genes[-1].modules[-1].domains:
                        if dom.type == 'A' or dom.type == 'AT':
                            substrate_smiles = None
                            custom = False
                            if dom.module.type == 'NRPS':
                                if specificity in nrps_name_to_smiles:
                                    substrate_smiles = nrps_name_to_smiles[specificity]
                                    substrate_name = specificity
                                else:
                                    substrate_smiles = specificity
                                    substrate_name = 'FA'
                                    custom = True
                            elif dom.module.type == 'PKS':
                                if specificity in pks_name_to_smiles:
                                    substrate_smiles = pks_name_to_smiles[specificity]
                                    substrate_name = specificity
                                else:
                                    substrate_smiles = specificity
                                    if substrate_smiles in pks_starter_reverse:
                                        substrate_name = pks_starter_reverse[substrate_smiles]
                                    else:
                                        substrate_name = 'start'
                                        custom = True

                            assert substrate_smiles

                            dom.substrate = Substrate(substrate_name, substrate_smiles, dom.module.type, custom=custom)

    return genes






