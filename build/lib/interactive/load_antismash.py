import os

from Bio import SeqIO
from interactive.gene import Gene

from interactive.substrate import Substrate
from interactive.parsers import parse_smiles
import interactive.flatfiles

FLATFILES = os.path.dirname(interactive.flatfiles.__file__)
PARAS_SMILES = os.path.join(FLATFILES, "PARAS_smiles.txt")


NRPS_NAME_TO_SMILES = parse_smiles(PARAS_SMILES)

PKS_NAME_TO_SMILES = {'malonylcoa': "CC(=O)S",
                      'methylmalonylcoa': "CCC(=O)S",
                      'methoxymalonylcoa': "COCC(=O)S",
                      'ethylmalonylcoa': "CCCC(=O)S",
                      'wildcard': "[*]CC(=O)S"}

NRPS_PKS_RULES = {'biosynthetic-additional (rule-based-clusters) AMP-binding',
                  'biosynthetic-additional (rule-based-clusters) PP-binding',
                  'biosynthetic (rule-based-clusters) NRPS: Condensation',
                  'biosynthetic (rule-based-clusters) T1PKS: PKS_AT',
                  'biosynthetic (rule-based-clusters) T1PKS: mod_KS'}

AS_TO_NRPS = {"ala": 'Alanine',
              "arg": "Arginine",
              "asn": "Asparagine",
              "asp": "AsparticAcid",
              "cys": "Cysteine",
              "gln": "Glutamine",
              "glu": "GlutamicAcid",
              "gly": "Glycine",
              "his": "Histidine",
              "ile": "Isoleucine",
              "leu": "Leucine",
              "lys": "Lysine",
              "met": "Methionine",
              "phe": "Phenylalanine",
              "pro": "Proline",
              "ser": "Serine",
              "thr": "Threonine",
              "trp": "Tryptophan",
              "tyr": "Tyrosine",
              "val": "Valine",
              "3-me-glu": "4-methylglutamicacid",
              "4ppro": "nrp",
              "aad": "2-Aminoadipicacid",
              'abu': "2-Aminobutyricacid",
              'aeo': '2-Amino-9,10-epoxy-8-oxodecanoidacid',
              'ala-b': 'beta-Alanine',
              'ala-d': 'D-Alanine',
              'allo-thr': "allo-Threonine",
              'b-ala': 'beta-Alanine',
              'beta-ala': 'beta-Alanine',
              'bmt': "4-butenyl-4-methylthreonine",
              'cap': "capreomycidine",
              'bht': 'nrp',
              'dab': "2,4-diaminobutyricacid",
              'dhb': "2,3-dihydroxybenzoicacid",
              'dhpg': "3,5-dihydroxyphenylglycine",
              'dht': "Dehydrobutyrine",
              'dpg': "3,5-dihydroxyphenylglycine",
              'hiv': "2-hydroxyisovalerate",
              'hiv-d': "D-2-hydroxyisovalerate",
              'hmp-d': 'nrp',
              'horn': 'nrp',
              'hpg': "HydroxyPhenylGlycine",
              'hyv': "4-Hydroxyvaline",
              'hyv-d': 'nrp',
              'iva': "isovaline",
              'lys-b': "beta-lysine",
              'orn': "Ornithine",
              'phg': "phenylglycine",
              'pip': "Pipecolicacid",
              'sal': "salicylicacid",
              'tcl': 'nrp',
              'vol': "valinol",
              'LDAP': 'nrp',
              'meval': "tert-Leu",
              'alle': "allo-Isoleucine",
              'alaninol': "Alaninol",
              'N-(1,1-dimethyl-1-allyl)Trp': 'nrp',
              'd-lyserg': "D-lysergicacid",
              'ser-thr': 'nrp',
              'mephe': 'nrp',
              'haorn': 'nrp',
              'hasn': 'nrp',
              'hforn': 'nrp',
              's-nmethoxy-trp': 'nrp',
              'alpha-hydroxy-isocaproic-acid': 'nrp',
              'MeHOval': 'nrp',
              '2-oxo-isovaleric-acid': "alpha-ketoisovalericacid",
              'aoda': 'nrp',
              'X': 'nrp'}

AS_TO_PKS = {'mal': 'malonylcoa',
             'mmal': 'methylmalonylcoa',
             'mxmal': 'methoxymalonylcoa',
             'emal': 'ethylmalonylcoa',
             }


def sort_genes(genes):
    genes.sort(key=lambda gene: gene.dna_coords[0])
    gene_groups = []
    gene_group = []
    previous_direction = None
    for i, gene in enumerate(genes):
        current_direction = gene.strand
        if current_direction != previous_direction:
            if gene_group:
                if previous_direction == -1:
                    gene_group.reverse()
                gene_groups += gene_group

            gene_group = [gene]
            previous_direction = current_direction
        else:
            gene_group.append(gene)
            if i == len(genes) - 1:
                if current_direction == -1:
                    gene_group.reverse()
                gene_groups += gene_group

    return gene_groups


def parse_antismash_modules(antismash_gbk, screen):
    gene_name_to_genes = {}
    module_nr = -1
    gene_name = None

    for record in SeqIO.parse(antismash_gbk, "genbank"):
        for feature in record.features:
            if feature.type == 'CDS':
                add_gene = False
                module_nr = -1
                qualifiers = feature.qualifiers

                if 'gene' in qualifiers:
                    gene_name = qualifiers['gene'][0]
                elif 'protein' in qualifiers:
                    gene_name = qualifiers['protein'][0]
                elif 'locus_tag' in qualifiers:
                    gene_name = qualifiers['locus_tag'][0]
                else:
                    gene_name = f'gene {len(gene_name_to_genes) + 1}'

                assert gene_name

                if 'gene_functions' in qualifiers:

                    for gene_function in qualifiers['gene_functions']:
                        if gene_function in NRPS_PKS_RULES:
                            add_gene = True

                    if add_gene:
                        gene = Gene(screen, len(gene_name_to_genes))
                        gene.name = gene_name
                        gene.dna_coords = (feature.location.start, feature.location.end)
                        gene.strand = feature.location.strand
                        gene_name_to_genes[gene_name] = gene

            elif feature.type == 'aSModule':

                module_type = feature.qualifiers['type'][0].upper()
                if gene_name in gene_name_to_genes:
                    gene = gene_name_to_genes[gene_name]
                    if module_type == 'NRPS' or module_type == 'PKS':
                        module_nr += 1
                        gene.add_module(module_nr, module_type)
                        for domain in gene.modules[module_nr].domains[:]:
                            gene.modules[module_nr].remove_domain(domain)

                        gene.modules[module_nr].dna_coords = (feature.location.start, feature.location.end)

                        for domain in feature.qualifiers['domains']:
                            if 'PKS_KS' in domain:
                                gene.modules[module_nr].add_domain('KS', domain_label=domain)
                            elif 'PKS_AT' in domain:
                                gene.modules[module_nr].add_domain('AT', domain_label=domain)
                            elif 'PKS_PP' in domain or 'ACP' in domain:
                                if gene.modules[module_nr].type == 'NRPS':
                                    gene.modules[module_nr].add_domain('PCP', domain_label=domain)
                                else:
                                    gene.modules[module_nr].add_domain('ACP', domain_label=domain)
                            elif 'PKS_KR' in domain:
                                gene.modules[module_nr].add_domain('KR', domain_label=domain)
                            elif 'PKS_DH' in domain:
                                gene.modules[module_nr].add_domain('DH', domain_label=domain)
                            elif 'PKS_ER' in domain:
                                gene.modules[module_nr].add_domain('ER', domain_label=domain)
                            elif 'Thioesterase' in domain:
                                gene.modules[module_nr].add_domain('TE', domain_label=domain)
                            elif 'AMP-binding' in domain:
                                gene.modules[module_nr].add_domain('A', domain_label=domain)
                            elif 'Epimerization' in domain:
                                gene.modules[module_nr].add_domain('E', domain_label=domain)
                            elif 'Condensation' in domain or 'Cglyc' in domain:
                                gene.modules[module_nr].add_domain('C', domain_label=domain)
                            elif 'nMT' in domain:
                                gene.modules[module_nr].add_domain('nMT', domain_label=domain)
                            elif 'PCP' in domain:
                                if gene.modules[module_nr].type == 'NRPS':
                                    gene.modules[module_nr].add_domain('PCP', domain_label=domain)
                                else:
                                    gene.modules[module_nr].add_domain('ACP', domain_label=domain)

    return gene_name_to_genes


def parse_antismash_domains(antismash_gbk, gene_name_to_genes):

    gene_name = None
    genes = []

    for record in SeqIO.parse(antismash_gbk, "genbank"):
        for feature in record.features:
            if feature.type == 'CDS':

                qualifiers = feature.qualifiers

                if 'gene' in qualifiers:

                    gene_name = qualifiers['gene'][0]

                elif 'protein' in qualifiers:
                    gene_name = qualifiers['protein'][0]
                elif 'locus_tag' in qualifiers:
                    gene_name = qualifiers['locus_tag'][0]
                else:
                    gene_name = f'gene {len(genes) + 1}'

                assert gene_name
                genes.append(gene_name)

            elif feature.type == 'aSDomain':
                if gene_name in gene_name_to_genes:
                    gene = gene_name_to_genes[gene_name]
                    domains = gene.get_domains()
                    domain_id = feature.qualifiers['domain_id'][0]

                    for domain in domains:
                        if domain.domain_label == domain_id:
                            if 'specificity' in feature.qualifiers:
                                for spec in feature.qualifiers['specificity']:
                                    if 'consensus:' in spec:
                                        specificity = spec.split('consensus:')[1].strip()
                                        if specificity in AS_TO_NRPS:
                                            substrate_name = AS_TO_NRPS[specificity]
                                            substrate_smiles = NRPS_NAME_TO_SMILES[substrate_name]
                                            domain.substrate = Substrate(substrate_name, substrate_smiles, 'NRPS')
                                        elif specificity in AS_TO_PKS:
                                            substrate_name = AS_TO_PKS[specificity]
                                            substrate_smiles = PKS_NAME_TO_SMILES[substrate_name]
                                            domain.substrate = Substrate(substrate_name, substrate_smiles, 'PKS')

                                    elif 'KR stereochemistry:' in spec:
                                        specificity = spec.split('KR stereochemistry:')[1].strip()
                                        if specificity != '(unknown)':
                                            domain.subtype = specificity


def genes_from_antismash(antismash_gbk, screen):
    gene_name_to_genes = parse_antismash_modules(antismash_gbk, screen)
    parse_antismash_domains(antismash_gbk, gene_name_to_genes)

    genes = []

    for gene in gene_name_to_genes.values():
        genes.append(gene)

    genes = sort_genes(genes)
    for i, gene in enumerate(genes):
        gene.gene_number = i
        gene.sort_modules()

    return genes

