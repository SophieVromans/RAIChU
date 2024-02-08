import random
import os
from sys import argv
from raichu.run_raichu import draw_cluster
from raichu.representations import ClusterRepresentation, ModuleRepresentation, DomainRepresentation
from paras.features import _METADATA
from raichu.substrate import PksStarterSubstrate, PksElongationSubstrate
from raichu.domain.domain_types import KRDomainSubtype, ERDomainSubtype
import traceback


PROTEINOGENIC_AA = ['alanine',
                    'cysteine',
                    'aspartic acid',
                    'glutamic acid',
                    'phenylalanine',
                    'glycine',
                    'histidine',
                    'isoleucine',
                    'lysine',
                    'leucine',
                    'methionine',
                    'asparagine',
                    'proline',
                    'glutamine',
                    'arginine',
                    'serine',
                    'threonine',
                    'valine',
                    'tryptophan',
                    'tyrosine']

AA_STARTER_CHOICES = []
AA_MODULE_CHOICES = []

for substrate_name, metadata in _METADATA.items():
    if metadata.type in ['amino_acid', 'beta_amino_acid']:
        AA_MODULE_CHOICES.append(substrate_name)
    if metadata.type in ['amino_acid', 'beta_amino_acid', 'acid']:
        AA_STARTER_CHOICES.append(substrate_name)

PKS_STARTER_SUBSTRATE_CHOICES = []
PKS_ELONGATION_SUBSTRATE_CHOICES = []

for pks_substrate in PksStarterSubstrate:
    PKS_STARTER_SUBSTRATE_CHOICES.append(pks_substrate.name)
for pks_substrate in PksElongationSubstrate:
    PKS_ELONGATION_SUBSTRATE_CHOICES.append(pks_substrate.name)

KR_DOMAIN_SUBTYPE_CHOICES = []
for subtype in KRDomainSubtype:
    KR_DOMAIN_SUBTYPE_CHOICES.append(subtype.name)

ER_DOMAIN_SUBTYPE_CHOICES = []
for subtype in ERDomainSubtype:
    ER_DOMAIN_SUBTYPE_CHOICES.append(subtype.name)

TERMINATION_DOMAIN_CHOICES = ["TD"] + 10 * ["TE"]


def choose(nr_true: int, nr_false: int):
    return random.choice([True] * nr_true + [False] * nr_false)


def generate_domain(gene_nr, domain_type, domain_subtype=None, name=None, active=True):

    go_to_next_gene = choose(1, 20)
    if go_to_next_gene:
        gene_nr += 1

    gene_name = f"gene {gene_nr}"

    domain = DomainRepresentation(gene_name, domain_type, domain_subtype, name=name, active=active)

    return domain, gene_nr


def generate_cis_pks_module(gene_nr: int, module_nr: int, terminal_module: bool = False) -> \
        tuple[ModuleRepresentation, int]:

    ks_domain = None
    kr_domain = None
    dh_domain = None
    er_domain = None
    te_domain = None

    unknown_domains = []

    if module_nr == 0:
        substrate = random.choice(PKS_STARTER_SUBSTRATE_CHOICES)
        has_ks = choose(1, 10)
        if has_ks:
            ks_domain, gene_nr = generate_domain(gene_nr, "KS")

    else:
        substrate = random.choice(PKS_ELONGATION_SUBSTRATE_CHOICES)
        ks_domain, gene_nr = generate_domain(gene_nr, "KS")

    has_unknown_domains = choose(1, 10)
    nr_unknowns = 0
    if has_unknown_domains:
        nr_unknowns = random.choice([1] * 15 + [2] * 4 + [3] * 1)

    has_kr = choose(1, 1)
    has_dh = choose(1, 9)
    has_er = choose(1, 19)

    if has_kr:

        has_dh = random.choice([True, False])
        if has_dh:
            has_er = random.choice([True, False])

    at_domain, gene_nr = generate_domain(gene_nr, "AT")

    if has_dh:
        dh_active = choose(19, 1)
        dh_domain, gene_nr = generate_domain(gene_nr, "DH", active=dh_active)

    if has_er:
        er_active = choose(19, 1)
        er_domain_subtype = random.choice(ER_DOMAIN_SUBTYPE_CHOICES)
        er_domain, gene_nr = generate_domain(gene_nr, 'ER', er_domain_subtype, active=er_active)

    if has_kr:
        kr_domain_subtype = random.choice(KR_DOMAIN_SUBTYPE_CHOICES)
        kr_active = choose(19, 1)
        kr_domain, gene_nr = generate_domain(gene_nr, "KR", kr_domain_subtype, active=kr_active)

    for i in range(nr_unknowns):
        domain, gene_nr = generate_domain(gene_nr, "UNKNOWN", name=f"unknown_{i + 1}")
        unknown_domains.append(domain)

    acp_domain, gene_nr = generate_domain(gene_nr, "ACP")

    if terminal_module:
        domain_type = random.choice(TERMINATION_DOMAIN_CHOICES)
        te_domain, gene_nr = generate_domain(gene_nr, domain_type)

    putative_domains = [ks_domain, at_domain, dh_domain, er_domain, kr_domain] + unknown_domains + \
                       [acp_domain, te_domain]

    domains: list[DomainRepresentation] = []
    for domain in putative_domains:
        if domain:
            domains.append(domain)

    module = ModuleRepresentation("PKS", "PKS_CIS", substrate, domains)

    return module, gene_nr


def generate_nrps_module(gene_nr, module_nr, terminal_module=False):

    has_proteinogenic_substrate = choose(3, 1)

    c_domain = None
    e_domain = None
    nmt_domain = None
    cyc_domain = None
    ox_domain = None
    te_domain = None

    unknown_domains = []

    if module_nr == 0:
        if has_proteinogenic_substrate:
            substrate = random.choice(PROTEINOGENIC_AA)
        else:
            substrate = random.choice(AA_STARTER_CHOICES)

    else:
        c_domain, gene_nr = generate_domain(gene_nr, "C")
        if has_proteinogenic_substrate:
            substrate = random.choice(PROTEINOGENIC_AA)
        else:
            substrate = random.choice(AA_MODULE_CHOICES)

    has_unknown_domains = choose(1, 10)
    nr_unknowns = 0
    if has_unknown_domains:
        nr_unknowns = random.choice([1] * 15 + [2] * 4 + [3] * 1)

    has_e = choose(1, 1)
    has_nmt = choose(1, 10)
    has_cyc = choose(1, 30)
    has_ox = choose(1, 30)

    if substrate in ['threonine', 'serine', 'cysteine']:
        has_cyc = choose(1, 6)

    if has_cyc:
        has_ox = choose(1, 1)

    if has_cyc:
        cyc_domain, gene_nr = generate_domain(gene_nr, "CYC")

    a_domain, gene_nr = generate_domain(gene_nr, "A")

    if has_ox:
        ox_domain, gene_nr = generate_domain(gene_nr, "OX")

    if has_nmt:
        nmt_domain, gene_nr = generate_domain(gene_nr, "nMT")

    for i in range(nr_unknowns):
        domain, gene_nr = generate_domain(gene_nr, "UNKNOWN", name=f"unknown_{i + 1}")
        unknown_domains.append(domain)

    pcp_domain, gene_nr = generate_domain(gene_nr, "PCP")
    if has_e:
        e_domain, gene_nr = generate_domain(gene_nr, "E")

    if terminal_module:
        domain_type = random.choice(TERMINATION_DOMAIN_CHOICES)
        te_domain, gene_nr = generate_domain(gene_nr, domain_type)

    putative_domains = [c_domain, cyc_domain, a_domain, ox_domain, nmt_domain] + unknown_domains + \
                       [pcp_domain, e_domain, te_domain]

    domains: list[DomainRepresentation] = []
    for domain in putative_domains:
        if domain:
            domains.append(domain)

    module = ModuleRepresentation("NRPS", None, substrate, domains)

    return module, gene_nr


def generate_trans_pks_module(gene_nr, module_nr, terminal_module=False):
    pass
    # return module, gene_nr


def generate_modular_cluster(nr_modules, drawing_dir, failed_dir, cluster_nr, cis_pks=True, nrps=True):

    gene_nr = 1
    choices = []
    if cis_pks:
        choices.append('cis-pks')
    if nrps:
        choices.append('nrps')

    terminal_module = False

    modules = []

    for i in range(nr_modules):
        if i != 0:
            go_to_next_gene = choose(1, 3)
            if go_to_next_gene:
                gene_nr += 1

        if i == nr_modules - 1:
            terminal_module = True

        module_type = random.choice(choices)
        if module_type == 'cis-pks':
            module, gene_nr = generate_cis_pks_module(gene_nr, i, terminal_module)
        elif module_type == 'nrps':
            module, gene_nr = generate_nrps_module(gene_nr, i, terminal_module)
        elif module_type == 'trans-pks':
            raise NotImplementedError
        else:
            raise ValueError(f"Module type must be 'cis-pks', 'nrps', 'trans-pks'. Got f{module_type}")

        modules.append(module)

    cluster = ClusterRepresentation(modules)
    try:
        drawing_out = os.path.join(drawing_dir, f"cluster_{cluster_nr}.svg")
        draw_cluster(cluster, drawing_out)

    except Exception:
        print(cluster.modules[0].substrate)
        print(traceback.format_exc())
        failed_out = os.path.join(failed_dir, f"cluster_{cluster_nr}")
        cluster.write_cluster(failed_out)


def generate_random_clusters(nr_clusters, out_folder, nrps=True, cis_pks=True):
    if not os.path.exists(out_folder):
        os.mkdir(out_folder)

    drawing_dir = os.path.join(out_folder, 'drawings')
    failed_dir = os.path.join(out_folder, 'failed')

    if not os.path.exists(drawing_dir):
        os.mkdir(drawing_dir)
    if not os.path.exists(failed_dir):
        os.mkdir(failed_dir)

    for i in range(nr_clusters):
        nr_modules = random.randint(2, 13)
        print(f"Drawing cluster {i + 1}")
        generate_modular_cluster(nr_modules, drawing_dir, failed_dir, i + 1, nrps=nrps, cis_pks=cis_pks)


if __name__ == "__main__":
    out_folder = argv[1]
    if not os.path.exists(out_folder):
        os.mkdir(out_folder)

    nrps_folder = os.path.join(out_folder, "nrps")
    pks_folder = os.path.join(out_folder, "pks")
    hybrid_folder = os.path.join(out_folder, "hybrid")

    generate_random_clusters(500, nrps_folder, nrps=True, cis_pks=False)
    generate_random_clusters(500, pks_folder, nrps=False, cis_pks=True)
    generate_random_clusters(500, hybrid_folder, nrps=True, cis_pks=True)




