from typing import List, Union
import os

from raichu.cluster import Cluster
from raichu.domain.domain import TailoringDomain, CarrierDomain, SynthesisDomain, RecognitionDomain, \
    TerminationDomain, UnknownDomain, Domain
from raichu.module import PKSModuleSubtype, NRPSModule, LinearPKSModule, IterativePKSModule, TransATPKSModule,\
    ModuleType
from raichu.domain.domain_types import TailoringDomainType, TerminationDomainType, CarrierDomainType, \
    SynthesisDomainType, RecognitionDomainType
from dataclasses import dataclass


DOMAIN_TO_SUPERTYPE = {}
for domain_name in TailoringDomainType.__members__:
    DOMAIN_TO_SUPERTYPE[domain_name] = TailoringDomain
for domain_name in CarrierDomainType.__members__:
    DOMAIN_TO_SUPERTYPE[domain_name] = CarrierDomain
for domain_name in SynthesisDomainType.__members__:
    DOMAIN_TO_SUPERTYPE[domain_name] = SynthesisDomain
for domain_name in RecognitionDomainType.__members__:
    DOMAIN_TO_SUPERTYPE[domain_name] = RecognitionDomain
for domain_name in TerminationDomainType.__members__:
    DOMAIN_TO_SUPERTYPE[domain_name] = TerminationDomain
DOMAIN_TO_SUPERTYPE["UNKNOWN"] = UnknownDomain


@dataclass
class DomainRepresentation:
    gene_name: Union[str, None]
    type: str
    subtype: Union[str, None]
    name: Union[str, None]
    active: bool
    used: bool


@dataclass
class ModuleRepresentation:
    type: str
    subtype: Union[str, None]
    substrate: str
    domains: List[DomainRepresentation]


@dataclass
class ClusterRepresentation:
    modules: List[ModuleRepresentation]


def make_domain(domain_repr: DomainRepresentation, substrate: str, strict: bool = True) -> Domain:
    domain_class = DOMAIN_TO_SUPERTYPE.get(domain_repr.type)
    if domain_class:
        if domain_class == RecognitionDomain:
            domain = domain_class(domain_repr.type, substrate, domain_subtype=domain_repr.subtype,
                                  active=domain_repr.active,
                                  used=domain_repr.used)
        elif domain_class == UnknownDomain:
            domain = UnknownDomain(domain_repr.name)
        else:
            domain = domain_class(domain_repr.type, domain_subtype=domain_repr.subtype, active=domain_repr.active,
                                  used=domain_repr.used)
    elif strict:
        raise ValueError(f"Unrecognised domain type: {domain_repr.type}")
    else:
        domain = UnknownDomain(domain_repr.name)

    return domain


def build_cluster(cluster_repr: ClusterRepresentation, strict: bool = True) -> Cluster:

    modules = []
    for i, module_repr in enumerate(cluster_repr.modules):
        if i == 0:
            starter = True
        else:
            starter = False

        if i == len(cluster_repr.modules) - 1:
            terminator = True
        else:
            terminator = False

        domains = []
        for domain_repr in module_repr.domains:
            domain = make_domain(domain_repr, module_repr.substrate,
                                 strict=strict)
            if domain_repr.gene_name is not None:
                domain.set_gene(domain_repr.gene_name)
            domains.append(domain)

        module_type = ModuleType.from_string(module_repr.type)
        if module_type.name == 'PKS':
            if module_repr.subtype is not None:
                module_subtype = PKSModuleSubtype.from_string(module_repr.subtype)

                if module_subtype.name == 'PKS_CIS':
                    module = LinearPKSModule(i, domains, starter=starter, terminator=terminator)
                elif module_subtype.name == 'PKS_TRANS':
                    module = TransATPKSModule(i, domains, starter=starter, terminator=terminator)
                elif module_subtype.name == 'PKS_ITER':
                    module = IterativePKSModule(i, domains, starter=starter, terminator=terminator)
                else:
                    raise ValueError(f"Unrecognised PKS module subtype: {module_subtype}.")
            else:
                raise ValueError("PKS module subtype must be specified.")

        elif module_type.name == 'NRPS':
            if module_repr.subtype is None:
                module = NRPSModule(i, domains, starter=starter, terminator=terminator)
            else:
                raise ValueError("NRPS module subtypes are currently not supported. Please pass None.")
        else:
            raise ValueError(f"Unrecognised module type: {module_repr.type}")

        modules.append(module)

    cluster = Cluster(modules)

    # return cluster_with_processed_trans_at_pks

    return cluster


def get_spaghettis(cluster_repr: ClusterRepresentation) -> List[str]:

    cluster = build_cluster(cluster_repr)
    # for module in cluster.modules:
    #     print(module.domains)
    cluster.compute_structures(compute_cyclic_products=False)
    cluster.draw_cluster()
    spaghettis = cluster.draw_spaghettis()

    return spaghettis


if __name__ == "__main__":
    cluster_repr = ClusterRepresentation([ModuleRepresentation("PKS", "PKS_CIS", "ACETYL_COA",
                                                               [DomainRepresentation("gene 1", 'AT', None, None, True,
                                                                                     True),
                                                                DomainRepresentation("gene 1", 'KR', 'A1', None, True,
                                                                                     True),
                                                                DomainRepresentation("gene 1", 'DH', None, None, True,
                                                                                     True),
                                                                DomainRepresentation("gene 1", 'ACP', None, None, True,
                                                                                     True)
                                                                ]),
                                          ModuleRepresentation("PKS", "PKS_TRANS", "METHYLMALONYL_COA",
                                                               [DomainRepresentation("gene 1", 'KS',
                                                                                     "TRANS_AT_PKS_BETA_OH", None, True,

                                                                                     True),
                                                                DomainRepresentation("gene 1", 'AT', None, None, True,
                                                                                     True),
                                                                DomainRepresentation("gene 1", 'KR', "A1", None, True,
                                                                                     True),
                                                                DomainRepresentation("gene 1", 'DH', None, None, True,
                                                                                     True),
                                                                DomainRepresentation("gene 1", 'ER', None, None, True,
                                                                                     True),
                                                                DomainRepresentation("gene 1", 'ACP', None, None, True,
                                                                                     True)
                                                                ]),
                                          ModuleRepresentation("PKS", "PKS_TRANS", "METHYLMALONYL_COA",
                                                               [DomainRepresentation("gene 1", 'KS',
                                                                                     "TRANS_AT_PKS_BETA_OH", None, True,
                                                                                     True),
                                                                DomainRepresentation("gene 1", 'AT', None, None, True,
                                                                                     True),
                                                                DomainRepresentation("gene 1", 'KR', "A1", None, True,
                                                                                     True),
                                                                DomainRepresentation("gene 1", 'DH', None, None, True,
                                                                                     True),
                                                                DomainRepresentation("gene 1", 'ER', None, None, True,
                                                                                     True),
                                                                DomainRepresentation("gene 1", 'ACP', None, None, True,
                                                                                     True)
                                                                ]),
                                          ModuleRepresentation("PKS", "PKS_TRANS", "METHYLMALONYL_COA",
                                                               [DomainRepresentation("gene 1", 'KS', "TRANS_AT_PKS_BETA_OH", None, True,
                                                                                     True),
                                                                DomainRepresentation("gene 1", 'AT', None, None, True,
                                                                                     True),
                                                                DomainRepresentation("gene 1", 'KR', "B1", None, True,
                                                                                     True),
                                                                DomainRepresentation("gene 1", 'ER', "S", None, True,
                                                                                     True),
                                                                DomainRepresentation("gene 1", 'DH', None, None, True,
                                                                                     True),
                                                                DomainRepresentation("gene 1", 'UNKNOWN', None, "MyDOM",
                                                                                     True, False),
                                                                DomainRepresentation("gene 1", 'ACP', None, None, True,
                                                                                     True),
                                                                DomainRepresentation("gene 1", 'ACP', None, None, True,
                                                                                     False)
                                                                ]),
                                          ModuleRepresentation("NRPS", None, "tyrosine",
                                                               [DomainRepresentation("gene 1", 'C', None, None, True,
                                                                                     True),
                                                                DomainRepresentation("gene 1", 'nMT', None, None, True,
                                                                                     True),
                                                                DomainRepresentation("gene 1", 'A', None, None, True,
                                                                                     True),
                                                                DomainRepresentation("gene 1", 'PCP', None, None, True,
                                                                                     True),
                                                                DomainRepresentation("gene 1", 'E', None, None, True,
                                                                                     True),
                                                                DomainRepresentation("gene 1", 'TE', None, None, True,
                                                                                     True)
                                                                ]),
                                          ])

    spaghettis = get_spaghettis(cluster_repr)
    out_folder = "starterpks"
    if not os.path.exists(out_folder):
        os.mkdir(out_folder)
    for i, spaghetti in enumerate(spaghettis):
        with open(os.path.join(out_folder, f"{i}.svg"), 'w') as out:
            out.write(spaghetti)
