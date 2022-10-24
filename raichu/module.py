from typing import List, Union

from pikachu.general import read_smiles
from pikachu.chem.structure import Structure

from raichu.domain.domain import Domain, TailoringDomain, RecognitionDomain, \
    SynthesisDomain, CarrierDomain, TerminationDomain
from raichu.central_chain_detection.label_central_chain import label_pk_central_chain, label_nrp_central_chain
from raichu.attach_to_domain import attach_to_domain_pk, attach_to_domain_nrp
from enum import Enum, unique


@unique
class ModuleType(Enum):
    NRPS = 1
    PKS = 2

    @staticmethod
    def from_string(label: str) -> "ModuleType":
        for value in ModuleType:
            if str(value.name) == label:
                return value
        raise ValueError(f"Unknown module type: {label}")


@unique
class PKSModuleSubtype(Enum):
    PKS_CIS = 1
    PKS_TRANS = 2
    PKS_ITER = 3

    @staticmethod
    def from_string(label: str) -> "PKSModuleSubtype":
        for value in PKSModuleSubtype:
            if str(value.name) == label:
                return value
        raise ValueError(f"Unknown PKS module subtype: {label}")


@unique
class NRPSDomainType(Enum):
    A = 1
    C = 2
    PCP = 3
    E = 4
    nMT = 5
    TE = 6
    TD = 7
    UNKNOWN = 8

    @staticmethod
    def from_string(label: str) -> "NRPSDomainType":
        for value in NRPSDomainType:
            if str(value.name) == label:
                return value
        raise ValueError(f"Unknown NRPS domain type: {label}")


@unique
class PKSDomainType(Enum):
    AT = 1
    KS = 2
    ACP = 3
    KR = 4
    DH = 5
    ER = 6
    TE = 7
    TD = 8
    UNKNOWN = 9
    TRANS_AT_PKS_KS_CLADE_1 = 10
    TRANS_AT_PKS_KS_CLADE_10 = 11
    TRANS_AT_PKS_KS_CLADE_100 = 12
    TRANS_AT_PKS_KS_CLADE_101 = 13
    TRANS_AT_PKS_KS_CLADE_102 = 14
    TRANS_AT_PKS_KS_CLADE_103 = 15
    TRANS_AT_PKS_KS_CLADE_104 = 16
    TRANS_AT_PKS_KS_CLADE_108 = 17
    TRANS_AT_PKS_KS_CLADE_109 = 18
    TRANS_AT_PKS_KS_CLADE_11 = 19
    TRANS_AT_PKS_KS_CLADE_111 = 20
    TRANS_AT_PKS_KS_CLADE_112 = 21
    TRANS_AT_PKS_KS_CLADE_113 = 22
    TRANS_AT_PKS_KS_CLADE_114 = 23
    TRANS_AT_PKS_KS_CLADE_115 = 24
    TRANS_AT_PKS_KS_CLADE_116 = 25
    TRANS_AT_PKS_KS_CLADE_118 = 26
    TRANS_AT_PKS_KS_CLADE_12 = 27
    TRANS_AT_PKS_KS_CLADE_120 = 28
    TRANS_AT_PKS_KS_CLADE_121 = 29
    TRANS_AT_PKS_KS_CLADE_122 = 30
    TRANS_AT_PKS_KS_CLADE_123 = 31
    TRANS_AT_PKS_KS_CLADE_124 = 32
    TRANS_AT_PKS_KS_CLADE_125 = 33
    TRANS_AT_PKS_KS_CLADE_126 = 34
    TRANS_AT_PKS_KS_CLADE_127 = 35
    TRANS_AT_PKS_KS_CLADE_128 = 36
    TRANS_AT_PKS_KS_CLADE_129 = 37
    TRANS_AT_PKS_KS_CLADE_13 = 38
    TRANS_AT_PKS_KS_CLADE_130 = 39
    TRANS_AT_PKS_KS_CLADE_131 = 40
    TRANS_AT_PKS_KS_CLADE_132 = 41
    TRANS_AT_PKS_KS_CLADE_133 = 42
    TRANS_AT_PKS_KS_CLADE_134 = 43
    TRANS_AT_PKS_KS_CLADE_135 = 44
    TRANS_AT_PKS_KS_CLADE_136 = 45
    TRANS_AT_PKS_KS_CLADE_137 = 46
    TRANS_AT_PKS_KS_CLADE_138 = 47
    TRANS_AT_PKS_KS_CLADE_139 = 48
    TRANS_AT_PKS_KS_CLADE_14 = 49
    TRANS_AT_PKS_KS_CLADE_140 = 50
    TRANS_AT_PKS_KS_CLADE_141 = 51
    TRANS_AT_PKS_KS_CLADE_142 = 52
    TRANS_AT_PKS_KS_CLADE_2 = 53
    TRANS_AT_PKS_KS_CLADE_21 = 54
    TRANS_AT_PKS_KS_CLADE_23 = 55
    TRANS_AT_PKS_KS_CLADE_25 = 56
    TRANS_AT_PKS_KS_CLADE_26 = 57
    TRANS_AT_PKS_KS_CLADE_27 = 58
    TRANS_AT_PKS_KS_CLADE_28 = 59
    TRANS_AT_PKS_KS_CLADE_30 = 60
    TRANS_AT_PKS_KS_CLADE_31 = 61
    TRANS_AT_PKS_KS_CLADE_32 = 62
    TRANS_AT_PKS_KS_CLADE_33 = 63
    TRANS_AT_PKS_KS_CLADE_34 = 64
    TRANS_AT_PKS_KS_CLADE_35 = 65
    TRANS_AT_PKS_KS_CLADE_36 = 66
    TRANS_AT_PKS_KS_CLADE_38 = 67
    TRANS_AT_PKS_KS_CLADE_39 = 68
    TRANS_AT_PKS_KS_CLADE_40 = 69
    TRANS_AT_PKS_KS_CLADE_42 = 70
    TRANS_AT_PKS_KS_CLADE_43 = 71
    TRANS_AT_PKS_KS_CLADE_44 = 72
    TRANS_AT_PKS_KS_CLADE_45 = 73
    TRANS_AT_PKS_KS_CLADE_46 = 74
    TRANS_AT_PKS_KS_CLADE_49 = 75
    TRANS_AT_PKS_KS_CLADE_5 = 76
    TRANS_AT_PKS_KS_CLADE_51 = 77
    TRANS_AT_PKS_KS_CLADE_52 = 78
    TRANS_AT_PKS_KS_CLADE_53 = 79
    TRANS_AT_PKS_KS_CLADE_55 = 80
    TRANS_AT_PKS_KS_CLADE_56 = 81
    TRANS_AT_PKS_KS_CLADE_57 = 82
    TRANS_AT_PKS_KS_CLADE_60 = 83
    TRANS_AT_PKS_KS_CLADE_61 = 84
    TRANS_AT_PKS_KS_CLADE_62 = 85
    TRANS_AT_PKS_KS_CLADE_64 = 86
    TRANS_AT_PKS_KS_CLADE_65 = 87
    TRANS_AT_PKS_KS_CLADE_66 = 88
    TRANS_AT_PKS_KS_CLADE_67 = 89
    TRANS_AT_PKS_KS_CLADE_68 = 90
    TRANS_AT_PKS_KS_CLADE_7 = 91
    TRANS_AT_PKS_KS_CLADE_70 = 92
    TRANS_AT_PKS_KS_CLADE_73 = 93
    TRANS_AT_PKS_KS_CLADE_74 = 94
    TRANS_AT_PKS_KS_CLADE_75 = 95
    TRANS_AT_PKS_KS_CLADE_76 = 96
    TRANS_AT_PKS_KS_CLADE_78 = 97
    TRANS_AT_PKS_KS_CLADE_79 = 98
    TRANS_AT_PKS_KS_CLADE_8 = 99
    TRANS_AT_PKS_KS_CLADE_80 = 100
    TRANS_AT_PKS_KS_CLADE_81 = 101
    TRANS_AT_PKS_KS_CLADE_82 = 102
    TRANS_AT_PKS_KS_CLADE_83 = 103
    TRANS_AT_PKS_KS_CLADE_84 = 104
    TRANS_AT_PKS_KS_CLADE_85 = 105
    TRANS_AT_PKS_KS_CLADE_86 = 106
    TRANS_AT_PKS_KS_CLADE_88 = 107
    TRANS_AT_PKS_KS_CLADE_89 = 108
    TRANS_AT_PKS_KS_CLADE_9 = 109
    TRANS_AT_PKS_KS_CLADE_90 = 110
    TRANS_AT_PKS_KS_CLADE_92 = 111
    TRANS_AT_PKS_KS_CLADE_93 = 112
    TRANS_AT_PKS_KS_CLADE_94 = 113
    TRANS_AT_PKS_KS_CLADE_95 = 114
    TRANS_AT_PKS_KS_CLADE_96 = 115
    TRANS_AT_PKS_KS_CLADE_97 = 116
    TRANS_AT_PKS_KS_CLADE_98 = 117
    TRANS_AT_PKS_KS_CLADE_99 = 118
    DUMMY_ALMT = 119 #Alpha-L-Methyltransferase
    DUMMY_AMT = 120 #Alpha-Methyltransferase
    DUMMY_SC = 121 #Smalles cyclase for creating pyran/furan rings
    DUMMY_ZDH = 122 #E-configured double bonds
    DUMMY_EDH = 123 #Z-configured double bonds
    DUMMY_AH = 124 #Alpha-hydroxylase
    DUMMY_GDH = 125 #Gamma-beta-dehydrogenase
    DUMMY_ZGDH = 126 #Z-Gamma-beta-dehydrogenase
    DUMMY_EGDH = 127 #E-Gamma-beta-dehydrogenase
    DUMMY_OMT = 128 #Beta-Hydroxymethyltransferase
    DUMMY_BMT = 129 #Beta-Methyltransferase

    @staticmethod
    def from_string(label: str) -> "PKSDomainType":
        for value in PKSDomainType:
            if str(value.name) == label:
                return value
        raise ValueError(f"Unknown PKS domain type: {label}")


class _Module:

    def __init__(self, nr: int, module_type: str, domains: List[Domain], *,
                 module_subtype: Union[str, None] = None,
                 starter: bool = False, terminator: bool = False) -> None:
        self.id = nr
        self.type = ModuleType.from_string(module_type)
        self.subtype = None

        if module_subtype is not None:
            if self.type.name == 'PKS':
                self.subtype = PKSModuleSubtype.from_string(module_subtype)
            else:
                raise ValueError(f"Module subtypes not supported for {self.type.name} modules.")

        self.domains = domains

        self.recognition_domain = None
        self.synthesis_domain = None
        self.carrier_domain = None
        self.termination_domain = None

        self.tailoring_domains = []

        self.is_starter_module = starter
        self.is_termination_module = terminator

        for domain in self.domains:
            if self.type.name == 'NRPS':
                assert domain.type.name in [v.name for v in NRPSDomainType]
            if self.type.name == 'PKS':
                assert domain.type.name in [v.name for v in PKSDomainType]

            if domain.used and domain.active:
                if isinstance(domain, TailoringDomain):
                    if domain.type.name not in [d.type.name for d in self.tailoring_domains]:
                        self.tailoring_domains.append(domain)
                    else:
                        raise ValueError(f"Cannot have two used tailoring domains of type {domain.type.name} in one \
module. Remove domain or set the 'used' or 'active' flag to False")

                if isinstance(domain, RecognitionDomain):
                    if not self.recognition_domain:
                        self.recognition_domain = domain
                    else:
                        raise ValueError(f"Cannot have more than one used recognition domain in one \
module. Remove a domain or set the 'used' or 'active' flag to False")

                if isinstance(domain, SynthesisDomain) and domain.is_elongating:
                    if not self.synthesis_domain:
                        self.synthesis_domain = domain
                    else:
                        raise ValueError("Cannot have more than one used and elongating synthesis domains \
in one module. Remove a domain, set the 'used' or 'active' flag to False, or set the 'is_elongating' flag to False.")

                if isinstance(domain, CarrierDomain):
                    if not self.carrier_domain:
                        self.carrier_domain = domain
                    else:
                        raise ValueError("Cannot have more than one used carrier domain in one \
module. Remove a domain or set the 'used' or 'active' flag to False")

                if isinstance(domain, TerminationDomain):
                    if not self.termination_domain:
                        self.termination_domain = domain
                    else:
                        raise ValueError("Cannot have more than one used termination domain in one \
module. Remove a domain or set the 'used' or 'active' flag to False")

        if not self.is_starter_module and not self.synthesis_domain:
            if self.type.name == 'NRPS':
                self.synthesis_domain = SynthesisDomain('C')
            elif self.type.name == 'PKS':
                self.synthesis_domain = SynthesisDomain('KS')

        if self.is_termination_module and not self.termination_domain:
            self.termination_domain = TerminationDomain("DUMMY_TE")

    def run_module(self, structure: Union[Structure, None] = None):
        raise NotImplementedError

    def get_tailoring_domain(self, domain_name: str) -> Union[TailoringDomain, None]:
        for domain in self.tailoring_domains:
            if domain.type.name == domain_name:
                return domain

        return

    def do_pks_tailoring(self, structure: Structure) -> Structure:
        kr_domain = self.get_tailoring_domain("KR")
        dh_domain = self.get_tailoring_domain("DH")
        er_domain = self.get_tailoring_domain("ER")

        if kr_domain and kr_domain.active:
            assert kr_domain.subtype is not None
            structure = kr_domain.do_tailoring(structure)
            if not kr_domain.subtype.name == 'C1' and not kr_domain.subtype.name == 'C2':
                if dh_domain and dh_domain.active:
                    structure = dh_domain.do_tailoring(structure)
                    if er_domain and er_domain.active:
                        structure = er_domain.do_tailoring(structure)

        return structure

    def do_nrps_tailoring(self, structure: Structure) -> Structure:
        e_domain = self.get_tailoring_domain('E')
        n_mt_domain = self.get_tailoring_domain('nMT')

        if e_domain and e_domain.active:
            structure = e_domain.do_tailoring(structure)
        if n_mt_domain and n_mt_domain.active:
            structure = n_mt_domain.do_tailoring(structure)

        return structure

    def release_chain(self, structure: Structure) -> Structure:
        assert structure
        assert self.termination_domain

        released_structure = self.termination_domain.release_chain(structure)

        return released_structure


class LinearPKSModule(_Module):

    def __init__(self, nr, domains, starter=False, terminator=False):
        super().__init__(nr, "PKS", domains, starter=starter, terminator=terminator, module_subtype="PKS_CIS")

    def run_module(self, structure: Union[Structure, None] = None) -> Structure:
        if structure is None:
            assert self.is_starter_module
            starter_unit = read_smiles(self.recognition_domain.substrate.smiles)
            label_pk_central_chain(starter_unit)
            structure = attach_to_domain_pk(starter_unit)
        else:
            structure = self.synthesis_domain.do_elongation(structure, self.recognition_domain.substrate)

        structure = self.do_pks_tailoring(structure)

        return structure


class IterativePKSModule(_Module):
    def __init__(self, nr, domains, starter=False, terminator=False, iterations: int = 1) -> None:
        super().__init__(nr, "PKS", domains, starter=starter, terminator=terminator, module_subtype="PKS_ITER")
        self.iterations = iterations

    def run_module(self, structure=None):
        self.do_pks_tailoring(structure)
        raise NotImplementedError


class TransATPKSModule(_Module):
    def __init__(self, nr, domains, substrate_name="MALONYL_COA", starter=False, terminator=False) -> None:
        super().__init__(nr, "PKS", domains, starter=starter, terminator=terminator, module_subtype="PKS_TRANS")
        if not self.recognition_domain:
            self.recognition_domain = RecognitionDomain("DUMMY_AT", substrate_name)

    def run_module(self, structure=None):
        self.do_pks_tailoring(structure)
        raise NotImplementedError


class NRPSModule(_Module):

    def __init__(self, nr, domains, starter=False, terminator=False):
        super().__init__(nr, "NRPS", domains, starter=starter, terminator=terminator)

    def run_module(self, structure=None) -> Structure:
        if structure is None:
            assert self.is_starter_module
            starter_unit = read_smiles(self.recognition_domain.substrate.smiles)
            label_nrp_central_chain(starter_unit)
            structure = attach_to_domain_nrp(starter_unit)
        else:
            structure = self.synthesis_domain.do_elongation(structure, self.recognition_domain.substrate)

        structure = self.do_nrps_tailoring(structure)

        return structure
