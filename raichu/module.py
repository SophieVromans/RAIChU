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
    TRANS_AT_PKS_KS_CLADE_1_ELONGATING = 10
    TRANS_AT_PKS_KS_CLADE_10_ELONGATING = 11
    TRANS_AT_PKS_KS_CLADE_100_ELONGATING = 12
    TRANS_AT_PKS_KS_CLADE_101_ELONGATING = 13
    TRANS_AT_PKS_KS_CLADE_102_NON_ELONGATING = 14
    TRANS_AT_PKS_KS_CLADE_103_ELONGATING = 15
    TRANS_AT_PKS_KS_CLADE_104_ELONGATING = 16
    TRANS_AT_PKS_KS_CLADE_108_ELONGATING = 17
    TRANS_AT_PKS_KS_CLADE_109_ELONGATING = 18
    TRANS_AT_PKS_KS_CLADE_11_ELONGATING = 19
    TRANS_AT_PKS_KS_CLADE_111_ELONGATING = 20
    TRANS_AT_PKS_KS_CLADE_112_ELONGATING = 21
    TRANS_AT_PKS_KS_CLADE_113_ELONGATING = 22
    TRANS_AT_PKS_KS_CLADE_114_ELONGATING = 23
    TRANS_AT_PKS_KS_CLADE_115_ELONGATING = 24
    TRANS_AT_PKS_KS_CLADE_116_ELONGATING = 25
    TRANS_AT_PKS_KS_CLADE_118_ELONGATING = 26
    TRANS_AT_PKS_KS_CLADE_12_ELONGATING = 27
    TRANS_AT_PKS_KS_CLADE_120_NON_ELONGATING = 28
    TRANS_AT_PKS_KS_CLADE_121_ELONGATING = 29
    TRANS_AT_PKS_KS_CLADE_122_ELONGATING = 30
    TRANS_AT_PKS_KS_CLADE_123_ELONGATING = 31
    TRANS_AT_PKS_KS_CLADE_124_ELONGATING = 32
    TRANS_AT_PKS_KS_CLADE_125_ELONGATING = 33
    TRANS_AT_PKS_KS_CLADE_126_ELONGATING = 34
    TRANS_AT_PKS_KS_CLADE_127_ELONGATING = 35
    TRANS_AT_PKS_KS_CLADE_128_ELONGATING = 36
    TRANS_AT_PKS_KS_CLADE_129_ELONGATING = 37
    TRANS_AT_PKS_KS_CLADE_13_ELONGATING = 38
    TRANS_AT_PKS_KS_CLADE_130_ELONGATING = 39
    TRANS_AT_PKS_KS_CLADE_131_ELONGATING = 40
    TRANS_AT_PKS_KS_CLADE_132_ELONGATING = 41
    TRANS_AT_PKS_KS_CLADE_133_ELONGATING = 42
    TRANS_AT_PKS_KS_CLADE_134_ELONGATING = 43
    TRANS_AT_PKS_KS_CLADE_135_ELONGATING = 44
    TRANS_AT_PKS_KS_CLADE_136_ELONGATING = 45
    TRANS_AT_PKS_KS_CLADE_137_ELONGATING = 46
    TRANS_AT_PKS_KS_CLADE_138_ELONGATING = 47
    TRANS_AT_PKS_KS_CLADE_139_NON_ELONGATING = 48
    TRANS_AT_PKS_KS_CLADE_14_ELONGATING = 49
    TRANS_AT_PKS_KS_CLADE_140_ELONGATING = 50
    TRANS_AT_PKS_KS_CLADE_141_ELONGATING = 51
    TRANS_AT_PKS_KS_CLADE_142_NON_ELONGATING = 52
    TRANS_AT_PKS_KS_CLADE_2_ELONGATING = 53
    TRANS_AT_PKS_KS_CLADE_21_ELONGATING = 54
    TRANS_AT_PKS_KS_CLADE_23_ELONGATING = 55
    TRANS_AT_PKS_KS_CLADE_25_ELONGATING = 56
    TRANS_AT_PKS_KS_CLADE_26_ELONGATING = 57
    TRANS_AT_PKS_KS_CLADE_27_ELONGATING = 58
    TRANS_AT_PKS_KS_CLADE_28_ELONGATING = 59
    TRANS_AT_PKS_KS_CLADE_30_NON_ELONGATING = 60
    TRANS_AT_PKS_KS_CLADE_31_NON_ELONGATING = 61
    TRANS_AT_PKS_KS_CLADE_32_ELONGATING = 62
    TRANS_AT_PKS_KS_CLADE_33_ELONGATING = 63
    TRANS_AT_PKS_KS_CLADE_34_ELONGATING = 64
    TRANS_AT_PKS_KS_CLADE_35_ELONGATING = 65
    TRANS_AT_PKS_KS_CLADE_36_NON_ELONGATING = 66
    TRANS_AT_PKS_KS_CLADE_38_NON_ELONGATING = 67
    TRANS_AT_PKS_KS_CLADE_39_ELONGATING = 68
    TRANS_AT_PKS_KS_CLADE_40_ELONGATING = 69
    TRANS_AT_PKS_KS_CLADE_42_ELONGATING = 70
    TRANS_AT_PKS_KS_CLADE_43_ELONGATING = 71
    TRANS_AT_PKS_KS_CLADE_44_NON_ELONGATING = 72
    TRANS_AT_PKS_KS_CLADE_45_NON_ELONGATING = 73
    TRANS_AT_PKS_KS_CLADE_46_NON_ELONGATING = 74
    TRANS_AT_PKS_KS_CLADE_49_ELONGATING = 75
    TRANS_AT_PKS_KS_CLADE_5_ELONGATING = 76
    TRANS_AT_PKS_KS_CLADE_51_ELONGATING = 77
    TRANS_AT_PKS_KS_CLADE_52_ELONGATING = 78
    TRANS_AT_PKS_KS_CLADE_53_ELONGATING = 79
    TRANS_AT_PKS_KS_CLADE_55_ELONGATING = 80
    TRANS_AT_PKS_KS_CLADE_56_ELONGATING = 81
    TRANS_AT_PKS_KS_CLADE_57_ELONGATING = 82
    TRANS_AT_PKS_KS_CLADE_60_ELONGATING = 83
    TRANS_AT_PKS_KS_CLADE_61_ELONGATING = 84
    TRANS_AT_PKS_KS_CLADE_62_ELONGATING = 85
    TRANS_AT_PKS_KS_CLADE_64_NON_ELONGATING = 86
    TRANS_AT_PKS_KS_CLADE_65_ELONGATING = 87
    TRANS_AT_PKS_KS_CLADE_66_ELONGATING = 88
    TRANS_AT_PKS_KS_CLADE_67_ELONGATING = 89
    TRANS_AT_PKS_KS_CLADE_68_ELONGATING = 90
    TRANS_AT_PKS_KS_CLADE_7_ELONGATING = 91
    TRANS_AT_PKS_KS_CLADE_70_ELONGATING = 92
    TRANS_AT_PKS_KS_CLADE_73_ELONGATING = 93
    TRANS_AT_PKS_KS_CLADE_74_ELONGATING = 94
    TRANS_AT_PKS_KS_CLADE_75_NON_ELONGATING = 95
    TRANS_AT_PKS_KS_CLADE_76_NON_ELONGATING = 96
    TRANS_AT_PKS_KS_CLADE_78_NON_ELONGATING = 97
    TRANS_AT_PKS_KS_CLADE_79_ELONGATING = 98
    TRANS_AT_PKS_KS_CLADE_8_ELONGATING = 99
    TRANS_AT_PKS_KS_CLADE_80_NON_ELONGATING = 100
    TRANS_AT_PKS_KS_CLADE_81_NON_ELONGATING = 101
    TRANS_AT_PKS_KS_CLADE_82_ELONGATING = 102
    TRANS_AT_PKS_KS_CLADE_83_NON_ELONGATING = 103
    TRANS_AT_PKS_KS_CLADE_84_ELONGATING = 104
    TRANS_AT_PKS_KS_CLADE_85_NON_ELONGATING = 105
    TRANS_AT_PKS_KS_CLADE_86_ELONGATING = 106
    TRANS_AT_PKS_KS_CLADE_88_NON_ELONGATING = 107
    TRANS_AT_PKS_KS_CLADE_89_ELONGATING = 108
    TRANS_AT_PKS_KS_CLADE_9_ELONGATING = 109
    TRANS_AT_PKS_KS_CLADE_90_ELONGATING = 110
    TRANS_AT_PKS_KS_CLADE_92_NON_ELONGATING = 111
    TRANS_AT_PKS_KS_CLADE_93_NON_ELONGATING = 112
    TRANS_AT_PKS_KS_CLADE_94_NON_ELONGATING = 113
    TRANS_AT_PKS_KS_CLADE_95_ELONGATING = 114
    TRANS_AT_PKS_KS_CLADE_96_ELONGATING = 115
    TRANS_AT_PKS_KS_CLADE_97_ELONGATING = 116
    TRANS_AT_PKS_KS_CLADE_98_ELONGATING = 117
    TRANS_AT_PKS_KS_CLADE_99_ELONGATING = 118
    ALMT = 119 #Alpha-L-Methyltransferase
    AMT = 120 #Alpha-Methyltransferase
    SC = 121 #Smalles cyclase for creating pyran/furan rings
    ZDH = 122 #E-configured double bonds
    EDH = 123 #Z-configured double bonds
    AH = 124 #Alpha-hydroxylase
    GDH = 125 #Gamma-beta-dehydrogenase
    ZGDH = 126 #Z-Gamma-beta-dehydrogenase
    EGDH = 127 #E-Gamma-beta-dehydrogenase
    OMT = 128 #Beta-Hydroxymethyltransferase
    BMT = 129 #Beta-Methyltransferase

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
