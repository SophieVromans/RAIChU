from typing import List, Union

from pikachu.general import read_smiles
from pikachu.chem.structure import Structure
from pikachu.chem.atom import Atom

from raichu.class_domain import CarrierDomain as ClassCarrierDomain
from raichu.domain.domain import Domain, TailoringDomain, RecognitionDomain, \
    SynthesisDomain, CarrierDomain, TerminationDomain
from raichu.central_chain_detection.label_central_chain import label_pk_central_chain, label_nrp_central_chain
from raichu.attach_to_domain import attach_to_domain_pk, attach_to_domain_nrp
from enum import Enum, unique
from raichu.substrate import PKSSubstrate


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
    CYC = 6
    TE = 7
    TD = 8
    UNKNOWN = 9
    OX = 10

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
    DUMMY_ALMT = 10  # Alpha-L-Methyltransferase
    DUMMY_AMT = 11  # Alpha-Methyltransferase
    DUMMY_SC = 12  # Smalles cyclase for creating pyran/furan rings
    DUMMY_ZDH = 13  # E-configured double bonds
    DUMMY_EDH = 14  # Z-configured double bonds
    DUMMY_AH = 15  # Alpha-hydroxylase
    DUMMY_GDH = 16  # Gamma-beta-dehydrogenase
    DUMMY_ZGDH = 17  # Z-Gamma-beta-dehydrogenase
    DUMMY_EGDH = 18  # E-Gamma-beta-dehydrogenase
    DUMMY_OMT = 19  # Beta-Hydroxymethyltransferase
    DUMMY_BMT = 20  # Beta-Methyltransferase

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

        self.is_broken = False
        self.is_starter_module = starter
        self.is_termination_module = terminator

        for domain in self.domains:
            if self.type.name == 'NRPS':
                if domain.type.name not in [v.name for v in NRPSDomainType]:
                    raise ValueError(f'Unknown NRPS domain type {domain.type.name}')
            if self.type.name == 'PKS':
                if domain.type.name not in [v.name for v in PKSDomainType]:
                    raise ValueError(
                        f'Unknown PKS domain type {domain.type.name}')

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
                        if not self.is_termination_module:
                            self.termination_domain.used = False
                    else:
                        raise ValueError("Cannot have more than one used termination domain in one \
module. Remove a domain or set the 'used' or 'active' flag to False")
        # I would rather implement a "dont do anything"- broken module
        if not self.is_starter_module and not self.synthesis_domain:
            if "CYC" in [domain.type.name for domain in self.tailoring_domains]:
                self.synthesis_domain = SynthesisDomain("DUMMY_C")
            else: 
                self.is_broken = True

        if self.is_termination_module and not self.termination_domain:
            self.termination_domain = TerminationDomain("DUMMY_TE")

    def run_module(self, structure: Union[Structure, None] = None):
        raise NotImplementedError

    def get_tailoring_domain(self, domain_name: str) -> Union[TailoringDomain, None]:
        for domain in self.tailoring_domains:
            if domain.type.name == domain_name:
                return domain

        return

    def add_module_label(self, structure):
        for atom in structure.graph:
            if type(atom) == Atom and not atom.annotations.has_annotation('module_nr'):
                atom.annotations.add_annotation('module_nr', f"module_{self.id + 1:02}")

            if type(atom) == ClassCarrierDomain:
                for neighbour in structure.graph[atom]:
                    if neighbour.type == 'S':
                        if not neighbour.annotations.has_annotation('module_nr'):
                            neighbour.annotations.add_annotation('module_nr', f"module_{self.id + 1:02}")
                        else:
                            neighbour.annotations.module_nr = f"module_{self.id + 1:02}"

    def do_pks_tailoring(self, structure: Structure) -> Structure:
        kr_domain = self.get_tailoring_domain("KR")
        dh_domain = self.get_tailoring_domain("DH")
        er_domain = self.get_tailoring_domain("ER")

        if kr_domain and kr_domain.active and kr_domain.used:
            assert kr_domain.subtype is not None
            structure, kr_tailored = kr_domain.do_tailoring(structure)
            if kr_tailored:
                if not kr_domain.subtype.name == 'C1' and not kr_domain.subtype.name == 'C2':
                    if dh_domain and dh_domain.active and dh_domain.used:
                        structure, dh_tailored = dh_domain.do_tailoring(structure)
                        if dh_tailored:
                            if er_domain and er_domain.active and er_domain.used:
                                structure, er_tailored = er_domain.do_tailoring(structure)
                                if not er_tailored:
                                    er_domain.used = False

                        else:
                            dh_domain.used = False
                            if er_domain:
                                er_domain.used = False
            else:
                kr_domain.used = False
                if dh_domain:
                    dh_domain.used = False
                if er_domain:
                    er_domain.used = False
        else:
            if dh_domain:
                dh_domain.used = False
            if er_domain:
                er_domain.used = False

        return structure

    def do_nrps_tailoring(self, structure: Structure) -> Structure:
        e_domain = self.get_tailoring_domain('E')
        n_mt_domain = self.get_tailoring_domain('nMT')
        cyc_domain = self.get_tailoring_domain('CYC')
        ox_domain = self.get_tailoring_domain('OX')

        if e_domain and e_domain.active and e_domain.used:
            structure, epimerized = e_domain.do_tailoring(structure)
            if not epimerized:
                e_domain.used = False
        if n_mt_domain and n_mt_domain.active and n_mt_domain.used:
            structure, methylated = n_mt_domain.do_tailoring(structure)
            if not methylated:
                n_mt_domain.used = False
        if cyc_domain and cyc_domain.active and cyc_domain.used:
            structure, cyclised = cyc_domain.do_tailoring(structure)
            if not cyclised:
                cyc_domain.used = False

            if ox_domain and ox_domain.active and ox_domain.used:
                structure, oxidated = ox_domain.do_tailoring(structure)
                if not oxidated:
                    ox_domain.used = False

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
        if self.is_broken:
            return structure
        else:
            if self.recognition_domain:
                if structure is None:
                    assert self.is_starter_module
                    structure = self.recognition_domain.substrate.starter_monomer.attach_to_acp()
                else:
                    structure = self.synthesis_domain.do_elongation(structure, self.recognition_domain.substrate)

                structure = self.do_pks_tailoring(structure)
            self.add_module_label(structure)
            return structure


class IterativePKSModule(_Module):
    def __init__(self, nr, domains, starter=False, terminator=False, iterations: int = 1) -> None:
        super().__init__(nr, "PKS", domains, starter=starter, terminator=terminator, module_subtype="PKS_ITER")
        self.iterations = iterations

    def run_module(self, structure=None):
        if self.is_broken:
            return structure
        else:
            for i in range(0, self.iterations):
                if self.recognition_domain:
                    # make new substrate to enable reuse of substrate
                    substrate = PKSSubstrate(self.recognition_domain.substrate_name)
                    if structure is None:
                        assert self.is_starter_module
                        structure = substrate.starter_monomer.attach_to_acp()
                    else:
                        structure = self.synthesis_domain.do_elongation(
                            structure, substrate)

                    structure = self.do_pks_tailoring(structure)
            self.add_module_label(structure)
            return structure


class TransATPKSModule(_Module):
    def __init__(self, nr, domains, substrate_name="MALONYL_COA", starter=False, terminator=False) -> None:
        super().__init__(nr, "PKS", domains, starter=starter, terminator=terminator, module_subtype="PKS_TRANS")
        if not self.recognition_domain:
            self.recognition_domain = RecognitionDomain("DUMMY_AT", substrate_name)

    def do_pks_tailoring(self, structure: Structure) -> Structure:
        
        kr_domain = self.get_tailoring_domain("KR")
        if not kr_domain:
            kr_domain = self.get_tailoring_domain("DUMMY_KR")
        dh_domain = self.get_tailoring_domain("DH")
        if not dh_domain:
            dh_domain = self.get_tailoring_domain("DUMMY_DH")
        edh_domain = self.get_tailoring_domain("EDH")
        if not edh_domain:
            edh_domain = self.get_tailoring_domain("DUMMY_EDH")
        zdh_domain = self.get_tailoring_domain("ZDH")
        if not zdh_domain:
            zdh_domain = self.get_tailoring_domain("DUMMY_ZDH")
        er_domain = self.get_tailoring_domain("ER")
        if not er_domain:
            er_domain = self.get_tailoring_domain("DUMMY_ER")
        almt_domain = self.get_tailoring_domain("ALMT")
        if not almt_domain:
            almt_domain = self.get_tailoring_domain("DUMMY_ALMT")
        amt_domain = self.get_tailoring_domain("AMT")
        if not amt_domain:
            amt_domain = self.get_tailoring_domain("DUMMY_AMT")
        emo_domain = self.get_tailoring_domain("EMO")
        if not amt_domain:
            emo_domain = self.get_tailoring_domain("DUMMY_EMO")
        sc_domain = self.get_tailoring_domain("SC")
        if not sc_domain:
            sc_domain = self.get_tailoring_domain("DUMMY_SC")
        ah_domain = self.get_tailoring_domain("AH")
        if not ah_domain:
            ah_domain = self.get_tailoring_domain("DUMMY_AH")
        gdh_domain = self.get_tailoring_domain("GDH")
        if not gdh_domain:
            gdh_domain = self.get_tailoring_domain("DUMMY_GDH")
        zgdh_domain = self.get_tailoring_domain("ZGDH")
        if not zgdh_domain:
            zgdh_domain = self.get_tailoring_domain("DUMMY_ZGDH")
        egdh_domain = self.get_tailoring_domain("EGDH")
        if not egdh_domain:
            egdh_domain = self.get_tailoring_domain("DUMMY_EGDH")
        omt_domain = self.get_tailoring_domain("OMT")
        if not omt_domain:
            omt_domain = self.get_tailoring_domain("DUMMY_OMT")
        bmt_domain = self.get_tailoring_domain("BMT")
        if not bmt_domain:
            bmt_domain = self.get_tailoring_domain("DUMMY_BMT")

        if amt_domain and amt_domain.active and amt_domain.used:
            structure, amt_tailored = amt_domain.do_tailoring(structure)
            if not amt_tailored:
                amt_domain.used = False
            if emo_domain and emo_domain.active and emo_domain.used:
                structure, emo_tailored = emo_domain.do_tailoring(structure)
                if not emo_tailored:
                    emo_domain.used = False
            structure.refresh_structure()
        if almt_domain and almt_domain.active and almt_domain.used:
            structure, almt_tailored = almt_domain.do_tailoring(structure)
            if not almt_tailored:
                almt_domain.used = False
            if emo_domain and emo_domain.active and emo_domain.used:
                structure, emo_tailored = emo_domain.do_tailoring(structure)
                if not emo_tailored:
                    emo_domain.used = False
        if ah_domain and ah_domain.active and ah_domain.used:
            structure, ah_tailored = ah_domain.do_tailoring(structure)
            if not ah_tailored:
                ah_domain.used = False
        if kr_domain and kr_domain.active and kr_domain.used:
            assert kr_domain.subtype is not None

            structure, kr_tailored = kr_domain.do_tailoring(structure)
            if not kr_tailored or kr_domain.subtype.name == 'C1' or kr_domain.subtype.name == 'C2':
                kr_domain.used = False
                if omt_domain:
                    omt_domain.used = False
                if dh_domain:
                    dh_domain.used = False
                if bmt_domain:
                    bmt_domain.used = False
                if er_domain:
                    er_domain.used = False
            if not kr_domain.subtype.name == 'C1' and not kr_domain.subtype.name == 'C2':
                if bmt_domain and bmt_domain.active and bmt_domain.used:
                    structure, bmt_tailored = bmt_domain.do_tailoring(structure)
                    if not bmt_tailored:
                        bmt_domain.used = False
                if omt_domain and omt_domain.active and omt_domain.used:
                    structure, omt_tailored = omt_domain.do_tailoring(structure)
                    if not omt_tailored:
                        omt_domain.used = False
                elif dh_domain and dh_domain.active and dh_domain.used:
                    structure, dh_tailored = dh_domain.do_tailoring(structure)
                    if not dh_tailored:
                        dh_domain.used = False
                        if er_domain:
                            er_domain.used = False
                    else:
                        if er_domain and er_domain.active and er_domain.used:
                            structure, er_tailored = er_domain.do_tailoring(structure)
                            if not er_tailored:
                                er_domain.used = False
                elif edh_domain and edh_domain.active and edh_domain.used:
                    structure, edh_tailored = edh_domain.do_tailoring(structure)
                    if not edh_tailored:
                        edh_domain.used = False

                        if er_domain:
                            er_domain.used = False
                        
                    else:
                        if er_domain and er_domain.active and er_domain.used:
                            structure, er_tailored = er_domain.do_tailoring(structure)
                            if not er_tailored:
                                er_domain.used = False
                elif zdh_domain and zdh_domain.active and zdh_domain.used:
                    structure, zdh_tailored = zdh_domain.do_tailoring(structure)
                    if not zdh_tailored:
                        zdh_domain.used = False
                        if er_domain:
                            er_domain.used = False
                    else:
                        if er_domain and er_domain.active and er_domain.used:
                            structure, er_tailored = er_domain.do_tailoring(structure)
                            if not er_tailored:
                                er_domain.used = False
                elif gdh_domain and gdh_domain.active and gdh_domain.used:
                    structure, gdh_domain_tailored = gdh_domain.do_tailoring(structure)
                    if not gdh_domain_tailored:
                        gdh_domain.used = False
                        if er_domain:
                            er_domain.used = False
                    else:
                        if er_domain and er_domain.active and er_domain.used:
                            structure, er_tailored = er_domain.do_tailoring(structure)
                            if not er_tailored:
                                er_domain.used = False
                elif egdh_domain and egdh_domain.active and egdh_domain.used:
                    structure, egdh_domain_tailored = egdh_domain.do_tailoring(structure)
                    if not egdh_domain_tailored:
                        egdh_domain.used = False
                        if er_domain:
                            er_domain.used = False
                    else:
                        if er_domain and er_domain.active and er_domain.used:
                            structure, er_tailored = er_domain.do_tailoring(structure)
                            if not er_tailored:
                                er_domain.used = False
                elif zgdh_domain and zgdh_domain.active and zgdh_domain.used:
                    structure, zgdh_domain_tailored = zgdh_domain.do_tailoring(structure)
                    if not zgdh_domain_tailored:
                        zgdh_domain.used = False
                        if er_domain:
                            er_domain.used = False
                    else:
                        if er_domain and er_domain.active and er_domain.used:
                            structure, er_tailored = er_domain.do_tailoring(structure)
                            if not er_tailored:
                                er_domain.used = False
        if sc_domain and sc_domain.active and sc_domain.used:
            structure, sc_tailored = sc_domain.do_tailoring(structure)
            if not sc_tailored:
                sc_domain.used = False

        structure.refresh_structure()
        return structure

    def run_module(self, structure: Union[Structure, None] = None) -> Structure:
        if self.is_broken:
            return structure
        else:
            if self.recognition_domain:
                if structure is None:
                    assert self.is_starter_module
                    structure = self.recognition_domain.substrate.starter_monomer.attach_to_acp()
                else:
                    structure = self.synthesis_domain.do_elongation(structure, self.recognition_domain.substrate)
                structure = self.do_pks_tailoring(structure)
            self.add_module_label(structure)
            return structure


class NRPSModule(_Module):

    def __init__(self, nr, domains, starter=False, terminator=False):
        super().__init__(nr, "NRPS", domains, starter=starter, terminator=terminator)

    def run_module(self, structure=None) -> Structure:
        if self.is_broken:
            return structure
        else:
            if self.recognition_domain:
                if structure is None:
                    assert self.is_starter_module
                    starter_unit = read_smiles(self.recognition_domain.substrate.smiles)
                    label_nrp_central_chain(starter_unit)
                    structure = attach_to_domain_nrp(starter_unit)
                else:
                    structure = self.synthesis_domain.do_elongation(structure, self.recognition_domain.substrate)

                structure = self.do_nrps_tailoring(structure)
            self.add_module_label(structure)
            return structure
