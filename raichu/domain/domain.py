from typing import Union, Tuple
from pikachu.chem.structure import Structure
from pikachu.general import read_smiles
from raichu.substrate import NRPSSubstrate, PKSSubstrate
from raichu.reactions.pks_tailoring_reactions import ketoreduction, enoylreduction, dehydration
from raichu.reactions.nrps_tailoring_reactions import epimerize, n_methylate
from raichu.reactions.pks_elongation_reactions import pks_elongation
from raichu.reactions.nrps_elongation_reactions import nrps_elongation
from raichu.reactions.chain_release import release_linear_reduction, release_linear_thioesterase
from raichu.domain.domain_types import DomainSuperClass, RecognitionDomainType, CarrierDomainType, \
    TailoringDomainType, SynthesisDomainType, TerminationDomainType, KSDomainSubtype, KRDomainSubtype, ERDomainSubtype

from dataclasses import dataclass    
    
 
@dataclass
class Domain:
    supertype: DomainSuperClass
    type: Union[TailoringDomainType, CarrierDomainType, SynthesisDomainType, RecognitionDomainType,
                TerminationDomainType]
    subtype: Union[None, KRDomainSubtype, KSDomainSubtype]
    domain_name: Union[str, None]
    active: bool = True
    gene: Union[str, None] = None
    used: bool = True

    def __post_init__(self):
        if self.domain_name is None:
            self.domain_name = self.type.name

    def __repr__(self):
        return self.domain_name

    def set_gene(self, gene_name: str) -> None:
        self.gene = gene_name


class UnknownDomain(Domain):
    def __init__(self, domain_name: str, active: bool = True) -> None:
        superclass = DomainSuperClass.from_string("TAILORING")
        domain_type = TailoringDomainType.from_string("UNKNOWN")
        assert domain_name
        super().__init__(superclass, domain_type, None, domain_name, active=active, used=False)


class TailoringDomain(Domain):
    def __init__(self, domain_type: str, domain_subtype: Union[str, None] = None, active: bool = True,
                 domain_name: Union[str, None] = None, used: bool = True) -> None:
        superclass = DomainSuperClass.from_string("TAILORING")
        domain_type = TailoringDomainType.from_string(domain_type)
        if domain_subtype is not None:
            if domain_type.name == 'KR' or domain_type.name == 'DUMMY_KR':
                domain_subtype = KRDomainSubtype.from_string(domain_subtype)
            elif domain_type.name == 'ER' or domain_type.name == 'DUMMY_ER':
                domain_subtype = ERDomainSubtype.from_string(domain_subtype)
            else:
                raise ValueError(f"RAIChU does not support domain subtypes for {domain_type.name}")
        else:
            if domain_type.name == 'KR' or domain_type.name == 'DUMMY_KR':
                domain_subtype = KRDomainSubtype.from_string("UNKNOWN")
            elif domain_type.name == 'ER' or domain_type.name == 'DUMMY_ER':
                domain_subtype = ERDomainSubtype.from_string("UNKNOWN")

        super().__init__(superclass, domain_type, domain_subtype, domain_name, active=active, used=used)

    def do_tailoring(self, structure) -> Tuple[Structure, bool]:
        """
        Performs tailoring reaction
        """
        if self.type.name == 'KR' or self.type.name == 'DUMMY_KR':
            return ketoreduction(structure, self.subtype)
        elif self.type.name == 'DH' or self.type.name == 'DUMMY_DH':
            return dehydration(structure)
        elif self.type.name == 'ER' or self.type.name == 'DUMMY_ER':
            return enoylreduction(structure, self.subtype)
        elif self.type.name == 'E':
            return epimerize(structure)
        elif self.type.name == 'nMT':
            return n_methylate(structure)
        else:
            raise Warning(f"Tailoring domain {self.domain_name} not recognised by RAIChU. Ignored.")


class SynthesisDomain(Domain):
    def __init__(self, domain_type: str, domain_subtype: Union[str, None] = None, active: bool = True,
                 domain_name: Union[str, None] = None, is_elongating: bool = True,
                 used: bool = True) -> None:
        superclass = DomainSuperClass.from_string("SYNTHESIS")
        domain_type = SynthesisDomainType.from_string(domain_type)

        if domain_subtype is not None:
            if domain_type.name == 'KS' or domain_type.name == "DUMMY_KS":
                domain_subtype = KSDomainSubtype.from_string(domain_subtype)
            else:
                raise ValueError(f"RAIChU does not support domain subtypes for {domain_type.name}")

        super().__init__(superclass, domain_type, domain_subtype, domain_name, active=active, used=used)
        self.is_elongating = is_elongating

    def do_elongation(self, structure, substrate) -> Structure:
        """
        Performs elongation reaction
        """
        if self.type.name == 'C' or self.type.name == "DUMMY_C":
            if self.is_elongating:
                building_block = read_smiles(substrate.smiles)
                return nrps_elongation(building_block, structure)
        elif self.type.name == 'KS' or self.type.name == "DUMMY_KS":
            if self.is_elongating:
                if self.subtype is None or self.subtype.name == 'CIS' or self.subtype.name == 'UNKNOWN':
                    building_block = substrate.elongation_monomer
                    return pks_elongation(structure, building_block)
                else:
                    raise ValueError(f"RAIChU does not support domain subtype {self.subtype.name}")
        else:
            raise ValueError(f"Unknown type of synthesis domain: {self.type}")


class RecognitionDomain(Domain):
    def __init__(self, domain_type: str, substrate_name: str, domain_subtype: Union[str, None] = None,
                 active: bool = True, domain_name: Union[str, None] = None, used: bool = True) -> None:

        superclass = DomainSuperClass.from_string("RECOGNITION")
        domain_type = RecognitionDomainType.from_string(domain_type)
 

        if domain_subtype is not None:
            raise ValueError(f"RAIChU does not support domain subtypes for {domain_type.name}")

        super().__init__(superclass, domain_type, domain_subtype, domain_name, active=active, used=used)

        if self.type.name == 'A' or self.type.name == 'DUMMY_A':
            self.substrate = NRPSSubstrate(substrate_name)
        elif self.type.name == 'AT' or self.type.name == 'DUMMY_AT':
            self.substrate = PKSSubstrate(substrate_name)
        else:
            raise ValueError(f"Unknown type of recognition domain: {self.type}")

    @property
    def smiles(self):
        return self.substrate.smiles

    @property
    def substrate_name(self):
        return self.substrate.name


class CarrierDomain(Domain):
    def __init__(self, domain_type: str, domain_subtype: Union[str, None] = None,
                 domain_name: Union[str, None] = None, active: bool = True,
                 used: bool = True) -> None:
        superclass = DomainSuperClass.from_string("CARRIER")
        domain_type = CarrierDomainType.from_string(domain_type)

        if domain_subtype is not None:
            raise ValueError(f"RAIChU does not support domain subtypes for {domain_type.name}")

        super().__init__(superclass, domain_type, domain_subtype, domain_name, active=active, used=used)


class TerminationDomain(Domain):
    def __init__(self, domain_type: str, domain_subtype: Union[str, None] = None, active: bool = True,
                 domain_name: Union[str, None] = None, used: bool = True) -> None:
        superclass = DomainSuperClass.from_string("TERMINATION")
        domain_type = TerminationDomainType.from_string(domain_type)

        if domain_subtype is not None:
            raise ValueError(f"RAIChU does not support domain subtypes for {domain_type.name}")

        super().__init__(superclass, domain_type, domain_subtype, domain_name, active=active, used=used)

    def release_chain(self, structure) -> Structure:
        """
        Performs chain release reaction and returns True if reaction was successful, False otherwise
        """
        if self.type.name == "TE" or self.type.name == "DUMMY_TE":
            return release_linear_thioesterase(structure)
        elif self.type.name == "TD" or self.type.name == "DUMMY_TD":
            return release_linear_reduction(structure)
        else:
            raise ValueError(f"Unknown type of termination domain: {self.type}")
