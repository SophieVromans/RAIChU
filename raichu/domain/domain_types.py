from enum import Enum, unique


@unique
class DomainSuperClass(Enum):
    """
    An Enum representing the superclass of domains supported by RAIChU
    """
    RECOGNITION = 1
    SYNTHESIS = 2
    CARRIER = 3
    TAILORING = 4
    TERMINATION = 5

    @staticmethod
    def from_string(label: str) -> "DomainSuperClass":
        for value in DomainSuperClass:
            if str(value.name) == label:
                return value
        raise ValueError(f"Unknown domain superclass: {label}")


@unique
class RecognitionDomainType(Enum):
    """
    An Enum representing the types of recognition domains supported by RAIChU
    """
    A = 1
    AT = 2
    DUMMY_A = 3
    DUMMY_AT = 4

    @staticmethod
    def from_string(label: str) -> "RecognitionDomainType":
        for value in RecognitionDomainType:
            if str(value.name) == label:
                return value
        raise ValueError(f"Unknown recognition domain type: {label}")


@unique
class TerminationDomainType(Enum):
    """
    An Enum representing the types of termination domains supported by RAIChU
    """
    TE = 1
    TD = 2
    DUMMY_TE = 3
    DUMMY_TD = 4

    @staticmethod
    def from_string(label: str) -> "TerminationDomainType":
        for value in TerminationDomainType:
            if str(value.name) == label:
                return value
        raise ValueError(f"Unknown termination domain type: {label}")


@unique
class SynthesisDomainType(Enum):
    """
    An Enum representing the types of synthesis domains supported by RAIChU
    """
    C = 1
    KS = 2
    DUMMY_C = 3
    DUMMY_KS = 4

    @staticmethod
    def from_string(label: str) -> "SynthesisDomainType":
        for value in SynthesisDomainType:
            if str(value.name) == label:
                return value
        raise ValueError(f"Unknown synthesis domain type: {label}")


@unique
class CarrierDomainType(Enum):
    """
    An Enum representing the types of carrier domains supported by RAIChU
    """
    PCP = 1
    ACP = 2
    DUMMY_PCP = 3
    DUMMY_ACP = 4

    @staticmethod
    def from_string(label: str) -> "CarrierDomainType":
        for value in CarrierDomainType:
            if str(value.name) == label:
                return value
        raise ValueError(f"Unknown carrier domain type: {label}")


@unique
class TailoringDomainType(Enum):
    """
    An Enum representing the types of tailoring domains supported by RAIChU
    """
    ER = 1
    KR = 2
    DH = 3
    E = 5
    nMT = 6
    UNKNOWN = 7
    DUMMY_ALMT = 8 #Alpha-L-Methyltransferase
    DUMMY_AMT = 9 #Alpha-Methyltransferase
    DUMMY_SC = 10 #Smalles cyclase for creating pyran/furan rings
    DUMMY_ZDH = 11 #E-configured double bonds
    DUMMY_EDH = 12 #Z-configured double bonds
    DUMMY_AH = 13 #Alpha-hydroxylase
    DUMMY_GDH = 14 #Gamma-beta-dehydrogenase
    DUMMY_ZGDH = 15 #Z-Gamma-beta-dehydrogenase
    DUMMY_EGDH = 16 #E-Gamma-beta-dehydrogenase
    DUMMY_OMT = 17 #Beta-Hydroxymethyltransferase
    DUMMY_BMT = 18 #Beta-Methyltransferase
    DUMMY_ER = 19
    DUMMY_KR = 20
    DUMMY_DH = 21
    DUMMY_E = 22
    DUMMY_nMT = 23
    @staticmethod
    def from_string(label: str) -> "TailoringDomainType":
        for value in TailoringDomainType:
            if str(value.name) == label:
                return value
        raise ValueError(f"Unknown tailoring domain type: {label}")


@unique
class KRDomainSubtype(Enum):
    """
    An Enum representing the subtypes of KR domain supported by RAIChU
    """
    A1 = 1
    A2 = 2
    B1 = 3
    B2 = 4
    C1 = 5
    C2 = 6
    UNKNOWN = 7

    @staticmethod
    def from_string(label: str) -> "KRDomainSubtype":
        for value in KRDomainSubtype:
            if str(value.name) == label:
                return value
        raise ValueError(f"Unknown KR domain subtype: {label}")


@unique
class KSDomainSubtype(Enum):
    """
    An Enum representing the subtypes of KS domain supported by RAIChU
    """
    CIS = 1
    UNKNOWN = 2
    TRANS_AT_PKS_NON_ELOGATING_BETA_L_OH = 3
    TRANS_AT_PKS_BETA_OH = 4
    TRANS_AT_PKS_NON_ELOGATING_ALPHAME_EDB = 5
    TRANS_AT_PKS_OXI = 6
    TRANS_AT_PKS_ST = 7
    TRANS_AT_PKS_BETA_OH_EDB = 8
    TRANS_AT_PKS_UNST = 9
    TRANS_AT_PKS_BETA_D_OH = 10
    TRANS_AT_PKS_NON_ELOGATING_BETA_OH = 11
    TRANS_AT_PKS_LACST = 12
    TRANS_AT_PKS_OUT = 13
    TRANS_AT_PKS_SHDB = 14
    TRANS_AT_PKS_BETA_OH_KETO = 15
    TRANS_AT_PKS_OXA = 16
    TRANS_AT_PKS_ACST = 17
    TRANS_AT_PKS_ALPHAME_BETAOH = 18
    TRANS_AT_PKS_NON_ELOGATING_OXA = 19
    TRANS_AT_PKS_BETA_D_OME = 20
    TRANS_AT_PKS_ARST = 21
    TRANS_AT_PKS_ALPHABETA_OH = 22
    TRANS_AT_PKS_RED_SHDB = 23
    TRANS_AT_PKS_BR = 24
    TRANS_AT_PKS_BETA_ME = 25
    TRANS_AT_PKS_EDB = 26
    TRANS_AT_PKS_BETA_MEDB = 27
    TRANS_AT_PKS_NON_ELOGATING_DB = 28
    TRANS_AT_PKS_ALPHAME_EDB = 29
    TRANS_AT_PKS_ALPHAME = 30
    TRANS_AT_PKS_ALPHAME_BETA_L_OH = 31
    TRANS_AT_PKS = 32
    TRANS_AT_PKS_AA = 33
    TRANS_AT_PKS_DB = 34
    TRANS_AT_PKS_RED  = 35
    TRANS_AT_PKS_PYR = 36
    TRANS_AT_PKS_ALPHAME_BETA_D_OH = 37
    TRANS_AT_PKS_NON_ELOGATING = 38
    TRANS_AT_PKS_MEOST = 39
    TRANS_AT_PKS_BETA_L_OH = 40
    TRANS_AT_PKS_ZDB = 41
    TRANS_AT_PKS_KETO = 42
    TRANS_AT_PKS_ALPHA_OH = 43
    TRANS_AT_PKS_NON_ELOGATING_PYR = 44

    @staticmethod
    def from_string(label: str) -> "KSDomainSubtype":
        for value in KSDomainSubtype:
            if str(value.name) == label:
                return value
        raise ValueError(f"Unknown KS domain subtype: {label}")

@unique
class ERDomainSubtype(Enum):
    """
    An Enum representing the subtypes of KS domain supported by RAIChU
    """
    S = 1
    R = 2
    UNKNOWN = 3

    @staticmethod
    def from_string(label: str) -> "ERDomainSubtype":
        for value in ERDomainSubtype:
            if str(value.name) == label:
                return value
        raise ValueError(f"Unknown ER domain subtype: {label}")
