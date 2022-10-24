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
    DUMMY_BMT = 28 #Beta-Methyltransferase
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
    TRANS_AT_PKS = 3
    TRANS_AT_PKS_AA = 4
    TRANS_AT_PKS_AA = 5
    TRANS_AT_PKS_AA = 6
    TRANS_AT_PKS_AA = 7
    TRANS_AT_PKS_ACST = 8
    TRANS_AT_PKS_ACST = 9
    TRANS_AT_PKS_ALPHA_OH = 10
    TRANS_AT_PKS_ALPHABETA_OH = 11
    TRANS_AT_PKS_ALPHABETA_OH = 12
    TRANS_AT_PKS_ALPHAME = 13
    TRANS_AT_PKS_ALPHAME = 14
    TRANS_AT_PKS_ALPHAME = 15
    TRANS_AT_PKS_ALPHAME = 16
    TRANS_AT_PKS_ALPHAME = 17
    TRANS_AT_PKS_ALPHAME_BETA_D_OH = 18
    TRANS_AT_PKS_ALPHAME_BETA_L_OH = 19
    TRANS_AT_PKS_ALPHAME_BETA_L_OH = 20
    TRANS_AT_PKS_ALPHAME_BETA_L_OH = 21
    TRANS_AT_PKS_ALPHAME_BETA_L_OH = 22
    TRANS_AT_PKS_ALPHAME_BETAOH = 23
    TRANS_AT_PKS_ALPHAME_EDB = 24
    TRANS_AT_PKS_ALPHAME_EDB = 25
    TRANS_AT_PKS_ALPHAME_EDB = 26
    TRANS_AT_PKS_ALPHAME_EDB = 27
    TRANS_AT_PKS_ALPHAME_EDB = 28
    TRANS_AT_PKS_ALPHAME_EDB = 29
    TRANS_AT_PKS_ALPHAME_EDB = 30
    TRANS_AT_PKS_ALPHAME_EDB = 31
    TRANS_AT_PKS_ALPHAME_EDB = 32
    TRANS_AT_PKS_ALPHAME_EDB = 33
    TRANS_AT_PKS_ARST = 34
    TRANS_AT_PKS_BETA_D_OH = 35
    TRANS_AT_PKS_BETA_D_OH = 36
    TRANS_AT_PKS_BETA_D_OH = 37
    TRANS_AT_PKS_BETA_D_OME = 38
    TRANS_AT_PKS_BETA_L_OH = 39
    TRANS_AT_PKS_BETA_L_OH = 40
    TRANS_AT_PKS_BETA_L_OH = 41
    TRANS_AT_PKS_BETA_L_OH = 42
    TRANS_AT_PKS_BETA_ME = 43
    TRANS_AT_PKS_BETA_ME = 44
    TRANS_AT_PKS_BETA_MEDB = 45
    TRANS_AT_PKS_BETA_OH = 46
    TRANS_AT_PKS_BETA_OH/EDB = 47
    TRANS_AT_PKS_BETA_OH/EDB = 48
    TRANS_AT_PKS_BETA_OH/KETO = 49
    TRANS_AT_PKS_BR = 50
    TRANS_AT_PKS_DB = 51
    TRANS_AT_PKS_DB = 52
    TRANS_AT_PKS_DB = 53
    TRANS_AT_PKS_DB = 54
    TRANS_AT_PKS_DB = 55
    TRANS_AT_PKS_DB = 56
    TRANS_AT_PKS_EDB = 57
    TRANS_AT_PKS_EDB = 58
    TRANS_AT_PKS_EDB = 59
    TRANS_AT_PKS_EDB = 60
    TRANS_AT_PKS_EDB = 61
    TRANS_AT_PKS_EDB = 62
    TRANS_AT_PKS_EDB = 63
    TRANS_AT_PKS_EDB = 64
    TRANS_AT_PKS_KETO = 65
    TRANS_AT_PKS_KETO = 66
    TRANS_AT_PKS_KETO = 67
    TRANS_AT_PKS_LACST = 68
    TRANS_AT_PKS_MEOST = 69
    TRANS_AT_PKS_NON_ELOGATING = 70
    TRANS_AT_PKS_NON_ELOGATING = 71
    TRANS_AT_PKS_NON_ELOGATING_ALPHAME_EDB = 72
    TRANS_AT_PKS_NON_ELOGATING_BETA_L_OH = 73
    TRANS_AT_PKS_NON_ELOGATING_BETA_L_OH = 74
    TRANS_AT_PKS_NON_ELOGATING_BETA_OH = 75
    TRANS_AT_PKS_NON_ELOGATING_BETA_OH = 76
    TRANS_AT_PKS_NON_ELOGATING_BETA_OH = 77
    TRANS_AT_PKS_NON_ELOGATING_BETA_OH = 78
    TRANS_AT_PKS_NON_ELOGATING_BETA_OH = 79
    TRANS_AT_PKS_NON_ELOGATING_BETA_OH = 80
    TRANS_AT_PKS_NON_ELOGATING_BETA_OH = 81
    TRANS_AT_PKS_NON_ELOGATING_BETA_OH = 82
    TRANS_AT_PKS_NON_ELOGATING_BETA_OH = 83
    TRANS_AT_PKS_NON_ELOGATING_BETA_OH = 84
    TRANS_AT_PKS_NON_ELOGATING_BETA_OH = 85
    TRANS_AT_PKS_NON_ELOGATING_DB = 86
    TRANS_AT_PKS_NON_ELOGATING_DB = 87
    TRANS_AT_PKS_NON_ELOGATING_DB = 88
    TRANS_AT_PKS_NON_ELOGATING_DB = 89
    TRANS_AT_PKS_NON_ELOGATING_OXA = 90
    TRANS_AT_PKS_NON_ELOGATING_OXA = 91
    TRANS_AT_PKS_NON_ELOGATING_PYR = 92
    TRANS_AT_PKS_OUT = 93
    TRANS_AT_PKS_OXA = 94
    TRANS_AT_PKS_OXI = 95
    TRANS_AT_PKS_PYR = 96
    TRANS_AT_PKS_PYR = 97
    TRANS_AT_PKS_RED  = 98
    TRANS_AT_PKS_RED_SHDB = 99
    TRANS_AT_PKS_RED_SHDB = 100
    TRANS_AT_PKS_RED_SHDB = 101
    TRANS_AT_PKS_RED_SHDB = 102
    TRANS_AT_PKS_RED_SHDB = 103
    TRANS_AT_PKS_SHDB = 104
    TRANS_AT_PKS_ST = 105
    TRANS_AT_PKS_ST = 106
    TRANS_AT_PKS_ST = 107
    TRANS_AT_PKS_UNST = 108
    TRANS_AT_PKS_ZDB = 109
    TRANS_AT_PKS_ZDB = 110
    TRANS_AT_PKS_ZDB = 111
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
