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
    TRANS_AT_PKS_KS_CLADE_1 = 3
    TRANS_AT_PKS_KS_CLADE_10 = 4
    TRANS_AT_PKS_KS_CLADE_100 = 5
    TRANS_AT_PKS_KS_CLADE_101 = 6
    TRANS_AT_PKS_KS_CLADE_102 = 7
    TRANS_AT_PKS_KS_CLADE_103 = 8
    TRANS_AT_PKS_KS_CLADE_104 = 9
    TRANS_AT_PKS_KS_CLADE_108 = 10
    TRANS_AT_PKS_KS_CLADE_109 = 11
    TRANS_AT_PKS_KS_CLADE_11 = 12
    TRANS_AT_PKS_KS_CLADE_111 = 13
    TRANS_AT_PKS_KS_CLADE_112 = 14
    TRANS_AT_PKS_KS_CLADE_113 = 15
    TRANS_AT_PKS_KS_CLADE_114 = 16
    TRANS_AT_PKS_KS_CLADE_115 = 17
    TRANS_AT_PKS_KS_CLADE_116 = 18
    TRANS_AT_PKS_KS_CLADE_118 = 19
    TRANS_AT_PKS_KS_CLADE_12 = 20
    TRANS_AT_PKS_KS_CLADE_120 = 21
    TRANS_AT_PKS_KS_CLADE_121 = 22
    TRANS_AT_PKS_KS_CLADE_122 = 23
    TRANS_AT_PKS_KS_CLADE_123 = 24
    TRANS_AT_PKS_KS_CLADE_124 = 25
    TRANS_AT_PKS_KS_CLADE_125 = 26
    TRANS_AT_PKS_KS_CLADE_126 = 27
    TRANS_AT_PKS_KS_CLADE_127 = 28
    TRANS_AT_PKS_KS_CLADE_128 = 29
    TRANS_AT_PKS_KS_CLADE_129 = 30
    TRANS_AT_PKS_KS_CLADE_13 = 31
    TRANS_AT_PKS_KS_CLADE_130 = 32
    TRANS_AT_PKS_KS_CLADE_131 = 33
    TRANS_AT_PKS_KS_CLADE_132 = 34
    TRANS_AT_PKS_KS_CLADE_133 = 35
    TRANS_AT_PKS_KS_CLADE_134 = 36
    TRANS_AT_PKS_KS_CLADE_135 = 37
    TRANS_AT_PKS_KS_CLADE_136 = 38
    TRANS_AT_PKS_KS_CLADE_137 = 39
    TRANS_AT_PKS_KS_CLADE_138 = 40
    TRANS_AT_PKS_KS_CLADE_139 = 41
    TRANS_AT_PKS_KS_CLADE_14 = 42
    TRANS_AT_PKS_KS_CLADE_140 = 43
    TRANS_AT_PKS_KS_CLADE_141 = 44
    TRANS_AT_PKS_KS_CLADE_142 = 45
    TRANS_AT_PKS_KS_CLADE_2 = 46
    TRANS_AT_PKS_KS_CLADE_21 = 47
    TRANS_AT_PKS_KS_CLADE_23 = 48
    TRANS_AT_PKS_KS_CLADE_25 = 49
    TRANS_AT_PKS_KS_CLADE_26 = 50
    TRANS_AT_PKS_KS_CLADE_27 = 51
    TRANS_AT_PKS_KS_CLADE_28 = 52
    TRANS_AT_PKS_KS_CLADE_30 = 53
    TRANS_AT_PKS_KS_CLADE_31 = 54
    TRANS_AT_PKS_KS_CLADE_32 = 55
    TRANS_AT_PKS_KS_CLADE_33 = 56
    TRANS_AT_PKS_KS_CLADE_34 = 57
    TRANS_AT_PKS_KS_CLADE_35 = 58
    TRANS_AT_PKS_KS_CLADE_36 = 59
    TRANS_AT_PKS_KS_CLADE_38 = 60
    TRANS_AT_PKS_KS_CLADE_39 = 61
    TRANS_AT_PKS_KS_CLADE_40 = 62
    TRANS_AT_PKS_KS_CLADE_42 = 63
    TRANS_AT_PKS_KS_CLADE_43 = 64
    TRANS_AT_PKS_KS_CLADE_44 = 65
    TRANS_AT_PKS_KS_CLADE_45 = 66
    TRANS_AT_PKS_KS_CLADE_46 = 67
    TRANS_AT_PKS_KS_CLADE_49 = 68
    TRANS_AT_PKS_KS_CLADE_5 = 69
    TRANS_AT_PKS_KS_CLADE_51 = 70
    TRANS_AT_PKS_KS_CLADE_52 = 71
    TRANS_AT_PKS_KS_CLADE_53 = 72
    TRANS_AT_PKS_KS_CLADE_55 = 73
    TRANS_AT_PKS_KS_CLADE_56 = 74
    TRANS_AT_PKS_KS_CLADE_57 = 75
    TRANS_AT_PKS_KS_CLADE_60 = 76
    TRANS_AT_PKS_KS_CLADE_61 = 77
    TRANS_AT_PKS_KS_CLADE_62 = 78
    TRANS_AT_PKS_KS_CLADE_64 = 79
    TRANS_AT_PKS_KS_CLADE_65 = 80
    TRANS_AT_PKS_KS_CLADE_66 = 81
    TRANS_AT_PKS_KS_CLADE_67 = 82
    TRANS_AT_PKS_KS_CLADE_68 = 83
    TRANS_AT_PKS_KS_CLADE_7 = 84
    TRANS_AT_PKS_KS_CLADE_70 = 85
    TRANS_AT_PKS_KS_CLADE_73 = 86
    TRANS_AT_PKS_KS_CLADE_74 = 87
    TRANS_AT_PKS_KS_CLADE_75 = 88
    TRANS_AT_PKS_KS_CLADE_76 = 89
    TRANS_AT_PKS_KS_CLADE_78 = 90
    TRANS_AT_PKS_KS_CLADE_79 = 91
    TRANS_AT_PKS_KS_CLADE_8 = 92
    TRANS_AT_PKS_KS_CLADE_80 = 93
    TRANS_AT_PKS_KS_CLADE_81 = 94
    TRANS_AT_PKS_KS_CLADE_82 = 95
    TRANS_AT_PKS_KS_CLADE_83 = 96
    TRANS_AT_PKS_KS_CLADE_84 = 97
    TRANS_AT_PKS_KS_CLADE_85 = 98
    TRANS_AT_PKS_KS_CLADE_86 = 99
    TRANS_AT_PKS_KS_CLADE_88 = 100
    TRANS_AT_PKS_KS_CLADE_89 = 101
    TRANS_AT_PKS_KS_CLADE_9 = 102
    TRANS_AT_PKS_KS_CLADE_90 = 103
    TRANS_AT_PKS_KS_CLADE_92 = 104
    TRANS_AT_PKS_KS_CLADE_93 = 105
    TRANS_AT_PKS_KS_CLADE_94 = 106
    TRANS_AT_PKS_KS_CLADE_95 = 107
    TRANS_AT_PKS_KS_CLADE_96 = 108
    TRANS_AT_PKS_KS_CLADE_97 = 109
    TRANS_AT_PKS_KS_CLADE_98 = 110
    TRANS_AT_PKS_KS_CLADE_99 = 111


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
