from typing import List

from pikachu.drawing.drawing import Drawer

from raichu.reactions.chain_release import cyclic_release
from raichu.drawing.drawer import RaichuDrawer
from raichu.module import _Module
from raichu.domain.domain import TailoringDomain, CarrierDomain, SynthesisDomain, RecognitionDomain, \
    TerminationDomain, UnknownDomain, Domain, make_domain, ModuleRepresentation, DomainRepresentation, ClusterRepresentation

class Cluster:
    def __init__(self, modules: List[_Module]) -> None:
        self.modules = modules
        self.chain_intermediate = None

        self.structure_intermediates = []
        self.linear_product = None
        self.cyclised_products = []
        self.module_mechanisms = []
        self.handle_transat()

    def handle_transat(self):
        transator_clade_to_tailoring_reactions={'TRANS_AT_PKS_NON_ELOGATING_DB': ['DUMMY_KR', 'DUMMY_DH'], 'TRANS_AT_PKS_BETA_OH_KETO': [], 'TRANS_AT_PKS_NON_ELOGATING_BETA_OH': ['DUMMY_KR', 'DUMMY_DH'], 'TRANS_AT_PKS_ALPHA_OH': ['DUMMY_AH'], 'TRANS_AT_PKS_EDB': ['DUMMY_KR', 'DUMMY_EDH'], 'TRANS_AT_PKS_ARST': [], 'TRANS_AT_PKS_PYR': ['DUMMY_SC'], 'TRANS_AT_PKS_ALPHAME_EDB': ['DUMMY_AMT', 'DUMMY_KR', 'DUMMY_EDH'], 'TRANS_AT_PKS_OXI': [], 'TRANS_AT_PKS_AA': [], 'TRANS_AT_PKS_ALPHAME_BETA_L_OH': ['DUMMY_AMT', 'DUMMY_KR_B1'], 'TRANS_AT_PKS_NON_ELOGATING_BETA_L_OH': ['DUMMY_KR_B1'], 'TRANS_AT_PKS_RED_SHDB': ['DUMMY_KR', 'DUMMY_GDH'], 'TRANS_AT_PKS_DB': ['DUMMY_KR', 'DUMMY_DH'], 'TRANS_AT_PKS_KETO': [], 'TRANS_AT_PKS_NON_ELOGATING': [], 'TRANS_AT_PKS_BETA_D_OH': ['DUMMY_KR'], 'TRANS_AT_PKS_BETA_L_OH': ['DUMMY_KR_B1'], 'TRANS_AT_PKS_BR': [], 'TRANS_AT_PKS_ALPHABETA_OH': ['DUMMY_AH', 'DUMMY_KR'], 'TRANS_AT_PKS_UNST': [], 'TRANS_AT_PKS_ST': [], 'TRANS_AT_PKS_MEOST': [], 'TRANS_AT_PKS_ACST': [], 'TRANS_AT_PKS_ALPHAME': ['DUMMY_KR', 'DUMMY_DH', 'DUMMY_ER', 'DUMMY_AMT'], 'TRANS_AT_PKS_BETA_D_OME': ['DUMMY_KR', 'DUMMY_OMT'], 'TRANS_AT_PKS_NON_ELOGATING_PYR': ['DUMMY_SC'], 'TRANS_AT_PKS_BETA_ME': ['DUMMY_KR', 'DUMMY_DH', 'DUMMY_ER', 'DUMMY_BMT'], 'TRANS_AT_PKS_BETA_OH_EDB': ['DUMMY_KR'], 'TRANS_AT_PKS_NON_ELOGATING_OXA': [], 'TRANS_AT_PKS_LACST': [], 'TRANS_AT_PKS_SHDB': ['DUMMY_KR', 'DUMMY_GDH'], 'TRANS_AT_PKS_OUT': [], 'TRANS_AT_PKS_BETA_OH': ['DUMMY_KR'], 'TRANS_AT_PKS_ZDB': ['DUMMY_KR', 'DUMMY_ZDH'], 'TRANS_AT_PKS_BETA_MEDB': ['DUMMY_KR', 'DUMMY_EDH', 'DUMMY_BMT'], 'TRANS_AT_PKS': [], 'TRANS_AT_PKS_OXA': [], 'TRANS_AT_PKS_ALPHAME_BETA_D_OH': ['DUMMY_ALMT', 'DUMMY_KR_A1 '], 'TRANS_AT_PKS_ALPHAME_BETAOH': ['DUMMY_AMT', 'DUMMY_KR'], 'TRANS_AT_PKS_NON_ELOGATING_ALPHAME_EDB': ['DUMMY_AMT', 'DUMMY_KR', 'DUMMY_EDH'], 'TRANS_AT_PKS_RED': ['DUMMY_SC']}
        transator_clade_to_starter_substrate={'TRANS_AT_PKS_NON_ELOGATING_DB': 'OC(=O)CC(S)=O', 'TRANS_AT_PKS_BETA_OH_KETO': 'OC(C(O)=O)C(S)=O OC(=O)C(=O)C(S)=O', 'TRANS_AT_PKS_NON_ELOGATING_BETA_OH': 'OC(=O)CC(S)=O', 'TRANS_AT_PKS_ALPHA_OH': 'OC(O)CC(O)S', 'TRANS_AT_PKS_EDB': '[H]\\C(O)=C/C(S)=O', 'TRANS_AT_PKS_ARST': 'SC(=O)C1=CC=CC=C3', 'TRANS_AT_PKS_PYR': 'OC(=O)CC(S)=O', 'TRANS_AT_PKS_ALPHAME_EDB': 'C\\C(O)=C/C(S)=O', 'TRANS_AT_PKS_OXI': 'OC(=O)CC(S)=O', 'TRANS_AT_PKS_AA': 'C(C(=O)O)N', 'TRANS_AT_PKS_ALPHAME_BETA_L_OH': 'CC(O)[C@@H](O)C(O)S', 'TRANS_AT_PKS_NON_ELOGATING_BETA_L_OH': 'O[C@H](C(O)=O)C(S)=O', 'TRANS_AT_PKS_RED_SHDB': 'OC(=O)CC(S)=O', 'TRANS_AT_PKS_DB': '[H]\\C(O)=C/C(S)=O', 'TRANS_AT_PKS_KETO': 'OC(=O)C(=O)C(S)=O', 'TRANS_AT_PKS_NON_ELOGATING': 'OC(=O)CC(S)=O', 'TRANS_AT_PKS_BETA_D_OH': 'OC(C(O)=O)C(S)=O', 'TRANS_AT_PKS_BETA_L_OH': 'O[C@H](C(O)=O)C(S)=O', 'TRANS_AT_PKS_BR': 'OC(=O)CC(S)=O', 'TRANS_AT_PKS_ALPHABETA_OH': 'OC(O)C(O)C(O)S', 'TRANS_AT_PKS_UNST': 'OC(=O)CC(S)=O', 'TRANS_AT_PKS_ST': 'CC(S)=O', 'TRANS_AT_PKS_MEOST': 'COC(S)=O', 'TRANS_AT_PKS_ACST': 'CCC(S)=O CC(S)=O', 'TRANS_AT_PKS_ALPHAME': 'CC(O)CC(S)=O CC(O)C(O)C(S)=O CC(O)C(=O)C(S)=O', 'TRANS_AT_PKS_BETA_D_OME': 'COC(C(O)=O)C(S)=O', 'TRANS_AT_PKS_NON_ELOGATING_PYR': 'OC(=O)CC(S)=O', 'TRANS_AT_PKS_BETA_ME': 'OC(=O)C(=C)C(S)=O CC(C(O)=O)C(S)=O', 'TRANS_AT_PKS_BETA_OH_EDB': 'OC(C(O)=O)C(S)=O O\\C=C\\C(S)=O', 'TRANS_AT_PKS_NON_ELOGATING_OXA': 'CC1=NC(C)=CO1 CC1=NC(C)=CS8', 'TRANS_AT_PKS_LACST': 'C[C@@H](O)C(S)=O', 'TRANS_AT_PKS_SHDB': 'OC(=O)CC(S)=O', 'TRANS_AT_PKS_OUT': 'OC(=O)CC(S)=O', 'TRANS_AT_PKS_BETA_OH': 'OC(C(O)=O)C(S)=O', 'TRANS_AT_PKS_ZDB': '[H]\\C(O)=C\\C(S)=O', 'TRANS_AT_PKS_BETA_MEDB': '[H]\\C(O)=C(/C)C(S)=O', 'TRANS_AT_PKS': 'OC(=O)CC(S)=O', 'TRANS_AT_PKS_OXA': 'OC(=O)CC(S)=O', 'TRANS_AT_PKS_ALPHAME_BETA_D_OH': 'C[C@H](O)[C@H](O)C(O)S', 'TRANS_AT_PKS_ALPHAME_BETAOH': 'CC(O)C(O)C(O)S', 'TRANS_AT_PKS_NON_ELOGATING_ALPHAME_EDB': 'C\\C(O)=C/C(S)=O', 'TRANS_AT_PKS_RED': 'OC(=O)CC(S)=O'}
        transator_clade_to_elongating={'TRANS_AT_PKS_NON_ELOGATING_DB': 'non-elongating', 'TRANS_AT_PKS_BETA_OH_KETO': 'elongating', 'TRANS_AT_PKS_NON_ELOGATING_BETA_OH': 'non-elongating', 'TRANS_AT_PKS_ALPHA_OH': 'elongating', 'TRANS_AT_PKS_EDB': 'elongating', 'TRANS_AT_PKS_ARST': 'elongating', 'TRANS_AT_PKS_PYR': 'elongating', 'TRANS_AT_PKS_ALPHAME_EDB': 'elongating', 'TRANS_AT_PKS_OXI': 'elongating', 'TRANS_AT_PKS_AA': 'elongating', 'TRANS_AT_PKS_ALPHAME_BETA_L_OH': 'elongating', 'TRANS_AT_PKS_NON_ELOGATING_BETA_L_OH': 'non-elongating', 'TRANS_AT_PKS_RED_SHDB': 'elongating', 'TRANS_AT_PKS_DB': 'elongating', 'TRANS_AT_PKS_KETO': 'elongating', 'TRANS_AT_PKS_NON_ELOGATING': 'non-elongating', 'TRANS_AT_PKS_BETA_D_OH': 'elongating', 'TRANS_AT_PKS_BETA_L_OH': 'elongating', 'TRANS_AT_PKS_BR': 'elongating', 'TRANS_AT_PKS_ALPHABETA_OH': 'elongating', 'TRANS_AT_PKS_UNST': 'elongating', 'TRANS_AT_PKS_ST': 'elongating', 'TRANS_AT_PKS_MEOST': 'elongating', 'TRANS_AT_PKS_ACST': 'elongating', 'TRANS_AT_PKS_ALPHAME': 'elongating', 'TRANS_AT_PKS_BETA_D_OME': 'elongating', 'TRANS_AT_PKS_NON_ELOGATING_PYR': 'non-elongating', 'TRANS_AT_PKS_BETA_ME': 'elongating', 'TRANS_AT_PKS_BETA_OH_EDB': 'elongating', 'TRANS_AT_PKS_NON_ELOGATING_OXA': 'non-elongating', 'TRANS_AT_PKS_LACST': 'elongating', 'TRANS_AT_PKS_SHDB': 'elongating', 'TRANS_AT_PKS_OUT': 'elongating', 'TRANS_AT_PKS_BETA_OH': 'elongating', 'TRANS_AT_PKS_ZDB': 'elongating', 'TRANS_AT_PKS_BETA_MEDB': 'elongating', 'TRANS_AT_PKS': 'elongating', 'TRANS_AT_PKS_OXA': 'elongating', 'TRANS_AT_PKS_ALPHAME_BETA_D_OH': 'elongating', 'TRANS_AT_PKS_ALPHAME_BETAOH': 'elongating', 'TRANS_AT_PKS_NON_ELOGATING_ALPHAME_EDB': 'non-elongating', 'TRANS_AT_PKS_RED': 'elongating'}
        for index, module in enumerate(self.modules):
            if module.type.name=="PKS":
                if module.subtype.name=="PKS_TRANS":
                    if index<len(self.modules)-1:
                        next_module=self.modules[index+1]
                        if next_module.type.name=="PKS":
                            if module.is_starter_module and next_module.subtype.name=="PKS_TRANS":
                                try:
                                    self.modules[index].recognition_domain.substrate.smiles=transator_clade_to_starter_substrate[next_module.synthesis_domain.subtype.name]
                                except:
                                    self.modules[index].recognition_domain.substrate="MALONYL_COA"
                            elif module.is_starter_module and not next_module.subtype.name=="PKS_TRANS":
                                self.modules[index].recognition_domain.substrate="MALONYL_COA"
                            elif not module.is_starter_module and not module.is_termination_module and next_module.subtype.name=="PKS_TRANS":
                                self.modules[index].recognition_domain.substrate="MALONYL_COA"
                                for domain in module.tailoring_domains:
                                    domain.used=False
                                for new_domain in transator_clade_to_tailoring_reactions[next_module.synthesis_domain.subtype.name]:
                                    if "KR" in new_domain:
                                        subtype=new_domain.split("_")[-1]
                                        if subtype=="KR":
                                            subtype="A1"
                                        self.modules[index].tailoring_domains+=[make_domain(DomainRepresentation("dummy_gene", "DUMMY_KR", subtype, None, True, True),"MALONYL_COA",True)]
                                    else:
                                        self.modules[index].tailoring_domains+=[make_domain(DomainRepresentation("dummy_gene", new_domain, None, None, True, True),"MALONYL_COA",True)]    
                            elif not module.is_starter_module and not module.is_termination_module and not next_module.subtype.name=="PKS_TRANS":
                                self.modules[index].recognition_domain.substrate="MALONYL_COA"
                        else:
                                module.recognition_domain.substrate="MALONYL_COA"
                    else:
                            module.recognition_domain.substrate="MALONYL_COA"

    def compute_structures(self, compute_cyclic_products=True):
        for module in self.modules:
            structure = module.run_module(self.chain_intermediate)
            self.structure_intermediates.append(structure.deepcopy())
            self.chain_intermediate = structure
            if module.is_termination_module:
                self.linear_product = module.release_chain(structure)
                break
        else:
            raise ValueError("Cluster must contain at least one termination module.")

        if compute_cyclic_products:
            self.cyclise_all()

    def cyclise(self, atom):
        cyclic_release(self.linear_product, atom)

    def cyclise_all(self):
        pass

    def draw_spaghettis(self):
        spaghetti_svgs = []
        for structure in self.structure_intermediates:
            drawing = RaichuDrawer(structure, dont_show=True)
            drawing.draw_structure()
            svg_string = drawing.save_svg_string()
            spaghetti_svgs.append(svg_string)

        linear_drawing = Drawer(self.linear_product)
        linear_svg = linear_drawing.save_svg_string()

        return spaghetti_svgs + [linear_svg]

    def draw_cluster(self):
        pass


class Mechanism:
    def __init__(self):
        self.structures = []

    def add_structure(self, structure):
        self.structures.append(structure.deepcopy())

    def draw_mechanism(self):
        pass

