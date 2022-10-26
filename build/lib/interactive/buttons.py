import pygame
import os

from interactive.style import BUTTON_TEXT_COLOUR, FONT, BUTTON_HIGHLIGHT_COLOUR, BLACK, BUTTON_COLOUR, \
    BUTTON_PANEL_COLOUR, HEIGHT, WIDTH, DOMAIN_BUTTON_SIZE, SUBSTRATE_GROUP_BUTTONS_PER_LINE, \
    SUBSTRATE_GROUP_BUTTON_PADDING, SUBSTRATE_GROUP_BUTTON_SIZE, SUBSTRATE_BUTTON_SIZE, KR_BUTTON_SIZE
import interactive.images.domains
import interactive.flatfiles
import interactive.images.kr_subtypes
from interactive.gene import Gene
from interactive.textbox import TextBox
from interactive.parsers import parse_smiles
from interactive.substrate import Substrate, SubstrateGroup, \
    PROTEINOGENIC_SUBSTRATES, NONPROTEINOGENIC_SUBSTRATES, FATTY_ACIDS, NON_AMINO_ACIDS

DOMAIN_IMAGE_DIR = os.path.dirname(interactive.images.domains.__file__)
KR_IMAGE_DIR = os.path.dirname(interactive.images.kr_subtypes.__file__)
FLATFILES = os.path.dirname(interactive.flatfiles.__file__)
PARAS_SMILES = os.path.join(FLATFILES, "PARAS_smiles.txt")
PKS_SMILES = os.path.join(FLATFILES, "at_specificities.txt")


class Button:
    def __init__(self, text, position, dimensions, font_size=None):

        self.text = text

        self.x, self.y = position
        self.width, self.height = dimensions
        self.rectangle = pygame.Rect(position, dimensions)
        self.font_size = int(self.height) - 10

        self.text_x = 0
        self.text_y = 0

        self.set_text_position()
        self.font = pygame.font.SysFont(FONT, self.font_size, bold=False)

        self.rendered_text = self.font.render(text, True, BUTTON_TEXT_COLOUR)
        self.selected = False

    def __hash__(self):
        return hash(self.text)

    def __eq__(self, other):
        return self.text == other.text

    def __repr__(self):
        return self.text

    def highlight(self, screen):
        pygame.draw.rect(screen, BUTTON_HIGHLIGHT_COLOUR, self.rectangle)
        pygame.draw.rect(screen, BLACK, self.rectangle, 2)
        screen.blit(self.rendered_text, (self.text_x, self.text_y))

    def draw(self, screen):
        pygame.draw.rect(screen, BUTTON_COLOUR, self.rectangle)
        pygame.draw.rect(screen, BLACK, self.rectangle, 2)
        screen.blit(self.rendered_text, (self.text_x, self.text_y))

    def erase_button(self, screen):
        pygame.draw.rect(screen, BUTTON_PANEL_COLOUR, self.rectangle)
        pygame.draw.rect(screen, BUTTON_PANEL_COLOUR, self.rectangle, 2)

    def set_text_position(self):
        self.text_x = self.x + 5
        self.text_y = self.y + int((self.height - self.font_size) / 2)

    def set_font(self):
        self.font = pygame.font.SysFont(FONT, self.font_size, bold=True)


class RenderClusterButton(Button):

    def __init__(self):
        position = (int(0.77 * WIDTH), int(0.92 * HEIGHT))
        dimensions = (int(0.2 * WIDTH), int(HEIGHT / 25))

        super().__init__("Render cluster", position, dimensions)

    def do_action(self, screen, active_buttons):
        position = (int(0.77 * WIDTH), int(0.92 * HEIGHT))
        dimensions = (int(0.2 * WIDTH), int(HEIGHT / 25))
        text_box = TextBox(position, dimensions)
        text_box.draw(screen)
        hide_buttons(ALL_BUTTONS, screen, active_buttons)
        show_button(SAVE_CLUSTER_BUTTON, screen, active_buttons)
        return text_box


class KRSubtypeButton(Button):
    def __init__(self, position, kr_subtype):
        self.kr_subtype = kr_subtype
        dimensions = (KR_BUTTON_SIZE, KR_BUTTON_SIZE)

        super().__init__(kr_subtype, position, dimensions)

        if kr_subtype:

            self.image = os.path.join(KR_IMAGE_DIR, f"KR_{kr_subtype}.png")
            self.highlight_image = os.path.join(KR_IMAGE_DIR, f"KR_{kr_subtype}_highlight.png")
        else:
            self.image = os.path.join(KR_IMAGE_DIR, f"KR.png")
            self.highlight_image = os.path.join(KR_IMAGE_DIR, f"KR_highlight.png")


    def draw(self, screen):
        kr_image = pygame.image.load(self.image)
        kr_image_scaled = pygame.transform.smoothscale(kr_image, (KR_BUTTON_SIZE, KR_BUTTON_SIZE))
        screen.blit(kr_image_scaled, self.rectangle)

    def highlight(self, screen):
        kr_image = pygame.image.load(self.highlight_image)
        kr_image_scaled = pygame.transform.smoothscale(kr_image, (KR_BUTTON_SIZE, KR_BUTTON_SIZE))
        screen.blit(kr_image_scaled, self.rectangle)

    def do_action(self, domain, mouse, screen, active_buttons):
        if self.kr_subtype:
            domain.subtype = self.kr_subtype
        else:
            domain.subtype = None

        domain.module.gene.erase()
        domain.module.gene.draw(mouse)
        reset_buttons(screen, active_buttons)


class DomainButton(Button):
    def __init__(self, position, text):
        dimensions = (DOMAIN_BUTTON_SIZE, DOMAIN_BUTTON_SIZE)
        super().__init__(text, position, dimensions)
        self.image = None
        self.highlight_image = None
        self.domain_type = text

        self.image = os.path.join(DOMAIN_IMAGE_DIR, f'{text.lower()}.png')
        self.highlight_image = os.path.join(DOMAIN_IMAGE_DIR, f'{text.lower()}_highlight.png')

    def draw(self, screen):
        domain_image = pygame.image.load(self.image)
        domain_image_scaled = pygame.transform.smoothscale(domain_image, (DOMAIN_BUTTON_SIZE, DOMAIN_BUTTON_SIZE))
        screen.blit(domain_image_scaled, self.rectangle)

    def highlight(self, screen):
        domain_image = pygame.image.load(self.highlight_image)
        domain_image_scaled = pygame.transform.smoothscale(domain_image, (DOMAIN_BUTTON_SIZE, DOMAIN_BUTTON_SIZE))
        screen.blit(domain_image_scaled, self.rectangle)

    def do_action(self, module, mouse, screen, active_buttons):
        module.gene.erase()
        module.add_domain(self.domain_type)
        module.gene.draw(mouse)
        reset_buttons(screen, active_buttons)


class AddGeneButton(Button):
    def __init__(self):
        position = (int(0.02 * WIDTH), int(0.57 * HEIGHT))
        dimensions = (int(0.2 * WIDTH), int(HEIGHT / 25))

        super().__init__("Add gene", position, dimensions)

    def do_action(self, screen, genes, text_box, active_buttons, mouse):
        gene_number = len(genes)
        gene = Gene(screen, gene_number)
        gene.name = text_box.text
        text_box.erase(screen)
        gene.draw(mouse)
        genes.append(gene)
        reset_buttons(screen, active_buttons)


class RemoveModuleButton(Button):

    def __init__(self):
        position = (int(0.02 * WIDTH), int(0.57 * HEIGHT))
        dimensions = (int(0.2 * WIDTH), int(HEIGHT / 25))

        super().__init__("Remove module", position, dimensions)

    def do_action(self, screen, current_module, active_buttons, mouse):
        gene = current_module.gene
        gene.remove_module(current_module.id)

        gene.erase()
        gene.draw(mouse)

        reset_buttons(screen, active_buttons)


class RemoveGeneButton(Button):
    def __init__(self):
        position = (int(0.02 * WIDTH), int(0.57 * HEIGHT))
        dimensions = (int(0.2 * WIDTH), int(HEIGHT / 25))

        super().__init__("Remove gene", position, dimensions)

    def do_action(self, screen, genes, current_gene, active_buttons, mouse):
        gene_number = current_gene.gene_number
        index_to_remove = None
        for i, gene in enumerate(genes):
            if gene.gene_number == gene_number:
                index_to_remove = i
            elif gene.gene_number > gene_number:
                gene.gene_number -= 1

        genes.pop(index_to_remove)
        current_gene.erase()

        for gene in genes:
            gene.erase()
            gene.draw(mouse)

        reset_buttons(screen, active_buttons)


class CreateGeneButton(Button):
    def __init__(self):
        position = (int(0.02 * WIDTH), int(0.52 * HEIGHT))
        dimensions = (int(0.2 * WIDTH), int(HEIGHT / 25))

        super().__init__("Create gene", position, dimensions)

    def do_action(self, screen, active_buttons):
        position = (int(0.02 * WIDTH), int(0.52 * HEIGHT))
        dimensions = (int(0.2 * WIDTH), int(HEIGHT / 25))
        text_box = TextBox(position, dimensions)
        text_box.draw(screen)
        hide_buttons(ALL_BUTTONS, screen, active_buttons)
        show_button(ADD_GENE_BUTTON, screen, active_buttons)
        return text_box


class SubstrateGroupButton(Button):
    def __init__(self, text, position, dimensions):
        super().__init__(text, position, dimensions)


class AddModuleButton(Button):
    def __init__(self):
        position = (int(0.02 * WIDTH), int(0.62 * HEIGHT))
        dimensions = (int(0.2 * WIDTH), int(HEIGHT / 25))

        super().__init__("Add module", position, dimensions)


class NRPSModuleButton(Button):
    def __init__(self):
        position = (int(0.02 * WIDTH), int(0.67 * HEIGHT))
        dimensions = (int(0.09 * WIDTH), int(HEIGHT / 25))

        super().__init__('NRPS', position, dimensions)

    def do_action(self, screen, current_gene, active_buttons, insertion_point, mouse):
        module_nr = insertion_point.insertion_point

        current_gene.add_module(module_nr, 'NRPS')
        current_gene.erase()
        current_gene.insertion_point = None
        current_gene.draw(mouse)
        reset_buttons(screen, active_buttons)


class PKSModuleButton(Button):
    def __init__(self):
        position = (int(0.13 * WIDTH), int(0.67 * HEIGHT))
        dimensions = (int(0.09 * WIDTH), int(HEIGHT / 25))

        super().__init__('PKS', position, dimensions)

    def do_action(self, screen, current_gene, active_buttons, insertion_point, mouse):
        module_nr = insertion_point.insertion_point

        current_gene.add_module(module_nr, 'PKS')
        current_gene.erase()
        current_gene.insertion_point = None
        current_gene.draw(mouse)
        reset_buttons(screen, active_buttons)


class AddDomainButton(Button):
    def __init__(self):
        position = (int(0.02 * WIDTH), int(0.62 * HEIGHT))
        dimensions = (int(0.2 * WIDTH), int(HEIGHT / 25))

        super().__init__("Add domain", position, dimensions)


class ACPDomainButton(DomainButton):
    def __init__(self):
        position = (int(0.35 * WIDTH), int(0.52 * HEIGHT))

        super().__init__(position, 'ACP')


class PCPDomainButton(DomainButton):
    def __init__(self):
        position = (int(0.35 * WIDTH), int(0.52 * HEIGHT))
        super().__init__(position, 'PCP')


class ADomainButton(DomainButton):
    def __init__(self):
        position = (int(0.30 * WIDTH), int(0.52 * HEIGHT))
        super().__init__(position, 'A')


class CDomainButton(DomainButton):
    def __init__(self):
        position = (int(0.25 * WIDTH), int(0.52 * HEIGHT))
        super().__init__(position, 'C')


class KSDomainButton(DomainButton):
    def __init__(self):
        position = (int(0.25 * WIDTH), int(0.52 * HEIGHT))
        super().__init__(position, 'KS')


class KRDomainButton(DomainButton):
    def __init__(self):
        position = (int(0.25 * WIDTH), int(0.62 * HEIGHT))
        super().__init__(position, 'KR')


class DHDomainButton(DomainButton):
    def __init__(self):
        position = (int(0.30 * WIDTH), int(0.62 * HEIGHT))
        super().__init__(position, 'DH')


class ERDomainButton(DomainButton):
    def __init__(self):
        position = (int(0.35 * WIDTH), int(0.62 * HEIGHT))
        super().__init__(position, 'ER')


class EDomainButton(DomainButton):
    def __init__(self):
        position = (int(0.25 * WIDTH), int(0.62 * HEIGHT))
        super().__init__(position, 'E')


class NMTDomainButton(DomainButton):
    def __init__(self):
        position = (int(0.30 * WIDTH), int(0.62 * HEIGHT))
        super().__init__(position, 'nMT')


class ATDomainButton(DomainButton):
    def __init__(self):
        position = (int(0.30 * WIDTH), int(0.52 * HEIGHT))
        super().__init__(position, 'AT')


class TEDomainButton(DomainButton):
    def __init__(self):
        position = (int(0.25 * WIDTH), int(0.72 * HEIGHT))
        super().__init__(position, 'TE')


class RemoveDomainButton(Button):
    def __init__(self):
        position = (int(0.02 * WIDTH), int(0.57 * HEIGHT))
        dimensions = (int(0.2 * WIDTH), int(HEIGHT / 25))

        super().__init__("Remove domain", position, dimensions)

    def do_action(self, current_domain, screen, active_buttons, mouse):
        gene = current_domain.module.gene
        current_domain.module.remove_domain(current_domain)
        gene.erase()
        gene.draw(mouse)
        reset_buttons(screen, active_buttons)


class SelectSubstrateButton(Button):
    def __init__(self):
        position = (int(0.75 * WIDTH), int(0.52 * HEIGHT))
        dimensions = (int(0.2 * WIDTH), int(HEIGHT / 25))

        super().__init__("Select substrate", position, dimensions)


class SelectDomainTypeButton(Button):
    def __init__(self):
        position = (int(0.75 * WIDTH), int(0.57 * HEIGHT))
        dimensions = (int(0.2 * WIDTH), int(HEIGHT / 25))

        super().__init__("Select subtype", position, dimensions)


class SetDomainInactiveButton(Button):

    def __init__(self):
        position = (int(0.75 * WIDTH), int(0.62 * HEIGHT))
        dimensions = (int(0.2 * WIDTH), int(HEIGHT / 25))

        super().__init__("Set to inactive", position, dimensions)


class SubstrateButton(Button):
    def __init__(self, substrate, position, dimensions):

        super().__init__(substrate.name, position, dimensions)
        self.substrate = substrate
        self.image = substrate.image
        self.highlight_image = substrate.highlight_image

    def draw(self, screen):
        substrate_image = pygame.image.load(self.image)
        substrate_image_scaled = pygame.transform.smoothscale(substrate_image, (SUBSTRATE_BUTTON_SIZE,
                                                                                SUBSTRATE_BUTTON_SIZE))
        screen.blit(substrate_image_scaled, self.rectangle)

    def highlight(self, screen):
        substrate_image = pygame.image.load(self.highlight_image)
        substrate_image_scaled = pygame.transform.smoothscale(substrate_image, (SUBSTRATE_BUTTON_SIZE,
                                                                                SUBSTRATE_BUTTON_SIZE))
        screen.blit(substrate_image_scaled, self.rectangle)

    def do_action(self):
        pass


class NRPSWildcardButton(SubstrateButton):
    def __init__(self):
        substrate = Substrate('nrp', 'NC(C(=O)O)[*]', 'NRPS')
        position = (int(0.52 * WIDTH), int(0.52 * HEIGHT))
        substrate.x = position[0]
        substrate.y = position[1]
        substrate.rectangle = pygame.Rect(substrate.x, substrate.y, substrate.width, substrate.height)
        dimensions = (substrate.width, substrate.height)

        super().__init__(substrate, position, dimensions)


class SubstrateSupergroupButton(Button):
    def __init__(self, text, position, group):
        dimensions = (int(0.2 * WIDTH), int(HEIGHT / 25))
        self.group = group

        super().__init__(text, position, dimensions)


class ProteinogenicButton(SubstrateSupergroupButton):
    def __init__(self):
        position = (int(0.25 * WIDTH), int(0.52 * HEIGHT))

        super().__init__("Proteinogenic", position,
                         PROTEINOGENIC_SUBSTRATES)


class NonProteinogenicButton(SubstrateSupergroupButton):
    def __init__(self):
        position = (int(0.25 * WIDTH), int(0.57 * HEIGHT))

        super().__init__("Non-proteinogenic", position,
                         NONPROTEINOGENIC_SUBSTRATES)


class NonAminoAcidButton(SubstrateSupergroupButton):
    def __init__(self):
        position = (int(0.25 * WIDTH), int(0.52 * HEIGHT))

        super().__init__("Non-amino acids", position,
                         NON_AMINO_ACIDS)


class FattyAcidButton(SubstrateSupergroupButton):
    def __init__(self):
        position = (int(0.25 * WIDTH), int(0.57 * HEIGHT))

        super().__init__("Fatty acids", position,
                         FATTY_ACIDS)


class RenderProductsButton(Button):
    def __init__(self):
        position = (int(0.52 * WIDTH), int(0.92 * HEIGHT))
        dimensions = (int(0.2 * WIDTH), int(HEIGHT / 25))

        super().__init__("Render products", position, dimensions)

    def do_action(self, screen, active_buttons):
        position = (int(0.02 * WIDTH), int(0.92 * HEIGHT))
        dimensions = (int(0.2 * WIDTH), int(HEIGHT / 25))
        text_box = TextBox(position, dimensions)
        text_box.draw(screen)
        hide_buttons(ALL_BUTTONS, screen, active_buttons)
        show_button(SAVE_PRODUCTS_BUTTON, screen, active_buttons)
        return text_box


class YesButton(Button):
    def __init__(self):
        position = (int(0.02 * WIDTH), int(0.92 * HEIGHT))
        dimensions = (int(0.05 * WIDTH), int(HEIGHT / 25))

        super().__init__("Yes", position, dimensions)


class NoButton(Button):
    def __init__(self):
        position = (int(0.08 * WIDTH), int(0.92 * HEIGHT))
        dimensions = (int(0.05 * WIDTH), int(HEIGHT / 25))

        super().__init__("No", position, dimensions)


class FattyAcidSuperOptionButton(Button):
    def __init__(self, text, position, dimensions):
        super().__init__(text, position, dimensions)


class FattyAcidOptionButton(Button):
    def __init__(self, text, position, dimensions):
        super().__init__(text, position, dimensions)


class SaveClusterButton(Button):
    def __init__(self):
        position = (int(0.02 * WIDTH), int(0.82 * HEIGHT))
        dimensions = (int(0.2 * WIDTH), int(HEIGHT / 25))

        super().__init__("Save cluster", position, dimensions)

    def do_action(self, screen, active_buttons):
        position = (int(0.02 * WIDTH), int(0.87 * HEIGHT))
        dimensions = (int(0.2 * WIDTH), int(HEIGHT / 25))
        text_box = TextBox(position, dimensions)
        text_box.draw(screen)
        hide_buttons(ALL_BUTTONS, screen, active_buttons)
        show_button(SAVE_TO_PNG_BUTTON, screen, active_buttons)
        return text_box


class ExportTabularButton(Button):
    def __init__(self):
        position = (int(0.27 * WIDTH), int(0.92 * HEIGHT))
        dimensions = (int(0.2 * WIDTH), int(HEIGHT / 25))

        super().__init__("Export to tabular", position, dimensions)

    def do_action(self, screen, active_buttons):
        position = (int(0.02 * WIDTH), int(0.87 * HEIGHT))
        dimensions = (int(0.2 * WIDTH), int(HEIGHT / 25))
        text_box = TextBox(position, dimensions)
        text_box.draw(screen)
        hide_buttons(ALL_BUTTONS, screen, active_buttons)
        show_button(SAVE_TO_TXT_BUTTON, screen, active_buttons)
        return text_box


class SaveToTxtButton(Button):
    def __init__(self):
        position = (int(0.02 * WIDTH), int(0.82 * HEIGHT))
        dimensions = (int(0.2 * WIDTH), int(HEIGHT / 25))

        super().__init__("Save to tabular", position, dimensions)


class SaveToPngButton(Button):
    def __init__(self):
        position = (int(0.02 * WIDTH), int(0.82 * HEIGHT))
        dimensions = (int(0.2 * WIDTH), int(HEIGHT / 25))

        super().__init__("Save to png", position, dimensions)


class CButton(Button):
    def __init__(self, nr, position, dimensions):
        self.nr = nr
        super().__init__(f"C{nr}", position, dimensions)


class SaveToFolderButton(Button):
    def __init__(self):
        position = (int(0.02 * WIDTH), int(0.82 * HEIGHT))
        dimensions = (int(0.2 * WIDTH), int(HEIGHT / 25))

        super().__init__("Save to folder", position, dimensions)


class SaveProductsButton(Button):
    def __init__(self):
        position = (int(0.02 * WIDTH), int(0.82 * HEIGHT))
        dimensions = (int(0.2 * WIDTH), int(HEIGHT / 25))

        super().__init__("Save products", position, dimensions)

    def do_action(self, screen, active_buttons):
        position = (int(0.02 * WIDTH), int(0.87 * HEIGHT))
        dimensions = (int(0.2 * WIDTH), int(HEIGHT / 25))
        text_box = TextBox(position, dimensions)
        text_box.draw(screen)
        hide_buttons(ALL_BUTTONS, screen, active_buttons)
        show_button(SAVE_TO_FOLDER_BUTTON, screen, active_buttons)
        return text_box


class StarterButton(Button):
    def __init__(self):
        position = (int(0.25 * WIDTH), int(0.52 * HEIGHT))
        dimensions = (int(0.2 * WIDTH), int(HEIGHT / 25))

        super().__init__("Starter substrates", position,
                         dimensions)


class ElongationButton(Button):
    def __init__(self):
        position = (int(0.25 * WIDTH), int(0.57 * HEIGHT))
        dimensions = (int(0.2 * WIDTH), int(HEIGHT / 25))

        super().__init__("Elongation substrates", position,
                         dimensions)


RENDER_CLUSTER_BUTTON = RenderClusterButton()
RENDER_PRODUCTS_BUTTON = RenderProductsButton()
SAVE_CLUSTER_BUTTON = SaveClusterButton()
SAVE_PRODUCTS_BUTTON = SaveProductsButton()
SAVE_TO_FOLDER_BUTTON = SaveToFolderButton()
SAVE_TO_PNG_BUTTON = SaveToPngButton()
SAVE_TO_TXT_BUTTON = SaveToTxtButton()
EXPORT_TABULAR_BUTTON = ExportTabularButton()

YES_BUTTON = YesButton()
NO_BUTTON = NoButton()

CREATE_GENE_BUTTON = CreateGeneButton()
ADD_GENE_BUTTON = AddGeneButton()
REMOVE_GENE_BUTTON = RemoveGeneButton()

ADD_MODULE_BUTTON = AddModuleButton()
REMOVE_MODULE_BUTTON = RemoveModuleButton()

NRPS_MODULE_BUTTON = NRPSModuleButton()
PKS_MODULE_BUTTON = PKSModuleButton()

ADD_DOMAIN_BUTTON = AddDomainButton()
REMOVE_DOMAIN_BUTTON = RemoveDomainButton()

SELECT_SUBSTRATE_BUTTON = SelectSubstrateButton()
SELECT_DOMAIN_TYPE_BUTTON = SelectDomainTypeButton()
SET_DOMAIN_INACTIVE_BUTTON = SetDomainInactiveButton()


KS_DOMAIN_BUTTON = KSDomainButton()
AT_DOMAIN_BUTTON = ATDomainButton()
ACP_DOMAIN_BUTTON = ACPDomainButton()
DH_DOMAIN_BUTTON = DHDomainButton()
ER_DOMAIN_BUTTON = ERDomainButton()
KR_DOMAIN_BUTTON = KRDomainButton()

C_DOMAIN_BUTTON = CDomainButton()
A_DOMAIN_BUTTON = ADomainButton()
PCP_DOMAIN_BUTTON = PCPDomainButton()
E_DOMAIN_BUTTON = EDomainButton()
NMT_DOMAIN_BUTTON = NMTDomainButton()
TE_DOMAIN_BUTTON = TEDomainButton()

STARTER_BUTTON = StarterButton()
ELONGATION_BUTTON = ElongationButton()

PKS_DOMAIN_TO_BUTTON = {'KS': KS_DOMAIN_BUTTON,
                        'AT': AT_DOMAIN_BUTTON,
                        'DH': DH_DOMAIN_BUTTON,
                        'ER': ER_DOMAIN_BUTTON,
                        'KR': KR_DOMAIN_BUTTON,
                        'ACP': ACP_DOMAIN_BUTTON,
                        'TE': TE_DOMAIN_BUTTON}

NRPS_DOMAIN_TO_BUTTON = {'C': C_DOMAIN_BUTTON,
                         'A': A_DOMAIN_BUTTON,
                         'nMT': NMT_DOMAIN_BUTTON,
                         'PCP': PCP_DOMAIN_BUTTON,
                         'E': E_DOMAIN_BUTTON,
                         'TE': TE_DOMAIN_BUTTON}

PROTEINOGENIC_BUTTON = ProteinogenicButton()
NON_PROTEINOGENIC_BUTTON = NonProteinogenicButton()
NON_AMINO_ACID_BUTTON = NonAminoAcidButton()
FATTY_ACID_BUTTON = FattyAcidButton()
NRPS_WILDCARD_BUTTON = NRPSWildcardButton()

NRPS_SUPERGROUP_BUTTONS = [PROTEINOGENIC_BUTTON,
                           NON_PROTEINOGENIC_BUTTON,
                           NON_AMINO_ACID_BUTTON,
                           FATTY_ACID_BUTTON,
                           NRPS_WILDCARD_BUTTON]

ALL_BUTTONS = [RENDER_CLUSTER_BUTTON,
               RENDER_PRODUCTS_BUTTON,
               CREATE_GENE_BUTTON,
               ADD_GENE_BUTTON,
               ADD_MODULE_BUTTON,
               NRPS_MODULE_BUTTON,
               PKS_MODULE_BUTTON,
               KS_DOMAIN_BUTTON,
               AT_DOMAIN_BUTTON,
               ACP_DOMAIN_BUTTON,
               DH_DOMAIN_BUTTON,
               ER_DOMAIN_BUTTON,
               KR_DOMAIN_BUTTON,
               C_DOMAIN_BUTTON,
               A_DOMAIN_BUTTON,
               PCP_DOMAIN_BUTTON,
               E_DOMAIN_BUTTON,
               TE_DOMAIN_BUTTON,
               NMT_DOMAIN_BUTTON,
               ADD_DOMAIN_BUTTON,
               REMOVE_DOMAIN_BUTTON,
               REMOVE_MODULE_BUTTON,
               REMOVE_GENE_BUTTON,
               SELECT_SUBSTRATE_BUTTON,
               SELECT_DOMAIN_TYPE_BUTTON,
               SET_DOMAIN_INACTIVE_BUTTON,
               PROTEINOGENIC_BUTTON,
               NON_PROTEINOGENIC_BUTTON,
               NON_AMINO_ACID_BUTTON,
               FATTY_ACID_BUTTON,
               SAVE_CLUSTER_BUTTON,
               SAVE_PRODUCTS_BUTTON,
               NRPS_WILDCARD_BUTTON,
               SAVE_TO_FOLDER_BUTTON,
               SAVE_TO_PNG_BUTTON,
               YES_BUTTON,
               NO_BUTTON,
               STARTER_BUTTON,
               ELONGATION_BUTTON,
               EXPORT_TABULAR_BUTTON,
               SAVE_TO_TXT_BUTTON]


def show_domain_buttons(module, screen, active_buttons):
    reset_buttons(screen, active_buttons)
    domain_types = set()

    for domain in module.domains:
        domain_types.add(domain.type)

    if module.type == 'NRPS':

        for domain, button in NRPS_DOMAIN_TO_BUTTON.items():
            if domain not in domain_types:
                show_button(button, screen, active_buttons)

    elif module.type == 'PKS':

        for domain, button in PKS_DOMAIN_TO_BUTTON.items():
            if domain not in domain_types:
                show_button(button, screen, active_buttons)


def hide_domain_buttons(screen, active_buttons):
    for button in active_buttons:
        if type(button) == DomainButton:
            hide_button(button, screen, active_buttons)


def highlight_buttons(active_buttons, screen, mouse):
    for button in active_buttons:
        if button.rectangle.collidepoint(mouse) or button.selected:
            button.highlight(screen)
        else:
            button.draw(screen)


def get_mouse_button(active_buttons, mouse):
    for button in active_buttons:
        if button.rectangle.collidepoint(mouse):
            return button

    return None


def make_buttons(screen):
    buttons = set()
    buttons.add(CREATE_GENE_BUTTON)
    CREATE_GENE_BUTTON.draw(screen)

    buttons.add(RENDER_CLUSTER_BUTTON)
    RENDER_CLUSTER_BUTTON.draw(screen)

    buttons.add(RENDER_PRODUCTS_BUTTON)
    RENDER_PRODUCTS_BUTTON.draw(screen)

    buttons.add(EXPORT_TABULAR_BUTTON)
    EXPORT_TABULAR_BUTTON.draw(screen)

    return buttons


def draw_buttons(active_buttons, screen, mouse):
    for button in active_buttons:
        if button.rectangle.collidepoint(mouse) or button.selected:
            button.highlight(screen)
        else:
            button.draw(screen)


def reset_buttons(screen, active_buttons):
    hide_buttons(ALL_BUTTONS, screen, active_buttons)
    show_button(CREATE_GENE_BUTTON, screen, active_buttons)
    show_button(RENDER_CLUSTER_BUTTON, screen, active_buttons)
    show_button(RENDER_PRODUCTS_BUTTON, screen, active_buttons)
    show_button(EXPORT_TABULAR_BUTTON, screen, active_buttons)


def hide_button(button, screen, active_buttons):
    try:
        active_buttons.remove(button)
    except KeyError:
        pass

    button.erase_button(screen)


def hide_buttons(buttons, screen, active_buttons):
    for button in buttons:
        hide_button(button, screen, active_buttons)


def show_button(button, screen, active_buttons):
    button.draw(screen)
    active_buttons.add(button)


def show_buttons(buttons, screen, active_buttons):
    for button in buttons:
        show_button(button, screen, active_buttons)


def make_raster(substrate_groups, buttons, buttons_per_line, button_size):
    x_coord = int(0.25 * WIDTH)
    y_coord = int(0.52 * HEIGHT)

    for i, group in enumerate(substrate_groups):
        if i % buttons_per_line == 0:
            if i != 0:
                y_coord += int(0.05 * HEIGHT)
            x_coord = int(0.25 * WIDTH)
        else:
            x_coord += button_size[0] + SUBSTRATE_GROUP_BUTTON_PADDING

        button = SubstrateGroupButton(group, (x_coord, y_coord), button_size)
        buttons.add(button)


def make_nrps_substrate_group_buttons():
    proteinogenic_buttons = set()
    nonproteinogenic_buttons = set()
    fatty_acid_buttons = set()
    non_amino_acid_buttons = set()

    make_raster(PROTEINOGENIC_SUBSTRATES, proteinogenic_buttons, SUBSTRATE_GROUP_BUTTONS_PER_LINE, SUBSTRATE_GROUP_BUTTON_SIZE)
    make_raster(NONPROTEINOGENIC_SUBSTRATES, nonproteinogenic_buttons, SUBSTRATE_GROUP_BUTTONS_PER_LINE, (int(WIDTH * 0.08), int(HEIGHT / 25)))
    make_raster(FATTY_ACIDS, fatty_acid_buttons, 1, (int(WIDTH * 0.15), int(HEIGHT / 25)))
    make_raster(NON_AMINO_ACIDS, non_amino_acid_buttons, 1, (int(WIDTH * 0.15), int(HEIGHT / 25)))

    return [proteinogenic_buttons, nonproteinogenic_buttons, fatty_acid_buttons, non_amino_acid_buttons]


def make_nrps_substrate_groups(supergroup, name_to_smiles, starter=False):
    group_to_buttons = {}
    for group, substrate_names in supergroup.items():
        group_to_buttons[group] = set()

        substrates = []

        for substrate_name in substrate_names:
            smiles = name_to_smiles[substrate_name]
            substrate = Substrate(substrate_name, smiles, 'NRPS', starter=starter)
            substrates.append(substrate)

        substrate_group = SubstrateGroup(group, substrates)
        for substrate in substrate_group.substrates:
            group_to_buttons[group].add(SubstrateButton(substrate,
                                                        (substrate.x, substrate.y),
                                                        (substrate.width, substrate.height)))

    return group_to_buttons


def make_nrps_substrate_buttons():

    name_to_smiles = parse_smiles(PARAS_SMILES)
    proteinogenic_group_to_buttons = make_nrps_substrate_groups(PROTEINOGENIC_SUBSTRATES,
                                                                name_to_smiles)

    nonproteinogenic_group_to_buttons = make_nrps_substrate_groups(NONPROTEINOGENIC_SUBSTRATES,
                                                                   name_to_smiles)
    fattyacid_group_to_buttons = make_nrps_substrate_groups(FATTY_ACIDS,
                                                            name_to_smiles, starter=True)
    nonaminoacid_group_to_buttons = make_nrps_substrate_groups(NON_AMINO_ACIDS,
                                                               name_to_smiles, starter=True)

    return proteinogenic_group_to_buttons, nonproteinogenic_group_to_buttons, fattyacid_group_to_buttons, \
           nonaminoacid_group_to_buttons


def make_fatty_acid_buttons():
    dimensions = (0.2 * WIDTH, HEIGHT / 25)

    fatty_acid_buttons = [FattyAcidSuperOptionButton('Set isoform', (0.25 * WIDTH, 0.52 * HEIGHT), dimensions),
                          FattyAcidSuperOptionButton('Add methyl group', (0.25 * WIDTH, 0.57 * HEIGHT), dimensions),
                          FattyAcidSuperOptionButton('Add amino group', (0.25 * WIDTH, 0.62 * HEIGHT), dimensions),
                          FattyAcidSuperOptionButton('Add -OH group', (0.25 * WIDTH, 0.67 * HEIGHT), dimensions),
                          FattyAcidSuperOptionButton('Add double bond', (0.25 * WIDTH, 0.72 * HEIGHT), dimensions),
                          FattyAcidSuperOptionButton('Set substrate', (0.25 * WIDTH, 0.82 * HEIGHT), dimensions)]

    small_dimensions = (0.1 * WIDTH, HEIGHT / 25)
    options_buttons = [FattyAcidOptionButton('iso', (0.47 * WIDTH, 0.52 * HEIGHT), small_dimensions),
                       FattyAcidOptionButton('anteiso', (0.47 * WIDTH, 0.57 * HEIGHT), dimensions),
                       FattyAcidOptionButton('cis', (0.47 * WIDTH, 0.52 * HEIGHT), small_dimensions),
                       FattyAcidOptionButton('trans', (0.47 * WIDTH, 0.57 * HEIGHT), dimensions),
                       FattyAcidOptionButton('undefined', (0.47 * WIDTH, 0.62 * HEIGHT), dimensions)]

    tiny_dimensions = (0.05 * WIDTH, HEIGHT / 25)

    c_buttons = [CButton(1, (0.47 * WIDTH, 0.52 * HEIGHT), tiny_dimensions),
                 CButton(2, (0.54 * WIDTH, 0.52 * HEIGHT), tiny_dimensions),
                 CButton(3, (0.61 * WIDTH, 0.52 * HEIGHT), tiny_dimensions),
                 CButton(4, (0.68 * WIDTH, 0.52 * HEIGHT), tiny_dimensions),
                 CButton(5, (0.75 * WIDTH, 0.52 * HEIGHT), tiny_dimensions),

                 CButton(6, (0.47 * WIDTH, 0.57 * HEIGHT), tiny_dimensions),
                 CButton(7, (0.54 * WIDTH, 0.57 * HEIGHT), tiny_dimensions),
                 CButton(8, (0.61 * WIDTH, 0.57 * HEIGHT), tiny_dimensions),
                 CButton(9, (0.68 * WIDTH, 0.57 * HEIGHT), tiny_dimensions),
                 CButton(10, (0.75 * WIDTH, 0.57 * HEIGHT), tiny_dimensions),

                 CButton(11, (0.47 * WIDTH, 0.62 * HEIGHT), tiny_dimensions),
                 CButton(12, (0.54 * WIDTH, 0.62 * HEIGHT), tiny_dimensions),
                 CButton(13, (0.61 * WIDTH, 0.62 * HEIGHT), tiny_dimensions),
                 CButton(14, (0.68 * WIDTH, 0.62 * HEIGHT), tiny_dimensions),
                 CButton(15, (0.75 * WIDTH, 0.62 * HEIGHT), tiny_dimensions),

                 CButton(16, (0.47 * WIDTH, 0.67 * HEIGHT), tiny_dimensions),
                 CButton(17, (0.54 * WIDTH, 0.67 * HEIGHT), tiny_dimensions),
                 CButton(18, (0.61 * WIDTH, 0.67 * HEIGHT), tiny_dimensions),
                 CButton(19, (0.68 * WIDTH, 0.67 * HEIGHT), tiny_dimensions),
                 CButton(20, (0.75 * WIDTH, 0.67 * HEIGHT), tiny_dimensions)]

    return fatty_acid_buttons, options_buttons, c_buttons


FATTY_ACID_SUPER_OPTIONS_BUTTONS, FATTY_ACID_OPTIONS_BUTTONS, C_BUTTONS = make_fatty_acid_buttons()

SET_ISOFORM_BUTTON, ADD_METHYL_GROUP_BUTTON, ADD_AMINO_GROUP_BUTTON, ADD_OH_GROUP_BUTTON, \
    ADD_DOUBLE_BOND_BUTTON, SET_SUBSTRATE_BUTTON = FATTY_ACID_SUPER_OPTIONS_BUTTONS

ISO_BUTTON, ANTEISO_BUTTON, CIS_BUTTON, TRANS_BUTTON, UNDEFINED_BUTTON = FATTY_ACID_OPTIONS_BUTTONS

ALL_BUTTONS += FATTY_ACID_SUPER_OPTIONS_BUTTONS
ALL_BUTTONS += FATTY_ACID_OPTIONS_BUTTONS
ALL_BUTTONS += C_BUTTONS


def show_carbon_nr_buttons(screen, active_buttons, mode='size', fatty_acid=None):
    if mode == 'size':
        show_buttons(C_BUTTONS[1:], screen, active_buttons)
    else:
        involved_positions = fatty_acid.count_involvement()
        if mode == 'double bond':
            c_number = fatty_acid.c_nr
            if fatty_acid.isoform:
                c_number -= 1

            double_bond_numbers = [bond[0] for bond in fatty_acid.double_bonds]
            for button in C_BUTTONS[:c_number - 2]:
                double_bond_possible = True
                for nr in double_bond_numbers:
                    if abs(button.nr - nr) <= 1:
                        double_bond_possible = False
                        break

                if button.nr in involved_positions and involved_positions[button.nr] >= 2 or \
                        (button.nr + 1) in involved_positions and involved_positions[button.nr + 1] >= 2:
                    double_bond_possible = False
                if double_bond_possible:
                    show_button(button, screen, active_buttons)

        else:
            c_number = fatty_acid.c_nr
            if fatty_acid.isoform:
                c_number -= 1

            for button in C_BUTTONS[:c_number - 1]:
                modification_possible = True

                if button.nr in involved_positions and involved_positions[button.nr] >= 2:
                    modification_possible = False

                if modification_possible:
                    show_button(button, screen, active_buttons)


def make_pks_substrate_starter_buttons():
    buttons = set()

    name_to_smiles = parse_smiles(PKS_SMILES)

    substrates = []

    for substrate_name, smiles in name_to_smiles.items():
        smiles = name_to_smiles[substrate_name]
        substrate = Substrate(substrate_name, smiles, 'PKS_starter', starter=True)
        substrates.append(substrate)

    substrate_group = SubstrateGroup('pks_starter', substrates)
    for substrate in substrate_group.substrates:
        buttons.add(SubstrateButton(substrate,
                                    (substrate.x, substrate.y),
                                    (substrate.width, substrate.height)))

    return buttons

PKS_STARTER_SUBSTRATE_BUTTONS = make_pks_substrate_starter_buttons()

def make_pks_substrate_buttons():
    buttons = set()

    name_to_smiles = {'malonylcoa': "CC(=O)S",
                      'methylmalonylcoa': "CCC(=O)S",
                      'methoxymalonylcoa': "COCC(=O)S",
                      'ethylmalonylcoa': "CCCC(=O)S",
                      'wildcard': "[*]CC(=O)S"}

    x = {'Malonyl-CoA': 'mal',
         'Methylmalonyl-CoA': 'mmal',
         'Methoxymalonyl-CoA': 'mxmal',
         'Ethylmalonyl-CoA': 'emal',
         'Isobutyryl-CoA': 'isobut',
         '2-Methylbutyryl-CoA': '2metbut',
         'trans-1,2-CPDA': 'trans-1,2-CPDA',
         'Acetyl-CoA': 'Acetyl-CoA',
         'Benzoyl-CoA': 'benz',
         'Propionyl-CoA': 'prop',
         '3-Methylbutyryl-CoA': '3metbut',
         'CE-Malonyl-CoA': 'cemal',
         '2-Rhyd-Malonyl-CoA': '2Rhydmal',
         'CHC-CoA': 'CHC-CoA',
         'inactive': 'inactive'}

    substrates = []

    for substrate_name, smiles in name_to_smiles.items():

        smiles = name_to_smiles[substrate_name]
        substrate = Substrate(substrate_name, smiles, 'PKS')
        substrates.append(substrate)

    substrate_group = SubstrateGroup('pks', substrates)
    for substrate in substrate_group.substrates:
        buttons.add(SubstrateButton(substrate,
                                    (substrate.x, substrate.y),
                                    (substrate.width, substrate.height)))

    return buttons

KR_SUBTYPES = ['', 'A1', 'A2', 'B1', 'B2', 'C1', 'C2']


def make_kr_buttons():
    buttons = set()
    x_coord = int(0.25 * WIDTH)
    y_coord = int(0.52 * HEIGHT)
    for kr_subtype in KR_SUBTYPES:

        buttons.add(KRSubtypeButton((x_coord, y_coord), kr_subtype))
        x_coord += 10 + KR_BUTTON_SIZE

    return buttons


PROTEINOGENIC_BUTTONS, NON_PROTEINOGENIC_BUTTONS, FATTY_ACID_BUTTONS, NON_AMINO_ACID_BUTTONS = \
    make_nrps_substrate_group_buttons()
PKS_SUBSTRATE_BUTTONS = make_pks_substrate_buttons()

PROTEINOGENIC_GROUP_TO_BUTTONS, NONPROTEINOGENIC_GROUP_TO_BUTTONS, FATTY_ACID_GROUP_TO_BUTTONS, \
           NON_AMINOACID_GROUP_TO_BUTTONS = make_nrps_substrate_buttons()

KR_BUTTONS = make_kr_buttons()

ALL_BUTTONS += PROTEINOGENIC_BUTTONS
ALL_BUTTONS += NON_PROTEINOGENIC_BUTTONS
ALL_BUTTONS += FATTY_ACID_BUTTONS
ALL_BUTTONS += NON_AMINO_ACID_BUTTONS
ALL_BUTTONS += PKS_SUBSTRATE_BUTTONS
ALL_BUTTONS += KR_BUTTONS
ALL_BUTTONS += PKS_STARTER_SUBSTRATE_BUTTONS

for buttons in PROTEINOGENIC_GROUP_TO_BUTTONS.values():
    ALL_BUTTONS += buttons

for buttons in NONPROTEINOGENIC_GROUP_TO_BUTTONS.values():
    ALL_BUTTONS += buttons

for buttons in FATTY_ACID_GROUP_TO_BUTTONS.values():
    ALL_BUTTONS += buttons

for buttons in NON_AMINOACID_GROUP_TO_BUTTONS.values():
    ALL_BUTTONS += buttons
