import os

import pygame

from interactive.buttons import make_buttons, get_mouse_button, Button, DomainButton, hide_button, \
    AddGeneButton, CreateGeneButton, show_buttons, ADD_MODULE_BUTTON, show_button, \
    AddModuleButton, NRPSModuleButton, NRPS_MODULE_BUTTON, PKSModuleButton, PKS_MODULE_BUTTON, \
    show_domain_buttons, reset_buttons, AddDomainButton, ADD_DOMAIN_BUTTON, \
    REMOVE_MODULE_BUTTON, REMOVE_DOMAIN_BUTTON, RemoveDomainButton, RemoveModuleButton, RemoveGeneButton, \
    REMOVE_GENE_BUTTON, SELECT_SUBSTRATE_BUTTON, SelectSubstrateButton, SELECT_DOMAIN_TYPE_BUTTON, \
    SelectDomainTypeButton, PKS_SUBSTRATE_BUTTONS, SubstrateSupergroupButton,\
    PROTEINOGENIC_BUTTONS, NON_PROTEINOGENIC_BUTTONS, FATTY_ACID_BUTTONS, NON_AMINO_ACID_BUTTONS, \
    ProteinogenicButton, NonProteinogenicButton, FattyAcidButton, NonAminoAcidButton, \
    PROTEINOGENIC_GROUP_TO_BUTTONS, NONPROTEINOGENIC_GROUP_TO_BUTTONS, FATTY_ACID_GROUP_TO_BUTTONS, \
    NON_AMINOACID_GROUP_TO_BUTTONS, SubstrateGroupButton, SubstrateButton, RenderClusterButton, \
    KR_BUTTONS, KRSubtypeButton, SaveClusterButton, SaveProductsButton, RenderProductsButton, \
    NRPSWildcardButton, SaveToPngButton, SaveToFolderButton, FATTY_ACID_SUPER_OPTIONS_BUTTONS, \
    hide_buttons, ALL_BUTTONS, SAVE_PRODUCTS_BUTTON, SAVE_CLUSTER_BUTTON, YES_BUTTON, NO_BUTTON, YesButton, NoButton, \
    FattyAcidOptionButton, FattyAcidSuperOptionButton, ISO_BUTTON, ANTEISO_BUTTON, show_carbon_nr_buttons, \
    CIS_BUTTON, TRANS_BUTTON, UNDEFINED_BUTTON, CButton, SET_ISOFORM_BUTTON, StarterButton, ElongationButton, \
    PKS_STARTER_SUBSTRATE_BUTTONS, PROTEINOGENIC_BUTTON, FATTY_ACID_BUTTON, NON_PROTEINOGENIC_BUTTON, \
    NON_AMINO_ACID_BUTTON, STARTER_BUTTON, ELONGATION_BUTTON, ExportTabularButton, SaveToTxtButton, NRPS_WILDCARD_BUTTON
from interactive.domain import Domain
from interactive.module import Module
from interactive.gene import Gene
from interactive.substrate import FattyAcid
from interactive.insertion_point import InsertionPoint
from interactive.render_cluster import render_cluster, render_products, export_tabular
from interactive.style import RENDER_WINDOW_SIZE, WHITE, REPLACE_TEXT_1, REPLACE_TEXT_2, HEIGHT, WIDTH, \
    SIZE, FATTY_ACID_IMAGE_SIZE, FAILED_TEXT_1, FAILED_TEXT_3, FAILED_TEXT_2


class RaichuManager:
    states = {'add gene',
              'remove_gene',
              'add module',
              'remove module',
              'add domain',
              'remove domain',
              'set domain subtype',

              }

    def __init__(self, screen):
        self.genes = []
        self.state = None
        self.active_buttons = make_buttons(screen)
        self.screen = screen

        self.selected_module = None
        self.selected_gene = None
        self.selected_domain = None
        self.selected_substrate_group = None

        self.text_box = None
        self.insertion_point = None
        self.cluster_image = None
        self.product_images = []
        self.replace_text = [REPLACE_TEXT_1, REPLACE_TEXT_2]
        self.fatty_acid_mode = 'size'
        self.fatty_acid = None
        self.fatty_acid_bond_type = None
        self.fatty_acid_modified = False

    def get_mouse_domain(self, mouse):
        for gene in self.genes:
            for module in gene.modules:
                for domain in module.domains:
                    if domain.rectangle.collidepoint(mouse):
                        return domain

        return None

    def get_mouse_module(self, mouse):
        for gene in self.genes:
            for module in gene.modules:
                if module.rectangle.collidepoint(mouse):
                    return module

        return None

    def get_mouse_gene(self, mouse):
        for gene in self.genes:
            if gene.rectangle.collidepoint(mouse):
                return gene

        return None

    def detect_click_location(self, mouse):

        button = get_mouse_button(self.active_buttons, mouse)
        domain = self.get_mouse_domain(mouse)
        module = self.get_mouse_module(mouse)
        gene = self.get_mouse_gene(mouse)

        if button:
            return button
        elif domain:
            self.reset_selections()
            domain.selected = True
            self.selected_domain = domain
            return domain
        elif module:
            self.reset_selections()
            module.selected = True
            self.selected_module = module
            return module
        elif gene:
            self.reset_selections()
            gene.selected = True
            self.selected_gene = gene
            return gene
        else:
            return None

    def reset_selections(self):
        if self.selected_gene:
            self.selected_gene.insertion_point = None
            self.selected_gene.selected = False
            self.selected_gene = None
        if self.selected_module:
            self.selected_module.selected = False
            self.selected_module = None
        if self.selected_domain:
            self.selected_domain.selected = False
            self.selected_domain = None
        self.selected_substrate_group = None

    def do_button_action(self, button, mouse):
        if type(button) == DomainButton or issubclass(type(button), DomainButton):
            button.do_action(self.selected_module, mouse, self.screen, self.active_buttons)
            self.selected_module.selected = False
            self.selected_module = None
            hide_button(button, self.screen, self.active_buttons)
        elif issubclass(type(button), SubstrateSupergroupButton):
            reset_buttons(self.screen, self.active_buttons)
            if type(button) == ProteinogenicButton:
                show_buttons(PROTEINOGENIC_BUTTONS, self.screen, self.active_buttons)
                self.selected_substrate_group = PROTEINOGENIC_GROUP_TO_BUTTONS
            elif type(button) == NonProteinogenicButton:
                show_buttons(NON_PROTEINOGENIC_BUTTONS, self.screen, self.active_buttons)
                self.selected_substrate_group = NONPROTEINOGENIC_GROUP_TO_BUTTONS
            elif type(button) == FattyAcidButton:
                show_buttons(FATTY_ACID_BUTTONS, self.screen, self.active_buttons)
                self.selected_substrate_group = FATTY_ACID_GROUP_TO_BUTTONS
            elif type(button) == NonAminoAcidButton:
                show_buttons(NON_AMINO_ACID_BUTTONS, self.screen, self.active_buttons)
                self.selected_substrate_group = NON_AMINOACID_GROUP_TO_BUTTONS
        elif type(button) == SubstrateGroupButton:
            reset_buttons(self.screen, self.active_buttons)
            if button.text != 'Custom':
                show_buttons(self.selected_substrate_group[button.text], self.screen, self.active_buttons)
                self.selected_substrate_group = None
            else:
                show_carbon_nr_buttons(self.screen, self.active_buttons)

        elif type(button) == StarterButton:
            reset_buttons(self.screen, self.active_buttons)
            if self.selected_domain.module.type == 'NRPS':
                show_buttons([FATTY_ACID_BUTTON, NON_AMINO_ACID_BUTTON, NRPS_WILDCARD_BUTTON], self.screen, self.active_buttons)
            elif self.selected_domain.module.type == 'PKS':
                show_buttons(PKS_STARTER_SUBSTRATE_BUTTONS, self.screen, self.active_buttons)

        elif type(button) == ElongationButton:
            reset_buttons(self.screen, self.active_buttons)
            if self.selected_domain.module.type == 'NRPS':
                show_buttons([PROTEINOGENIC_BUTTON, NON_PROTEINOGENIC_BUTTON, NRPS_WILDCARD_BUTTON], self.screen, self.active_buttons)
            elif self.selected_domain.module.type == 'PKS':
                show_buttons(PKS_SUBSTRATE_BUTTONS, self.screen, self.active_buttons)

        elif type(button) == FattyAcidSuperOptionButton:

            reset_buttons(self.screen, self.active_buttons)

            if button.text == 'Set substrate':
                self.fatty_acid.set_substrate(self.selected_domain)
                self.selected_domain.module.gene.erase()
                self.selected_domain.module.gene.draw(mouse)
                reset_buttons(self.screen, self.active_buttons)
                self.reset_selections()
                self.fatty_acid = None
                self.fatty_acid_mode = 'size'
                self.fatty_acid_bond_type = None
                self.fatty_acid_modified = False
                pygame.draw.rect(self.screen, WHITE,
                                 pygame.Rect(SIZE[0] - (FATTY_ACID_IMAGE_SIZE[0] + 20), SIZE[1] / 2 + 20,
                                             FATTY_ACID_IMAGE_SIZE[0], FATTY_ACID_IMAGE_SIZE[1]))
                pygame.draw.rect(self.screen, WHITE,
                                 pygame.Rect(SIZE[0] - (FATTY_ACID_IMAGE_SIZE[0] + 20), SIZE[1] / 2 + 20,
                                             FATTY_ACID_IMAGE_SIZE[0], FATTY_ACID_IMAGE_SIZE[1]), 2)

            else:

                if button.text == 'Add double bond':
                    self.fatty_acid_mode = 'double bond'
                    show_buttons([CIS_BUTTON, TRANS_BUTTON, UNDEFINED_BUTTON], self.screen, self.active_buttons)

                elif button.text == 'Set isoform':
                    if self.fatty_acid.c_nr > 4:
                        show_buttons([ISO_BUTTON, ANTEISO_BUTTON], self.screen, self.active_buttons)
                    elif self.fatty_acid.c_nr == 4:
                        show_button(ISO_BUTTON, self.screen, self.active_buttons)

                else:
                    if button.text == 'Add amino group':
                        self.fatty_acid_mode = 'amino'
                    elif button.text == 'Add methyl group':
                        self.fatty_acid_mode = 'methyl'
                    elif button.text == 'Add -OH group':
                        self.fatty_acid_mode = 'hydroxyl'

                    show_carbon_nr_buttons(self.screen, self.active_buttons, mode=self.fatty_acid_mode, fatty_acid=self.fatty_acid)
                show_buttons(FATTY_ACID_SUPER_OPTIONS_BUTTONS, self.screen, self.active_buttons)
                if self.fatty_acid_modified:
                    hide_button(SET_ISOFORM_BUTTON, self.screen, self.active_buttons)

        elif type(button) == FattyAcidOptionButton:
            reset_buttons(self.screen, self.active_buttons)
            if button.text == 'cis':
                self.fatty_acid_bond_type = 'cis'
            elif button.text == 'trans':
                self.fatty_acid_bond_type = 'trans'
            elif button.text == 'iso':
                self.fatty_acid.set_isoform('iso')
            elif button.text == 'anteiso':
                self.fatty_acid.set_isoform('anteiso')

            if button.text in {'cis', 'trans', 'undefined'}:
                self.fatty_acid_mode = 'double bond'
                show_carbon_nr_buttons(self.screen, self.active_buttons, mode=self.fatty_acid_mode,
                                       fatty_acid=self.fatty_acid)

            show_buttons(FATTY_ACID_SUPER_OPTIONS_BUTTONS, self.screen, self.active_buttons)
            if self.fatty_acid_modified:
                hide_button(SET_ISOFORM_BUTTON, self.screen, self.active_buttons)

        elif type(button) == CButton:
            if self.fatty_acid_mode == 'double bond':
                self.fatty_acid.add_double_bond(button.nr, self.fatty_acid_bond_type)
                self.fatty_acid_bond_type = None
                self.fatty_acid_modified = True
            elif self.fatty_acid_mode == 'size':
                self.fatty_acid = FattyAcid(button.nr)
            elif self.fatty_acid_mode == 'amino':
                self.fatty_acid.add_amino_group(button.nr)
                self.fatty_acid_modified = True
            elif self.fatty_acid_mode == 'methyl':
                self.fatty_acid.add_methyl_group(button.nr)
                self.fatty_acid_modified = True
            elif self.fatty_acid_mode == 'hydroxyl':
                self.fatty_acid.add_hydroxyl_group(button.nr)
                self.fatty_acid_modified = True
            self.fatty_acid_mode = 'size'

            reset_buttons(self.screen, self.active_buttons)
            show_buttons(FATTY_ACID_SUPER_OPTIONS_BUTTONS, self.screen, self.active_buttons)
            if self.fatty_acid_modified:
                hide_button(SET_ISOFORM_BUTTON, self.screen, self.active_buttons)
            elif self.fatty_acid and self.fatty_acid.c_nr < 4:
                hide_button(SET_ISOFORM_BUTTON, self.screen, self.active_buttons)


        elif type(button) == AddGeneButton:
            button.do_action(self.screen, self.genes, self.text_box, self.active_buttons, mouse)
            self.text_box = None
        elif type(button) == CreateGeneButton:
            self.reset_selections()
            self.text_box = button.do_action(self.screen, self.active_buttons)
        elif type(button) == AddModuleButton:
            show_buttons([PKS_MODULE_BUTTON, NRPS_MODULE_BUTTON], self.screen, self.active_buttons)
        elif type(button) == PKSModuleButton or type(button) == NRPSModuleButton:
            button.do_action(self.screen, self.selected_gene, self.active_buttons, self.insertion_point, mouse)
            self.reset_selections()
            self.insertion_point = None
        elif type(button) == RemoveGeneButton:
            button.do_action(self.screen, self.genes, self.selected_gene, self.active_buttons, mouse)
        elif type(button) == AddDomainButton:
            reset_buttons(self.screen, self.active_buttons)
            show_domain_buttons(self.selected_module, self.screen, self.active_buttons)
        elif type(button) == RemoveModuleButton:
            button.do_action(self.screen, self.selected_module, self.active_buttons, mouse)
            self.reset_selections()
        elif type(button) == RemoveDomainButton:
            button.do_action(self.selected_domain, self.screen, self.active_buttons, mouse)
            self.reset_selections()
        elif type(button) == RenderProductsButton:

            try:

                render_products(self.genes, os.path.join(os.getcwd(), 'tmp_out_raichu'))
                product_dir = os.path.join(os.getcwd(), 'tmp_out_raichu')

                product_images = []

                for product_name in os.listdir(product_dir):
                    if product_name.endswith('.png'):
                        product_file = os.path.join(os.path.join(product_dir, product_name))
                        product_image = pygame.image.load(product_file)
                        product_images.append(product_image)

                product_image_height = RENDER_WINDOW_SIZE[1] / 2
                product_image_width = RENDER_WINDOW_SIZE[0] / 3

                for i, image in enumerate(product_images):
                    image_width = image.get_width()
                    image_height = image.get_height()

                    if image_width <= product_image_width and image_height > product_image_height:
                        ratio = product_image_height / image_height
                        image_width *= ratio
                        image_height = product_image_height
                    elif image_width > product_image_width and image_height <= product_image_height:
                        ratio = product_image_width / image_width
                        image_height *= ratio
                        image_width = product_image_width
                    elif image_width > product_image_width and image_height > product_image_height:
                        ratio_1 = product_image_width / image_width
                        ratio_2 = product_image_height / image_height
                        if ratio_1 <= ratio_2:
                            ratio = ratio_1
                        else:
                            ratio = ratio_2

                        image_height *= ratio
                        image_width *= ratio

                    render_window_size = (int(image_width), int(image_height))

                    self.product_images.append(pygame.transform.smoothscale(image, render_window_size))
                hide_buttons(ALL_BUTTONS, self.screen, self.active_buttons)
                show_button(SAVE_PRODUCTS_BUTTON, self.screen, self.active_buttons)

            except Exception as e:
                print(e)
                self.screen.blit(FAILED_TEXT_3, (0.42 * WIDTH, 0.52 * HEIGHT))
                self.screen.blit(FAILED_TEXT_2, (0.42 * WIDTH, 0.57 * HEIGHT))

            self.reset_selections()
        elif type(button) == RenderClusterButton:

            try:

                render_cluster(self.genes, os.path.join(os.getcwd(), 'tmp_out_raichu.png'))

                cluster_image = pygame.image.load(os.path.join(os.getcwd(), 'tmp_out_raichu.png'))

                image_width = cluster_image.get_width()
                image_height = cluster_image.get_height()

                if image_width <= RENDER_WINDOW_SIZE[0] and image_height > RENDER_WINDOW_SIZE[1]:
                    ratio = RENDER_WINDOW_SIZE[1] / image_height
                    image_width *= ratio
                    image_height = RENDER_WINDOW_SIZE[1]
                elif image_width > RENDER_WINDOW_SIZE[0] and image_height <= RENDER_WINDOW_SIZE[1]:
                    ratio = RENDER_WINDOW_SIZE[0] / image_width
                    image_height *= ratio
                    image_width = RENDER_WINDOW_SIZE[0]
                elif image_width > RENDER_WINDOW_SIZE[0] and image_height > RENDER_WINDOW_SIZE[1]:
                    ratio_1 = RENDER_WINDOW_SIZE[0] / image_width
                    ratio_2 = RENDER_WINDOW_SIZE[1] / image_height
                    if ratio_1 <= ratio_2:
                        ratio = ratio_1
                    else:
                        ratio = ratio_2

                    image_height *= ratio
                    image_width *= ratio

                render_window_size = (int(image_width), int(image_height))

                self.cluster_image = pygame.transform.smoothscale(cluster_image, render_window_size)
                hide_buttons(ALL_BUTTONS, self.screen, self.active_buttons)
                show_button(SAVE_CLUSTER_BUTTON, self.screen, self.active_buttons)
            except Exception as e:
                print(e)
                self.screen.blit(FAILED_TEXT_1, (0.42 * WIDTH, 0.52 * HEIGHT))
                self.screen.blit(FAILED_TEXT_2, (0.42 * WIDTH, 0.57 * HEIGHT))

            self.reset_selections()

        elif type(button) == SaveProductsButton or type(button) == SaveClusterButton:
            self.reset_selections()
            self.text_box = button.do_action(self.screen, self.active_buttons)
        elif type(button) == SaveToTxtButton:
            if self.text_box.text.endswith('.txt'):
                new_file = os.path.join(os.getcwd(), self.text_box.text)
            else:
                new_file = os.path.join(os.getcwd(), self.text_box.text + '.txt')

            if os.path.exists(new_file):
                for i, replace_text in enumerate(self.replace_text):
                    self.screen.blit(replace_text, (0.02 * WIDTH, 0.52 * HEIGHT + 0.05 * HEIGHT * i))
                show_buttons([YES_BUTTON, NO_BUTTON], self.screen, self.active_buttons)
            else:
                export_tabular(self.genes, new_file)
                self.random_click(mouse)

        elif type(button) == ExportTabularButton:
            self.reset_selections()
            self.text_box = button.do_action(self.screen, self.active_buttons)
        elif type(button) == SaveToPngButton:
            png_dir = os.path.join(os.getcwd(), 'tmp_out_raichu.png')
            if self.text_box.text.endswith('.png'):
                new_file = os.path.join(os.getcwd(), self.text_box.text)
            else:
                new_file = os.path.join(os.getcwd(), self.text_box.text + '.png')
            new_png_dir = os.path.join(os.getcwd(), new_file)
            if os.path.exists(new_png_dir):
                for i, replace_text in enumerate(self.replace_text):
                    self.screen.blit(replace_text, (0.02 * WIDTH, 0.52 * HEIGHT + 0.05 * HEIGHT * i))
                show_buttons([YES_BUTTON, NO_BUTTON], self.screen, self.active_buttons)
            else:
                os.replace(png_dir, new_png_dir)
                self.random_click(mouse)
        elif type(button) == SaveToFolderButton:
            folder = os.path.join(os.getcwd(), 'tmp_out_raichu')
            new_folder = os.path.join(os.getcwd(), self.text_box.text)
            if os.path.isdir(new_folder):
                for i, replace_text in enumerate(self.replace_text):
                    self.screen.blit(replace_text, (0.02 * WIDTH, 0.52 * HEIGHT + 0.05 * HEIGHT * i))
                show_buttons([YES_BUTTON, NO_BUTTON], self.screen, self.active_buttons)
            else:
                os.replace(folder, new_folder)
                self.random_click(mouse)

        elif type(button) == YesButton:
            if self.product_images:
                folder = os.path.join(os.getcwd(), 'tmp_out_raichu')
                new_folder = os.path.join(os.getcwd(), self.text_box.text)
                for file_name in os.listdir(new_folder):
                    file_path = os.path.join(new_folder, file_name)
                    os.remove(file_path)
                os.rmdir(new_folder)
                os.replace(folder, new_folder)
                self.random_click(mouse)
            elif self.cluster_image:
                png_dir = os.path.join(os.getcwd(), 'tmp_out_raichu.png')
                if self.text_box.text.endswith('.png'):
                    new_file = os.path.join(os.getcwd(), self.text_box.text)
                else:
                    new_file = os.path.join(os.getcwd(), self.text_box.text + '.png')

                os.replace(png_dir, new_file)
                self.random_click(mouse)
            else:
                if self.text_box.text.endswith('.txt'):
                    directory = os.path.join(os.getcwd(), self.text_box.text)
                else:
                    directory = os.path.join(os.getcwd(), self.text_box.text + '.txt')
                export_tabular(self.genes, directory)
                self.random_click(mouse)

        elif type(button) == NoButton:
            self.random_click(mouse)

        elif type(button) == SelectDomainTypeButton:
            reset_buttons(self.screen, self.active_buttons)
            show_buttons(KR_BUTTONS, self.screen, self.active_buttons)
        elif type(button) == SelectSubstrateButton:
            if self.selected_domain.type == 'A':
                reset_buttons(self.screen, self.active_buttons)
                show_buttons([STARTER_BUTTON, ELONGATION_BUTTON], self.screen, self.active_buttons)

            elif self.selected_domain.type == 'AT':
                reset_buttons(self.screen, self.active_buttons)
                show_buttons([STARTER_BUTTON, ELONGATION_BUTTON], self.screen, self.active_buttons)
        elif type(button) == KRSubtypeButton:
            button.do_action(self.selected_domain, mouse, self.screen, self.active_buttons)
            self.reset_selections()

        elif type(button) == SubstrateButton or type(button) == NRPSWildcardButton:
            self.selected_domain.substrate = button.substrate
            self.selected_domain.module.gene.erase()
            self.selected_domain.module.gene.draw(mouse)
            reset_buttons(self.screen, self.active_buttons)
            self.reset_selections()

    def do_click_action(self, mouse):

        selected_entity = self.detect_click_location(mouse)
        if selected_entity:
            for gene in self.genes:
                gene.draw(mouse)
            if type(selected_entity) == Button or issubclass(type(selected_entity), Button):
                self.do_button_action(selected_entity, mouse)
            elif type(selected_entity) == Gene:
                reset_buttons(self.screen, self.active_buttons)
                self.selected_gene = selected_entity
                self.insertion_point = InsertionPoint(selected_entity, mouse)
                selected_entity.insertion_point = self.insertion_point
                selected_entity.draw(mouse)
                show_buttons([ADD_MODULE_BUTTON, REMOVE_GENE_BUTTON], self.screen, self.active_buttons)
            elif type(selected_entity) == Module:
                self.selected_module = selected_entity
                reset_buttons(self.screen, self.active_buttons)
                show_buttons([ADD_DOMAIN_BUTTON, REMOVE_MODULE_BUTTON], self.screen,
                             self.active_buttons)
            elif type(selected_entity) == Domain:
                self.selected_domain = selected_entity
                reset_buttons(self.screen, self.active_buttons)
                show_buttons([REMOVE_DOMAIN_BUTTON], self.screen, self.active_buttons)

                if selected_entity.type == 'A' or selected_entity.type == 'AT':
                    show_button(SELECT_SUBSTRATE_BUTTON, self.screen, self.active_buttons)
                elif selected_entity.type == 'KR':
                    show_button(SELECT_DOMAIN_TYPE_BUTTON, self.screen, self.active_buttons)

        else:
            self.random_click(mouse)

    def random_click(self, mouse):
        self.screen.fill(WHITE)

        reset_buttons(self.screen, self.active_buttons)
        for gene in self.genes:
            gene.draw(mouse)

        for gene in self.genes:
            gene.selected = False
            gene.insertion_point = None
            for module in gene.modules:
                module.selected = False
                for domain in module.domains:
                    domain.selected = False

        self.selected_gene = None
        self.selected_module = None
        self.selected_domain = None
        self.selected_substrate_group = None
        self.cluster_image = False
        self.product_images = []
        self.text_box = None
        self.fatty_acid = None
        self.fatty_acid_mode = 'size'
        self.fatty_acid_bond_type = None
        self.fatty_acid_modified = False
        pygame.draw.rect(self.screen, WHITE, pygame.Rect(SIZE[0] - (FATTY_ACID_IMAGE_SIZE[0] + 20), SIZE[1] / 2 + 20,
                                                         FATTY_ACID_IMAGE_SIZE[0], FATTY_ACID_IMAGE_SIZE[1]))
        pygame.draw.rect(self.screen, WHITE, pygame.Rect(SIZE[0] - (FATTY_ACID_IMAGE_SIZE[0] + 20), SIZE[1] / 2 + 20,
                                                         FATTY_ACID_IMAGE_SIZE[0], FATTY_ACID_IMAGE_SIZE[1]), 2)
        raichu_dir = os.path.join(os.getcwd(), 'tmp_out_raichu')
        raichu_png_dir = os.path.join(os.getcwd(), 'tmp_out_raichu.png')
        if os.path.exists(raichu_png_dir):
            os.remove(raichu_png_dir)
        if os.path.exists(raichu_dir):
            for file_name in os.listdir(raichu_dir):
                file_path = os.path.join(raichu_dir, file_name)
                os.remove(file_path)
            os.rmdir(raichu_dir)







