import pygame
from interactive.style import BLACK, GENE_COLOUR, GENE_LABEL_SIZE, GENE_SPACING, GENE_PADDING, \
    GENE_HEIGHT, BACKGROUND_COLOUR, FONT, MODULE_PADDING, DOMAIN_SIZE, GENE_HIGHLIGHT_COLOUR
from interactive.module import Module

GENE_X_START = 30


class Gene:
    def __init__(self, screen, gene_number):
        self.name = None
        self.gene_number = gene_number
        self.screen = screen
        self.modules = []
        self.gene_length = 0
        self.x = GENE_PADDING
        self.y = GENE_PADDING + self.gene_number * GENE_SPACING
        self.width = 0
        self.height = GENE_HEIGHT
        self.dna_coords = None
        self.strand = None

        self.rectangle = None
        self.set_rectangle()

        self.font_size = self.height - 10

        self.text_x = self.x + 2
        self.text_y = self.y + int((self.height - self.font_size) / 2)

        self.font = pygame.font.SysFont(FONT, self.font_size, bold=False)

        self.selected = False

        self.insertion_point = None

    def sort_modules(self):
        self.modules.sort(key=lambda module: module.dna_coords[0])
        if self.strand == -1:
            self.modules.reverse()

        for i, module in enumerate(self.modules):
            module.id = i

    def set_rectangle(self):
        self.y = GENE_PADDING + self.gene_number * GENE_SPACING
        self.width = GENE_LABEL_SIZE
        for module in self.modules:
            self.width += module.width
            self.width += MODULE_PADDING

        self.width += MODULE_PADDING

        self.rectangle = pygame.Rect(self.x, self.y, self.width, self.height)

    def draw(self, mouse):

        for module in self.modules:
            module.set_rectangle()

        self.set_rectangle()

        self.text_x = self.x + 2
        self.text_y = self.y + int((self.height - self.font_size) / 2)

        colour = GENE_COLOUR
        if self.rectangle.collidepoint(mouse) or self.selected:
            for module in self.modules:
                if module.rectangle.collidepoint(mouse) or module.selected:
                    break
            colour = GENE_HIGHLIGHT_COLOUR

        pygame.draw.rect(self.screen, colour, self.rectangle)
        pygame.draw.rect(self.screen, BLACK, self.rectangle, 1)

        rendered_text = self.font.render(self.name, True, BLACK)
        self.screen.blit(rendered_text, (self.text_x, self.text_y))

        for module in self.modules:
            module.draw(mouse)

        if self.insertion_point:
            self.insertion_point.draw(self.screen)

    def erase(self):
        height_difference = DOMAIN_SIZE - self.height
        rectangle = pygame.Rect(self.x, self.y - height_difference / 2 - 20, self.width, self.height + height_difference + 20)
        pygame.draw.rect(self.screen, BACKGROUND_COLOUR, rectangle)
        #pygame.draw.rect(self.screen, BACKGROUND_COLOUR, rectangle, 5)

    def adjust_module_indices(self, insertion_point):
        for module in self.modules:
            if module.id >= insertion_point:
                module.id += 1
                module.set_rectangle()

    def remove_module(self, module_nr):
        self.modules.pop(module_nr)

        for module in self.modules:
            if module.id > module_nr:
                module.id -= 1

    def get_domains(self):
        domains = []
        for module in self.modules:
            for domain in module.domains:
                domains.append(domain)

        return domains

    def add_module(self, insertion_point, module_type):

        self.adjust_module_indices(insertion_point)
        new_module = Module(self.screen, insertion_point, module_type, self)

        starter_module = False
        if insertion_point == 0 and self.gene_number == 0:
            starter_module = True

        domains = None

        if module_type == 'NRPS' and starter_module:
            domains = ['A', 'PCP']
        elif module_type == 'NRPS' and not starter_module:
            domains = ['C', 'A', 'PCP']
        elif module_type == 'PKS' and starter_module:
            domains = ['AT', 'ACP']
        elif module_type == 'PKS' and not starter_module:
            domains = ['KS', 'AT', 'ACP']

        assert domains

        for domain in domains:
            new_module.add_domain(domain)

        self.modules.insert(insertion_point, new_module)




