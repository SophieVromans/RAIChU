import pygame

from interactive.domain import Domain
from interactive.style import DOMAIN_SIZE, GENE_LABEL_SIZE, MODULE_SPACING, MODULE_COLOUR, \
    MODULE_LINE_COLOUR, MODULE_PADDING, MODULE_HIGHLIGHT_COLOUR


class Module:
    module_to_order = {'NRPS': ['C', 'A', 'nMT', 'PCP', 'E', 'TE'],
                       'PKS': ['KS', 'AT', 'DH', 'ER', 'KR', 'ACP', 'TE']}

    domain_to_subsequent_domains = {'C': ['A', 'nMT', 'PCP', 'E', 'TE'],
                                    'A': ['nMT', 'PCP', 'E', 'TE'],
                                    'nMT': ['PCP', 'E', 'TE'],
                                    'PCP': ['E', 'TE'],
                                    'E': ['TE'],
                                    'KS': ['AT', 'DH', 'ER', 'KR', 'ACP', 'TE'],
                                    'AT': ['DH', 'ER', 'KR', 'ACP', 'TE'],
                                    'DH': ['ER', 'KR', 'ACP', 'TE'],
                                    'ER': ['KR', 'ACP', 'TE'],
                                    'KR': ['ACP', 'TE'],
                                    'ACP': [],
                                    'TE': []}

    def __init__(self, screen, module_index, module_type, gene):
        assert module_type in ['NRPS', 'PKS']

        self.screen = screen
        self.id = module_index
        self.gene = gene
        self.type = module_type
        self.domains = []
        self.rectangle = None
        self.dna_coords = None

        self.x = 0
        self.y = 0

        self.width = 0
        self.height = DOMAIN_SIZE

        self.set_rectangle()

        self.selected = False

    def set_rectangle(self):
        self.y = self.gene.y - 5
        self.x = self.gene.x + GENE_LABEL_SIZE
        for module in self.gene.modules:
            if module.id != self.id:
                self.x += len(module.domains) * DOMAIN_SIZE
                self.x += MODULE_SPACING
                self.x += MODULE_PADDING
            else:
                break
        self.x += MODULE_PADDING

        self.width = len(self.domains) * DOMAIN_SIZE + MODULE_SPACING
        self.rectangle = pygame.Rect(self.x, self.y, self.width, self.height)

    def adjust_domain_indices(self, insertion_point):
        for domain in self.domains:
            if domain.id >= insertion_point:
                domain.id += 1
                domain.set_rectangle()

    def draw(self, mouse):
        self.set_rectangle()

        colour = MODULE_COLOUR
        if self.rectangle.collidepoint(mouse) or self.selected:
            for domain in self.domains:
                if domain.rectangle.collidepoint(mouse) or domain.selected:
                    break
            colour = MODULE_HIGHLIGHT_COLOUR

        pygame.draw.rect(self.screen, colour, self.rectangle)
        pygame.draw.rect(self.screen, MODULE_LINE_COLOUR, self.rectangle, 1)

        for domain in self.domains:
            domain.draw(mouse)

    def remove_domain(self, domain):
        domain_index = domain.id
        self.domains.pop(domain_index)

        for domain in self.domains:
            if domain.id > domain_index:
                domain.id -= 1

    def add_domain(self, domain_type, domain_label=None):
        assert domain_type not in [domain.type for domain in self.domains]
        assert domain_type in self.module_to_order[self.type]

        insertion_point = 0

        for i, domain in enumerate(self.domains):
            if domain.type not in self.domain_to_subsequent_domains[domain_type]:
                insertion_point += 1
            else:
                break

        self.adjust_domain_indices(insertion_point)

        domain = Domain(self.screen, domain_type, self, insertion_point, domain_label=domain_label)

        self.domains.insert(insertion_point, domain)









