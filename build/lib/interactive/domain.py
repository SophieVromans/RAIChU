import os
import interactive.images.domains
from interactive.style import DOMAIN_SIZE, MODULE_SPACING, FONT, BLACK, RED

import pygame

DOMAIN_IMAGE_DIR = os.path.dirname(interactive.images.domains.__file__)


class Domain:
    def __init__(self, screen, domain_type, module, domain_id, domain_label=None):
        self.screen = screen
        self.module = module
        self.type = domain_type
        self.image = os.path.join(DOMAIN_IMAGE_DIR, f"{self.type.lower()}.png")
        self.highlight_image = os.path.join(DOMAIN_IMAGE_DIR, f"{self.type.lower()}_highlight.png")
        self.subtype = None
        self.id = domain_id
        self.rectangle = None
        self.substrate_rectangle = None
        self.substrate = None
        self.domain_label = domain_label

        self.x = 0
        self.y = 0

        self.width = DOMAIN_SIZE
        self.height = DOMAIN_SIZE

        self.font_size = self.height - 16
        self.font = pygame.font.SysFont(FONT, self.font_size, bold=False)

        self.set_rectangle()

        self.selected = False

    def set_rectangle(self):

        self.x = self.module.x + MODULE_SPACING / 2 + self.id * DOMAIN_SIZE
        self.y = self.module.y
        self.rectangle = pygame.Rect(self.x, self.y, self.width, self.height)

    def draw_domain_subtype(self):
        rendered_text = self.font.render(self.subtype, False, BLACK)
        text_rectangle = rendered_text.get_rect(center=(self.x + self.width / 2, self.y - 10))
        self.screen.blit(rendered_text, text_rectangle)

    def draw_substrate_text(self):
        colour = BLACK
        if self.module.gene.gene_number != 0 or self.module.id != 0:
            if self.substrate.starter:
                colour = RED

        rendered_text = self.font.render(self.substrate.abbr, False, colour)
        text_rectangle = rendered_text.get_rect(center=(self.x + self.width/2, self.y - 10))
        self.screen.blit(rendered_text, text_rectangle)

    def set_domain_subtype(self, subtype):
        self.subtype = subtype

    def draw(self, mouse):
        self.set_rectangle()

        if not self.selected and not self.rectangle.collidepoint(mouse):
            domain_image = pygame.image.load(self.image)
        else:
            domain_image = pygame.image.load(self.highlight_image)

        domain_image_scaled = pygame.transform.smoothscale(domain_image, (DOMAIN_SIZE, DOMAIN_SIZE))
        self.screen.blit(domain_image_scaled, self.rectangle)
        if self.substrate:
            self.draw_substrate_text()
        if self.subtype:
            self.draw_domain_subtype()



