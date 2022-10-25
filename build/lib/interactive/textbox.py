import pygame
from interactive.style import BLACK, BUTTON_PANEL_COLOUR, FONT, WHITE


class TextBox:

    def __init__(self, position, dimensions):
        self.x, self.y = position
        self.width, self.height = dimensions
        self.rectangle = pygame.Rect(position, dimensions)
        self.text = ''

        self.font_size = self.height - 10

        self.text_x = self.x + 2
        self.text_y = self.y + int((self.height - self.font_size) / 2)

        self.font = pygame.font.SysFont(FONT, self.font_size, bold=False)

    def add_character(self, character):
        self.text += character

    def remove_character(self):
        self.text = self.text[:-1]

    def draw(self, screen):
        pygame.draw.rect(screen, WHITE, self.rectangle)
        pygame.draw.rect(screen, BLACK, self.rectangle, 1)
        rendered_text = self.font.render(self.text, True, BLACK)
        screen.blit(rendered_text, (self.text_x, self.text_y))

    def erase(self, screen):
        pygame.draw.rect(screen, BUTTON_PANEL_COLOUR, self.rectangle)
        pygame.draw.rect(screen, BUTTON_PANEL_COLOUR, self.rectangle, 1)
