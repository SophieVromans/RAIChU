import pygame

from interactive.style import BUTTON_PANEL_COLOUR, BUTTON_PANEL, BLACK, SIZE, RENDER_WINDOW_SIZE


class ButtonPanel:
    def __init__(self, screen):
        self.screen = screen

    def draw(self):
        # pygame.draw.rect(self.screen, BUTTON_PANEL_COLOUR, BUTTON_PANEL)
        pygame.draw.line(self.screen, BLACK, (BUTTON_PANEL.x, BUTTON_PANEL.y),
                         (BUTTON_PANEL.x + BUTTON_PANEL.width, BUTTON_PANEL.y))


class ProductRaster:
    def __init__(self, screen):
        self.screen = screen
        self.rectangles = self.get_rectangles()

    def get_rectangles(self):
        panel_width = RENDER_WINDOW_SIZE[0] / 3
        panel_height = RENDER_WINDOW_SIZE[1] / 2
        x_1 = SIZE[0] - RENDER_WINDOW_SIZE[0]
        x_2 = x_1 + panel_width
        x_3 = x_2 + panel_width
        y_1 = 0
        y_2 = panel_height

        rectangles = [pygame.Rect(x_1, y_1, panel_width, panel_height),
                      pygame.Rect(x_2, y_1, panel_width, panel_height),
                      pygame.Rect(x_3, y_1, panel_width, panel_height),
                      pygame.Rect(x_1, y_2, panel_width, panel_height),
                      pygame.Rect(x_2, y_2, panel_width, panel_height),
                      pygame.Rect(x_3, y_2, panel_width, panel_height)]

        return rectangles

    def draw_images(self, images):
        for i, image in enumerate(images):
            if i <= 5:
                center = self.rectangles[i].center
                rectangle = image.get_rect(center=center)
                self.screen.blit(image, rectangle)
                pygame.draw.rect(self.screen, BLACK, rectangle, 2)




