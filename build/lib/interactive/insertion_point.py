import pygame
from interactive.style import BLACK, MODULE_PADDING, GENE_LABEL_SIZE


class InsertionPoint:
    def __init__(self, gene, mouse):
        self.gene = gene
        self.insertion_point, self.insertion_x = self.get_insertion_point(mouse)

    def get_insertion_point(self, mouse):

        insertion_point = 0
        insertion_x = self.gene.x + GENE_LABEL_SIZE + MODULE_PADDING / 2

        for i, module in enumerate(self.gene.modules):
            if module.x < mouse[0]:
                insertion_point = i + 1
                insertion_x = module.x + module.width + MODULE_PADDING / 2
            else:
                break

        return insertion_point, insertion_x

    def draw(self, screen):
        length = self.gene.height
        y_start = self.gene.y
        y_end = y_start + length
        pygame.draw.line(screen, BLACK, (self.insertion_x, y_start), (self.insertion_x, y_end - 1))
