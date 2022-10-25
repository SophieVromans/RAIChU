from typing import List

from pikachu.drawing.drawing import Drawer
from matplotlib import pyplot as plt

from raichu.reactions.chain_release import cyclic_release
from raichu.drawing.drawer import RaichuDrawer
from raichu.module import _Module
from raichu.drawing.bubbles import draw_bubbles


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
        pass

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
        svg = draw_bubbles(self)
        with open('testy.svg', 'w') as svg_out:
            svg_out.write(svg)


class Mechanism:
    def __init__(self):
        self.structures = []

    def add_structure(self, structure):
        self.structures.append(structure.deepcopy())

    def draw_mechanism(self):
        pass
