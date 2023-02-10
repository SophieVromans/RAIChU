from typing import List

from pikachu.drawing.drawing import Drawer

from raichu.reactions.chain_release import cyclic_release
from raichu.drawing.drawer import RaichuDrawer
from raichu.module import _Module

from raichu.drawing.bubbles import draw_bubbles

from raichu.tailoring_enzymes import TailoringEnzyme
from raichu.substrate import PKSSubstrate
from raichu.data.trans_at import TRANSATOR_CLADE_TO_TAILORING_REACTIONS, TRANSATOR_CLADE_TO_STARTER_SUBSTRATE, \
    TRANSATOR_CLADE_TO_ELONGATING
from raichu.domain.domain import TailoringDomain


class Cluster:
    def __init__(self, modules: List[_Module], tailoring_enzymes_representation = None) -> None:
        self.modules = modules
        self.tailoring_enzymes_representation = tailoring_enzymes_representation
        self.chain_intermediate = None

        self.structure_intermediates = []
        self.linear_product = None
        self.cyclised_product = None
        self.cyclised_products = []
        self.module_mechanisms = []
        self.tailoring_enzymes = []
        self.handle_transat()

    def handle_transat(self):
        for i, module in enumerate(self.modules):
            if module.type.name == "PKS" and module.subtype.name == "PKS_TRANS":
                substrate = PKSSubstrate("MALONYL_COA")
                if i < len(self.modules) - 1:
                    j = i + 1
                    next_module = self.modules[j]
                    if next_module.is_broken:
                        while next_module.is_broken:
                            next_module = self.modules[j]
                            j += 1
                    if next_module.type.name == "PKS":
                        if module.is_starter_module and next_module.subtype.name == "PKS_TRANS":
                            if next_module.synthesis_domain:
                                if next_module.synthesis_domain.subtype and next_module.synthesis_domain.subtype.name != "UNKNOWN":
                                    substrate_name = TRANSATOR_CLADE_TO_STARTER_SUBSTRATE.get(next_module.synthesis_domain.subtype.name)
                                    if substrate_name is not None:
                                        substrate = PKSSubstrate(substrate_name)
                                    else:
                                        substrate = PKSSubstrate("ACETYL_COA")
                                else:
                                    substrate = PKSSubstrate("ACETYL_COA")
                            else:
                                substrate = PKSSubstrate("ACETYL_COA")
                                    # TODO: Transfer names of substrates to a dictionary in raichu.substrate
                        elif module.is_starter_module:
                            substrate = PKSSubstrate("ACETYL_COA")
                        if next_module.synthesis_domain:
                            if not module.is_termination_module and next_module.subtype.name == "PKS_TRANS" and next_module.synthesis_domain.subtype and next_module.synthesis_domain.subtype.name != "UNKNOWN":
                                # Ignore tailoring domains in the module itself
                                module.tailoring_domains = []
                                for dummy_domain_type, dummy_domain_subtype in TRANSATOR_CLADE_TO_TAILORING_REACTIONS[next_module.synthesis_domain.subtype.name]:
                                    self.modules[i].tailoring_domains.append(TailoringDomain(dummy_domain_type,
                                                                                             dummy_domain_subtype))
                if module.synthesis_domain:
                    if module.synthesis_domain.subtype and module.synthesis_domain.subtype.name != "UNKNOWN":
                        self.modules[i].synthesis_domain.is_elongating = TRANSATOR_CLADE_TO_ELONGATING[module.synthesis_domain.subtype.name]
                self.modules[i].recognition_domain.substrate = substrate

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
        self.cyclised_product = cyclic_release(self.linear_product, atom)
        return self.cyclised_product

    def cyclise_all(self):
        pass

    def initialize_tailoring_enzymes_on_structure(self):
        if self.tailoring_enzymes_representation:
            for tailoring_enzyme_representation in self.tailoring_enzymes_representation:
                modification_sites = []
                for atoms_for_reaction in tailoring_enzyme_representation.modification_sites:
                    atoms_for_reaction_with_numbers = map(
                        lambda atom: [int(''.join(filter(str.isdigit, atom))), atom], atoms_for_reaction)
                    atoms_in_structure = list(map(
                        str, self.chain_intermediate.atoms.values()))
                    atoms_for_reaction_initialized = [self.chain_intermediate.atoms[atom[0]] if atom[1] in atoms_in_structure else print(
                        f"Non-existing atom for tailoring {str(atom[1])}. RAIChU will skip this tailoring reaction.") for atom in atoms_for_reaction_with_numbers]
                    atoms_for_reaction_initialized = list(
                        filter(lambda atom: atom is not None, atoms_for_reaction_initialized))
                    modification_sites += [atoms_for_reaction_initialized]
                self.tailoring_enzymes += [TailoringEnzyme(
                    tailoring_enzyme_representation.gene_name, tailoring_enzyme_representation.type, modification_sites, tailoring_enzyme_representation.substrate)]

    def do_tailoring(self):
        self.initialize_tailoring_enzymes_on_structure()
        for tailoring_enzyme in self.tailoring_enzymes:
            self.linear_product = tailoring_enzyme.do_tailoring(self.linear_product)
            self.chain_intermediate = self.linear_product

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

    def get_drawings(self, whitespace=30):

        drawings = []
        widths = []

        for i, structure in enumerate(self.structure_intermediates):

            drawing = RaichuDrawer(structure, dont_show=True)
            drawing.flip_y_axis()
            drawing.move_to_positive_coords()
            drawing.convert_to_int()

            carrier_domain_pos = None

            for atom in drawing.structure.graph:
                if atom.annotations.domain_type:
                    carrier_domain_pos = atom.draw.position
                    atom.draw.positioned = False

            assert carrier_domain_pos

            min_x = 100000000
            max_x = -100000000

            for atom in drawing.structure.graph:
                if atom.draw.positioned:
                    if atom.draw.position.x < min_x:
                        min_x = atom.draw.position.x
                    if atom.draw.position.x > max_x:
                        max_x = atom.draw.position.x

            width = (carrier_domain_pos.x - min_x + 0.5 * whitespace,
                     max_x - carrier_domain_pos.x + 0.5 * whitespace)
            widths.append(width)
            drawings.append(drawing)

        return drawings, widths

    def draw_cluster(self, as_string=True, out_file=None):
        drawings, widths = self.get_drawings()
        bubble_svg, bubble_positions, last_domain_coord = draw_bubbles(self, widths)
        min_x = 100000000
        max_x = -100000000
        min_y = 100000000
        max_y = -100000000

        svg_strings = []
        squiggly_svgs = []
        padding = None
        svg_style = None

        for i, drawing in enumerate(drawings):
            drawing.set_structure_id(f"s{i}")
            svg_style = drawing.svg_style

            padding = drawing.options.padding

            carrier_domain_pos = None

            for atom in drawing.structure.graph:
                if atom.annotations.domain_type:
                    carrier_domain_pos = atom.draw.position
                    atom.draw.positioned = False

            assert carrier_domain_pos
            bubble_x, bubble_y = bubble_positions[i]
            x_translation = bubble_x - carrier_domain_pos.x
            y_translation = bubble_y - carrier_domain_pos.y

            drawing.move_structure(x_translation, y_translation + 15)
            svg = drawing.draw_svg()
            svg_strings.append(svg)

            sulphur_pos = None
            carrier_pos = None

            for atom in drawing.structure.graph:
                if atom.draw.positioned:
                    if atom.draw.position.x < min_x:
                        min_x = atom.draw.position.x
                    if atom.draw.position.y < min_y:
                        min_y = atom.draw.position.y
                    if atom.draw.position.x > max_x:
                        max_x = atom.draw.position.x
                    if atom.draw.position.y > max_y:
                        max_y = atom.draw.position.y
                if atom.annotations.domain_type:
                    sulphur_pos = atom.get_neighbour('S').draw.position
                    carrier_pos = atom.draw.position

            squiggly_svg = f'<path d="M {sulphur_pos.x} {sulphur_pos.y - 5} Q {sulphur_pos.x - 5} {sulphur_pos.y - (sulphur_pos.y - 5 - carrier_pos.y)/2}, {carrier_pos.x} {sulphur_pos.y - 5 - (sulphur_pos.y - 5 - carrier_pos.y)/2} T {carrier_pos.x} {carrier_pos.y}" stroke="grey" fill="white"/>'
            squiggly_svgs.append(squiggly_svg)
        assert padding is not None

        x1 = 0
        x2 = max([max_x + padding, last_domain_coord + padding])
        y1 = 0
        y2 = max_y + padding

        width = x2
        height = y2

        svg_string = f"""<svg width="{width}" height="{height}" viewBox="{x1} {y1} {x2} {y2}" xmlns="http://www.w3.org/2000/svg">\n"""
        if svg_style:
            svg_string += f"{svg_style}\n"
        svg_string += bubble_svg
        for string in svg_strings:
            svg_string += string
        for squiggly_svg in squiggly_svgs:
            svg_string += squiggly_svg
        svg_string += "</svg>"

        if as_string:
            return svg_string
        else:
            if out_file is None:
                raise ValueError("Must provide output svg directory if 'as_string' is set to False.")
            else:
                with open(out_file, 'w') as svg_out:
                    svg_out.write(svg_string)

    def draw_product(self, as_string=True, out_file=None):
            assert self.chain_intermediate
            drawing = RaichuDrawer(self.chain_intermediate, dont_show=True, add_url=True, draw_Cs_in_pink=False)
            drawing.draw_structure()
            svg_string = drawing.save_svg_string()
            if as_string:
                return svg_string
            else:
                if out_file is None:
                    raise ValueError("Must provide output svg directory if 'as_string' is set to False.")
                else:
                    with open(out_file, 'w') as svg_out:
                        svg_out.write(svg_string)

class Mechanism:
    def __init__(self):
        self.structures = []

    def add_structure(self, structure):
        self.structures.append(structure.deepcopy())

    def draw_mechanism(self):
        pass
