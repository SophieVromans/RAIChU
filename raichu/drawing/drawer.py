from pikachu.drawing.drawing import *
from pikachu.math_functions import *
from pikachu.smiles.smiles import Smiles

from raichu.central_chain_detection.find_central_chain import find_central_chain


class RaichuDrawer(Drawer):
    def __init__(self, structure, options=None, save_png=None, dont_show=False,
                 coords_only=True, dpi=100, save_svg=None):
        self.dont_show = dont_show
        self.dpi = dpi
        if options is None:
            self.options = Options()
        else:
            self.options = options
        if save_png is None:
            self.save_png = None
        else:
            # Check if filename is valid
            assert save_png.endswith('.png')
            self.save_png = save_png
        if save_svg is None:
            self.save_svg = None
        else:
            # Check if filename is valid
            assert save_svg.endswith('.svg')
            self.save_svg = save_svg
        super().__init__(structure, options=None, coords_only=coords_only)

    def center_on_carrier_domain(self):
        pass

    def export_to_json(self):
        json_structure = {"atoms": [],
                          "bonds": [],
                          "modules": []}

        ring_centers_x = []
        ring_centers_y = []

        for ring in self.rings:
            self.set_ring_center(ring)

            ring_centers_x.append(ring.center.x)
            ring_centers_y.append(ring.center.y)

        for bond_nr, bond in self.structure.bonds.items():
            if bond.atom_1.draw.positioned and bond.atom_2.draw.positioned:
                json_bond = {"id": bond.nr,
                             "atom_1": bond.atom_1.nr,
                             "atom_2": bond.atom_2.nr,
                             "type": bond.type,
                             "chiral": False,
                             "lines": []}

                if bond in self.chiral_bonds:

                    orientation, chiral_center = self.chiral_bond_to_orientation[bond]
                    json_bond["chiral"] = True
                    json_bond["wedge_orientation"] = orientation
                    json_bond["wedge_origin"] = chiral_center.nr
                elif (bond.atom_1.type == 'S' and bond.atom_2.annotations.domain_type or \
                      (bond.atom_2.type == 'S' and bond.atom_1.annotations.domain_type)):
                    json_bond["type"] = "squiggle"


        # If the starter unit contains an unknown moiety, number this R1
        for atom in self.structure.graph:
            if atom.type == '*' and not atom.annotations.unknown_index:
                atom.annotations.unknown_index = 1
        for atom in self.structure.graph:
            if atom.draw.positioned:
                json_atom = {"id": atom.nr,
                             "type": atom.type,
                             "text": '',
                             "atom_x": atom.draw.position.x,
                             "atom_y": atom.draw.position.y}

                if atom.type != 'C' and atom.draw.positioned:
                    text_h = ''
                    text_h_pos = None

                    if atom.type != 'C' or atom.draw.draw_explicit:
                        text = atom.type
                    else:
                        text = ''

                    if atom.annotations.domain_type:
                        text = atom.annotations.domain_type

                    horizontal_alignment = 'center'

                    orientation = self.get_hydrogen_text_orientation(atom)
                    if orientation == 'H_above_atom':
                        text_h_pos = Vector(atom.draw.position.x, atom.draw.position.y + 6)
                    if orientation == 'H_below_atom':
                        text_h_pos = Vector(atom.draw.position.x, atom.draw.position.y - 6)

                    atom_draw_position = Vector(atom.draw.position.x, atom.draw.position.y)

                    if atom.type == '*':
                        neighbouring_c = atom.get_neighbour('C')
                        assert neighbouring_c
                        # In order to not let the number of the sidechain overlap
                        # with the bond, move Rx symbol along hor. axis depending on
                        # if group is added to the left or right of the chain
                        delta_x_r = self.get_delta_x_sidechain(atom, neighbouring_c)
                        atom_draw_position.x += delta_x_r
                        text = fr'$R_{atom.annotations.unknown_index}$'

                    if not atom.charge and (atom.type != 'C' or atom.draw.draw_explicit):

                        if atom.draw.has_hydrogen:
                            hydrogen_count = 0
                            for neighbour in atom.neighbours:
                                if neighbour.type == 'H' and not neighbour.draw.is_drawn:
                                    hydrogen_count += 1

                            if hydrogen_count:

                                if hydrogen_count > 1:
                                    if orientation == 'H_before_atom':
                                        text = r'$H_{hydrogens}{atom_type}$'.format(hydrogens=hydrogen_count,
                                                                                    atom_type=atom.type)
                                        horizontal_alignment = 'right'
                                        atom_draw_position.x += 3
                                    elif orientation == 'H_below_atom' or orientation == 'H_above_atom':
                                        text = atom.type
                                        text_h = r'$H_{hydrogens}$'.format(hydrogens=hydrogen_count)

                                    else:
                                        text = r'${atom_type}H_{hydrogens}$'.format(hydrogens=hydrogen_count,
                                                                                    atom_type=atom.type)
                                        horizontal_alignment = 'left'
                                        atom_draw_position.x -= 3
                                elif hydrogen_count == 1:
                                    if orientation == 'H_before_atom':
                                        text = f'H{atom.type}'
                                        horizontal_alignment = 'right'
                                        atom_draw_position.x += 3
                                    elif orientation == 'H_below_atom' or orientation == 'H_above_atom':
                                        text = atom.type
                                        text_h = 'H'
                                    else:
                                        text = f'{atom.type}H'
                                        horizontal_alignment = 'left'
                                        atom_draw_position.x -= 3

                    elif atom.charge:
                        if atom.charge > 0:
                            charge_symbol = '+'
                        else:
                            charge_symbol = '-'

                        hydrogen_count = 0
                        for neighbour in atom.neighbours:
                            if neighbour.type == 'H' and not neighbour.draw.is_drawn:
                                hydrogen_count += 1

                        if not hydrogen_count:

                            if abs(atom.charge) > 1:

                                text = r'${atom_type}^{charge}{charge_symbol}$'.format(charge=atom.charge,
                                                                                       atom_type=atom.type,
                                                                                       charge_symbol=charge_symbol)
                            elif abs(atom.charge) == 1:
                                text = r'${atom_type}^{charge_symbol}$'.format(atom_type=atom.type,
                                                                               charge_symbol=charge_symbol)

                            horizontal_alignment = 'left'
                            atom_draw_position.x -= 3
                        else:

                            if hydrogen_count > 1:
                                if orientation == 'H_before_atom':
                                    if abs(atom.charge) > 1:
                                        text = r'$H_{hydrogens}{atom_type}^{charge}{charge_symbol}$'.format(
                                            hydrogens=hydrogen_count,
                                            atom_type=atom.type,
                                            charge=atom.charge,
                                            charge_symbol=charge_symbol)
                                    elif abs(atom.charge) == 1:
                                        text = r'$H_{hydrogens}{atom_type}^{charge_symbol}$'.format(
                                            hydrogens=hydrogen_count,
                                            atom_type=atom.type,
                                            charge_symbol=charge_symbol)

                                    horizontal_alignment = 'right'
                                    atom_draw_position.x += 3
                                elif orientation == 'H_above_atom' or orientation == 'H_below_atom':
                                    text_h = r'$H_{hydrogens}$'.format(hydrogens=hydrogen_count)
                                    if abs(atom.charge) > 1:
                                        text = r'${atom_type}^{charge}{charge_symbol}$'.format(atom_type=atom.type,
                                                                                               charge=atom.charge,
                                                                                               charge_symbol=charge_symbol)
                                    elif abs(atom.charge) == 1:
                                        text = r'${atom_type}^{charge_symbol}$'.format(atom_type=atom.type,
                                                                                       charge_symbol=charge_symbol)
                                else:
                                    if abs(atom.charge) > 1:
                                        text = r'${atom_type}H_{hydrogens}^{charge}{charge_symbol}$'.format(
                                            hydrogens=hydrogen_count,
                                            atom_type=atom.type,
                                            charge=atom.charge,
                                            charge_symbol=charge_symbol)
                                    elif abs(atom.charge) == 1:
                                        text = r'${atom_type}H_{hydrogens}^{charge_symbol}$'.format(
                                            hydrogens=hydrogen_count,
                                            atom_type=atom.type,
                                            charge_symbol=charge_symbol)

                                    horizontal_alignment = 'left'
                                    atom_draw_position.x -= 3
                            elif hydrogen_count == 1:
                                if orientation == 'H_before_atom':
                                    if abs(atom.charge) > 1:

                                        text = r'$H{atom_type}^{charge}{charge_symbol}$'.format(atom_type=atom.type,
                                                                                                charge=atom.charge,
                                                                                                charge_symbol=charge_symbol)
                                    elif abs(atom.charge) == 1:
                                        text = r'$H{atom_type}^{charge_symbol}$'.format(atom_type=atom.type,
                                                                                        charge_symbol=charge_symbol)
                                    horizontal_alignment = 'right'
                                    atom_draw_position.x += 3
                                elif orientation == 'H_above_atom' or orientation == 'H_below_atom':
                                    text_h = 'H'
                                    if abs(atom.charge) > 1:

                                        text = r'${atom_type}^{charge}{charge_symbol}$'.format(atom_type=atom.type,
                                                                                               charge=atom.charge,
                                                                                               charge_symbol=charge_symbol)
                                    elif abs(atom.charge) == 1:
                                        text = r'${atom_type}^{charge_symbol}$'.format(atom_type=atom.type,
                                                                                       charge_symbol=charge_symbol)

                                else:
                                    if abs(atom.charge) > 1:
                                        text = r'${atom_type}H^{charge}{charge_symbol}$'.format(atom_type=atom.type,
                                                                                                charge=atom.charge,
                                                                                                charge_symbol=charge_symbol)

                                    elif abs(atom.charge) == 1:
                                        text = r'${atom_type}H^{charge_symbol}$'.format(atom_type=atom.type,
                                                                                        charge_symbol=charge_symbol)
                                    horizontal_alignment = 'left'
                                    atom_draw_position.x -= 3

                    if text:
                        plt.text(atom_draw_position.x, atom_draw_position.y,
                                 text,
                                 horizontalalignment=horizontal_alignment,
                                 verticalalignment='center',
                                 color=atom.draw.colour)
                    if text_h:
                        plt.text(text_h_pos.x, text_h_pos.y,
                                 text_h,
                                 horizontalalignment='center',
                                 verticalalignment='center',
                                 color=atom.draw.colour)

        # If a png filename is included in the initialization of the
        # Raichu_drawer object, don't show the structure, but do save it as a
        # png image to the provided filename
        if self.dont_show:
            return self

        elif self.save_png is None and not self.dont_show and self.save_svg is None:
            plt.show()

        else:
            if self.save_png:
                plt.savefig(self.save_png)
            elif self.save_svg:
                plt.savefig(self.save_svg)
            plt.clf()
            plt.close()

    def finetune_overlap_resolution(self, masked_bonds=None, highest_atom=None):

        if not masked_bonds:
            masked_bonds = set()

        else:
            masked_bonds = set(masked_bonds)

        if self.total_overlap_score > self.options.overlap_sensitivity:
            clashing_atoms = self.find_clashing_atoms()

            best_bonds = []
            for atom_1, atom_2 in clashing_atoms:
                shortest_path = self.find_shortest_path(atom_1, atom_2)
                rotatable_bonds = []
                distances = []

                for i, bond in enumerate(shortest_path):
                    distance_1 = i
                    distance_2 = len(shortest_path) - i

                    average_distance = len(shortest_path) / 2

                    distance_metric = abs(average_distance - distance_1) + abs(average_distance - distance_2)

                    if self.bond_is_rotatable(bond) and bond not in masked_bonds:
                        rotatable_bonds.append(bond)
                        distances.append(distance_metric)

                best_bond = None
                optimal_distance = float('inf')
                for i, distance in enumerate(distances):
                    if distance < optimal_distance:
                        if highest_atom:
                            if not atom_1.draw.position.x > highest_atom.draw.position.x and atom_2.draw.position.x > highest_atom.draw.position.x:
                                best_bond = rotatable_bonds[i]
                                optimal_distance = distance
                        else:
                            best_bond = rotatable_bonds[i]
                            optimal_distance = distance

                if best_bond:
                    best_bonds.append(best_bond)

            best_bonds = list(set(best_bonds))

            for best_bond in best_bonds:
                if self.total_overlap_score > self.options.overlap_sensitivity:

                    subtree_size_1 = self.get_subgraph_size(best_bond.atom_1, {best_bond.atom_2})
                    subtree_size_2 = self.get_subgraph_size(best_bond.atom_2, {best_bond.atom_1})

                    if subtree_size_1 < subtree_size_2:
                        rotating_atom = best_bond.atom_1
                        parent_atom = best_bond.atom_2
                    else:
                        rotating_atom = best_bond.atom_2
                        parent_atom = best_bond.atom_1

                    overlap_score, _, _ = self.get_overlap_score()

                    scores = [overlap_score]

                    for i in range(12):
                        self.rotate_subtree(rotating_atom, parent_atom, math.radians(30), parent_atom.draw.position)
                        new_overlap_score, _, _ = self.get_overlap_score()
                        scores.append(new_overlap_score)

                    assert len(scores) == 13

                    scores = scores[:12]

                    best_i = 0
                    best_score = scores[0]

                    for i, score in enumerate(scores):
                        if score < best_score:
                            best_score = score
                            best_i = i

                    self.total_overlap_score = best_score

                    self.rotate_subtree(rotating_atom, parent_atom, math.radians(30 * best_i + 1), parent_atom.draw.position)

    def draw_structure(self):
        min_x = 100000000
        max_x = -100000000
        min_y = 100000000
        max_y = -100000000

        for atom in self.structure.graph:
            if atom.draw.positioned:
                if atom.draw.position.x < min_x:
                    min_x = atom.draw.position.x
                if atom.draw.position.y < min_y:
                    min_y = atom.draw.position.y
                if atom.draw.position.x > max_x:
                    max_x = atom.draw.position.x
                if atom.draw.position.y > max_y:
                    max_y = atom.draw.position.y

        height = max_y - min_y
        width = max_x - min_x
        self.line_width = 2

        fig, ax = plt.subplots(figsize=((width + 2 * self.options.padding) /
        50.0, (height + 2 * self.options.padding) / 50.0), dpi=self.dpi)

        ax.set_aspect('equal', adjustable='box')
        ax.axis('off')

        ax.set_xlim(
            [min_x - self.options.padding, max_x + self.options.padding])
        ax.set_ylim(
            [min_y - self.options.padding, max_y + self.options.padding])
        plt.subplots_adjust(left=0, bottom=0, right=1, top=1, wspace=0,
                            hspace=0)

        params = {'mathtext.default': 'regular', }
        plt.rcParams.update(params)

        ring_centers_x = []
        ring_centers_y = []

        for ring in self.rings:
            self.set_ring_center(ring)

            ring_centers_x.append(ring.center.x)
            ring_centers_y.append(ring.center.y)

        for bond_nr, bond in self.structure.bonds.items():
            if bond.atom_1.draw.positioned and bond.atom_2.draw.positioned:
                line = Line(bond.atom_1.draw.position,
                            bond.atom_2.draw.position, bond.atom_1,
                            bond.atom_2)
                midpoint = line.get_midpoint()
                truncated_line = line.get_truncated_line(
                    self.options.short_bond_length)
                if bond.type == 'single':
                    if bond in self.chiral_bonds:
                        orientation, chiral_center = self.chiral_bond_to_orientation[bond]
                        self.plot_chiral_bond(orientation, chiral_center, line,
                                              ax, midpoint)
                    else:
                        if (bond.atom_1.type == 'S' and
                                bond.atom_2.annotations.domain_type or
                                (bond.atom_2.type == 'S' and
                                 bond.atom_1.annotations.domain_type)):
                            self.plot_halflines_s_domain(line, ax, midpoint)
                        else:
                            self.plot_halflines(line, ax, midpoint)
                elif bond.type == 'double':
                    if not self.is_terminal(
                            bond.atom_1) and not self.is_terminal(bond.atom_2):
                        self.plot_halflines(line, ax, midpoint)

                        common_ring_numbers = self.get_common_rings(
                            bond.atom_1, bond.atom_2)

                        if common_ring_numbers:
                            common_rings = []
                            for ring_nr in common_ring_numbers:
                                common_rings.append(self.get_ring(ring_nr))

                            common_rings.sort(key=lambda x: len(x.members))
                            common_ring = common_rings[0]
                            ring_centre = common_ring.center
                            second_line = line.double_line_towards_center(
                                ring_centre, self.options.bond_spacing,
                                self.options.double_bond_length)
                            second_line_midpoint = second_line.get_midpoint()
                            self.plot_halflines_double(second_line, ax,
                                                       second_line_midpoint)

                        else:
                            bond_neighbours = bond.atom_1.drawn_neighbours +\
                                              bond.atom_2.drawn_neighbours
                            if bond_neighbours:
                                vectors = [atom.draw.position for atom in
                                           bond_neighbours]
                                gravitational_point = Vector.get_average(
                                    vectors)
                                second_line = line.double_line_towards_center(
                                    gravitational_point,
                                    self.options.bond_spacing,
                                    self.options.double_bond_length)
                                second_line_midpoint = second_line.get_midpoint()
                                self.plot_halflines_double(second_line, ax,
                                                           second_line_midpoint)
                            else:
                                print("Shouldn't happen!")
                    else:
                        if self.is_terminal(bond.atom_1) and self.is_terminal(
                                bond.atom_2):
                            dummy_1 = Vector(bond.atom_1.draw.position.x + 1,
                                             bond.atom_1.draw.position.y + 1)
                            dummy_2 = Vector(bond.atom_1.draw.position.x - 1,
                                             bond.atom_1.draw.position.y - 1)
                            double_bond_line_1 = line.double_line_towards_center(
                                dummy_1,
                                self.options.bond_spacing / 2.0,
                                self.options.double_bond_length)
                            double_bond_line_1_midpoint = double_bond_line_1.get_midpoint()
                            double_bond_line_2 = line.double_line_towards_center(
                                dummy_2,
                                self.options.bond_spacing / 2.0,
                                self.options.double_bond_length)
                            double_bond_line_2_midpoint = double_bond_line_2.get_midpoint()

                            self.plot_halflines_double(double_bond_line_1, ax,
                                                       double_bond_line_1_midpoint)
                            self.plot_halflines_double(double_bond_line_2, ax,
                                                       double_bond_line_2_midpoint)

                        else:

                            if self.is_terminal(bond.atom_1):
                                terminal_atom = bond.atom_1
                                branched_atom = bond.atom_2
                            else:
                                terminal_atom = bond.atom_2
                                branched_atom = bond.atom_1

                            if len(branched_atom.drawn_neighbours) >= 3:
                                closest_two = self.get_sorted_distances_from_list(
                                    terminal_atom,
                                    branched_atom.drawn_neighbours)
                                closest_atom_1 = closest_two[0][1]
                                closest_atom_2 = closest_two[1][1]

                                line = Line(terminal_atom.draw.position,
                                            branched_atom.draw.position,
                                            terminal_atom, branched_atom)

                                bond_1_line = Line(branched_atom.draw.position,
                                                   closest_atom_1.draw.position,
                                                   branched_atom,
                                                   closest_atom_1)
                                bond_2_line = Line(branched_atom.draw.position,
                                                   closest_atom_2.draw.position,
                                                   branched_atom,
                                                   closest_atom_1)

                                double_bond_line_1 = line.double_line_towards_center(
                                    closest_atom_1.draw.position,
                                    self.options.bond_spacing / 2.0,
                                    self.options.double_bond_length)
                                double_bond_line_2 = line.double_line_towards_center(
                                    closest_atom_2.draw.position,
                                    self.options.bond_spacing / 2.0,
                                    self.options.double_bond_length)

                                double_bond_line_1_midpoint = \
                                    double_bond_line_1.get_midpoint()
                                double_bond_line_2_midpoint = \
                                    double_bond_line_2.get_midpoint()

                                intersection_1 = double_bond_line_1.find_intersection(
                                    bond_1_line)
                                intersection_2 = double_bond_line_2.find_intersection(
                                    bond_2_line)

                                if terminal_atom.draw.position.x >\
                                        branched_atom.draw.position.x:
                                    double_bond_line_1.point_1 = intersection_1
                                    double_bond_line_2.point_1 = intersection_2
                                else:
                                    double_bond_line_1.point_2 = intersection_1
                                    double_bond_line_2.point_2 = intersection_2

                                self.plot_halflines(double_bond_line_1, ax,
                                                    double_bond_line_1_midpoint)
                                self.plot_halflines(double_bond_line_2, ax,
                                                    double_bond_line_2_midpoint)

                            else:
                                pass

        # If the starter unit contains an unknown moiety, number this R1
        for atom in self.structure.graph:
            if atom.type == '*' and not atom.annotations.unknown_index:
                atom.annotations.unknown_index = 1
        for atom in self.structure.graph:
            if atom.type != 'C' and atom.draw.positioned:
                text_h = ''
                text_h_pos = None

                if atom.type != 'C' or atom.draw.draw_explicit:
                    text = atom.type
                else:
                    text = ''

                if atom.annotations.domain_type:
                    text = atom.annotations.domain_type

                horizontal_alignment = 'center'

                orientation = self.get_hydrogen_text_orientation(atom)
                if orientation == 'H_above_atom':
                    text_h_pos = Vector(atom.draw.position.x, atom.draw.position.y + 6)
                if orientation == 'H_below_atom':
                    text_h_pos = Vector(atom.draw.position.x, atom.draw.position.y - 6)

                atom_draw_position = Vector(atom.draw.position.x, atom.draw.position.y)

                if atom.type == '*':
                    neighbouring_c = atom.get_neighbour('C')
                    assert neighbouring_c
                    # In order to not let the number of the sidechain overlap
                    # with the bond, move Rx symbol along hor. axis depending on
                    # if group is added to the left or right of the chain
                    delta_x_r = self.get_delta_x_sidechain(atom, neighbouring_c)
                    atom_draw_position.x += delta_x_r
                    text = fr'$R_{atom.annotations.unknown_index}$'

                if not atom.charge and (atom.type != 'C' or atom.draw.draw_explicit):

                    if atom.draw.has_hydrogen:
                        hydrogen_count = 0
                        for neighbour in atom.neighbours:
                            if neighbour.type == 'H' and not neighbour.draw.is_drawn:
                                hydrogen_count += 1

                        if hydrogen_count:

                            if hydrogen_count > 1:
                                if orientation == 'H_before_atom':
                                    text = r'$H_{hydrogens}{atom_type}$'.format(hydrogens=hydrogen_count,
                                                                                atom_type=atom.type)
                                    horizontal_alignment = 'right'
                                    atom_draw_position.x += 3
                                elif orientation == 'H_below_atom' or orientation == 'H_above_atom':
                                    text = atom.type
                                    text_h = r'$H_{hydrogens}$'.format(hydrogens=hydrogen_count)

                                else:
                                    text = r'${atom_type}H_{hydrogens}$'.format(hydrogens=hydrogen_count,
                                                                                atom_type=atom.type)
                                    horizontal_alignment = 'left'
                                    atom_draw_position.x -= 3
                            elif hydrogen_count == 1:
                                if orientation == 'H_before_atom':
                                    text = f'H{atom.type}'
                                    horizontal_alignment = 'right'
                                    atom_draw_position.x += 3
                                elif orientation == 'H_below_atom' or orientation == 'H_above_atom':
                                    text = atom.type
                                    text_h = 'H'
                                else:
                                    text = f'{atom.type}H'
                                    horizontal_alignment = 'left'
                                    atom_draw_position.x -= 3

                elif atom.charge:
                    if atom.charge > 0:
                        charge_symbol = '+'
                    else:
                        charge_symbol = '-'

                    hydrogen_count = 0
                    for neighbour in atom.neighbours:
                        if neighbour.type == 'H' and not neighbour.draw.is_drawn:
                            hydrogen_count += 1

                    if not hydrogen_count:

                        if abs(atom.charge) > 1:

                            text = r'${atom_type}^{charge}{charge_symbol}$'.format(charge=atom.charge,
                                                                                   atom_type=atom.type,
                                                                                   charge_symbol=charge_symbol)
                        elif abs(atom.charge) == 1:
                            text = r'${atom_type}^{charge_symbol}$'.format(atom_type=atom.type,
                                                                           charge_symbol=charge_symbol)

                        horizontal_alignment = 'left'
                        atom_draw_position.x -= 3
                    else:

                        if hydrogen_count > 1:
                            if orientation == 'H_before_atom':
                                if abs(atom.charge) > 1:
                                    text = r'$H_{hydrogens}{atom_type}^{charge}{charge_symbol}$'.format(
                                        hydrogens=hydrogen_count,
                                        atom_type=atom.type,
                                        charge=atom.charge,
                                        charge_symbol=charge_symbol)
                                elif abs(atom.charge) == 1:
                                    text = r'$H_{hydrogens}{atom_type}^{charge_symbol}$'.format(
                                        hydrogens=hydrogen_count,
                                        atom_type=atom.type,
                                        charge_symbol=charge_symbol)

                                horizontal_alignment = 'right'
                                atom_draw_position.x += 3
                            elif orientation == 'H_above_atom' or orientation == 'H_below_atom':
                                text_h = r'$H_{hydrogens}$'.format(hydrogens=hydrogen_count)
                                if abs(atom.charge) > 1:
                                    text = r'${atom_type}^{charge}{charge_symbol}$'.format(atom_type=atom.type,
                                                                                           charge=atom.charge,
                                                                                           charge_symbol=charge_symbol)
                                elif abs(atom.charge) == 1:
                                    text = r'${atom_type}^{charge_symbol}$'.format(atom_type=atom.type,
                                                                                   charge_symbol=charge_symbol)
                            else:
                                if abs(atom.charge) > 1:
                                    text = r'${atom_type}H_{hydrogens}^{charge}{charge_symbol}$'.format(
                                        hydrogens=hydrogen_count,
                                        atom_type=atom.type,
                                        charge=atom.charge,
                                        charge_symbol=charge_symbol)
                                elif abs(atom.charge) == 1:
                                    text = r'${atom_type}H_{hydrogens}^{charge_symbol}$'.format(
                                        hydrogens=hydrogen_count,
                                        atom_type=atom.type,
                                        charge_symbol=charge_symbol)

                                horizontal_alignment = 'left'
                                atom_draw_position.x -= 3
                        elif hydrogen_count == 1:
                            if orientation == 'H_before_atom':
                                if abs(atom.charge) > 1:

                                    text = r'$H{atom_type}^{charge}{charge_symbol}$'.format(atom_type=atom.type,
                                                                                            charge=atom.charge,
                                                                                            charge_symbol=charge_symbol)
                                elif abs(atom.charge) == 1:
                                    text = r'$H{atom_type}^{charge_symbol}$'.format(atom_type=atom.type,
                                                                                    charge_symbol=charge_symbol)
                                horizontal_alignment = 'right'
                                atom_draw_position.x += 3
                            elif orientation == 'H_above_atom' or orientation == 'H_below_atom':
                                text_h = 'H'
                                if abs(atom.charge) > 1:

                                    text = r'${atom_type}^{charge}{charge_symbol}$'.format(atom_type=atom.type,
                                                                                           charge=atom.charge,
                                                                                           charge_symbol=charge_symbol)
                                elif abs(atom.charge) == 1:
                                    text = r'${atom_type}^{charge_symbol}$'.format(atom_type=atom.type,
                                                                                   charge_symbol=charge_symbol)

                            else:
                                if abs(atom.charge) > 1:
                                    text = r'${atom_type}H^{charge}{charge_symbol}$'.format(atom_type=atom.type,
                                                                                            charge=atom.charge,
                                                                                            charge_symbol=charge_symbol)

                                elif abs(atom.charge) == 1:
                                    text = r'${atom_type}H^{charge_symbol}$'.format(atom_type=atom.type,
                                                                                    charge_symbol=charge_symbol)
                                horizontal_alignment = 'left'
                                atom_draw_position.x -= 3

                if text:
                    plt.text(atom_draw_position.x, atom_draw_position.y,
                             text,
                             horizontalalignment=horizontal_alignment,
                             verticalalignment='center',
                             color=atom.draw.colour)
                if text_h:
                    plt.text(text_h_pos.x, text_h_pos.y,
                             text_h,
                             horizontalalignment='center',
                             verticalalignment='center',
                             color=atom.draw.colour)

        # If a png filename is included in the initialization of the
        # Raichu_drawer object, don't show the structure, but do save it as a
        # png image to the provided filename
        if self.dont_show:
            return self

        elif self.save_png is None and not self.dont_show and self.save_svg is None:
            plt.show()

        else:
            if self.save_png:
                plt.savefig(self.save_png)
            elif self.save_svg:
                plt.savefig(self.save_svg)
            plt.clf()
            plt.close()


    def process_structure(self):
        self.position()
        self.structure.refresh_structure()
        self.restore_ring_information()
        self.resolve_primary_overlaps()
        self.total_overlap_score, sorted_overlap_scores, atom_to_scores = self.get_overlap_score()

        for i in range(self.options.overlap_resolution_iterations):
            for bond in self.drawn_bonds:
                if self.bond_is_rotatable(bond):

                    tree_depth_1 = self.get_subgraph_size(bond.atom_1,
                                                          {bond.atom_2})
                    tree_depth_2 = self.get_subgraph_size(bond.atom_2,
                                                          {bond.atom_1})

                    atom_1 = bond.atom_2
                    atom_2 = bond.atom_1

                    if tree_depth_1 > tree_depth_2:
                        atom_1 = bond.atom_1
                        atom_2 = bond.atom_2

                    subtree_overlap_score, _ = self.get_subtree_overlap_score(
                        atom_2, atom_1, atom_to_scores)

                    if subtree_overlap_score > self.options.overlap_sensitivity:
                        neighbours_2 = atom_2.drawn_neighbours[:]
                        neighbours_2.remove(atom_1)

                        if len(neighbours_2) == 1:
                            neighbour = neighbours_2[0]
                            angle = neighbour.draw.position.get_rotation_away_from_vector(
                                atom_1.draw.position, atom_2.draw.position,
                                math.radians(120))

                            self.rotate_subtree(neighbour, atom_2, angle,
                                                atom_2.draw.position)

                            new_overlap_score, _, _ = self.get_overlap_score()
                            if new_overlap_score > self.total_overlap_score:
                                self.rotate_subtree(neighbour, atom_2,
                                                    -angle,
                                                    atom_2.draw.position)
                            else:
                                self.total_overlap_score = new_overlap_score

                        elif len(neighbours_2) == 2:
                            if atom_2.draw.rings and atom_1.draw.rings:
                                continue

                            neighbour_1 = neighbours_2[0]
                            neighbour_2 = neighbours_2[1]
                            if len(neighbour_1.draw.rings) == 1 and len(
                                    neighbour_2.draw.rings) == 1:
                                # If the neighbours are in different rings,
                                # or in rings at all, do nothing
                                if neighbour_1.draw.rings[0] != \
                                        neighbour_2.draw.rings[0]:
                                    continue
                                elif neighbour_1.draw.rings or neighbour_2.draw.rings:
                                    continue
                                else:
                                    angle_1 = neighbour_1.draw.position.get_rotation_away_from_vector(
                                        atom_1.position, atom_2.position,
                                        math.radians(120))
                                    angle_2 = neighbour_2.draw.position.get_rotation_away_from_vector(
                                        atom_1.position, atom_2.position,
                                        math.radians(120))

                                    self.rotate_subtree(neighbour_1,
                                                        atom_2,
                                                        angle_1,
                                                        atom_2.position)
                                    self.rotate_subtree(neighbour_2,
                                                        atom_2,
                                                        angle_2,
                                                        atom_2.position)

                                    new_overlap_score, _, _ = self.get_overlap_score()

                                    if new_overlap_score > self.total_overlap_score:
                                        self.rotate_subtree(neighbour_1,
                                                            atom_2,
                                                            -angle_1,
                                                            atom_2.position)
                                        self.rotate_subtree(neighbour_2,
                                                            atom_2,
                                                            -angle_2,
                                                            atom_2.position)
                                    else:
                                        self.total_overlap_score = new_overlap_score
                        self.total_overlap_score, sorted_overlap_scores, atom_to_scores = self.get_overlap_score()


        # substructure search to see if structure is a polyketide or NRP

        is_polyketide = False
        is_nrp = False
        if self.structure.find_substructures(
                Smiles('SC(=O)').smiles_to_structure()):
            if self.structure.find_substructures(
                    Smiles('SC(CN)=O').smiles_to_structure()):
                is_nrp = True
            else:
                is_polyketide = True

        attached_to_domain = False
        for atom in self.structure.graph:
            if atom.type == 'S':
                sulphur = atom
                for neighbour in sulphur.neighbours:
                    if neighbour.annotations.domain_type:
                        attached_to_domain = True
                        domain = neighbour

        self.structure.refresh_structure()

        ### NRPS + PK code: Force pk/peptide backbone to be drawn straight:
        # If struct=PK/NRP, find central chain and attached domain
        if attached_to_domain and (is_nrp or is_polyketide):
            backbone_atoms = find_central_chain(self.structure)
            pcp = None
            for atom in self.structure.graph:
                if atom.annotations.domain_type:
                    pcp = atom

            assert pcp

            # Fix position of the S and C atom and the S-C angle
            sulphur = backbone_atoms[0]
            first_carbon = backbone_atoms[1]

            sulphur.draw.position.x = first_carbon.draw.position.x
            sulphur.draw.position.y = first_carbon.draw.position.y + 15

            angle = get_angle(first_carbon.draw.position,
                              sulphur.draw.position)
            angle_degrees = round(math.degrees(angle), 3)

            correct_angle_deg = -120
            delta_angle_deg = correct_angle_deg - angle_degrees
            delta_angle_rad = math.radians(delta_angle_deg)

            self.rotate_subtree(sulphur, first_carbon, delta_angle_rad,
                                first_carbon.draw.position)

            # Fix position domain straight above sulphur atom
            pcp.draw.position.x = sulphur.draw.position.x
            pcp.draw.position.y = sulphur.draw.position.y + 15

            # Rotate all other bonds in peptide backbone of NRP
            i = 0
            fixed_atoms = set()
            while i < (len(backbone_atoms) - 1):
                atom1 = backbone_atoms[i]
                atom2 = backbone_atoms[i + 1]
                angle = get_angle(atom1.draw.position, atom2.draw.position)
                angle_degrees = round(math.degrees(angle), 3)

                # Save angle last two atoms, needed later
                if i == (len(backbone_atoms) - 2):
                    last_angle_degrees = angle_degrees

                # Special case for amino acids with cyclic backbone (proline)
                if atom1.inside_ring and atom2.inside_ring:
                    if angle_degrees != 90.0:
                        correct_angle_deg = 90
                        delta_angle_deg = correct_angle_deg - angle_degrees
                        delta_angle_rad = math.radians(delta_angle_deg)
                        self.rotate_subtree(atom2, atom1, delta_angle_rad,
                                            atom1.draw.position)
                    i += 1

                    # Flip cyclic backbone amino acid if necessary
                    if atom1.draw.position.x > backbone_atoms[i - 2].draw.position.x:
                        first_atom_cycle = None
                        for next_atom in atom1.neighbours:
                            if next_atom.type != 'H' and next_atom not in backbone_atoms:
                                first_atom_cycle = next_atom
                        assert first_atom_cycle
                        if atom1.draw.position.x > first_atom_cycle.draw.position.x:
                            masked = {atom1, atom2}
                            for atom in self.traverse_substructure(
                                    first_atom_cycle, masked):
                                delta_x = 2 * (atom1.draw.position.x - atom.draw.position.x)
                                atom.draw.position.x += delta_x

                    if atom1.draw.position.x < backbone_atoms[i - 2].draw.position.x:
                        first_atom_cycle = None
                        for next_atom in atom1.neighbours:
                            if next_atom.type != 'H' and next_atom not in backbone_atoms:
                                first_atom_cycle = next_atom
                        assert first_atom_cycle
                        if atom1.draw.position.x < first_atom_cycle.draw.position.x:
                            masked = {atom1, atom2}
                            for atom in self.traverse_substructure(first_atom_cycle, masked):
                                delta_x = 2 * (atom.draw.position.x - atom1.draw.position.x)
                                atom.draw.position.x -= delta_x

                # Fix bond angle backbone atoms if second backbone atom (one
                # with larger position.x) is inside ring, and the first is not
                elif atom2 != backbone_atoms[-1] and atom2.inside_ring and not atom1.inside_ring and \
                        backbone_atoms[i + 2].inside_ring:
                    if atom2 not in fixed_atoms:
                        if round(math.degrees(get_angle(
                                backbone_atoms[i - 1].draw.position,
                                backbone_atoms[i].draw.position)),
                                 3) == 120.0:
                            correct_angle_deg = 60.0
                            first_angle_cyclic = 60.0
                        elif round(math.degrees(get_angle(
                                backbone_atoms[i - 1].draw.position,
                                backbone_atoms[i].draw.position)),
                                   3) == 60.0:
                            correct_angle_deg = 120.0
                            first_angle_cyclic = 120.0
                        else:
                            print('should not happen!!!')
                        delta_angle_deg = correct_angle_deg - angle_degrees
                        delta_angle_rad = math.radians(delta_angle_deg)
                        self.rotate_subtree(atom2, atom1, delta_angle_rad,
                                            atom1.draw.position)

                        fixed_atoms.add(atom2)
                        i = 0
                    else:
                        if angle_degrees == 120.0:
                            first_angle_cyclic = 120.0
                        elif angle_degrees == 60.0:
                            first_angle_cyclic = 60.0
                        i += 1

                # Other way around... (first one in ring, after cycle)
                elif atom1.inside_ring and not atom2.inside_ring and \
                        backbone_atoms[i - 1].inside_ring:
                    if first_angle_cyclic == 60.0:
                        correct_angle_deg = 120.0
                    elif first_angle_cyclic == 120.0:
                        correct_angle_deg = 60.0
                    if angle_degrees != correct_angle_deg:
                        delta_angle_deg = correct_angle_deg - angle_degrees
                        delta_angle_rad = math.radians(delta_angle_deg)
                        self.rotate_subtree(atom2, atom1, delta_angle_rad,
                                            atom1.draw.position)
                        i = 0
                    else:
                        i += 1
                else:
                    if angle_degrees != 120.0 and angle_degrees != 60:
                        if round(math.degrees(get_angle(
                                backbone_atoms[i - 1].draw.position,
                                backbone_atoms[i].draw.position)),
                                 3) == 120.0:
                            correct_angle_deg = 60.0
                        elif round(math.degrees(get_angle(
                                backbone_atoms[i - 1].draw.position,
                                backbone_atoms[i].draw.position)),
                                   3) == 60.0:
                            correct_angle_deg = 120.0
                        delta_angle_deg = correct_angle_deg - angle_degrees
                        delta_angle_rad = math.radians(delta_angle_deg)
                        self.rotate_subtree(atom2, atom1, delta_angle_rad,
                                            atom1.draw.position)
                        i = 0
                    else:
                        if i >= 1:
                            if angle_degrees == 120.0:
                                if round(math.degrees(get_angle(
                                        backbone_atoms[
                                            i - 1].draw.position,
                                        backbone_atoms[i].draw.position)),
                                         3) == 120:
                                    correct_angle_deg = 60.0
                                    delta_angle_deg = correct_angle_deg - angle_degrees
                                    delta_angle_rad = math.radians(
                                        delta_angle_deg)
                                    self.rotate_subtree(atom2, atom1,
                                                        delta_angle_rad,
                                                        atom1.draw.position)
                                    i = 0
                                else:
                                    i += 1
                            elif angle_degrees == 60.0:
                                if round(math.degrees(get_angle(
                                        backbone_atoms[
                                            i - 1].draw.position,
                                        backbone_atoms[i].draw.position)),
                                         3) == 60:
                                    correct_angle_deg = 120.0
                                    delta_angle_deg = correct_angle_deg - angle_degrees
                                    delta_angle_rad = math.radians(
                                        delta_angle_deg)
                                    self.rotate_subtree(atom2, atom1,
                                                        delta_angle_rad,
                                                        atom1.draw.position)
                                    i = 0
                                else:
                                    i += 1
                        else:
                            i += 1

            # Force pk/amino acid sidechains to stick out straight from each side
            i = 1
            while i < (len(backbone_atoms)):
                terminal_carboxylic_acid = False
                atom = backbone_atoms[i]
                atom_neighbours = []
                atom_neighbour_types = []
                connected_to_sidechain = False
                for neighbour in atom.neighbours:
                    atom_neighbours.append(neighbour)
                    atom_neighbour_types.append(neighbour.type)

                    if neighbour not in backbone_atoms and \
                            neighbour.type != 'H' and neighbour.type != 'S' and\
                            neighbour.type == 'O' and \
                            self.structure.bond_lookup[neighbour][atom].type == 'double':
                        first_atom_sidechain = neighbour
                        connected_to_sidechain = True
                        if backbone_atoms[i - 1].draw.position.x < \
                                backbone_atoms[i].draw.position.x:
                            sidechain_orientation = 'right'
                        elif backbone_atoms[i - 1].draw.position.x > \
                                backbone_atoms[i].draw.position.x:
                            sidechain_orientation = 'left'
                    elif neighbour not in backbone_atoms and \
                            neighbour.type != 'H' and neighbour.type != 'S' and\
                            not neighbour.inside_ring:
                        first_atom_sidechain = neighbour
                        connected_to_sidechain = True
                        if backbone_atoms[i - 1].draw.position.x < \
                                backbone_atoms[i].draw.position.x:
                            sidechain_orientation = 'right'
                        elif backbone_atoms[i - 1].draw.position.x > \
                                backbone_atoms[i].draw.position.x:
                            sidechain_orientation = 'left'

                    elif neighbour not in backbone_atoms and \
                            neighbour.type != 'H' and neighbour.type != 'S' and \
                            neighbour.inside_ring and any(
                            bond.type == 'double' for bond in
                            neighbour.bonds) and atom == backbone_atoms[-1]:
                        first_atom_sidechain = neighbour
                        connected_to_sidechain = True
                        if backbone_atoms[i - 1].draw.position.x < \
                                backbone_atoms[i].draw.position.x:
                            sidechain_orientation = 'diagonal_right'
                        elif backbone_atoms[i - 1].draw.position.x > \
                                backbone_atoms[i].draw.position.x:
                            sidechain_orientation = 'diagonal_left'

                    elif neighbour not in backbone_atoms and \
                            neighbour.type != 'H' and neighbour.type != 'S' and \
                            neighbour.in_ring(self.structure) and not atom.in_ring(self.structure):
                        first_atom_sidechain = neighbour
                        connected_to_sidechain = True
                        if backbone_atoms[i - 1].draw.position.x < \
                                backbone_atoms[i].draw.position.x:
                            sidechain_orientation = 'right'
                        elif backbone_atoms[i - 1].draw.position.x > \
                                backbone_atoms[i].draw.position.x:
                            sidechain_orientation = 'left'

                    # If bond is directly connected to sidechain, check&rotate
                    if connected_to_sidechain:
                        angle = get_angle(atom.draw.position,
                                          first_atom_sidechain.draw.position)
                        angle_degrees = round(math.degrees(angle), 3)
                        # Check if it is a terminal carboxylic acid group
                        if atom_neighbour_types.count('O') == 2 and len(
                                atom_neighbour_types) == 3:
                            for next_atom in atom_neighbours:
                                if next_atom.type == 'O':
                                    if \
                                    self.structure.bond_lookup[next_atom][atom].type == 'single':
                                        hydroxyl = next_atom
                                        angle2 = get_angle(
                                            atom.draw.position,
                                            hydroxyl.draw.position)
                                        angle_degrees2 = round(
                                            math.degrees(angle2),
                                            3)
                                        if last_angle_degrees == 120.0:
                                            correct_angle_deg = 60.0
                                        elif last_angle_degrees == 60.0:
                                            correct_angle_deg = 120.0
                                        delta_angle_deg = correct_angle_deg - angle_degrees2
                                        delta_angle_rad = math.radians(
                                            delta_angle_deg)
                                        self.rotate_subtree(hydroxyl,
                                                            atom,
                                                            delta_angle_rad,
                                                            atom.draw.position)
                                    elif \
                                    self.structure.bond_lookup[next_atom][
                                        atom].type == 'double':
                                        carbonyl = next_atom
                                        angle2 = get_angle(
                                            atom.draw.position,
                                            carbonyl.draw.position)
                                        angle_degrees2 = round(
                                            math.degrees(angle2),
                                            3)
                                        correct_angle_deg = 0
                                        delta_angle_deg = correct_angle_deg - angle_degrees2
                                        delta_angle_rad = math.radians(
                                            delta_angle_deg)
                                        self.rotate_subtree(carbonyl, atom,
                                                            delta_angle_rad,
                                                            atom.draw.position)
                                    terminal_carboxylic_acid = True
                        if not terminal_carboxylic_acid:
                            if sidechain_orientation == 'right':
                                correct_angle_deg = 180
                            elif sidechain_orientation == 'left':
                                correct_angle_deg = 0
                            elif sidechain_orientation == 'diagonal_left':
                                correct_angle_deg = 120
                            elif sidechain_orientation == 'diagonal_right':
                                correct_angle_deg = 60
                            delta_angle_deg = correct_angle_deg - angle_degrees
                            delta_angle_rad = math.radians(delta_angle_deg)
                            self.rotate_subtree(first_atom_sidechain, atom,
                                                delta_angle_rad,
                                                atom.draw.position)
                i += 1

            # If the drawer rotated the entire structure, correct this
            angle = get_angle(pcp.draw.position,
                              sulphur.draw.position)
            angle_degrees = round(math.degrees(angle), 3)
            if angle_degrees != 90.0:
                correct_angle_deg = 90
                delta_angle_deg = correct_angle_deg - angle_degrees
                delta_angle_rad = math.radians(delta_angle_deg)
                self.rotate_subtree(sulphur, pcp, delta_angle_rad,
                                    pcp.draw.position)

            # Fix rotation bulky sidechains so carboxyl groups dont need to move
            i = 1
            while i < (len(backbone_atoms)):
                atom = backbone_atoms[i]
                atom_neighbours = []
                atom_neighbour_types = []
                for neighbour in atom.neighbours:
                    atom_neighbours.append(neighbour)
                    atom_neighbour_types.append(neighbour.type)
                    if neighbour not in backbone_atoms and \
                            neighbour.type != 'H' and neighbour.type != 'S' and\
                            not neighbour.inside_ring:
                        first_atom_sidechain = neighbour
                        for further_atom in first_atom_sidechain.neighbours:
                            types = []
                            for further_atom_neighbour in further_atom.neighbours:
                                types.append(further_atom_neighbour.type)
                            if further_atom.type == 'C' and further_atom not in backbone_atoms and (
                                    (types.count(
                                            'C') == 3) or (
                                            types.count('O') == 2 and
                                            types.count('H') == 0) or
                                    (types.count('N') == 1 and
                                     types.count('O') == 1 and
                                     types.count('H') == 0) or
                                    (types.count('C') == 2 and
                                     types.count('N') == 1 and
                                     types.count('H') == 0)):
                                angle_bulky_sidechain = get_angle(
                                    first_atom_sidechain.draw.position,
                                    further_atom.draw.position)
                                angle_bulky_sidechain = round(
                                    math.degrees(angle_bulky_sidechain), 3)
                                if atom.draw.position.x < first_atom_sidechain.draw.position.x:
                                    correct_angle_deg = 160
                                elif further_atom.draw.position.x < first_atom_sidechain.draw.position.x:
                                    correct_angle_deg = 20
                                delta_angle_deg = correct_angle_deg - angle_bulky_sidechain
                                delta_angle_rad = math.radians(
                                    delta_angle_deg)
                                self.rotate_subtree(further_atom,
                                                    first_atom_sidechain,
                                                    delta_angle_rad,
                                                    first_atom_sidechain.draw.position)
                    # Fix incorrect rotation carbonyl groups
                    if neighbour.type == 'O' and \
                            self.structure.bond_lookup[atom][
                                neighbour].type == 'double':
                        correctly_rotated = True
                        if neighbour.draw.position.x > atom.draw.position.x and atom.draw.position.x < \
                                backbone_atoms[i - 1].draw.position.x:
                            correct_angle_deg = 0
                            correctly_rotated = False
                        elif neighbour.draw.position.y != atom.draw.position.y and atom.draw.position.x < \
                                backbone_atoms[i - 1].draw.position.x:
                            correct_angle_deg = 0
                            correctly_rotated = False
                        elif neighbour.draw.position.x < atom.draw.position.x and atom.draw.position.x > \
                                backbone_atoms[i - 1].draw.position.x:
                            correct_angle_deg = 180
                            correctly_rotated = False
                        elif neighbour.draw.position.y != atom.draw.position.y and atom.draw.position.x > \
                                backbone_atoms[i - 1].draw.position.x:
                            correct_angle_deg = 180
                            correctly_rotated = False
                        if not correctly_rotated:
                            angle = get_angle(atom.draw.position,
                                              neighbour.draw.position)
                            angle_degrees = round(math.degrees(angle), 3)
                            delta_angle_deg = correct_angle_deg - angle_degrees
                            delta_angle_rad = math.radians(delta_angle_deg)
                            self.rotate_subtree(neighbour,
                                                atom, delta_angle_rad,
                                                atom.draw.position)
                i += 1

            # self.resolve_primary_overlaps()
            self.total_overlap_score, sorted_overlap_scores, atom_to_scores = self.get_overlap_score()
            central_chain_bonds = set()
            for bond in self.structure.bonds.values():
                if bond.atom_1.annotations.in_central_chain or\
                        bond.atom_2.annotations.in_central_chain:
                    central_chain_bonds.add(bond)

            self.finetune_overlap_resolution(masked_bonds=central_chain_bonds, highest_atom=sulphur)
            self.resolve_secondary_overlaps(sorted_overlap_scores)

    # End NRPS rotation code

    def plot_halflines_s_domain(self, line, ax, midpoint):
        truncated_line = line.get_truncated_line(self.options.short_bond_length)
        self.plot_line_dashed(truncated_line, ax, color='#a6a6a6')

    def plot_line_dashed(self, line, ax, color='grey'):
        with matplotlib.rc_context({'path.sketch': (5, 10, 1)}):
            ax.plot([line.point_1.x, line.point_2.x],
                    [line.point_1.y, line.point_2.y], color=color,
                    linewidth=self.line_width/1.5)

    def get_delta_x_sidechain(self, sidechain_atom, neighbouring_carbon):
        """Returns the change tot the x-position of the sidegroup symbol in
        the structure drawing as float
        sidechain_atom: Atom object of the unknown sidechain group
        neighoburin_carbon: Atom object of the carbon in the central chain
        bound to the unknown sidechain group
        """
        sidechainatom_x = sidechain_atom.draw.position.x
        carbon_x = neighbouring_carbon.draw.position.x
        if (carbon_x - sidechainatom_x) < -10:
            delta_x = 1.5
        elif (carbon_x - sidechainatom_x) > 10:
            delta_x = -1.5
        else:
            delta_x = 0.0

        return delta_x

def get_angle(vector1, vector2):
    difference = Vector.subtract_vectors(vector1, vector2)
    difference_angle = difference.angle()
    return difference_angle


