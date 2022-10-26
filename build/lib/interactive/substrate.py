import os
import pygame

import interactive.flatfiles
import interactive.images.nrps_substrates
import interactive.images.pks_substrates
import interactive.images.pks_starter_substrates

from pikachu.general import png_from_smiles

from interactive.style import SUBSTRATE_BUTTON_SIZE, SUBSTRATE_BUTTON_PADDING, \
    SUBSTRATE_BUTTONS_PER_LINE, WIDTH, HEIGHT, FATTY_ACID_IMAGE_SIZE

SUBSTRATE_DIR = os.path.dirname(interactive.flatfiles.__file__)
NRPS_SUBSTRATE_IMG = os.path.dirname(interactive.images.nrps_substrates.__file__)
PKS_SUBSTRATE_IMG = os.path.dirname(interactive.images.pks_substrates.__file__)
PKS_STARTER_SUBSTRATE_IMG = os.path.dirname(interactive.images.pks_starter_substrates.__file__)
SUBSTRATE_TO_ABBREVIATION_FILE = os.path.join(SUBSTRATE_DIR, 'substrate_to_abbreviation.txt')

PROTEINOGENIC_SUBSTRATES = {"Ala": ['Alanine',
                                    'D-Alanine',
                                    'dehydroalanine',
                                    'beta-Alanine',
                                    '2-methyl-beta-alanine'],
                            "Arg": ["Arginine",
                                    "D-Arginine",
                                    "4,5-dehydroarginine",
                                    "homoarginine"],
                            "Asn": ["Asparagine",
                                    "D-Asparagine"],
                            "Asp": ["AsparticAcid",
                                    "D-AsparticAcid"],
                            "Cys": ["Cysteine",
                                    "D-Cysteine",
                                    "dehydro-cysteine",
                                    "cysteicacid",
                                    "D-cysteicacid",
                                    "alpha-methylcysteine",
                                    "D-alpha-methylcysteine"],
                            "Gln": ["Glutamine",
                                    "D-Glutamine",
                                    "3-Hydroxyglutamine"],
                            "Glu": ["GlutamicAcid",
                                    "D-GlutamicAcid",
                                    "4-methylglutamicacid"],
                            "Gly": ["Glycine"],
                            "His": ["Histidine",
                                    "D-Histidine"],
                            "Ile": ["Isoleucine",
                                    "D-Isoleucine",
                                    "allo-Isoleucine",
                                    "D-allo-Isoleucine",
                                    "homoisoleucine"],
                            "Leu": ["Leucine",
                                    "D-Leucine",
                                    "tert-Leu",
                                    "D-tert-Leucine"],
                            "Lys": ["Lysine",
                                    "D-Lysine",
                                    "beta-lysine"],
                            "Met": ["Methionine",
                                    "D-Methionine"],
                            "Phe": ["Phenylalanine",
                                    "D-Phenylalanine",
                                    "Homophenylalanine",
                                    "N-acetylphenylalanine",
                                    "alpha-amino-phenyl-valericacid",
                                    "D-alpha-amino-phenyl-valericacid"],
                            "Pro": ["Proline",
                                    "D-Proline",
                                    "4-methylproline"],
                            "Ser": ["Serine",
                                    "D-Serine",
                                    "Homoserine",
                                    "D-Homoserine",
                                    "2-Methylserine"],
                            "Thr": ["Threonine",
                                    "D-Threonine",
                                    "allo-Threonine",
                                    "D-allo-Threonine",
                                    "4-butenyl-4-methylthreonine"],
                            "Trp": ["Tryptophan",
                                    "D-Tryptophan",
                                    "2,3-Dehydro-Tryptophan"],
                            "Tyr": ["Tyrosine",
                                    "D-Tyrosine",
                                    "dehydrotyrosine",
                                    "(R)beta-tyrosine",
                                    "Homotyrosine",
                                    "D-Homotyrosine",
                                    "3-methyltyrosine",
                                    "4,5-Dihydroxyhomotyrosine",
                                    "3-hydroxy-4-O-methyl-5-methyltyrosine"],
                            "Val": ["Valine",
                                    "D-Valine",
                                    "4-Hydroxyvaline",
                                    "isovaline",
                                    "D-isovaline",
                                    "alpha-ethylnorvaline",
                                    "D-alpha-ethylnorvaline",
                                    "Norvaline"]}

NONPROTEINOGENIC_SUBSTRATES = {"Aad": ["2-Aminoadipicacid",
                                       "D-2-Aminoadipicacid"],
                               "Abu": ["2-Aminobutyricacid",
                                       "D-2-Aminobutyricacid",
                                       "2-Aminoisobutyricacid",
                                       "Dehydrobutyrine"],
                               "Cit": ["Citrulline",
                                       "D-Citrulline"],
                               "Cor": ["coronamicacid",
                                       "R-coronamicacid",
                                       "norcoronamicacid",
                                       "R-norcoronamicacid"],
                               "Dab": ["2,4-diaminobutyricacid",
                                       "D-2,4-diaminobutyricacid",
                                       "(2S)-2,3-diaminobutyricacid",
                                       "(2S,3S)-2,3-diaminobutyricacid",
                                       "(2R,3R)-2,3-diaminobutyricacid",
                                       "2,3-diamino-3-methyl-propanoicacid",
                                       "D-2,3-Diaminopropionicacid"],
                               "End": ["enduracididine",
                                       "D-enduracididine"],
                               "Kyn": ["kynurenine",
                                       "D-kynurenine"],
                               "Orn": ["Ornithine",
                                       "D-Ornithine",
                                       "4,5-Dihydroxyornithine"],
                               "PheGly": ["phenylglycine",
                                          "D-PhenylGlycine",
                                          "HydroxyPhenylGlycine",
                                          "D-HydroxyPhenylGlycine",
                                          "3,5-dihydroxyphenylglycine",
                                          "D-3,5-dihydroxyphenylglycine"],
                               "Pip": ["Pipecolicacid",
                                       "D-Pipecolicacid",
                                       "L-piperazicacid"],
                               "Other": ["Anticapsin",
                                         "Beta-hydroxycyclohex-2-enylalanine",
                                         "capreomycidine",
                                         "phosphinothricin",
                                         "D-phosphinothricin",
                                         "(2S,3R)-2-amino-3-hydroxy-4-(4-nitrophenyl)butanoate",
                                         "3-(3-Pyridyl)-Alanine",
                                         "3-[(1R,2R)-2-Nitrocyclopropyl]-L-alanine"]}

NON_AMINO_ACIDS = {'Aromatic': ["2-carboxyquinoxaline",
                                "Quinoline-2-carboxylicacid",
                                "2,3-dihydroxybenzoicacid",
                                "4-Hydroxyindole-3-carboxylicacid",
                                "dihydroxyphenylthiazolgroup",
                                "D-dihydroxyphenylthiazolgroup",
                                "D-lysergicacid",
                                "hydrocinnamicacid",
                                "hydroxypicolinicacid",
                                "phenylaceticacid",
                                "salicylicacid"],
                   # 'Cyclic': ["cycloOrnithine",
                   #            "D-cycloOrnithine",
                   #            "homoserinelactone",
                   #            "D-homoserinelactone"],
                   'Other': ["2-hydroxyisovalerate",
                             "D-2-hydroxyisovalerate",
                             # "3-aminopropanimidamide",
                             # "Alaninol",
                             "alpha-ketoisocaproicacid",
                             # "alpha-ketoisovalericacid",
                             "Glycocyamine",
                             "Glycolicacid"]
                             # "valinol"]
                   }

# FATTY_ACIDS = {'Amino group': ['2-Amino-9,10-epoxy-8-oxodecanoidacid',
#                                'D-methyl-2-aminooctanoicacid',
#                                'methyl-2-aminooctanoicacid'],
#
#                'No amino group': [],
#                'Custom': []}

FATTY_ACIDS = {'Amino group': ['methyl-2-aminooctanoicacid',
                               'D-methyl-2-aminooctanoicacid'],
               'Custom': []}


def parse_substrate_to_abbr():
    substrate_to_abbreviation = {}
    with open(SUBSTRATE_TO_ABBREVIATION_FILE, 'r') as abbreviation_file:
        for line in abbreviation_file:
            line = line.strip()
            if line:
                substrate, abbreviation = line.split()
                substrate_to_abbreviation[substrate] = abbreviation

    return substrate_to_abbreviation


SUBSTRATE_TO_ABBREVIATION = parse_substrate_to_abbr()


class SubstrateGroup:
    def __init__(self, name, substrates):
        self.name = name
        self.substrates = substrates
        self.coords = []
        self.create_raster()
        self.set_rectangles()

    def create_raster(self):

        x_coord = int(0.23 * WIDTH)
        y_coord = int(0.52 * HEIGHT)
        column_nr = 0

        for i, substrate in enumerate(self.substrates):
            if i % SUBSTRATE_BUTTONS_PER_LINE == 0:
                x_coord = int(0.23 * WIDTH)
                column_nr += 1
                if i != 0 and column_nr <= 2:
                    y_coord += SUBSTRATE_BUTTON_SIZE + SUBSTRATE_BUTTON_PADDING
                elif i != 0:
                    x_coord = int(0.23 * WIDTH - SUBSTRATE_BUTTON_SIZE - SUBSTRATE_BUTTON_PADDING)

            else:
                x_coord += SUBSTRATE_BUTTON_SIZE + SUBSTRATE_BUTTON_PADDING

            self.coords.append((x_coord, y_coord))

    def set_rectangles(self):

        for i, substrate in enumerate(self.substrates):
            x, y = self.coords[i]
            substrate.x = x
            substrate.y = y
            substrate.rectangle = pygame.Rect(x, y, SUBSTRATE_BUTTON_SIZE, SUBSTRATE_BUTTON_SIZE)


class Substrate:
    def __init__(self, name, smiles, module_type, set_images=True, custom=False, starter=False):
        self.name = name
        self.abbr = self.name
        if self.name in SUBSTRATE_TO_ABBREVIATION:
            self.abbr = SUBSTRATE_TO_ABBREVIATION[self.name]
        self.smiles = smiles
        self.image = None
        self.highlight_image = None
        self.custom = custom
        self.starter = starter

        if set_images:
            if module_type == 'NRPS':
                self.image = os.path.join(NRPS_SUBSTRATE_IMG, f'{name}.png')
                self.highlight_image = os.path.join(NRPS_SUBSTRATE_IMG, f'{name}_highlight.png')

            if module_type == 'PKS':
                self.image = os.path.join(PKS_SUBSTRATE_IMG, f'{name}.png')
                self.highlight_image = os.path.join(PKS_SUBSTRATE_IMG, f'{name}_highlight.png')

            if module_type == 'PKS_starter':
                self.image = os.path.join(PKS_STARTER_SUBSTRATE_IMG, f'{name}.png')
                self.highlight_image = os.path.join(PKS_STARTER_SUBSTRATE_IMG, f'{name}_highlight.png')

        self.rectangle = None
        self.x = 0
        self.y = 0
        self.width = SUBSTRATE_BUTTON_SIZE
        self.height = SUBSTRATE_BUTTON_SIZE


class FattyAcid:
    def __init__(self, c_nr):
        self.c_nr = c_nr
        self.isoform = None
        self.double_bonds = []
        self.tailoring_groups = []
        self.smiles = ''
        self.to_smiles()
        self.image = None
        self.scaled_image = None
        self.draw()
        self.name = None
        self.set_name()

    def set_name(self):
        name = f"C{self.c_nr}"
        if self.isoform:
            name = f"{self.isoform}-{name}"
        if self.double_bonds:
            name = f"{name}:{len(self.double_bonds)}"

        self.name = name

    def set_isoform(self, isoform):
        self.isoform = isoform
        self.draw()

    def set_substrate(self, domain):
        self.to_smiles()
        self.set_name()
        substrate = Substrate(self.name, self.smiles, 'NRPS', set_images=False, custom=True, starter=True)
        domain.substrate = substrate

    def add_methyl_group(self, position):
        self.tailoring_groups.append(('C', position))
        self.draw()

    def add_amino_group(self, position):
        self.tailoring_groups.append(('N', position))
        self.draw()

    def add_hydroxyl_group(self, position):
        self.tailoring_groups.append(('O', position))
        self.draw()

    def add_double_bond(self, position, stereochem):
        self.double_bonds.append((position, stereochem))
        self.draw()

    def count_involvement(self):
        positions = {}
        for pos, cistrans in self.double_bonds:
            if pos not in positions:
                positions[pos] = 0
            positions[pos] += 1

            if pos + 1 not in positions:
                positions[pos + 1] = 0
            positions[pos + 1] += 1

        for group, pos in self.tailoring_groups:
            if pos not in positions:
                positions[pos] = 0

            positions[pos] += 1

        if self.isoform == 'iso':
            if 2 not in positions:
                positions[2] = 0
            positions[2] += 1

        elif self.isoform == 'anteiso':
            if 3 not in positions:
                positions[3] = 0
            positions[3] += 1

        return positions

    def draw(self):
        self.to_smiles()
        path = os.path.join(os.getcwd(), 'temp_fatty_acid_raichu.png')
        png_from_smiles(self.smiles, path)
        self.image = pygame.image.load(path)

        image_width = self.image.get_width()
        image_height = self.image.get_height()

        fatty_acid_image_width, fatty_acid_image_height = FATTY_ACID_IMAGE_SIZE

        if image_width <= fatty_acid_image_width and image_height > fatty_acid_image_height:
            ratio = fatty_acid_image_height / image_height
            image_width *= ratio
            image_height = fatty_acid_image_height
        elif image_width > fatty_acid_image_width and image_height <= fatty_acid_image_height:
            ratio = fatty_acid_image_width / image_width
            image_height *= ratio
            image_width = fatty_acid_image_width
        elif image_width > fatty_acid_image_width and image_height > fatty_acid_image_height:
            ratio_1 = fatty_acid_image_width / image_width
            ratio_2 = fatty_acid_image_height / image_height
            if ratio_1 <= ratio_2:
                ratio = ratio_1
            else:
                ratio = ratio_2

            image_height *= ratio
            image_width *= ratio

        render_window_size = (int(image_width), int(image_height))
        self.scaled_image = pygame.transform.smoothscale(self.image, render_window_size)

    def to_smiles(self):
        smiles = self.c_nr * ['C'] + ['(=O)O']

        if self.isoform == 'iso':
            smiles.pop(0)
            smiles[1] += '(C)'
        elif self.isoform == 'anteiso':
            smiles.pop(0)
            smiles[2] += '(C)'

        self.double_bonds.sort(reverse=True, key=lambda x: x[0])
        for double_bond, stereochem in self.double_bonds:
            smiles[double_bond - 1] += '='

            if stereochem and double_bond != 1:
                carbon_2 = smiles[double_bond]

                if carbon_2.endswith('/') and stereochem == 'trans':
                    smiles[double_bond - 2] += '/'
                elif carbon_2.endswith('/') and stereochem == 'cis':
                    smiles[double_bond - 2] += '\\'
                elif carbon_2.endswith('\\') and stereochem == 'trans':
                    smiles[double_bond - 2] += '\\'
                elif carbon_2.endswith('\\') and stereochem == 'cis':
                    smiles[double_bond - 2] += '/'
                elif stereochem == 'trans':
                    smiles[double_bond] += '/'
                    smiles[double_bond - 2] += '/'
                elif stereochem == 'cis':
                    smiles[double_bond] += '/'
                    smiles[double_bond - 2] += '\\'

        self.tailoring_groups.sort(reverse=True, key=lambda x: x[1])

        for group, pos in self.tailoring_groups:
            carbon = smiles[pos - 1]
            if not carbon.endswith('C'):
                carbon_pos = None
                for i, symbol in enumerate(carbon):
                    if symbol == 'C':
                        carbon_pos = i
                        break
                assert carbon_pos != None
                smiles[pos - 1] = carbon[:carbon_pos + 1] + f'({group})' + carbon[carbon_pos + 1:]

            else:
                smiles[pos - 1] += f'({group})'

        self.smiles = ''.join(smiles)
