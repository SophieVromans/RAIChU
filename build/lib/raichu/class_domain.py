from pikachu.chem.structure import *
from pikachu.chem.atom import AtomAnnotations, AtomDrawProperties, Atom
from raichu.data.attributes import ATTRIBUTES


class CarrierDomain(Atom):
    def __init__(self, atom_type, atom_nr, chiral, charge, aromatic):
        super().__init__(atom_type, atom_nr, chiral, charge, aromatic)
        self.type = atom_type
        self.domain_type = None
        self.domain_nr = atom_nr
        self.annotations = AtomAnnotations()
        for annotation in ATTRIBUTES:
            self.annotations.add_annotation(annotation, default=False)

    def __repr__(self):
        return f'{self.annotations.domain_type}_{self.domain_nr}'

    def __hash__(self):
        return self.nr

    def __eq__(self, domain):
        if type(domain) == CarrierDomain:
            return self.nr == domain.nr
        else:
            return False

    def print_domaintype(self):
        print(self.domain_type)

    def copy(self):
        domain_copy = CarrierDomain(self.type, self.nr, self.chiral, self.charge, self.aromatic)
        domain_copy.hybridisation = self.hybridisation
        domain_copy.pyrrole = self.pyrrole
        domain_copy.furan = self.furan
        domain_copy.thiophene = self.thiophene
        domain_copy.draw = AtomDrawProperties()
        domain_copy.draw.colour = self.draw.colour
        domain_copy.annotations = self.annotations.copy()
        domain_copy.domain_nr = domain_copy.nr
        domain_copy.domain_type = self.domain_type
        connectivity = []

        for connection in self.connectivity:
            connectivity.append(connection)

        domain_copy.connectivity = tuple(connectivity)

        for neighbour in self.neighbours:
            domain_copy.neighbours.append(neighbour)

        # Note: when copying a structure, neighbours need to be refreshed again
        # after all atoms have been changed!

        return domain_copy

    def set_domain_type(self, domain_type):
        assert domain_type == 'ACP' or domain_type == 'AT' or\
               domain_type == 'KS' or domain_type == 'KR' or\
               domain_type == 'DH' or domain_type == 'ER' or\
               domain_type == 'TE' or domain_type == 'PCP'
        self.domain_type = domain_type
        self.annotations.set_annotation('domain_type', domain_type)


def make_domain(domain_type, next_atom_nr):
    """
    Returns a Domain of the input type containing the appropriate nr in the
    Structure object

    domain_type: Str, domain type
    next_atom_nr: Int, 'atom nr' of the domain in the structure
    """
    domain = CarrierDomain('I', next_atom_nr, None, 0, False)
    domain.set_domain_type(domain_type)
    return domain


def make_scaffold_domain(domain_type):
    """
    Returns a Domain of the input type containing the appropriate nr in the
    Structure object

    domain_type: Str, domain type
    next_atom_nr: Int, 'atom nr' of the domain in the structure
    """
    structure = Structure()

    domain = CarrierDomain('I', 0, None, 0, False)
    domain.set_domain_type(domain_type)

    structure.add_disconnected_atom(domain)

    sulphur = Atom('S', 1, None, 0, False)
    structure.add_bond(domain, sulphur, 'single', 0)
    structure.refine_structure()

    return structure




