from raichu.data.attributes import ATTRIBUTES


def label_rest_groups(subunit, chain_intermediate):
    atom_types = []

    for atom in chain_intermediate.graph:
        atom_types.append(atom.type)

    nr_unknown_atoms = atom_types.count('*')

    counter = 0
    for atom in subunit.graph:
        if atom.type == '*':
            counter += 1
            # TODO: replace with atom attribute
            atom.unknown_index = nr_unknown_atoms + counter


def initialise_atom_attributes(structure):
    for atom in structure.graph:
        if not hasattr(atom.annotations, 'in_central_chain'):

            for attribute in ATTRIBUTES:
                atom.annotations.add_annotation(attribute, False)


def reset_nrp_annotations(structure):
    for atom in structure.graph:
        atom.annotations.chiral_c_ep = False
        atom.annotations.n_atom_nmeth = False
        atom.annotations.c2_acid = False
        atom.annotations.leaving_oh_h = False
        atom.annotations.leaving_oh_o = False
