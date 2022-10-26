
from pikachu.smiles.smiles import *
from pikachu.general import read_smiles
from raichu.central_chain_detection.label_central_chain import label_pk_central_chain


def find_central_chain(pks_nrps_attached):
    """Identifies the central chain atoms in the input structure from the
    in_central_chain Atom attribute, and returns these Atom objects as a list

    pks_nrps_attached: PIKAChU Structure object, input (hybrid) PK/NRP structure
    """
    # If the structure is a PK starter unit, find the central carbon chain

    if not any(atom.annotations.in_central_chain for atom in pks_nrps_attached.graph) and\
            len(pks_nrps_attached.find_substructures(read_smiles('C(=O)S'))) > 0:
        pks_nrps_attached = label_pk_central_chain(pks_nrps_attached)

    # Find atoms in the structure inside a cycle
    pks_nrps_attached.find_cycles()
    for atom in pks_nrps_attached.graph:
        if atom.in_ring(pks_nrps_attached):
            atom.inside_ring = True
        else:
            atom.inside_ring = False

    # Identify starting point central chain attached NRP/PK
    sulphur = None
    for atom in pks_nrps_attached.graph:
        if atom.type == 'S' and any(neighbour.annotations.domain_type
                                    for neighbour in atom.neighbours):
            sulphur = atom
        if not atom.annotations.in_central_chain:
            atom.annotations.in_central_chain = False
    assert sulphur

    central_chain = [sulphur]
    visited = [sulphur]
    atom_central_chain = sulphur
    end_atom = False

    # Identify complete central chain from in_central_chain Atom attributes
    while not end_atom:
        for neighbour in atom_central_chain.neighbours:
            if neighbour.annotations.in_central_chain and neighbour not in visited:
                central_chain.append(neighbour)
                visited.append(neighbour)
                atom_central_chain = neighbour
            elif not neighbour.annotations.in_central_chain:
                neighbours = []
                for next_atom in atom_central_chain.neighbours:
                    neighbours.append(next_atom)
                if all(atom in visited for atom in neighbours):
                    end_atom = True
                else:
                    visited.append(neighbour)

    return central_chain

