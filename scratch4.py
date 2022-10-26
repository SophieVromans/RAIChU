from raichu.visualize_cluster import *
from pikachu.general import  draw_structure, structure_to_smiles, draw_smiles

daptomycin_cluster = [['module_1', 'starter_module_nrps', 'tryptophan'],
                      ['module_2', 'elongation_module_nrps', 'd-asparagine', []],
                      ['module_3', 'elongation_module_nrps', 'asparticacid', []],
                      ['module_4', 'elongation_module_nrps', 'threonine', []],
                      ['module_5', 'elongation_module_nrps', 'glycine', []],
                      ['module_6', 'elongation_module_nrps', 'ornithine', []],
                      ['module_7', 'elongation_module_nrps', 'asparticacid', []],
                      ['module_8', 'elongation_module_nrps', 'd-alanine', []],
                      ['module_9', 'elongation_module_nrps', 'asparticacid', []],
                      ['module_10', 'elongation_module_nrps', 'glycine', []],
                      ['module_11', 'elongation_module_nrps', 'd-serine', []],
                      ['module_12', 'elongation_module_nrps', '4-methylglutamicacid', []],
                      ['module_13', 'terminator_module_nrps', 'kynurenine', []]]

clust = [['module 1', 'starter_module_pks', 'SC(=O)CC'],
         ['module 2', 'elongation_module_pks', 'pk', []],
         ['module 3', 'elongation_module_nrps', 'nrp', []]]
# product = cluster_to_structure(clust, attach_to_acp=True)
# lin_te = thioesterase_linear_product(product)
# print(structure_to_smiles(lin_te))
# RaichuDrawer(lin_te, save_svg='test.svg')
# draw_structure(lin_te)
# clust = [['module 1', 'starter_module_pks', 'SC(=O)CC'],
#          ['module 2', 'elongation_module_pks', 'ethylmalonylcoa', ['KR']],
#          ['module 3', 'elongation_module_pks', 'methylmalonylcoa', ['KR']],
#          ['module 4', 'elongation_module_pks', 'methylmalonylcoa', ['KR']],
#          ['module 5', 'elongation_module_pks', 'pk', []]]

clust = [['module 1', 'starter_module_pks', 'SC(C(C(O)=O)OC)=O', ''],
 ['module 2', 'elongation_module_pks', 'ethylmalonylcoa', ['KR', 'DH', 'ER'], ''],
 ['module 3', 'elongation_module_pks', 'malonylcoa', ['KR', 'DH', 'ER'], ''],
 ['module 4', 'elongation_module_pks', 'methylmalonylcoa', ['KR', 'DH'], ''],
 ['module 5', 'elongation_module_pks', 'methylmalonylcoa', ['KR'], ''],
 ['module 6', 'terminator_module_pks', 'pk', ['KR'], '']]

draw_cluster(clust)
clust2 = [['module 1', 'starter_module_nrps', 'valine', ['nMT', 'E']],
          ['module 2', 'terminator_module_nrps', 'serine', ['E']],
          ['module 3', 'terminator_module_nrps', 'asparagine', ['nMT', 'E']],
          ['module 4', 'terminator_module_nrps', 'valine', []]]

clust_wrong = [['module 1', 'starter_module_pks', 'SC(=O)CC'],
            ['module 2', 'elongation_module_nrps', 'valine', ['E', 'nMT']],
               ['module 2', 'elongation_module_nrps', 'proline', ['E', 'nMT']],
               ['module 3', 'terminator_module_pks', 'methoxymalonylacp', ['KR', 'DH', 'ER']]]

clust_wrong = [['module 1', 'starter_module_pks', 'SC(=O)CC'],
            ['module 2', 'elongation_module_pks', 'methylmalonylcoa', ['KR']],
            ['module 3', 'elongation_module_pks', 'malonylcoa', ['KR_B1']]]

# draw_cluster(clust_wrong)
# draw_cluster(clust2, interactive=True)
# product = cluster_to_structure(clust2)
# for atom in product.graph:
#     print(atom, atom.annotations.chiral_c_ep)

# ep_product = nrps_epimerization(product)
# RaichuDrawer(ep_product)
# draw_cluster(clust)
#
# # asparagine = read_smiles('CC(C)[C@@H](C(=O)O)N')
# asparagine = read_smiles(('C([C@@H](C(=O)O)N)C(=O)N'))
# asparagine.add_attributes(ATTRIBUTES, boolean=True)
# product = condensation_nrps(asparagine, intermediate)
#
#
# product = attach_to_domain_nrp(product, 'PCP')
#
# RaichuDrawer(product)


####
# cysteine = read_smiles('C([C@@H](C(=O)O)N)S')
# isoleucine = read_smiles('CC[C@H](C)[C@@H](C(=O)O)N')
# valine = read_smiles('CC(C)[C@@H](C(=O)O)N')
# valine.add_attributes(ATTRIBUTES, boolean=True)
# valine = condensation_nrps(isoleucine, valine)
# valine = condensation_nrps(cysteine, valine)
# #
# valine = attach_to_domain_nrp(valine, 'PCP')
# # prod = nrps_epimerization(valine)
# # prod = nrps_methylation(prod)
# RaichuDrawer(valine)

# clustje = [['module 1', 'starter_module_pks', 'SC(C(C(O)=O)CC)=O', 'gene1'],
#          ['module 2', 'elongation_module_pks', 'malonylcoa', [], 'gene1'],
#          ['module 3', 'elongation_module_nrps', 'valine', ['E', 'nMT'], 'gene2'],
#          ['module 4', 'elongation_module_pks', 'methoxymalonylacp', ['KR'], 'gene2'],
#            ['module 3', 'elongation_module_nrps', 'valine', ['E', 'nMT'],
#             'gene2'],
#            ['module 3', 'elongation_module_nrps', 'valine', ['E', 'nMT'],
#             'gene2'],
#            ['module 3', 'elongation_module_nrps', 'valine', ['E', 'nMT'],
#             'gene2'],
#          ['module 5', 'terminator_module_pks', 'methoxymalonylacp', ['KR','DH','ER'], 'gene3']]
# draw_cluster(clustje, save_fig='TEST.png')