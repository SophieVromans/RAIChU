"""
This script contains examples of the general functionalities of RAIChU.
"""
from visualize_cluster import draw_cluster


if __name__ == "__main__":

    #Calculate structure directly from input PKS cluster
    erythromycin_cluster = [['module_1', 'starter_module_pks', 'SC(=O)CC'],
                           ['module_2', 'elongation_module_pks', 'methylmalonylcoa', ['KR_B2']],
                           ['module_3', 'elongation_module_pks', 'methylmalonylcoa', ['KR_A1']],
                           ['module_4', 'elongation_module_pks', 'methylmalonylcoa', ['KR_C2']],
                           ['module_5', 'elongation_module_pks', 'methylmalonylcoa', ['KR', 'DH', 'ER']],
                           ['module_6', 'elongation_module_pks', 'methylmalonylcoa', ['KR_A1']],
                           ['module_7', 'terminator_module_pks', 'methylmalonylcoa', ['KR_A1']]]

    draw_cluster(erythromycin_cluster)

    #Visualise PKS cluster (interactive mode)
    draw_cluster(erythromycin_cluster, interactive=True)

    #Visualise NRPS cluster
    nrps_cluster = [['NRPS module 1', 'starter_module_nrps', 'd-threonine'],
                    ['NRPS module 2', 'elongation_module_nrps', 'valine', []],
                    ['NRPS module 3', 'elongation_module_nrps', 'serine', ['E', 'nMT']],
                    ['NRPS module 4', 'elongation_module_nrps', 'cysteine', ['E']],
                    ['NRPS module 5', 'elongation_module_nrps', 'glutamicacid', []],
                    ['NRPS module 6', 'elongation_module_nrps', 'alanine', ['E']],
                    ['NRPS module 7', 'terminator_module_nrps', 'valine', []]]
    draw_cluster(nrps_cluster)

    #Visualise hybrid PKS/NRPS cluster
    hybrid =       [['PKS module 1', 'starter_module_pks', 'SC(=O)CC'],
                    ['NRPS module 2', 'elongation_module_nrps', 'valine', []],
                    ['NRPS module 3', 'elongation_module_nrps', 'serine', []],
                    ['PKS module 4', 'elongation_module_pks', 'malonylcoa', ['KR','DH']],
                    ['NRPS module 5', 'elongation_module_nrps', 'glutamicacid', ['E', 'nMT']],
                    ['PKS module 6', 'elongation_module_pks', 'methoxymalonylacp', ['KR','DH','ER']],
                    ['NRPS module 7', 'terminator_module_nrps', 'valine', []]]
    draw_cluster(hybrid, interactive=True)