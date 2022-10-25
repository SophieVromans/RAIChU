from validation.create_random_clusters import *
import time

if __name__ == "__main__":
    df1 = open('speed_assessment_NRPS.txt', 'a')
    for i in range(1, 501):
        cluster = generate_random_nrps_cluster()
        print(cluster)
        start = time.time()
        RaichuDrawer(cluster, save_fig=f'NRPS_cluster_{i}.png')
        end = time.time()
        delta_time = end - start
        print(len(cluster), ' modules,', delta_time, ' seconds')
        df1.write(f'{delta_time:.5f}\t{len(cluster)}')
        df1.write('\n')
        plt.close('all')
    df1.close()

    df2 = open('speed_assessment_PKS.txt', 'a')
    for i in range(1, 501):
        generated_cluster = generate_random_pks_cluster()
        print(generated_cluster)
        start = time.time()
        RaichuDrawer(generated_cluster, save_fig=f'PKS_cluster_{i}.png')
        end = time.time()
        delta_time = end - start
        print(len(generated_cluster), ' modules,', delta_time, ' seconds')
        df2.write(f'{delta_time:.5f}\t{len(generated_cluster)}')
        df2.write('\n')
        plt.close('all')
    df2.close()

    df3 = open('speed_assessment_hybrid_PKS_NRPS.txt', 'a')
    for i in range(1, 501):
        generated_cluster = generate_random_hybrid_cluster()
        print(generated_cluster)
        start = time.time()
        RaichuDrawer(generated_cluster,
                         save_fig=f'hybrid_PKS_NRPS_cluster_{i}.png')
        end = time.time()
        delta_time = end - start
        print(len(generated_cluster), ' modules,', delta_time, ' seconds')
        df3.write(f'{delta_time:.5f}\t{len(generated_cluster)}')
        df3.write('\n')
        plt.close('all')
    df3.close()

