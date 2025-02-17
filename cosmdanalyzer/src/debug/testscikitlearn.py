"""scikit_learnが3次元空間特徴を利用したクラスタリングが可能か調べる"""
#import time
#from typing import Any
#import sklearn.cluster
#import scipy.sparse
#
#
#def debug_scikit_learn():
#    """scikit_learnが3次元空間特徴を利用したクラスタリングが可能か調べる"""
#    data, neg = generate_input(10)
#    print('start', flush=True)
#    start_time = time.perf_counter()
#
#    clustering = sklearn.cluster.AgglomerativeClustering(
#        n_clusters=2,
#        metric='euclidean',
#        # connectivity=neg,
#        linkage='ward').fit(data)
#
#    end_time = time.perf_counter()
#    print(end_time - start_time)
#    print(clustering.labels_)
#
#
#def generate_input(size: int) -> tuple[list[tuple[float, float, float]], Any]:
#    data = []
#    for z in range(size):
#        for y in range(size):
#            for x in range(size):
#                data.append((x, y, z))
#    neg_mat_data = []
#    row = []
#    col = []
#    for i in range(size**3):
#        for j in range(max(0, i - 3), min(size, i + 3)):
#            neg_mat_data.append(1)
#            row.append(i)
#            col.append(j)
#    return (data,
#            scipy.sparse.csr_matrix((neg_mat_data, (row, col)),
#                                    shape=(size**3, size**3)))
