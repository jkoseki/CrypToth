"""3次元空間の特徴を利用したsingle-linkageクラスタリング"""
from collections.abc import Callable, Iterable, Set, Sequence
import operator
from typing import TypeVar
from . import clusterdata
from .. import index

_EL = TypeVar('_EL')


def single_linkage_3d(
        data_getter: Callable[[tuple[int, int, int]], _EL],
        indices: Set[tuple[int, int, int]] | Iterable[tuple[int, int, int]],
        distance_func: Callable[[tuple[int, int, int], _EL,
                                 tuple[int, int, int], _EL], float],
        distance_3d_threshold: float, distance_threshold: float,
        n_clusters: int,
        ) -> Sequence[Sequence[tuple[int, int, int]]]:
    """3次元空間の特徴を利用したsingle-linkageクラスタリング

    Args:
        data_getter: 3次元空間のインデックスを与えると対応するデータを返す関数.
                     範囲外のインデックスを参照した場合はNoneを返す.
        indices: data_getterが意味のある要素を返すインデックスの集合
        distance_func: 要素間距離を規定する関数,
                       2つのインデックスと要素の値を入力として距離を返す.
        distance_3d_threshold: 3次元空間上のL2距離が指定距離より離れている
                               要素間の評価距離を無限大にする.
        distance_threshold: 要素間距離がしきい値以上の場合無限大にする.
        n_clusters: 最低クラスタ数
    Returns:
        所属インデックスの配列で表されるクラスタの配列
    """
    neg_idxs = list(index.half_sphere_grid_index_iterator(
        distance_3d_threshold))
    cluster_table = clusterdata.ClusterTable[tuple[int, int, int]]()
    cluster_num = 0
    distance_list = []
    indices_set = indices if isinstance(indices, Set) else set(indices)
    for data_i in indices_set:
        data1 = data_getter(data_i)
        for neg_i in neg_idxs:
            neg_i = (neg_i[0] + data_i[0], neg_i[1] + data_i[1],
                     neg_i[2] + data_i[2])
            if neg_i not in indices_set:
                continue
            data2 = data_getter(neg_i)
            if data2 is not None:
                d = distance_func(data_i, data1, neg_i, data2)
                if d < distance_threshold:
                    distance_list.append((data_i, neg_i, d))
        cluster_table.add_element(data_i)
        cluster_num += 1
    distance_list.sort(key=operator.itemgetter(2))
    for ele1, ele2, d in distance_list:
        done_concat = cluster_table.concat_cluster(ele1, ele2)
        if done_concat:
            cluster_num -= 1
            if cluster_num <= n_clusters:
                break
    return cluster_table.create_cluster_to_element_sequence()
