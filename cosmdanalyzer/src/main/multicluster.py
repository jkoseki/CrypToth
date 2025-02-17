"""1つ以上のプローブでのクラスタリング処理"""
from collections.abc import (
    Callable, Collection, Hashable, Iterable, Iterator, Sequence)
import itertools
import math
from typing import TypeVar
from .. import clustering

_ID = TypeVar('_ID', bound=Hashable)
_EL = TypeVar('_EL', bound=Hashable)


def multi_mean_shift(
    multi_elements: Iterable[
        tuple[Iterable[_EL], Callable[[_EL], float], _ID]],
    distance_func: Callable[[_EL, _EL], float],
    bandwidth: float,
    sum_func: Callable[[_EL, _EL], _EL],
    mul_func: Callable[[_EL, float], _EL],
) -> Iterator[tuple[Iterator[_EL], set[_ID]]]:
    """
    Args:
        multi_elements: (要素集合,要素から重みを返す関数, ID)の集合
        distance_func: 要素間の距離関数
        bandwidth: 重心を計算する球の半径
        sum_func: 要素同士の加算
        mul_func: 要素と重みの乗算
    Returns:
        クラスタ毎に, クラスタの要素とクラスタに含まれる元要素のID集合を返す
    """
    el_sum = sum_elements(multi_elements, sum_func)
    clusters = clustering.weighted_mean_shift_detail(
        el_sum.keys(), (lambda el: el_sum[el][0]),
        distance_func, bandwidth, sum_func, mul_func)
    for cluster, _, conv in clusters:
        yield (cluster, el_sum[conv][1])


def sum_elements(
        multi_elements: Iterable[
            tuple[Iterable[_EL], Callable[[_EL], float], _ID]],
        sum_func: Callable[[_EL, _EL], _EL],
) -> dict[_EL, tuple[float, set[_ID]]]:
    """重みの和を返す

    Args:
        multi_elements: (要素集合,要素から重みを返す関数, ID)の集合
        sum_func: 要素同士の加算関数
    Returns:
    """
    el_sum: dict[_EL, tuple[float, set[_ID]]] = dict()
    for elements, weight_getter, eid in multi_elements:
        for el in elements:
            if el in el_sum:
                v, s = el_sum[el]
                s.add(eid)
                el_sum[el] = (weight_getter(el) + v, s)
            else:
                el_sum[el] = (weight_getter(el), {eid, })
    return el_sum


def marge_clusters(clusters: Sequence[tuple[Collection[_EL], _ID]],
                   th_rate: float
                   ) -> Iterator[tuple[Iterator[_EL], set[_ID]]]:
    """重複しているクラスタ結合する.
    結合基準はクラスタiの要素がクラスタjに指定率以上含まれる場合
    クラスタiとクラスタjを結合する.
    結合後のクラスタ対他のクラスタの重複率判定は行わない.

    Args:
        clusters: (クラスタの要素集合, クラスタ識別子)の集合
        th_rate: [0.0, 1.0]で表されこの率以上重複しているクラスタ同士を結合する
    Returns:
        (結合後のクラスタの要素集合, 元クラスタの識別子集合)を列挙する
    """
    all_table: dict[_EL, set[int]] = dict()
    for i, (cluster, _) in enumerate(clusters):
        for el in cluster:
            if el in all_table:
                all_table[el].add(i)
            else:
                all_table[el] = {i, }
    cl_table = clustering.ClusterTable(range(len(clusters)))
    for i, (cluster, _) in enumerate(clusters):
        other_count: dict[int, int] = dict()
        n_el = 0
        for el in cluster:
            n_el += 1
            for other in filter(lambda j: j != i, all_table[el]):
                if other in other_count:
                    other_count[other] += 1
                else:
                    other_count[other] = 1
        th = math.ceil(n_el * th_rate)
        for k, v in other_count.items():
            if v >= th:
                cl_table.concat_cluster(i, k)
    new_clusters = cl_table.create_cluster_to_element_sequence()
    for cl in new_clusters:
        src_set: set[int] = set()
        for src_cluster_idx in cl:
            src_set.add(clusters[src_cluster_idx][1])
        yield (itertools.chain.from_iterable(
            map(lambda i: clusters[i][0], cl)), src_set)
