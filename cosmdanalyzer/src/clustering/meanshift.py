"""mean_shiftクラスタリング"""
from collections.abc import Callable, Collection, Iterator
import itertools
from typing import TypeVar
from ..neighbors import vptree
from . import neighbor


_EL = TypeVar('_EL')


def weighted_mean_shift(
        elements: Collection[_EL],
        weight_getter: Callable[[_EL], float],
        distance_func: Callable[[_EL, _EL], float],
        bandwidth: float,
        sum_func: Callable[[_EL, _EL], _EL],
        mul_func: Callable[[_EL, float], _EL],
) -> Iterator[Iterator[_EL]]:
    """重み付きの要素に対するmean_shiftクラスタリング

    Args:
        elements: 要素集合
        weight_getter: 要素から重みを返す関数
        distance_func: 要素間の距離関数
        bandwidth: 重心を計算する球の半径
        sum_func: 要素同士の加算
        mul_func: 要素と重みの乗算
    Returns:
        クラスター毎の構成要素のイテレータ
    """
    conv = bandwidth * 0.001
    update_func = _GravityUpdate(
        vptree.VpTree(elements, distance_func),
        weight_getter, bandwidth, sum_func, mul_func)
    el_conv = tuple(map(lambda i: (i, _one_idx_update(
        i, update_func, distance_func, conv)), elements))
    neg_clusters = neighbor.neighbor(
        el_conv, lambda vl, vr: distance_func(vl[1], vr[1]), conv)
    for cl in neg_clusters:
        yield map(lambda v: v[0], cl)


def weighted_mean_shift_detail(
        elements: Collection[_EL],
        weight_getter: Callable[[_EL], float],
        distance_func: Callable[[_EL, _EL], float],
        bandwidth: float,
        sum_func: Callable[[_EL, _EL], _EL],
        mul_func: Callable[[_EL, float], _EL],
) -> Iterator[tuple[Iterator[_EL], _EL, _EL]]:
    """重み付きの要素に対するmean_shiftクラスタリング
    詳細情報として各クラスタの収束先も返す

    Args:
        elements: 要素集合
        weight_getter: 要素から重みを返す関数
        distance_func: 要素間の距離関数
        bandwidth: 重心を計算する球の半径
        sum_func: 要素同士の加算
        mul_func: 要素と重みの乗算
    Returns:
        クラスター毎の構成要素のイテレータと収束値,収束値の最近傍要素
    """
    conv = bandwidth * 0.001
    tree = vptree.VpTree(elements, distance_func)
    update_func = _GravityUpdate(
        tree, weight_getter, bandwidth, sum_func, mul_func)
    el_conv = tuple(map(lambda i: (i, _one_idx_update(
        i, update_func, distance_func, conv)), elements))
    neg_clusters = neighbor.neighbor(
        el_conv, lambda vl, vr: distance_func(vl[1], vr[1]), conv)
    for cl in neg_clusters:
        first = next(cl)
        near_el = tree.nearest_neighbor(first[1])
        yield (itertools.chain((first[0], ), map(lambda v: v[0], cl)),
               first[1], near_el[1])


def _one_idx_update(idx: _EL,
                    update_func: Callable[[_EL], _EL],
                    distance_func: Callable[[_EL, _EL], float],
                    threshold: float) -> _EL:
    """1つの要素を収束するまでupdate_funcで更新する.

    Args:
        idx: 1つの要素
        update_func: 更新関数
        distance_func: 要素間の距離関数
        threshold: 更新距離がしきい値未満になった場合に更新を終了する.
    Returns:
        入力要素の収束値
    """
    cur_idx = idx
    while True:
        next_idx = update_func(cur_idx)
        diff = distance_func(cur_idx, next_idx)
        cur_idx = next_idx
        if diff < threshold:
            return cur_idx


class _GravityUpdate:
    """範囲球内の重心を返す."""

    def __init__(self,
                 tree: vptree.VpTree,
                 weight_getter: Callable[[_EL], float],
                 bandwidth: float,
                 sum_func: Callable[[_EL, _EL], _EL],
                 mul_func: Callable[[_EL, float], _EL],
                 ):
        """

        Args:
            weight_getter: インデックスから値を返す関数
            bandwidth: 範囲球の半径
        """
        self._tree = tree
        self._weight_getter = weight_getter
        self._bandwidth = bandwidth
        self._sum_func = sum_func
        self._mul_func = mul_func

    def __call__(self, query: _EL) -> _EL:
        """範囲球内の重心を返す.

        Args:
            query: 球の中心値
        Returns:
            範囲球内の重心
        """
        itr = map((lambda i: (i, self._weight_getter(i))),
                  map(lambda a: a[1],
                      self._tree.neighbors(query, self._bandwidth)))
        try:
            sum_idx, sum_w = next(itr)
            sum_idx = self._mul_func(sum_idx, sum_w)
        except StopIteration:
            raise RuntimeError()
        for idx, w in itr:
            sum_w += w
            sum_idx = self._sum_func(sum_idx, self._mul_func(idx, w))
        return self._mul_func(sum_idx, 1.0 / sum_w)
