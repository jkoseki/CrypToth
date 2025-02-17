"""線形近傍探索
遅いので比較用
"""
from collections.abc import Callable, Iterable
from sys import float_info
from typing import TypeVar


_V = TypeVar("_V")


def search_nearest_neighbor(values: Iterable[_V],
                            distance_func: Callable[[_V, _V], float],
                            query: _V) -> tuple[float, _V]:
    """線形近傍探索
    遅いので比較用

    Args:
        values: 最近傍要素の候補集合
        distance_func: 距離を定義する関数
        query: この要素の最近傍要素をVP木から探す
    Returns:
        (最近傍要素の距離, 最近傍要素)
    """
    min_d = float_info.max
    min_v = None
    for v in values:
        d = distance_func(query, v)
        if min_d > d:
            min_d = d
            min_v = v
    return (min_d, min_v)
