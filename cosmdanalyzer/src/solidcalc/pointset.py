"""点集合に対する計算"""
from collections.abc import Iterable
from .typehint import Vector3f
import math
from sys import float_info
from ..neighbors import vptree

_min_growth = math.sqrt(float_info.min)


def calc_point_set_box_overwraped(
        l_set: Iterable[Vector3f],
        r_set: Iterable[Vector3f],
        r_box_size: float) -> float:
    """点集合l_setがr_setに重なっている割合を返す.
    重なっている条件は、それぞれの点集合に含まれる2点間の各座標の距離の
    最大のものがr_box_size/2以下である場合とする.

    Args:
        l_set: 点集合
        r_set: 点集合
        r_box_size: r_setの点が表すBOXの1辺の長さ
    """
    r_tree = vptree.VpTree(r_set, _ax_max_distance)
    overwrap_count = 0
    all_count = 0
    for l_val in l_set:
        all_count += 1
        d, v = r_tree.nearest_neighbor(l_val, (r_box_size * 0.5) * _min_growth)
        if v is not None:
            overwrap_count += 1
    return overwrap_count / all_count


def _ax_max_distance(lhs: Vector3f, rhs: Vector3f):
    """各軸について差の絶対値をとり,その最大値を距離とする"""
    return max(abs(lhs[i] - rhs[i]) for i in range(3))
