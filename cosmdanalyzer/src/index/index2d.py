"""3次元のインデックス操作"""
from collections.abc import Iterator
import math
from . import const


def circle_grid_index_iterator(r: float) -> Iterator[tuple[int, int]]:
    """原点を中心とする円内部のグリッドのインデックスを列挙する.
    円の境界部を含む.
    列挙は1軸->0軸の順でインクリメントされる.

    Args:
        r: 円の半径
    Returns:
        (i0, i1)の順に並んだインデックスを返すイテレータ
    """
    max_0 = math.floor(r)
    if (r - max_0 < const._MIN_FLOAT):
        r *= 1.0 + const._MIN_FLOAT
    for i0 in range(-max_0, max_0 + 1):
        max_1 = math.floor(math.sqrt(r**2 - i0**2))
        for i1 in range(-max_1, max_1 + 1):
            yield (i0, i1)


def circle_and_rect_grid_index_iterator(
        r: float, rect: tuple[tuple[int, int], tuple[int, int]]
        ) -> Iterator[tuple[int, int]]:
    """原点を中心とする円内部かつ指定矩形内部のグリッドのインデックスを列挙する.
    矩形、円ともに境界部を含む.
    列挙は1軸->0軸の順でインクリメントされる.

    Args:
        r: 円の半径
        rect: 開始点と大きさで表される矩形
    Returns:
        (i0, i1)の順に並んだインデックスを返すイテレータ
    """
    max_0 = math.floor(r)
    if (r - max_0 < const._MIN_FLOAT):
        r *= 1.0 + const._MIN_FLOAT
    rect_end = (rect[0][0] + rect[1][0], rect[0][1] + rect[1][1])
    for i0 in range(max(-max_0, rect[0][0]), min(max_0, rect_end[0]) + 1):
        max_1 = math.floor(math.sqrt(r**2 - i0**2))
        for i1 in range(max(-max_1, rect[0][1]), min(max_1, rect_end[1]) + 1):
            yield (i0, i1)


def half_circle_grid_index_iterator(r: float) -> Iterator[tuple[int, int]]:
    """原点を中心とする円内部のグリッドのうち,
    原点より多きいインデックスを列挙する.
    インデックスの大きさは, 1軸,0軸の順に優先的に比較して
    値が大きい方が大きいとする.
    列挙は1軸->0軸の順でインクリメントされる.

    Args:
        r: 円の半径
    Returns:
        (i0, i1)の順に並んだインデックスを返すイテレータ
    """
    max_r = math.floor(r)
    for i1 in range(1, max_r + 1):
        yield (0, i1)
    for i0 in range(1, max_r + 1):
        max_1 = math.floor(math.sqrt(r**2 - i0**2))
        for i1 in range(-max_1, max_1 + 1):
            yield (i0, i1)


def offset_circle_grid_index_iterator(
        pos: tuple[float, float], r: float) -> Iterator[tuple[int, int]]:
    """中心と半径を指定して円内部のグリッドのインデックスを列挙する.
    列挙は後ろの次元から順にインクリメントされる.

    Args:
        pos: 円の中心座標
        r: 円の半径
    Returns:
        2次元インデックスを返すイテレータ
    """
    min0 = math.ceil(pos[0] - r)
    max0 = math.floor(pos[0] + r)
    r_2 = r**2
    for idx0 in range(min0, max0 + 1):
        r1 = math.sqrt(abs(r_2 - (idx0 - pos[0])**2))
        min1 = math.ceil(pos[1] - r1)
        max1 = math.floor(pos[1] + r1)
        for idx1 in range(min1, max1 + 1):
            yield (idx0, idx1)
