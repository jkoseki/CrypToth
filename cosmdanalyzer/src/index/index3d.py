"""3次元のインデックス操作"""
from collections.abc import Iterable, Iterator, Set
import math
import sys
from . import index2d
from . import const


def add(idx0: tuple[int, int, int], idx1: tuple[int, int, int]
        ) -> tuple[int, int, int]:
    """2つの3次元インデックスの加算

    Args:
        idx0: 3次元インデックス
        idx1: 3次元インデックス
    Returns:
        idx0 + idx1
    """
    return (idx0[0] + idx1[0], idx0[1] + idx1[1], idx0[2] + idx1[2])


def cmp(idx0: tuple[int, int, int], idx1: tuple[int, int, int]
        ) -> int:
    """2つの3次元インデックスの加算

    Args:
        idx0: 3次元インデックス
        idx1: 3次元インデックス
    Returns:
        0軸->1軸->2軸の順に優先的に比較して
        idx0が小さいなら-1, 同じなら0, 大きいなら1
    """
    if idx0[0] < idx1[0]:
        return -1
    elif idx0[0] > idx1[0]:
        return 1
    else:
        if idx0[1] < idx1[1]:
            return -1
        elif idx0[1] > idx1[1]:
            return 1
        else:
            if idx0[2] < idx1[2]:
                return -1
            elif idx0[2] > idx1[2]:
                return 1
            return 0


def dence_matrix_3d_indices(size_0: int, size_1: int, size_2: int) \
        -> Iterator[tuple[int, int, int]]:
    """密な3次元行列の3次元インデックスを列挙する
    列挙は2軸->1軸->0軸の順でインクリメントされる.

    Args:
        size_0: 0次元目の要素数
        size_1: 1次元目の要素数
        size_2: 2次元目の要素数
    Returns:
        3次元インデックスを返すイテレータ
    """
    for v0 in range(size_0):
        for v1 in range(size_1):
            for v2 in range(size_2):
                yield (v0, v1, v2)


def sphere_grid_index_iterator(r: float) -> Iterator[tuple[int, int, int]]:
    """原点を中心とする球内部のグリッドのインデックスを列挙する.
    列挙は2軸->1軸->0軸の順でインクリメントされる.

    Args:
        r: 球の半径
    Returns:
        (i0, i1, i2)の順に並んだインデックスを返すイテレータ
    """
    max_r = math.floor(r)
    if (r - max_r < const._MIN_FLOAT):
        r *= 1.0 + const._MIN_FLOAT
    for i0 in range(-max_r, max_r + 1):
        r_2d = math.sqrt(r**2 - i0**2)
        for i1, i2 in index2d.circle_grid_index_iterator(r_2d):
            yield (i0, i1, i2)


def sphere_and_box_grid_index_iterator(
        r: float, box: tuple[tuple[int, int, int], tuple[int, int, int]]
        ) -> Iterator[tuple[int, int, int]]:
    """原点を中心とする球内部のグリッドのインデックスを列挙する.
    列挙は2軸->1軸->0軸の順でインクリメントされる.

    Args:
        r: 球の半径
        box: (開始点,大きさ)で表されるボックス
    Returns:
        (i0, i1, i2)の順に並んだインデックスを返すイテレータ
    """
    rect = (box[0][1:], box[1][1:])
    max_r = math.floor(r)
    if (r - max_r < const._MIN_FLOAT):
        r *= 1.0 + const._MIN_FLOAT
    for i0 in range(max(-max_r, box[0][0]),
                    min(max_r, box[0][0] + box[1][0]) + 1):
        r_2d = math.sqrt(r**2 - i0**2)
        for i1, i2 in index2d.circle_and_rect_grid_index_iterator(r_2d, rect):
            yield (i0, i1, i2)


def half_sphere_grid_index_iterator(r: float
                                    ) -> Iterator[tuple[int, int, int]]:
    """原点を中心とする球内部のグリッドのうち,
    原点より多きいインデックスを列挙する.
    インデックスの大きさは, 軸0,軸1,軸2の順に優先的に比較して
    値が大きい方が大きいとする.
    列挙は2軸->1軸->0軸の順でインクリメントされる.

    Args:
        r: 球の半径
    Returns:
        (i0, i1, i2)の順に並んだインデックスを返すイテレータ
    """
    for i1, i2 in index2d.half_circle_grid_index_iterator(r):
        yield (0, i1, i2)
    max_r = math.floor(r)
    for i0 in range(1, max_r + 1):
        r_2d = math.sqrt(r**2 - i0**2)
        for i1, i2 in index2d.circle_grid_index_iterator(r_2d):
            yield (i0, i1, i2)


def offset_sphere_grid_index_iterator(
        pos: tuple[float, float, float], r: float
        ) -> Iterator[tuple[int, int, int]]:
    """中心と半径を指定して球内部のグリッドのインデックスを列挙する.
    列挙は後ろの次元から順にインクリメントされる.

    Args:
        pos: 円の中心座標
        r: 円の半径
    Returns:
        3次元インデックスを返すイテレータ
    """
    min0 = math.ceil(pos[0] - r)
    max0 = math.floor(pos[0] + r)
    r_2 = r**2
    for idx0 in range(min0, max0 + 1):
        r0 = math.sqrt(abs(r_2 - (idx0 - pos[0])**2))
        for idx1, idx2 in index2d.offset_circle_grid_index_iterator(
                (pos[1], pos[2]), r0):
            yield (idx0, idx1, idx2)


def convert_3d_index_to_1d(idx_3d: tuple[int, int, int],
                           size_3d: tuple[int, int, int]) -> int:
    """3次元のインデックスを1次元のインデックスに変換する."""
    return (idx_3d[0] * size_3d[1] * size_3d[2]
            + idx_3d[1] * size_3d[2] + idx_3d[2])


def convert_1d_index_to_3d(idx_1d: int, size_3d: tuple[int, int, int]
                           ) -> tuple[int, int, int]:
    """1次元のインデックスを3次元のインデックスに変換する."""
    idx = math.floor(idx_1d / size_3d[2])
    idx2 = idx_1d - idx * size_3d[2]
    idx0 = math.floor(idx / size_3d[1])
    idx1 = idx - idx0 * size_3d[1]
    return (idx0, idx1, idx2)


def expand_idxs_float(
        idxs: Iterable[tuple[int, int, int]] | Set[tuple[int, int, int]],
        margin: float,
        box: tuple[tuple[int, int, int], tuple[int, int, int]] | None = None,
        ) -> set[tuple[int, int, int]]:
    """3次元インデックス集合を指定距離分拡大する.

    Args:
        idxs: 3次元インデックス集合
        margin: 拡大距離, 単位はインデック座標
        box: インデックスの許容範囲, (開始点, 大きさ)で表され境界を含む
    Returns:
        拡大後のインデックス集合
    """
    if not isinstance(idxs, Set):
        idxs = set(idxs)
    new_idx_set = set()
    if box is None:
        box = _MAX_BOX
    neg_diff_idxs = (((1, 0, 0), 0, box[1][0] + box[0][0]),
                     ((0, 1, 0), 1, box[1][1] + box[0][1]),
                     ((0, 0, 1), 2, box[1][2] + box[0][2]),
                     ((-1, 0, 0), 0, box[0][0]),
                     ((0, -1, 0), 1, box[0][1]),
                     ((0, 0, -1), 2, box[0][2]))
    for idx in idxs:
        for neg_idx, ax, edge in map(
                lambda n: (add(idx, n[0]), *n[1:]), neg_diff_idxs):
            if not ((neg_idx in idxs) or (neg_idx in new_idx_set)
                    or (idx[ax] == edge)):
                b = ((box[0][0] - idx[0],
                      box[0][1] - idx[1],
                      box[0][2] - idx[2]),
                     box[1])
                for add_idx in sphere_and_box_grid_index_iterator(margin, b):
                    add_idx = add(idx, add_idx)
                    if not (add_idx in idxs):
                        new_idx_set.add(add_idx)
                break
    return idxs | new_idx_set


_MAX_BOX = ((int(sys.maxsize / 4), ) * 3, (int(sys.maxsize / 2), ) * 3)
