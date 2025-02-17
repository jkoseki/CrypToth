"""compactnessのスコア計算"""
from collections.abc import Callable, Iterable
from collections import deque
import itertools
from sys import float_info
from ..neighbors import vptree
from ..solidcalc import vector3f
from ..solidcalc.typehint import Vector3f


def calc_patch_compactness(
        residues: Iterable[int],
        res_to_atoms: Callable[[int], Iterable[int]],
        atom_to_pos: Callable[[int], Vector3f],
        is_exposed_atom: Callable[[int], bool],
) -> float:
    """パッチのcompactnessを計算する.

    Args:
        residues: 残基IDの集合
        res_to_atoms: 残基IDを構成原子のID集合に変換する関数
        atom_to_pos: 原子IDを原子座標に変換する関数
        is_exposed_atom: 原子IDが表面原子を表している場合はTrueを返す関数
    Returns:
        compactness
        表面原子を内包する残基が2つ未満の場合は0.0
    """
    res_to_surface_atoms = (
        lambda res_id: filter(is_exposed_atom, res_to_atoms(res_id)))
    sum_count = 0
    comp = 0.0
    resides_queue = deque(residues)
    while len(resides_queue) > 0:
        i_atom = res_to_surface_atoms(resides_queue.popleft())
        i_atom_first = next(i_atom, None)
        if i_atom_first is None:
            continue
        i_tree = vptree.VpTree(
            map(atom_to_pos, itertools.chain((i_atom_first, ), i_atom)),
            _euclidean_distance)
        for j in resides_queue:
            min_d = float_info.max
            for pos in map(atom_to_pos, res_to_surface_atoms(j)):
                min_d = i_tree.nearest_neighbor(pos, thresthold=min_d)[0]
            if min_d < float_info.max:
                comp += min_d
                sum_count += 1
    if sum_count > 0:
        return comp / sum_count
    return 0.0


def _euclidean_distance(v0: Vector3f, v1: Vector3f):
    return vector3f.norm(vector3f.sub(v0, v1))
