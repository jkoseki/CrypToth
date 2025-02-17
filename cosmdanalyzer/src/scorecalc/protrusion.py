"""protrusionスコア計算"""
from collections.abc import Callable, Iterable
from ..solidcalc.typehint import Sphere, Vector3f


def calc_patch_protrusion(
        residues: Iterable[int],
        atom_to_pos: Callable[[int], Vector3f],
        res_to_atom: Callable[[int], Iterable[int]],
        d_atom_in_sphere: Callable[[Sphere], Iterable[float]],
) -> float | None:
    """1つのパッチのprotrusionを計算する.

    Args:
        residues: 残基IDの集合
        atom_to_pos: 原子IDを原子座標に変換する関数
        res_to_atoms: 残基IDを構成原子のID集合に変換する関数
        d_atom_in_sphere: 球内にある原子の球中心からの距離を返す関数
    Returns:
        1つのパッチのprotrusion, パッチが空集合の場合はNone
    """
    all_count = 0
    protrude_count = 0
    for res_idx in residues:
        if _calc_residue_protrusion(
                map(atom_to_pos, res_to_atom(res_idx)), d_atom_in_sphere):
            protrude_count += 1
        all_count += 1
    if all_count > 0:
        return protrude_count / all_count
    else:
        return None


def _calc_residue_protrusion(
        positions: Iterable[Vector3f],
        d_atom_in_sphere: Callable[[Sphere], Iterable[float]],
) -> bool:
    """1つの残基が突き出ているか判定する.
    残基が突き出ているとは構成原子の1/2以上が突き出ていると定義する.

    Args:
        d_atom_in_sphere: 球内にある原子の球中心からの距離を返す関数
    Returns:
        突き出ている場合はTrue, いない場合はFalse
    """
    all_count = 0
    protrude_count = 0
    for pos in positions:
        if _calc_atom_protrusion(pos, d_atom_in_sphere):
            protrude_count += 1
        all_count += 1
    return (protrude_count * 2) >= all_count


def _calc_atom_protrusion(
        pos: Vector3f,
        d_atom_in_sphere: Callable[[Sphere], Iterable[float]],
) -> bool:
    """1原子が突き出ているか判定する.
    原子が突き出ているとは原子中心からl2ノルムが8以上12未満の
    他の原子中心の数が120未満であると定義する.

    Args:
        pos: 原子中心の座標
        d_atom_in_sphere: 球内にある原子の球中心からの距離を返す関数
    Returns:
        突き出ている場合はTrue, いない場合はFalse
    """
    count = 0
    for d in d_atom_in_sphere((pos, 12.0)):
        if d >= 8.0:
            count += 1
            if count >= 120:
                return False
    return True
