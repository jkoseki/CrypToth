"""convexityのスコア計算"""
from collections.abc import Callable, Collection, Iterable, Iterator
import statistics
from .. import common
from ..solidcalc import vector3f
from ..solidcalc.typehint import Sphere, Vector3f


def calc_patch_convexity(res_ids: Collection[int],
                         res_to_atoms: Callable[[int], Iterable[int]],
                         res_to_ca: Callable[[int], int],
                         atom_to_res: Callable[[int], int],
                         atom_to_pos: Callable[[int], Vector3f],
                         atom_to_weight: Callable[[int], float],
                         atom_in_sphere: Callable[[Sphere], Iterable[int]],
                         is_atom_exposed: Callable[[int], bool],
                         distance: float,
                         ) -> float:
    """パッチのconvexityを計算する.

    Args:
        res_ids: 残基IDの集合
        res_to_atoms: 残基IDから原子ID集合を返す関数
        res_to_ca: 残基IDからCα原子のIDを返す関数
        atom_to_res:
        atom_to_pos:
        atom_to_weight:
        atom_in_sphere: 指定球内にある原子中心を原子IDで列挙する関数
        is_atom_exposed:
        distance:
    Returns:
        パッチのconvexity, 露出残基が存在しない場合は-1.0
    """
    neg_res: dict[int, tuple[int, ...]] = dict()
    neg_point = common.BufferdFunction(
        (lambda i:
         calc_residue_solvent_point(
             res_to_atoms(i), is_atom_exposed, atom_to_pos, atom_to_weight,
             atom_to_pos(res_to_ca(i))) if res_to_ca(i) is not None
         else None))
    for res_id in res_ids:
        neg_res[res_id] = tuple(detect_neighbor_residues(
            res_to_atoms(res_id),
            atom_to_pos, atom_to_res, atom_in_sphere, distance,
            (res_id, )))
    try:
        return statistics.mean(filter(
            lambda v: v is not None,
            map((lambda res_id: calc_residue_convexity(
                res_id,
                neg_res[res_id],
                neg_point)
            ),
                res_ids
            )
        ))
    except statistics.StatisticsError:
        return -1.0


def calc_residue_convexity(
        res_id: int,
        neg_res_ids: Iterable[int],
        res_to_point: Callable[[int], tuple[Vector3f, Vector3f] | None],
) -> float | None:
    """残基のconvexityを計算する.

    Args:
        res_id: 計算対象の残基ID
        neg_res_ids: 計算対象の残基に隣接する残基のID集合
        res_to_point: 残基IDから(solvent_point, exposed_centroid)を返す関数
                      値が存在しない場合はNoneを返す
    Returns:
        残基のconvexity, res_to_pointが対象残基でNoneを返す場合はNone
    """
    p = res_to_point(res_id)
    if p is None:
        return None
    solv, expo = p
    sum_conv = 0.0
    count_conv = 0
    for neg_solv, neg_expo in filter(
            lambda x: x is not None, map(res_to_point, neg_res_ids)):
        sum_conv += (vector3f.norm(vector3f.sub(neg_solv, solv))
                     / vector3f.norm(vector3f.sub(neg_expo, expo)) - 1)
        count_conv += 1
    if count_conv > 0:
        return sum_conv / count_conv
    return 0.0


def calc_residue_solvent_point(atom_ids: Iterable[int],
                               is_atom_exposed: Callable[[int], bool],
                               atom_to_pos: Callable[[int], Vector3f],
                               atom_to_weight: Callable[[int], float],
                               default_burried_cent: Vector3f,
                               ) -> tuple[Vector3f, Vector3f] | None:
    """残基のsolvent pointとexposed centroidを求める.

    Args:
        atom_ids: 原子ID集合
        is_atom_exposed: 原子IDを受け取り露出している場合はTrueを返す関数
        burried_atoms: 露出していない原子のID集合
        atom_to_pos: 原子IDから座標への変換関数
        atom_to_weight: 原子IDから重さへの変換関数
        default_burried_cent: 埋没原子が存在しない場合のburried_centの代用値
    Returns:
        (solvent_point, exposed_centroid)
        露出原子が存在しない場合はNone
    """
    exposed_cent = (0.0, 0.0, 0.0)
    burried_cent = (0.0, 0.0, 0.0)
    exposed_count = 0
    burried_count = 0

    for atom_id in atom_ids:
        wp = vector3f.mul(atom_to_pos(atom_id), atom_to_weight(atom_id))
        if is_atom_exposed(atom_id):
            exposed_cent = vector3f.add(exposed_cent, wp)
            exposed_count += 1
        else:
            burried_cent = vector3f.add(burried_cent, wp)
            burried_count += 1
    if exposed_count == 0:
        return None
    exposed_cent = vector3f.mul(exposed_cent, 1 / exposed_count)
    burried_cent = (vector3f.mul(burried_cent, 1 / burried_count)
                    if burried_count > 0 else default_burried_cent)
    solvent_point = vector3f.add(
        exposed_cent,
        vector3f.mul(
            vector3f.unit(vector3f.sub(exposed_cent, burried_cent)), 3)
    )
    return (solvent_point, exposed_cent)


def detect_neighbor_residues(atom_ids: Iterable[int],
                             atom_to_pos: Callable[[int], Vector3f],
                             atom_to_res: Callable[[int], int],
                             atom_in_sphere: Callable[[Sphere], Iterable[int]],
                             distance: float,
                             ignore_residues: Iterable[int] = tuple(),
                             ) -> Iterator[int]:
    """入力原子に隣接する残基を走査する.
    出力対象の残基には入力原子が属する残基も含まれる.

    Args:
        atom_ids: 残基に所属する原子のID集合
        atom_to_pos: 原子のIDを中心座標に変換する関数
        atom_to_res: 原子のIDを所属する残基のIDに変換する関数
        atom_in_sphere: 指定球内にある原子中心を原子IDで列挙する関数
        distance: 指定距離以下に原子間距離が接近している場合隣接残基とみなす
        ignore_residues: 出力に含めない残基のID集合
    Returns:
        入力原子に隣接する残基, 入力原子が属する残基も含まれる.
    """
    searched_res = set(ignore_residues)
    for atom_id in atom_ids:
        pos = atom_to_pos(atom_id)
        for neg_atom_idx in atom_in_sphere((pos, distance)):
            neg_res_idx = atom_to_res(neg_atom_idx)
            if not (neg_res_idx in searched_res):
                yield neg_res_idx
                searched_res.add(neg_res_idx)
