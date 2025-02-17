"""スポットと対応するパッチを見つける."""
from collections.abc import Callable, Hashable, Iterable, Iterator
import itertools
import math
from typing import NamedTuple, TypeVar
from .. import clustering
from .. import index
from ..neighbors import vptree
from ..solidcalc import vector3f
from ..solidcalc.typehint import Vector3f
from . import multicluster


_ID = TypeVar('_ID', bound=Hashable)


class SingleLinkageInput(NamedTuple):
    threshold: float


class DbscanInput(NamedTuple):
    epsilon: float
    min_pts: int


class MeanShiftInput(NamedTuple):
    bandwidth: float


def input_to_index_unit(
    input_data: (SingleLinkageInput | DbscanInput | MeanShiftInput),
    width_one_index: float,
) -> (SingleLinkageInput | DbscanInput | MeanShiftInput):
    """クラスタリングの入力パラメータをインデックス単位の距離に変換する."""
    w = width_one_index
    if isinstance(input_data, SingleLinkageInput):
        return SingleLinkageInput(to_index_unit(input_data.threshold, w))
    elif isinstance(input_data, DbscanInput):
        return DbscanInput(to_index_unit(input_data.epsilon, w),
                           input_data.min_pts)
    elif isinstance(input_data, MeanShiftInput):
        return MeanShiftInput(to_index_unit(input_data.bandwidth, w))
    raise TypeError


def to_index_unit(val: float, width_one_index: float):
    """valをindex座標系の値に変換する.
    結果が0以外の整数の場合はindex点が境界にならないよう微小量増加させる.
    """
    ret = val / width_one_index
    if ret.is_integer():
        if ret > 0:
            ret += 1.0e-8
        elif ret < 0:
            ret -= 1.0e-8
    return ret


def detect_hotspots(
        voxel_indicies: Iterable[tuple[int, int, int]],
        voxel_to_value: Callable[[tuple[int, int, int]], float],
        voxel_threshold: float,
        voxel_to_pos: Callable[[tuple[int, int, int]], Vector3f],
        exposed_atoms_itr: Iterator[Iterator[int]],
        atom_to_pos_itr: Iterator[Callable[[int], Vector3f]],
        pos_threshold: float,
        clustering_input: SingleLinkageInput | DbscanInput | MeanShiftInput,
        expand: float,
        idxs_box: tuple[tuple[int, int, int], tuple[int, int, int]],
        voxel_width: float,
) -> Iterator[tuple[tuple[int, int, int], ...]]:
    """ホットスポットを捜査する.

    Args:
        voxel_indicies: ボクセルの3次元インデックス集合
        voxel_to_value: インデックスに対応する値を返す関数
        voxel_threshold: インデックスの値が指定値以上の場合のみ使用する
        voxel_to_pos: ボクセルのインデックスに対応する座標を返す関数
        exposed_atoms_itr: フレーム毎の溶媒露出原子のID集合
        atom_to_pos_itr: フレーム毎の原子IDから原子座標を返す関数
        pos_threshold: 溶媒露出原子からの距離が指定値以下のボクセルのみ使う
        clustering_input: クラスタリングアルゴリズムへの入力パラメータ
                          単位はインデックス座標
        expand: ホットスポットを指定距離分拡大する
        idxs_box: ホットスポット拡大時のインデックスの許容範囲
                  (開始点, 大きさ)で表され境界を含む
        voxel_width: ボクセル1つあたりの幅
    Returns:
        ホットスポット毎のボクセルインデックス集合
    """
    all_atoms_pos = itertools.chain.from_iterable(
        map(to_pos, atoms)
        for to_pos, atoms in zip(atom_to_pos_itr, exposed_atoms_itr)
    )
    voxel_indicies = filter(
        create_voxel_filter(
            voxel_to_value, voxel_threshold, voxel_to_pos,
            all_atoms_pos, pos_threshold),
        voxel_indicies
    )
    clusters = clustering_voxels(
        voxel_indicies, voxel_to_value,
        input_to_index_unit(clustering_input, voxel_width))
    if expand > 0.0:
        clusters = map(
            lambda cl: index.expand_idxs_float(
                cl, to_index_unit(expand, voxel_width), idxs_box),
            clusters)
    return map(lambda cl: sorted(cl), clusters)


def detect_multi_hotspots(
        multi_voxels: Iterable[tuple[Iterable[tuple[int, int, int]],
                                     Callable[[tuple[int, int, int]], float],
                                     float,
                                     Iterable[Vector3f],
                                     _ID]],
        voxel_to_pos: Callable[[tuple[int, int, int]], Vector3f],
        pos_threshold: float,
        clustering_input: SingleLinkageInput | DbscanInput | MeanShiftInput,
        expand: float,
        idxs_box: tuple[tuple[int, int, int], tuple[int, int, int]],
        voxel_width: float,
        marge_rate: float,
) -> Iterator[tuple[tuple[tuple[int, int, int]], set[_ID]]]:
    """複数のボクセル集合を統合したホットスポットを捜査する.

    Args:
        multi_voxels: ボクセル集合毎に
                      (ボクセルのインデックス集合,
                       ボクセルインデックスから値を返す関数,
                       指定しきい値以上のボクセルのみ使用する,
                       全フレームの露出原子の座標集合, 識別子)
        voxel_to_pos: ボクセルのインデックスに対応する座標を返す関数
        pos_threshold: 溶媒露出原子からの距離が指定値以下のボクセルのみ使う
        clustering_input: クラスタリングアルゴリズムへの入力パラメータ
                          単位はインデックス座標
        expand: ホットスポットを指定距離分拡大する
        idxs_box: ホットスポット拡大時のインデックスの許容範囲
                  (開始点, 大きさ)で表され境界を含む
        voxel_width: ボクセル1つあたりの幅
        marge_rate: [0.0, 1.0]で表されこの率以上重複している
                    クラスタ同士を結合する
    Returns:
        (ホットスポット毎のボクセルインデックス集合, 元クラスタのID集合)
    """
    multi_voxels = (
        (filter(create_voxel_filter(to_v, threshold, voxel_to_pos,
                                    atoms_pos, pos_threshold),
                voxels),
         to_v, i)
        for voxels, to_v, threshold, atoms_pos, i in multi_voxels)
    clusters = multi_clustering_voxels(
        multi_voxels, input_to_index_unit(clustering_input, voxel_width),
        marge_rate)
    if expand > 0.0:
        clusters = (
            (index.expand_idxs_float(
                cl, to_index_unit(expand, voxel_width), idxs_box),
             ids)
            for cl, ids in clusters
        )
    return ((sorted(cl), ids) for cl, ids in clusters)


def create_voxel_filter(
        voxel_to_value: Callable[[tuple[int, int, int]], float],
        voxel_threshold: float,
        voxel_to_pos: Callable[[tuple[int, int, int]], Vector3f],
        exposed_atoms_pos_itr: Iterator[Vector3f],
        pos_threshold: float,
) -> Callable[[tuple[int, int, int]], bool]:
    """ボクセルインデックスを入力として有効な場合はTrueを返す関数を生成する.

    Args:
        voxel_to_value: インデックスに対応する値を返す関数
        voxel_threshold: インデックスの値が指定値以上の場合のみ使用する
        voxel_to_pos: ボクセルのインデックスに対応する座標を返す関数
        exposed_atoms_pos_itr: 全フレームの溶媒露出原子の座標集合
        pos_threshold: 溶媒露出原子からの距離が指定値以下のボクセルのみ使う
    Returns:
            ボクセルインデックスを入力として有効な場合はTrueを返す関数
    """
    tree = vptree.VpTree(exposed_atoms_pos_itr,
                         lambda vl, vr: vector3f.norm(vector3f.sub(vl, vr)))
    return (lambda i:
            (voxel_to_value(i) >= voxel_threshold)
            and tree.exists_neighbor(voxel_to_pos(i), pos_threshold)
            )


def multi_clustering_voxels(
        multi_voxels: Iterable[tuple[Iterable[tuple[int, int, int]],
                                     Callable[[tuple[int, int, int]], float],
                                     _ID]],
        clustering_input: SingleLinkageInput | DbscanInput | MeanShiftInput,
        marge_rate: float,
) -> Iterator[tuple[Iterator[tuple[int, int, int]], set[_ID]]]:
    """複数のボクセルデータからクラスタリングを行う

    Args:
        multi_voxels:
        clustering_input:
        marge_rate: [0.0, 1.0]で表されこの率以上重複している
                    クラスタ同士を結合する
    Returns:
        (結合後のクラスタの要素集合, 元クラスタの識別子集合)を列挙する
    """
    if isinstance(clustering_input, MeanShiftInput):
        clusters = multicluster.multi_mean_shift(
            multi_voxels, _distance_func,
            clustering_input.bandwidth, vector3f.add, vector3f.mul)
    else:
        if isinstance(clustering_input, SingleLinkageInput):
            cluster_func = (
                lambda v, v_to_val:
                clustering.single_linkage_3d(
                    v_to_val, v, _indexed_distance_func,
                    clustering_input.threshold,
                    clustering_input.threshold, 1
                )
            )
        elif isinstance(clustering_input, DbscanInput):
            cluster_func = (
                lambda v, _:
                clustering.dbscan(
                    v, _distance_func, clustering_input.epsilon,
                    clustering_input.min_pts
                )
            )
        cluster_itr = (zip(map(lambda c: tuple(c), cluster_func(v, v_to_val)),
                           itertools.repeat(i))
                       for v, v_to_val, i in multi_voxels)
        clusters = multicluster.marge_clusters(
            tuple(itertools.chain.from_iterable(cluster_itr)), marge_rate)
    return clusters


def clustering_voxels(
        voxel_indicies: Iterable[tuple[int, int, int]],
        voxel_to_value: Callable[[tuple[int, int, int]], float],
        clustering_input: SingleLinkageInput | DbscanInput | MeanShiftInput,
) -> Iterator[Iterable[tuple[int, int, int], ...]]:
    """ボクセルを指定アルゴリズムでクラスタリングする.

    Args:
        voxel_indicies: 操作対象の3次元インデックス集合
        voxel_to_value: インデックスに対応する値を取得する関数
        clustering_input: クラスタリングアルゴリズムへの入力パラメータ
                          単位はインデックス座標
    Returns:
        クラスター毎のボクセルインデックス集合
    """
    if isinstance(clustering_input, SingleLinkageInput):
        clusters = iter(clustering.single_linkage_3d(
            voxel_to_value, voxel_indicies, _indexed_distance_func,
            clustering_input.threshold, clustering_input.threshold, 1))
    elif isinstance(clustering_input, DbscanInput):
        clusters = clustering.dbscan(
            voxel_indicies, _distance_func, clustering_input.epsilon,
            clustering_input.min_pts)
    elif isinstance(clustering_input, MeanShiftInput):
        clusters = clustering.weighted_mean_shift(
            tuple(voxel_indicies),
            voxel_to_value, _distance_func, clustering_input.bandwidth,
            vector3f.add, vector3f.mul)
    return clusters


def detect_frame_union_patches(
        hotspots: Iterable[Iterable[Vector3f]],
        atom_to_res: Callable[[int], int],
        res_to_atoms: Callable[[int], Iterable[int]],
        atom_to_pos_itr: Iterable[Callable[[int], Vector3f]],
        exposed_atoms_itr: Iterable[Iterable[int]],
) -> Iterator[set[int]]:
    """すべてのフレームのパッチの和集合を求める.

    Args:
        hotspots: ホットスポット毎のボクセル座標集合
        atom_to_res: 原子IDから残基IDへの変換関数
        res_to_atoms: 残基IDから構成原子ID集合を返す関数
        atom_to_pos_itr: フレーム毎の原子IDから原子座標への変換関数
        exposed_atoms_itr: フレーム毎の溶媒接触可能な原子のID集合を返す関数
    Returns:
        ホットスポット毎に対応するスポットの残基ID集合
    """
    res_surface_atoms: dict[int, list[Vector3f]] = dict()
    for atom_to_pos, exposed_atoms in zip(atom_to_pos_itr, exposed_atoms_itr):
        searched_res: set[int] = set()
        surface_res_itr = map(
            atom_to_res,
            filter(lambda a: atom_to_res(a) not in searched_res, exposed_atoms)
        )
        for res_id in surface_res_itr:
            searched_res.add(res_id)
            if not (res_id in res_surface_atoms):
                res_surface_atoms[res_id] = list()
            res_surface_atoms[res_id].extend(
                map(atom_to_pos, res_to_atoms(res_id)))
    res_trees: list[tuple[int, vptree.VpTree]] = list()
    for res_id, res_atoms_pos in res_surface_atoms.items():
        tree = vptree.VpTree(
            res_atoms_pos,
            lambda vl, vr: vector3f.norm(vector3f.sub(vl, vr)))
        res_trees.append((res_id, tree))
    distance = 5.0
    for hotspot in hotspots:
        patch: set[int] = set()
        for res_id, tree in res_trees:
            for hotp in hotspot:
                if tree.exists_neighbor(hotp, distance):
                    patch.add(res_id)
                    break
        yield patch


def _distance_func(idx1: tuple[float, float, float],
                   idx2: tuple[float, float, float]):
    return math.sqrt((idx1[0] - idx2[0])**2
                     + (idx1[1] - idx2[1])**2
                     + (idx1[2] - idx2[2])**2)


def _indexed_distance_func(idx1: tuple[int, int, int],
                           data1: float,
                           idx2: tuple[int, int, int],
                           data2: float) -> float:
    d = 0.0
    for i1, i2 in zip(idx1, idx2):
        d += (i1 - i2)**2
    # if data1 > 0:
    #     data1 = data1 * 0.5 + 0.5
    # if data2 > 0:
    #     data2 = data2 * 0.5 + 0.5
    # d += (data1 - data2)**2
    return math.sqrt(d)
