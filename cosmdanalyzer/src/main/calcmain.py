"""計算部分のメインルーチン"""
from collections.abc import (
    Callable, Collection, Iterable, Iterator, MutableSequence, Sequence
)
import math
import os
import pathlib
from typing import NamedTuple
from gridData import Grid
from .. import chem
from .. import common
# from .. import visualization
from .. import index
from ..neighbors import vptree
from .. import solidcalc
from ..solidcalc import vector3f
from ..solidcalc.typehint import Sphere, Vector3f
from .. import fpocket
from . import calccharge
from . import fpocketscore
from . import output
from . import spot
from . import scoretype
from . import input
from ..scorecalc import gfe, rmsf


class MyGrid(NamedTuple):
    to_value: Callable[[tuple[int, int, int]], float]
    to_pos: Callable[[tuple[int, int, int]], Vector3f]
    shape: tuple[int, int, int]
    size: float


class SingleSystem(NamedTuple):
    mol: chem.Mol
    n_probe_heavy_atoms: int
    exposed_atom_set: tuple[set[int]]
    grid: MyGrid
    basename: str
    fpocket_pdb: str | bytes | os.PathLike
    fpocket_info: str | bytes | os.PathLike


def calc_main(src_system_infos: Iterable[input.SystemInfo],
              out_dir_path: str | bytes | os.PathLike,
              occupancy_threashold: float,
              clustering_input: (spot.SingleLinkageInput
                                 | spot.DbscanInput
                                 | spot.MeanShiftInput),
              hotspot_extend: float,
              fpocket_threthold: float,
              hydrophobicity_path: str | bytes | os.PathLike,
              charge_path: str | bytes | os.PathLike,
              temperature: float,
              solvent_radius: float,
              score_weight: Sequence[float],
              spot_marge_rate: float,
              resolution: int,
              output_detail: bool,
              verbose: bool):
    """計算部分のメインルーチン

    Args:
        temperature: トラジェクトリ作成時の温度(K)
        solvent_radius: 溶媒半径
        score_weight: 各スコアの重み
        spot_marge_rate: [0.0, 1.0]で表されこの率以上重複している
                         異なるプローブのスポットを結合する
        resolution: 球面を多面体で近似するときの頂点数
        output_detail: スコアの詳細を出力する場合はTrue, しない場合はFalse
        verbose: 標準出力に詳細な処理情報を表示する場合はTrue, しない場合はFalse
    """
    src_systems = tuple(init_single_system(
        info, solvent_radius, resolution)
        for info in src_system_infos)
    grid_idx_to_pos = src_systems[0].grid[1]
    grid_shape = src_systems[0].grid[2]
    grid_size = src_systems[0].grid[3]
    hotspot_idx_list = tuple(spot.detect_multi_hotspots(
        map(lambda v: to_detect_hotspot(
            v.mol, v.n_probe_heavy_atoms, v.exposed_atom_set, v.grid,
            v.basename, occupancy_threashold),
            src_systems),
        grid_idx_to_pos,
        5.0,
        clustering_input,
        hotspot_extend,
        ((0, 0, 0), (grid_shape[0] - 1, grid_shape[1] - 1, grid_shape[2] - 1)),
        grid_size,
        spot_marge_rate,
    ))
    hotspot_id_list = tuple(id for _, id in hotspot_idx_list)
    hotspot_idx_list = tuple(idx for idx, _ in hotspot_idx_list)
    if verbose:
        print('n_hotspot = {}'.format(len(hotspot_idx_list)))

    hotspot_list = tuple(
        tuple(map(grid_idx_to_pos, hotspot_idxs))
        for hotspot_idxs in hotspot_idx_list
    )

    def weight_func(i: int):
        return lambda s: (score_weight[i], s)
    mean_scores = [[0.0, ] * len(hotspot_list)
                   for _ in range((len(score_weight) + 1))]
    n_all_frames = 0
    for src_system in src_systems:
        if verbose:
            print('calc system {}, n_frame = {}'.format(
                src_system.basename, src_system.mol.get_num_conformers()))
        mol = src_system.mol
        protein_idxs = tuple(mol.get_atom_idxs())
        n_probe_heavy_atoms = src_system.n_probe_heavy_atoms
        exposed_atom_set = src_system.exposed_atom_set
        grid_idx_to_val = src_system.grid.to_value

        res_atom_idxs = mol.divide_to_residue(protein_idxs)
        res_to_atoms = (lambda r: iter(res_atom_idxs[r]))
        patch_list = tuple(spot.detect_frame_union_patches(
            hotspot_list, mol.atom_to_residue,
            res_to_atoms,
            ((lambda a: mol.atom_to_position(a, i))
             for i in range(mol.get_num_conformers())),
            exposed_atom_set,
        ))
        (score_gfe, score_fpocket, score_hydrophobicity,
         score_size, score_protrusion, score_convexity, score_compactness,
         score_charge_density, score_rmsf) = calc_scores(
            mol, protein_idxs, res_to_atoms, hotspot_idx_list,
            hotspot_list, patch_list, exposed_atom_set,
            grid_shape, grid_size, grid_idx_to_val, n_probe_heavy_atoms,
            solvent_radius, temperature,
            src_system.fpocket_info, src_system.fpocket_pdb,
            fpocket_threthold,
            hydrophobicity_path, charge_path,
            output_detail, resolution, verbose)
        sum_score = tuple(
            weighted_sum(s)
            for s in zip(
                map(weight_func(0), score_gfe),
                map(weight_func(1), score_size.get_result()),
                map(weight_func(2), score_protrusion.get_result()),
                map(weight_func(3), score_convexity.get_result()),
                map(weight_func(4), score_compactness.get_result()),
                map(weight_func(5), score_hydrophobicity),
                map(weight_func(6), score_charge_density.get_result()),
                map(weight_func(7), score_rmsf.get_result()),
                map(weight_func(8), score_fpocket),
            ))
        n_frame = mol.get_num_conformers()
        n_all_frames += n_frame
        mul_frame = (lambda v: v * n_frame)
        add_to_sequence(mean_scores[0], map(mul_frame, sum_score))
        add_to_sequence(mean_scores[1], map(mul_frame, score_gfe))
        add_to_sequence(mean_scores[2], map(
            mul_frame, score_size.get_result()))
        add_to_sequence(mean_scores[3], map(
            mul_frame, score_protrusion.get_result()))
        add_to_sequence(mean_scores[4], map(
            mul_frame, score_convexity.get_result()))
        add_to_sequence(mean_scores[5], map(
            mul_frame, score_compactness.get_result()))
        add_to_sequence(mean_scores[6], map(mul_frame, score_hydrophobicity))
        add_to_sequence(mean_scores[7], map(
            mul_frame, score_charge_density.get_result()))
        add_to_sequence(mean_scores[8], map(
            mul_frame, score_rmsf.get_result()))
        add_to_sequence(mean_scores[9], map(mul_frame, score_fpocket))
        # output
        write_pymol_src_wrapper(
            out_dir_path, src_system.basename, mol,
            hotspot_list, patch_list, res_to_atoms,
            sum_score, score_gfe, score_size, score_protrusion,
            score_convexity, score_compactness, score_hydrophobicity,
            score_charge_density, score_rmsf, score_fpocket,
            output_detail,
        )
    for mean_score in mean_scores:
        mul_scaler_to_sequence(mean_score, 1.0 / n_all_frames)
    write_mean_score_info_file(
        pathlib.Path(out_dir_path) / 'all_info.txt', mean_scores)
    write_hotspot_probe_file(pathlib.Path(out_dir_path) / 'spot_probe.toml',
                             hotspot_id_list)


def write_hotspot_probe_file(
        out_file: str | bytes | os.PathLike,
        hotspot_id_list: Iterable[Iterable[str]]) -> None:
    """ホットスポットに対応するプローブを出力する"""
    with open(out_file, 'w') as f:
        for i, hotspot_id in enumerate(hotspot_id_list):
            try:
                itr = iter(sorted(hotspot_id))
                f.write('Patch {} = {}'.format(i, next(itr)))
                for id in itr:
                    f.write(', {}'.format(id))
                f.write('\n')
            except StopIteration:
                pass


def write_mean_score_info_file(
        out_file: str | bytes | os.PathLike,
        mean_scores: Iterable[Iterable[float]]) -> None:
    mean_scores = iter(mean_scores)
    rot_mean_socres = []
    scores = next(mean_scores)
    for s in scores:
        rot_mean_socres.append([s, ])
    for scores in mean_scores:
        for s, r in zip(scores, rot_mean_socres):
            r.append(s)
    score_names = ('score', 'gfe', 'size', 'protrusion', 'convexity',
                   'compactness', 'hydrophobicity', 'charge_density',
                   'flexibility', 'fpocket')
    patch_scores = (zip(score_names, s) for s in rot_mean_socres)
    chem.write_score_info_file(out_file, patch_scores)


def write_pymol_src_wrapper(
        out_dir_root: str | bytes | os.PathLike,
        basename: str,
        mol: chem.Mol,
        hotspot_list: Iterable[Iterable[Vector3f]],
        patch_list: Iterable[Iterable[int]],
        res_to_atoms: Callable[[int], Iterable[int]],
        sum_score: Iterable[float],
        score_gfe: Iterable[float],
        score_size: scoretype.ScoreSize,
        score_protrusion: scoretype.ScoreProtrusion,
        score_convexity: scoretype.ScoreConvexity,
        score_compactness: scoretype.ScoreCompactness,
        score_hydrophobicity: Iterable[float],
        score_charge_density: scoretype.ScoreChargeDensity,
        score_rmsf: Iterable[float],
        score_fpocket: Iterable[float],
        output_detail: bool):
    def add_name(name: str):
        return lambda s: (name, s)
    scores = zip(
        map(add_name('score'), sum_score),
        map(add_name('gfe'), score_gfe),
        map(add_name('size'), score_size.get_result()),
        map(add_name('protrusion'), score_protrusion.get_result()),
        map(add_name('convexity'), score_convexity.get_result()),
        map(add_name('compactness'), score_compactness.get_result()),
        map(add_name('hydrophobicity'), score_hydrophobicity),
        map(add_name('charge_density'), score_charge_density.get_result()),
        map(add_name('flexibility'), score_rmsf.get_result()),
        map(add_name('fpocket'), score_fpocket),
    )
    all_patch_info: list[tuple[Iterable[int],
                               Iterable[Vector3f],
                               chem.PatchScores]] = []
    for patch_res_idxs, hotspot, scores in zip(
            patch_list, hotspot_list, scores):
        atom_idxs = res_to_atom_iterator(patch_res_idxs, res_to_atoms)
        all_patch_info.append((hotspot, atom_idxs, scores))
    out_dir_path = pathlib.Path(out_dir_root) / basename
    chem.write_pymol_src(basename, out_dir_path, mol.get_rdkit_mol(),
                         all_patch_info)
    if output_detail:
        out_file = os.path.join(out_dir_path, 'score_detail.csv')
        with open(out_file, mode='w') as out:
            output.output_detail_score_csv(
                out,
                zip(score_size.get_detail_result(),
                    score_protrusion.get_detail_result(),
                    score_convexity.get_detail_result(),
                    score_compactness.get_detail_result(),
                    score_charge_density.get_detail_result()),
                ('size', 'protrusion', 'convexity',
                 'compactness', 'charge_density'),
            )


def add_to_sequence(dst: MutableSequence[float], values: Iterable[float]
                    ) -> None:
    for i, v in enumerate(values):
        dst[i] += v


def mul_scaler_to_sequence(dst: MutableSequence[float], val: float
                           ) -> None:
    for i in range(len(dst)):
        dst[i] *= val


def calc_scores(
        mol: chem.Mol,
        protein_idxs: Collection[int],
        res_to_atoms: Callable[[int], Iterable[int]],
        hotspot_idx_list: Iterable[Iterable[tuple[int, int, int]]],
        hotspot_list: Iterable[Iterable[Vector3f]],
        patch_list: Collection[set[int]],
        exposed_atom_set: Sequence[set[int]],
        grid_shape: tuple[int, int, int],
        grid_size: float,
        grid_idx_to_val: Callable[[tuple[int, int, int]], float],
        n_probe_heavy_atoms: int,
        solvent_radius: float,
        temperature: float,
        fpocket_info: str | bytes | os.PathLike | None,
        fpocket_pdb: str | bytes | os.PathLike | None,
        fpocket_threthold: float,
        hydrophobicity_path: str | bytes | os.PathLike,
        charge_path: str | bytes | os.PathLike,
        output_detail: bool,
        resolution: float,
        verbose: bool,
) -> tuple[Sequence[float, ...], Sequence[float, ...], Sequence[float, ...],
           scoretype.ScoreSize,
           scoretype.ScoreProtrusion,
           scoretype.ScoreConvexity,
           scoretype.ScoreCompactness,
           scoretype.ScoreChargeDensity,
           rmsf.AllPatchRmsfCalc]:
    return (*calc_non_frame_scores(
        mol, protein_idxs, res_to_atoms,
        hotspot_idx_list, hotspot_list, patch_list, grid_shape,
        grid_size, grid_idx_to_val, n_probe_heavy_atoms,
        solvent_radius, temperature, fpocket_info, fpocket_pdb,
        fpocket_threthold,
        hydrophobicity_path, resolution),
        *calc_frame_scores(
        mol, protein_idxs, res_to_atoms,
        patch_list, exposed_atom_set, solvent_radius, output_detail,
        resolution, charge_path, verbose)
    )


def calc_frame_scores(
        mol: chem.Mol,
        protein_idxs: Collection[int],
        res_to_atoms: Callable[[int], Iterable[int]],
        patch_list: Collection[set[int]],
        exposed_atom_set: Sequence[set[int]],
        solvent_radius: float,
        output_detail: bool,
        resolution: int,
        charge_path: str | bytes | os.PathLike,
        verbose: bool,
) -> tuple[scoretype.ScoreSize,
           scoretype.ScoreProtrusion,
           scoretype.ScoreConvexity,
           scoretype.ScoreCompactness,
           scoretype.ScoreChargeDensity,
           rmsf.AllPatchRmsfCalc]:
    all_res_idxs: set[int] = set()
    for p in patch_list:
        all_res_idxs.update(p)
    atom_to_charge = calccharge.calc_atoms_charge_from_rtp_file(
        all_res_idxs,
        res_to_atoms,
        mol.atom_to_atomic_number,
        mol.get_neighbor_atoms,
        mol.atom_to_residue_symbol,
        charge_path)
    res_to_ca = common.BufferdFunction[int, int](
        lambda res_id: residue_to_ca_index(
            res_id, res_to_atoms, mol.atom_to_name))
    score_size = scoretype.ScoreSize(
        patch_list, res_to_atoms, mol.atom_to_vdw_radius, resolution,
        calc_detail=output_detail)
    score_protrusion = scoretype.ScoreProtrusion(
        patch_list, res_to_atoms, calc_detail=output_detail)
    score_convexity = scoretype.ScoreConvexity(
        patch_list, res_to_atoms, res_to_ca,
        mol.atom_to_residue, mol.atom_to_weight,
        4.0,
        calc_detail=output_detail)
    score_compactness = scoretype.ScoreCompactness(
        patch_list, res_to_atoms, calc_detail=output_detail)
    score_charge_density = scoretype.ScoreChargeDensity(
        patch_list, res_to_atoms,
        mol.atom_to_vdw_radius, solvent_radius,
        (lambda a: atom_to_charge[a]),
        resolution=resolution,
        calc_detail=output_detail)
    score_rmsf = rmsf.AllPatchRmsfCalc(
        res_to_atoms, patch_list, mol.atom_to_weight)
    for frame_idx in range(mol.get_num_conformers()):
        if verbose:
            print('.', end='')
        atom_to_pos = (lambda a: mol.atom_to_position(a, frame_idx))
        tree = vptree.VpTree[tuple[int, Vector3f]](
            map(lambda i: (i, atom_to_pos(i)), protein_idxs),
            lambda vl, vr: vector3f.norm(vector3f.sub(vl[1], vr[1])))
        d_atom_in_sphere = (lambda s: map(lambda v: v[0],
                                          tree.neighbors((0, s[0]), s[1])))
        atom_in_sphere = (lambda s: map(lambda v: v[1][0],
                                        tree.neighbors((0, s[0]), s[1])))
        atom_to_sphere = (lambda i: (atom_to_pos(i),
                                     mol.atom_to_vdw_radius(i)))
        vdw_col_sphere = gen_col_sphere(protein_idxs, atom_to_sphere)
        atom_to_as_sphere = (
            lambda i: (atom_to_pos(i),
                       mol.atom_to_vdw_radius(i) + solvent_radius))
        as_col_sphere = gen_col_sphere(protein_idxs, atom_to_as_sphere)
        is_exposed_atom = (lambda a: a in exposed_atom_set[frame_idx])
        score_rmsf.add_frame(atom_to_pos)
        score_size.add_frame(atom_to_pos, vdw_col_sphere)
        score_protrusion.add_frame(atom_to_pos, d_atom_in_sphere)
        score_convexity.add_frame(atom_to_pos, atom_in_sphere,
                                  is_exposed_atom)
        score_compactness.add_frame(atom_to_pos, is_exposed_atom)
        score_charge_density.add_frame(atom_to_pos, is_exposed_atom,
                                       as_col_sphere)
    if verbose:
        print()
    return (score_size, score_protrusion, score_convexity, score_compactness,
            score_charge_density, score_rmsf)


def calc_non_frame_scores(
        mol: chem.Mol,
        protein_idxs: Collection[int],
        res_to_atoms: Callable[[int], Iterable[int]],
        hotspot_idx_list: Iterable[Iterable[tuple[int, int, int]]],
        hotspot_list: Iterable[Iterable[Vector3f]],
        patch_list: Iterator[set[int]],
        grid_shape: tuple[int, int, int],
        grid_size: float,
        grid_idx_to_val: Callable[[tuple[int, int, int]], float],
        n_probe_heavy_atoms: int,
        solvent_radius: float,
        temperature: float,
        fpocket_info: str | bytes | os.PathLike | None,
        fpocket_pdb: str | bytes | os.PathLike | None,
        fpocket_threthold: float,
        hydrophobicity_path: str | bytes | os.PathLike,
        resolution: int,
) -> tuple[Sequence[float, ...], Sequence[float, ...], Sequence[float, ...]]:
    """

    Args:
        fpocket_threthold: [0,1]hotspotに含まれるfpocketの座標点の割合が
                                                        指定値未満の場合は無視する.

    """
    score_gfe = tuple(gfe.calc_all_gfe(
        protein_idxs,
        (lambda i: (mol.atom_to_position(i, mol.get_num_conformers() - 1))),
        mol.atom_to_vdw_radius,
        common.deep2_map(grid_idx_to_val, hotspot_idx_list),
        grid_shape[0] * grid_shape[1] * grid_shape[2] * grid_size**3,
        n_probe_heavy_atoms,
        temperature,
        solvent_radius,
        resolution,
    ))
    n_hotspot = len(score_gfe)
    # fpocket
    score_fpocket: list[float | None] = [0.0, ] * n_hotspot
    if (fpocket_pdb is not None) and (fpocket_info is not None):
        with open(fpocket_pdb, 'r') as pdb,\
                open(fpocket_info, 'r') as info:
            src_fpocket = tuple(fpocket.parse_fpocket(pdb, info))
        for i, hotspot in enumerate(hotspot_list):
            score_fpocket[i] = fpocketscore.fpocket_score(
                hotspot, grid_size, src_fpocket, fpocket_threthold)
    # hydrophobicity
    score_hydrophobicity: list[float] = []
    hydrophobicity_table = chem.ResidueHydrophobicity(hydrophobicity_path)
    for patch in patch_list:
        count = 0
        hydro_sum = 0.0
        for res in patch:
            try:
                hydro_sum += hydrophobicity_table(mol.atom_to_residue_symbol(
                    next(iter(res_to_atoms(res)))))
                count += 1
            except KeyError:
                pass
        score_hydrophobicity.append(hydro_sum / count
                                    if count > 0 else 0.0)
    return (score_gfe, score_fpocket, score_hydrophobicity)


def res_to_atom_iterator(
        res_idxs: Iterable[int], res_to_atom: Callable[[int], Iterable[int]]
) -> Iterator[int]:
    """残基のインデックス集合から原子のインデックス集合に変換する."""
    for res_idx in res_idxs:
        for atom_idx in res_to_atom(res_idx):
            yield atom_idx


def get_grid_access(grid: Grid) -> MyGrid:
    """OpenDX形式のグリッドオブジェクトから使用するアクセス情報を取得する.

    Args:
        grid: OpenDX形式のグリッドを読み出したオブジェクト
    Returns:
        (3次元インデックスから対応グリッドの値を返す関数,
         3次元インデックスから対応グリッドの座標を返す関数,
         グリッドの各軸の数,
         1グリッドの幅,
        )
    """
    return MyGrid(
        to_value=(lambda idx: grid.grid[idx[0], idx[1], idx[2]]),
        to_pos=(lambda idx:
                (grid.edges[0][idx[0]],
                 grid.edges[1][idx[1]],
                 grid.edges[2][idx[2]])),
        shape=grid.grid.shape,
        size=grid.edges[0][1] - grid.edges[0][0])


def gen_col_sphere(sphere_ids: Iterable[int],
                   sphere_getter: Callable[[int], Sphere]
                   ) -> Callable[[int], Iterator[int]]:
    """球のIDから衝突している球のIDを返す関数を生成する.

    Args:
        sphere_ids: 球のID集合
        sphere_getter: IDから対応する球を返す関数
    Returns:
        球のIDから衝突している球のIDを返す関数
    """
    col_dict = solidcalc.create_strict_collided_dict(
        sphere_ids, sphere_getter)
    return (lambda a: iter(col_dict.get(a, tuple())))


def residue_to_ca_index(res_id: int,
                        res_to_atoms: Callable[[int], Iterable[int]],
                        atom_to_name: Callable[[int], str]) -> int:
    """残基IDからCα原子のIDを返す."""
    for atom in res_to_atoms(res_id):
        if atom_to_name(atom) == ' CA ':
            return atom


def weighted_sum(w_v: Iterable[tuple[float, float]]) -> float:
    return sum(map(lambda v: v[0] * v[1], w_v))


def calc_exposed_atoms_set_all_frame(
        atom_ids: Collection[int],
        atom_to_pos: Callable[[int, int], Vector3f],
        atom_to_vdw_radius: Callable[[int], float],
        solvent_radius: float,
        frame_indicies: Iterable[int],
        resolution: int,
) -> Iterator[set[int]]:
    """すべてのフレームについて溶媒露出原子を求める.

    Args:
        atom_ids: 対象の原子ID集合
        atom_to_pos: 原子IDから原子座標を返す関数
        atom_to_vdw_radius: 原子IDからファンデルワールス半径を返す関数
        solvent_radius: 溶媒半径
        resolution: 球面を多面体で近似するときの頂点数
    Returns:
        フレーム毎の溶媒露出原子のID集合
    """
    def atom_to_as_sphere(frame_idx: int):
        return (lambda i: (atom_to_pos(i, frame_idx),
                           atom_to_vdw_radius(i) + solvent_radius))
    return (
        set(solidcalc.search_surface_spheres(
            atom_ids,
            atom_to_as_sphere(frame_idx),
            resolution,
        ))
        for frame_idx in frame_indicies
    )


def init_single_system(
        info: input.SystemInfo,
        solvent_radius: float,
        resolution: float,
) -> SingleSystem:
    """1プローブのトラジェクトリの初期処理を行う."""
    pdb_str, n_probe_heavy_atoms = input.trajectory_pdb_files_filter(info.pdbs)
    mol = chem.create_mol_from_pdb_str(pdb_str)
    protein_idxs = tuple(mol.get_atom_idxs())
    exposed_atom_set = tuple(
        calc_exposed_atoms_set_all_frame(
            protein_idxs, mol.atom_to_position, mol.atom_to_vdw_radius,
            solvent_radius, range(mol.get_num_conformers()), resolution
        )
    )
    return SingleSystem(
            mol=mol,
            n_probe_heavy_atoms=n_probe_heavy_atoms,
            exposed_atom_set=exposed_atom_set,
            grid=get_grid_access(Grid(common.path_to_str(info.dx))),
            basename=info.basename,
            fpocket_pdb=info.fpocket_pdb,
            fpocket_info=info.fpocket_info,
            )


def to_detect_hotspot(
        mol: chem.Mol,
        n_probe_heavy_atoms: int,
        exposed_atom_all_frame: Iterable[Iterable[int]],
        grid: MyGrid,
        id: str,
        occupancy_threashold: float,
) -> tuple[Iterable[tuple[int, int, int]],
           Callable[[tuple[int, int, int]], float],
           float,
           Iterable[Vector3f],
           int]:
    return (index.dence_matrix_3d_indices(*grid.shape),
            grid.to_value,
            occupancy_threashold * n_probe_heavy_atoms,
            to_all_atoms_pos(mol, exposed_atom_all_frame),
            id,
            )


def to_all_atoms_pos(mol: chem.Mol,
                     atom_all_frame: Iterable[Iterable[int]]
                     ) -> Iterator[Vector3f]:
    """フレーム毎に指定した原子集合のすべての座標を
    1次元のイテレータとして返す"""
    for i, exposed_atoms in enumerate(atom_all_frame):
        for a in exposed_atoms:
            yield mol.atom_to_position(a, i)
