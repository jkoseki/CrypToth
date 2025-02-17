"""GFE"""
from collections.abc import Callable, Collection, Iterable, Iterator
import math
import statistics
from .. import solidcalc
from ..solidcalc.typehint import Vector3f


def calc_all_gfe(protein_atom_ids: Collection[int],
                 atom_to_position: Callable[[int], Vector3f],
                 atom_to_vdw_radius: Callable[[int], float],
                 hotspot_values_list: Iterable[Iterable[float]],
                 system_volume: float,
                 n_probe_heavy_atoms: int,
                 temperature: float,
                 solvent_radius: float,
                 resolution: int,
                 ) -> Iterator[float]:
    """すべてのホットスポットのGFEを計算する.

    Args:
        protein_atom_ids: タンパク質の原子ID集合
        atom_to_position: 原子IDから原子座標に変換する関数
        atom_to_vdw_radius: 原子IDから原子半径に変換する関数
        hotspot_values_list: ホットスポット集合
        system_volume: 系全体の体積
        n_probe_heavy_atoms: 系のすべてのプローブ分子の重原子数の合計
        temperature: 系の温度(K)
        solvent_radius: 溶媒半径
        resolution: 原子表面を多面体で近似するときの頂点数
    Returns:
        各hostpotに対応するGFE
    """
    atom_to_last_as_sphere = (
        lambda i: (atom_to_position(i),
                   atom_to_vdw_radius(i) + solvent_radius))
    volume = solidcalc.calc_multi_sphere_volume(
            protein_atom_ids, atom_to_last_as_sphere, resolution)
    solvent_volume = system_volume - volume
    gfe_func = create_gfe_func(n_probe_heavy_atoms, temperature,
                               solvent_volume)
    for hotspot_values in hotspot_values_list:
        yield statistics.mean(map(gfe_func, hotspot_values))


_R = 8.31446261815324
_CALORIE = 4.184
_GFE_CONST = -(_R / (_CALORIE * 1000))


def create_gfe_func(
        n_probe_heavy_atoms: int,
        temperature: float,
        solvent_volume: float) -> Callable[[float], float]:
    """1つの共溶媒占有率から対応するGFEを返す関数を作成する.

    Args:
        n_probe_heavy_atoms: 系のすべてのプローブ分子の重原子数の合計
        temperature: voxel_occupancies計算時の温度(K)
        solvent_volume: 溶媒領域の体積
    Returns:
        1つの共溶媒占有率から対応するGFEを返す関数
    """
    mean_proba = n_probe_heavy_atoms / solvent_volume
    threshold = mean_proba * (math.e**(3 / (_GFE_CONST * temperature)))
    log_mean_proba = math.log(n_probe_heavy_atoms / solvent_volume)
    return (lambda o: (_GFE_CONST * temperature
                       * (math.log(o) - log_mean_proba))
            if o > threshold else 3.0)
