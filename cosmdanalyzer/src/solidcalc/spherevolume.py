"""球で構成された図形の体積を計算する."""
import math
from collections.abc import Callable, Collection, Hashable, Iterable
from typing import TypeVar
from . import commonutil
from . import spherepoint
from . import sweepprune
from .typehint import Sphere, Vector3f


_ID = TypeVar("_ID", bound=Hashable)


def calc_multi_sphere_volume(
        sphere_ids: Collection[_ID],
        sphere_getter: Callable[[_ID], Sphere],
        resolution: int,
) -> float:
    """複数の球で構成される図形の体積を計算する

    Args:
        sphere_ids: 球の識別子の集合
        sphere_getter: 識別子から[中心座標,半径]で表される球を返す関数
        resolution: 球面を多面体で近似するときの頂点数
    Returns:
        体積
    """
    col_dict = sweepprune.create_strict_collided_dict(
            sphere_ids, sphere_getter)
    col_func = (lambda a: iter(col_dict.get(a, tuple())))
    vol_func = create_one_volume_in_multi_spheres_func(
            sphere_getter, col_func, resolution)
    return sum(map(vol_func, sphere_ids))


def create_one_volume_in_multi_spheres_func(
        sphere_getter: Callable[[_ID], Sphere],
        collided_sphere_getter: Callable[[_ID], Iterable[_ID]],
        resolution: int,
) -> Callable[[_ID], float]:
    """複数の球からなる立体の中の1つの球の占める体積を計算する関数を返す

    Args:
        sphere_getter: 識別子から[中心座標,半径]で表される球を返す関数
        collided_sphere_getter: 球の識別子から
                                衝突している球の識別子集合を返す関数
        resolution: 球面を多面体で近似するときの頂点数
    Returns:
        球の識別子から対応する球の占める体積を返す関数
    """
    gen_point = spherepoint.SpherePointGenerator()
    return (lambda i: calc_one_sphere_volume(
            sphere_getter(i),
            map(sphere_getter, collided_sphere_getter(i)),
            resolution, lambda n: gen_point.normalized_sphere_points(n)))


def calc_one_sphere_volume(
        sp: Sphere,
        other_spheres: Iterable[Sphere],
        resolution: int,
        sphere_point_generator: Callable[[int], Iterable[Vector3f]],
) -> float:
    """他の球との干渉部分を省いた球1つの体積を返す.

    Args:
        sp: 球
        other_spheres: 干渉する球を(中心座標, 半径)で列挙する
        resolution: 球面を多面体で近似するときの頂点数
        sphere_point_generator: 原点中心半径1の球面上に指定数の点を生成する.
    Returns:
        体積
    """
    planes = tuple(commonutil.calc_sphere_div_plane(sp, s)
                   for s in other_spheres)
    return calc_plane_divided_one_sphere_volume(sp, planes, resolution,
                                                sphere_point_generator)


def calc_plane_divided_one_sphere_volume(
        sp: Sphere,
        planes: Collection[tuple[Vector3f, Vector3f]],
        resolution: int,
        sphere_point_generator: Callable[[int], Iterable[Vector3f]],
) -> float:
    """球を複数の平面で切断した立体の体積を計算する.
    球を分割された球面を底とする錐体の集合で近似して体積を求める.

    Args:
        sp: 球
        planes: 平面(平面に含まれる点, 平面に垂直な単位ベクトル)
        resolution: 球表面の分割数(大きいほど計算精度が高い)
        sphere_point_generator: 原点中心半径1の球面上に指定数の点を生成する.
    Returns:
        体積
    """
    sum_d3 = 0.0
    for v in sphere_point_generator(resolution):
        d = sp[1]
        for plane in planes:
            plane_d = commonutil.calc_intersection_line_plane(
                    (sp[0], v), plane)
            if (plane_d > 0) and (plane_d < d):
                d = plane_d
        sum_d3 += d**3
    return (4 * math.pi / (resolution * 3)) * sum_d3
