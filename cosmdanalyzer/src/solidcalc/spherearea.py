"""複数の球で構成される図形の表面積を計算する"""
from collections.abc import Callable, Collection, Hashable, Iterable, Iterator
import math
from typing import TypeVar
from . import spherepoint
from . import sweepprune
from .typehint import Vector3f, Sphere
from . import vector3f
from .. import common


_ID = TypeVar('_ID', bound=Hashable)


def calc_multi_sphere_area(
        sphere_ids: Iterable[_ID],
        sphere_getter: Callable[[_ID], Sphere],
        resolution: int,
) -> float:
    """複数の球で構成される図形の表面積を計算する.

    Args:
        sphere_ids: 表面積を計算する球の識別子の集合
        sphere_getter: 識別子から(中心座標,半径)で表される球を返す関数
        resolution: 球面を多面体で近似するときの頂点数
    Returns:
        表面積
    """
    col_dict = sweepprune.create_strict_collided_dict(
            sphere_ids, sphere_getter)
    col_func = (lambda a: iter(col_dict.get(a, tuple())))
    one_func = create_one_area_in_multi_spheres_func(
            sphere_getter, col_func, resolution)
    return sum(map(one_func, sphere_ids))


def create_one_area_in_multi_spheres_func(
        sphere_getter: Callable[[_ID], Sphere],
        collided_sphere_getter: Callable[[_ID], Iterable[_ID]],
        resolution: int,
) -> Callable[[_ID], float]:
    """複数の球からなる立体の中の1つの球の占める表面積を計算する関数を返す

    Args:
        sphere_getter: 識別子から[中心座標,半径]で表される球を返す関数
        collided_sphere_getter: 指定した識別子の球と
                                衝突している球の識別子の集合を返す関数
        resolution: 球面を多面体で近似するときの頂点数
    Returns:
        球IDから1つの球の占める表面積を返す関数
    """
    plot_func = create_one_plot_in_multi_spheres_func(
        sphere_getter, collided_sphere_getter, resolution)

    def _one_area_in_multi_spheres(i: _ID):
        itr, area = plot_func(i)
        return common.len_iterator(itr) * area
    return _one_area_in_multi_spheres


def plot_multi_sphere_surface(
        sphere_ids: Iterable[_ID],
        sphere_getter: Callable[[_ID], Sphere],
        resolution: int,
) -> Iterator[tuple[_ID, Iterator[Vector3f], float]]:
    """複数の球で構成される図形の表面を覆う点を生成する.

    Args:
        sphere_ids: 球の識別子の集合
        sphere_getter: 識別子から[中心座標,半径]で表される球を返す関数
        resolution: 球面を多面体で近似するときの頂点数
    Returns:
        対応する球毎に[識別子, 表面の点の座標集合, 1点あたりの表面積]を返す.
    """
    col_dict = sweepprune.create_strict_collided_dict(
            sphere_ids, sphere_getter)
    col_func = (lambda a: iter(col_dict.get(a, tuple())))
    one_func = create_one_plot_in_multi_spheres_func(
            sphere_getter, col_func, resolution)
    return map(lambda i: (i, *one_func(i)), sphere_ids)


def create_one_plot_in_multi_spheres_func(
        sphere_getter: Callable[[_ID], Sphere],
        collided_sphere_getter: Callable[[_ID], Iterable[_ID]],
        resolution: int,
) -> Callable[[_ID], tuple[Iterator[Vector3f], float]]:
    """複数の球からなる立体の中の1つの球表面を覆う点を計算する関数を返す

    Args:
        sphere_getter: 識別子から[中心座標,半径]で表される球を返す関数
        collided_sphere_getter: 指定した識別子の球と
                                衝突している球の識別子の集合を返す関数
        resolution: 球面を多面体で近似するときの頂点数
    Returns:
        球IDから[表面の点の座標集合, 1点あたりの表面積]を返す関数
    """
    gen_point = spherepoint.SpherePointGenerator()

    def _one_plot_in_multi_spheres(i: _ID):
        pos, r = sphere_getter(i)
        one_area = 4 * math.pi * (r**2) / resolution
        points = _remove_in_sphere_points(
            gen_point.sphere_points(resolution, pos, r),
            tuple(sphere_getter(i) for i in collided_sphere_getter(i)))
        return (points, one_area)
    return _one_plot_in_multi_spheres


def search_surface_spheres(
        sphere_ids: Iterable[_ID],
        sphere_getter: Callable[[_ID], Sphere],
        resolution: int,
) -> Iterator[_ID]:
    """複数の球で構成される図形の表面に存在する球を列挙する.

    Args:
        sphere_ids: 球の識別子の集合
        sphere_getter: 識別子から[中心座標,半径]で表される球を返す関数
        resolution: 球面を多面体で近似するときの頂点数
    Returns:
        表面に存在する球の識別子のイテレータ
    """
    col_dict = sweepprune.create_strict_collided_dict(
            sphere_ids, sphere_getter)
    col_func = (lambda a: iter(col_dict.get(a, tuple())))
    one_func = create_is_surface_in_multi_spheres_func(
            sphere_getter, col_func, resolution)
    return filter(one_func, sphere_ids)


def create_is_surface_in_multi_spheres_func(
        sphere_getter: Callable[[_ID], Sphere],
        collided_sphere_getter: Callable[[_ID], Iterable[_ID]],
        resolution: int,
) -> Callable[[_ID], bool]:
    """指定された球が複数の球からなる立体の表面に存在する場合Trueを返す関数
    を返す

    Args:
        sphere_getter: 識別子から[中心座標,半径]で表される球を返す関数
        collided_sphere_getter: 指定した識別子の球と
                                衝突している球の識別子の集合を返す関数
        resolution: 球面を多面体で近似するときの頂点数
    Returns:
        球IDから立体の表面に存在する場合Trueを返す関数
    """
    plot_func = create_one_plot_in_multi_spheres_func(
        sphere_getter, collided_sphere_getter, resolution)

    def _is_surface_in_multi_spheres(i: _ID):
        try:
            next(plot_func(i)[0])
            return True
        except StopIteration:
            return False
    return _is_surface_in_multi_spheres


def _remove_in_sphere_points(
        points: Iterable[Vector3f], spheres: Collection[Sphere]
) -> Iterator[Vector3f]:
    """入力座標イテレータから球の内部にある座標を削除したイテレータを返す.

    Args:
        points: 削除対象の座標リスト, 直接編集される.
        pos: 球の中心座標
        r: 球の半径
    Returns:
        級の内部に無い点を返すイテレータ
    """
    r2_iterable = tuple(r**2 for _, r in spheres)
    for pos in points:
        for (center, _), r2 in zip(spheres, r2_iterable):
            if vector3f.norm2(vector3f.sub(pos, center)) < r2:
                break
        else:
            yield pos
