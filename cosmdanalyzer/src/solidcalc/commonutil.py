"""名前空間外からは呼び出さないユーティリティー関数群"""
import sys
from .typehint import Vector3f
from . import Sphere, vector3f


_SMALL_FLOAT = sys.float_info.radix**int(sys.float_info.min_exp / 2)
_BIG_FLOAT = sys.float_info.radix**int(sys.float_info.max_exp / 2)


def calc_sphere_div_plane(sphere0: Sphere, sphere1: Sphere
                          ) -> tuple[Vector3f, Vector3f]:
    """2つの球を分割する平面を計算する.

    Args:
        sphere0: 球0
        sphere1: 球1
    Returns:
        (平面に含まれる1点の座標, 球0->1方向の平面に垂直なベクトル)
    """
    vec12 = vector3f.sub(sphere1[0], sphere0[0])
    norm12 = vector3f.norm(vec12)
    normal_vec12 = vector3f.mul(vec12, 1.0 / norm12)
    d = (norm12 + (sphere0[1]**2 - sphere1[1]**2) / norm12) * 0.5
    pos = vector3f.add(sphere1[0], vector3f.mul(normal_vec12, d))
    return (pos, normal_vec12)


def is_div_plane_innner(pos: Vector3f,
                        plane: tuple[Vector3f, Vector3f]) -> bool:
    """点が平面の垂直ベクトル側にあるか判定する.

    Args:
        pos: 点の座標
        plane: (平面に含まれる1点の座標, 平面に垂直なベクトル)
    Returns:
        点が平面の垂直ベクトル側にある場合はTrue, それ以外はFalse
    """
    return (vector3f.dot(vector3f.sub(pos, plane[0]), plane[1]) > 0)


def calc_intersection_line_plane(
        line: tuple[Vector3f, Vector3f],
        planes: tuple[Vector3f, Vector3f]) -> float:
    """半直線と平面の交点を開始点からの線上の距離で返す.

    Args:
        line: 半直線(開始座標, 方向ベクトル)
        line_vec: 半直線の方向ベクトル
        planes: 平面(平面に含まれる1点, 平面に垂直なベクトル)
    Return:
        半直線と平面の交点の開始点からの線上の距離
    """
    plane_pos = vector3f.sub(planes[0], line[0])
    distance_plane = vector3f.dot(plane_pos, planes[1])
    len_line_vec = vector3f.dot(line[1], planes[1])
    if abs(len_line_vec) > _SMALL_FLOAT:
        return distance_plane / len_line_vec
    else:
        return _BIG_FLOAT
