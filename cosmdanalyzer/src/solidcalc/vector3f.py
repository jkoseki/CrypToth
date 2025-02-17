"""tuple[float, float, float]で表された3次元ベクトルの演算"""
import math
from .typehint import Vector3f


def neg(v: Vector3f) -> Vector3f:
    """-v"""
    return (-v[0], -v[1], -v[2])


def add(v0: Vector3f, v1: Vector3f) -> Vector3f:
    """3次元ベクトルの同士の加算"""
    return ((v0[0] + v1[0]), (v0[1] + v1[1]), (v0[2] + v1[2]))


def sub(v0: Vector3f, v1: Vector3f) -> Vector3f:
    """3次元ベクトルの同士の加算"""
    return ((v0[0] - v1[0]), (v0[1] - v1[1]), (v0[2] - v1[2]))


def dot(v0: Vector3f, v1: Vector3f) -> float:
    """3次元ベクトルの同士の内積"""
    return (v0[0] * v1[0]) + (v0[1] * v1[1]) + (v0[2] * v1[2])


def mul(vec: Vector3f, scalar: float) -> Vector3f:
    """3次元ベクトルとスカラーの積"""
    return (vec[0] * scalar, vec[1] * scalar, vec[2] * scalar)


def norm(v: Vector3f) -> float:
    """ユーグリッドノルム"""
    return math.sqrt(v[0]**2 + v[1]**2 + v[2]**2)


def norm2(v: Vector3f) -> float:
    """ユーグリッドノルムの2乗"""
    return v[0]**2 + v[1]**2 + v[2]**2


def unit(v: Vector3f) -> Vector3f:
    """単位ベクトル"""
    return mul(v, 1 / norm(v))


def add_scalar(vec: Vector3f, scalar: float) -> Vector3f:
    """スカラーとの加算"""
    return (vec[0] + scalar, vec[1] + scalar, vec[2] + scalar)
