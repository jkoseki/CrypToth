"""型のエイリアス"""
from typing import TypeAlias

"""3次元ベクトル"""
Vector3f: TypeAlias = tuple[float, float, float]

"""球 (中心座標, 半径) """
Sphere: TypeAlias = tuple[Vector3f, float]

"""axis-aligned bounding box (最小点, 各軸の幅)"""
Box3d: TypeAlias = tuple[Vector3f, Vector3f]
