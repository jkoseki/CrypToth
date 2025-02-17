"""球表面を覆う均一な点"""
import math
from collections.abc import Iterator
from .typehint import Vector3f


class SpherePointGenerator:
    """球の表面に均一に発生させた点の座標を生成する."""

    def __init__(self):
        self._normalized_cache = dict()

    def normalized_sphere_points(self, plot_num: int) -> Iterator[Vector3f]:
        """半径1, 中心原点の球の表面を均一に覆う点の座標を列挙する.

        Args:
            plot_num: 球表面の分割数
        Returns:
            半径1, 中心原点の球の表面を均一に覆う点の座標のイテレータ
        """
        if not (plot_num in self._normalized_cache):
            self._normalized_cache[plot_num] = tuple(
                iterate_normalized_sphere_points(plot_num))
        return iter(self._normalized_cache[plot_num])

    def sphere_points(self, plot_num: int, pos: Vector3f, r: float
                      ) -> Iterator[Vector3f]:
        """半径, 中心座標を指定して球の表面を均一に覆う点の座標を列挙する.

        Args:
            plot_num: 球表面の分割数
            pos: 球の中心座標
            r: 球の半径
        Returns:
            球の表面を均一に覆う点の座標を列挙する.
        """
        for x, y, z in self.normalized_sphere_points(plot_num):
            yield (x * r + pos[0], y * r + pos[1], z * r + pos[2])


def iterate_normalized_sphere_points(plot_num: int) -> Iterator[Vector3f]:
    """半径1, 中心原点の球の表面を均一に覆う点の座標を列挙する.

    Args:
        plot_num: 球表面の分割数
    Returns:
        半径1, 中心原点の球の表面を均一に覆う点の座標を列挙する.
    """
    spiral_itr = _iterate_normalized_sphere_spiral_points(plot_num)
    if plot_num >= 14:
        buf = list(next(spiral_itr) for _ in range(6))
        mean = buf[0]
        for i in (1, 3, 4, 5):
            b = buf[i]
            mean = (mean[0] + b[0], mean[1] + b[1], mean[2] + b[2])
        yield (mean[0] / 5, mean[1] / 5, mean[2] / 5)
        for b in buf:
            yield b
        for _ in range(plot_num-14):
            yield next(spiral_itr)
        buf = list(next(spiral_itr) for _ in range(6))
        for b in buf:
            yield b
        mean = buf[0]
        for i in (1, 2, 4, 5):
            b = buf[i]
            mean = (mean[0] + b[0], mean[1] + b[1], mean[2] + b[2])
        yield (mean[0] / 5, mean[1] / 5, mean[2] / 5)
    else:
        return spiral_itr


def _iterate_normalized_sphere_spiral_points(plot_num: int
                                             ) -> Iterator[Vector3f]:
    """半径1, 中心原点の球の表面を螺旋状に覆う点の座標を列挙する.
    極は除外する.

    Args:
        plot_num: 球表面の分割数
    Returns:
        半径1, 中心原点の球の表面を螺旋状に覆う点の座標を列挙する.
        極は除外するため、(plot_num - 2)の点が列挙される.
    """
    c = 3.6 / math.sqrt(plot_num)
    phi = 0.0
    cos_theta = -1.0
    cos_theta_inc = 2 / (plot_num - 1)
    for i in range(1, plot_num - 1):
        cos_theta += cos_theta_inc
        sin_theta = math.sqrt(1 - cos_theta**2)
        phi = phi + c / sin_theta
        yield (sin_theta * math.cos(phi),
               sin_theta * math.sin(phi),
               cos_theta)


def normalized_spherical_coordinate_to_xyz(theta: float, phi: float
                                           ) -> Vector3f:
    """半径1の球面座標をX,Y,Z軸に変換する.
    x = sinθ cosφ
    y = sinθ sinφ
    z = cosθ

    Args:
        theta: θ
        phi: φ
    Returns:
        X,Y,Z軸で表される座標系
    """
    sin_theta = math.sin(theta)
    return (sin_theta * math.cos(phi),
            sin_theta * math.sin(phi),
            math.cos(theta))
