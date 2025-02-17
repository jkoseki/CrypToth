"""vptreeの実行時間計測"""
import random
import time
from ..neighbors import vptree
from ..solidcalc import vector3f
from ..solidcalc.typehint import Vector3f


def speed_vp_tree():
    n_points_list = [100000, 1000000]
    n_query = 10000
    for n_points in n_points_list:
        _speed_n_vptree(n_points, n_query)


def _speed_n_vptree(n_points: int, n_query: int):
    print("# vptree speed test. n_points = {}, n_query = {}.".format(
        n_points, n_query))
    size = 128.0
    points = _create_point_data(n_points, size)
    start_time = time.perf_counter()
    tree = vptree.VpTree(iter(points), _distance_func)
    pass_time = time.perf_counter() - start_time
    print("create time {}".format(pass_time))
    querys = _create_point_data(n_query, size)
    start_time = time.perf_counter()
    for q in querys:
        _, _ = tree.nearest_neighbor(q)
    pass_time = time.perf_counter() - start_time
    print("query time {}".format(pass_time))


def _create_point_data(n_points: int, size: float
                       ) -> tuple[tuple[int, int, int], ...]:
    return tuple(tuple(random.uniform(-size, size) for _ in range(3))
                 for _ in range(n_points))

def _distance_func(v0: Vector3f, v1: Vector3f) -> float:
    return vector3f.norm(vector3f.sub(v0, v1))
