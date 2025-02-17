import math
from typing import Any
from .. import visualization
from .. import clustering


def test_single_linkage():
    size = 50
    data = create_data(size)
    distance_3d_threshold = size * 0.1
    cluster, n_clusters = clustering.single_linkage_3d(
            lambda x, y, z: data.get((x, y, z), None),
            data.keys(),
            (lambda i1, d1, i2, d2: distance_func(
                i1, d1, i2, d2, 1.0 / distance_3d_threshold)),
            distance_3d_threshold,
            5)
    cluster_itr = ((p, cluster(*p)) for p in data.keys())
    visualization.plot_identified_3d_points(cluster_itr)


def distance_func(idx1: tuple[int, int, int],
                  data1: tuple[float, float, float, float, float],
                  idx2: tuple[int, int, int],
                  data2: tuple[float, float, float, float, float],
                  idx_weight: float) -> float:
    d = 0.0
    for i1, i2 in zip(idx1, idx2):
        d += ((i1 - i2) * idx_weight)**2
    for d1, d2 in zip(data1, data2):
        if d1 > 0:
            d1 = d1 * 0.5 + 0.5
        if d2 > 0:
            d2 = d2 * 0.5 + 0.5
        d += (d1 - d2)**2
    return math.sqrt(d)


def _plot_exists_points(data: dict[tuple[int, int, int], Any]):
    """データが存在するインデックスを可視化する."""
    visualization.plot_3d_points(data.keys())


def create_data(size: int) -> dict[tuple[int, int, int],
                                   tuple[float, float, float, float, float]]:
    points = []
    points.append(((0.2, 0.2, 0.2), 0.1))
    points.append(((0.3, 0.2, 0.2), 0.15))
    points.append(((0.8, 0.8, 0.5), 0.2))
    points.append(((0.5, 0.5, 0.5), 0.2))
    points.append(((0.8, 0.8, 0.8), 0.2))

    # points.append(((0.5, 0.1, 0.5), 0.05))
    # points.append(((0.5, 0.3, 0.5), 0.05))
    # points.append(((0.5, 0.5, 0.5), 0.05))
    # points.append(((0.5, 0.7, 0.5), 0.05))
    # points.append(((0.5, 0.9, 0.5), 0.05))
    for i in range(len(points)):
        points[i] = ((math.floor(points[i][0][0] * size),
                      math.floor(points[i][0][1] * size),
                      math.floor(points[i][0][2] * size)),
                     points[i][1] * size)
    data = dict()
    for z in range(size):
        for y in range(size):
            for x in range(size):
                approach = []
                exists = False
                for p, r in points:
                    d2 = (x - p[0])**2 + (y - p[1])**2 + (z - p[2])**2
                    r2 = r**2
                    if d2 < r2:
                        approach.append((r2 - d2) / r2)
                        exists = True
                    else:
                        approach.append(0.0)
                if exists:
                    data[(x, y, z)] = tuple(approach)
    return data
