import random
import unittest
from ..neighbors import liner, vptree
from ..solidcalc import vector3f
from ..solidcalc.typehint import Vector3f


class TestVpTree(unittest.TestCase):

    def test_create_vp_tree(self):
        n_points = 64
        size = 128.0
        points = tuple(tuple(random.uniform(-size, size) for _ in range(3))
                       for _ in range(n_points))
        tree = vptree.VpTree(iter(points), distance_func)
        n_query = 64
        for _ in range(n_query):
            query = tuple(tuple(random.uniform(-size, size) for _ in range(3)))
            td, tv = tree.nearest_neighbor(query)
            ad, av = liner.search_nearest_neighbor(
                iter(points), distance_func, query)
            self.assertEqual(td, ad)
            self.assertEqual(tv, av)
        radius = size / 8.0
        for _ in range(n_query):
            query = tuple(tuple(random.uniform(-size, size) for _ in range(3)))
            neg = set(map(lambda v: v[1], iter(tree.neighbors(query, radius))))
            for p in points:
                d = distance_func(query, p)
                if p in neg:
                    self.assertTrue(d < radius)
                else:
                    self.assertTrue(d >= radius)


def distance_func(v0: Vector3f, v1: Vector3f) -> float:
    return vector3f.norm(vector3f.sub(v0, v1))
