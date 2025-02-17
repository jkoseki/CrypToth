"""CollisionGrid関係"""
import random
import unittest
from src import solidcalc
from src.solidcalc import vector3f


class Collision(unittest.TestCase):

    def test_sphere(self):
        size = 128
        spheres = tuple((tuple(random.uniform(-size, size) for _ in range(3)),
                         random.uniform(size / 64, size / 4))
                        for _ in range(64))
        collided = set()
        for i0 in range(len(spheres)):
            pos0, r0 = spheres[i0]
            for i1 in range(i0 + 1, len(spheres)):
                pos1, r1 = spheres[i1]
                if vector3f.norm2(vector3f.sub(pos0, pos1)) < (r0 + r1)**2:
                    collided.add((i0, i1))
        grid = solidcalc.SphereCollisionGrid(
                (lambda i: spheres[i]), size / 32,
                ((float(-size), float(-size), float(-size)),
                 (float(size), float(size), float(size))))
        grid.add_spheres(range(len(spheres)))
        col_count = dict()
        for i0 in range(len(spheres)):
            for i1 in grid.spheres_in_sphere(spheres[i0], ignore=iter((i0,))):
                idx = (min(i0, i1), max(i0, i1))
                self.assertTrue(idx in collided)
                if idx in col_count:
                    col_count[idx] += 1
                else:
                    col_count[idx] = 1
        # 検出した衝突数のチェック
        all_count = 0
        for _, v in col_count.items():
            self.assertEqual(v, 2)
            all_count += v
        self.assertEqual(len(collided) * 2, all_count)
