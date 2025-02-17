import math
import random
import unittest
from src import index


class Index(unittest.TestCase):

    def test_dence_matrix_3d_indices(self):
        for _ in range(8):
            size = tuple(random.randrange(8, 128) for _ in range(3))
            itr = index.dence_matrix_3d_indices(*size)
            for i0 in range(size[0]):
                for i1 in range(size[1]):
                    for i2 in range(size[2]):
                        self.assertEqual(next(itr), (i0, i1, i2))

    def test_sphere_grid_index(self):
        for _ in range(8):
            r = random.randrange(16, 32)
            r2 = r**2
            idx_set = set()
            for idx in index.sphere_grid_index_iterator(r):
                self.assertTrue(idx[0]**2 + idx[1]**2 + idx[2]**2 <= r2)
                idx_set.add(idx)
            self.assertTrue(len(idx_set) >= (math.floor(r) * 6 + 1))
            for idx in index.sphere_grid_index_iterator(r + 2):
                if not (idx in idx_set):
                    self.assertTrue(
                            idx[0]**2 + idx[1]**2 + idx[2]**2 > r2,
                            msg="({},{},{}), r={}".format(*idx, r))

    def test_sphere_and_box_grid_index(self):
        for _ in range(8):
            r = random.randrange(16, 32)
            r2 = r**2
            idx_set = set()
            box_start = (-int(r * 0.5), -int(r * 0.8), -int(r * 0.9))
            box = (box_start, (int(r * 0.9) - box_start[0],
                               int(r * 1.0) - box_start[1],
                               int(r * 0.7) - box_start[2]))
            for idx in index.sphere_and_box_grid_index_iterator(r, box):
                self.assertTrue(idx[0]**2 + idx[1]**2 + idx[2]**2 <= r2)
                idx_set.add(idx)
            self.assertTrue(len(idx_set) >= (math.floor(r) * 6 + 1))
            for idx in index.sphere_and_box_grid_index_iterator(r + 2, box):
                if not (idx in idx_set):
                    self.assertTrue(
                            idx[0]**2 + idx[1]**2 + idx[2]**2 > r2,
                            msg="({},{},{}), r={}".format(*idx, r))
            i0 = box[0][0] - 1
            for i1 in range(box[0][1], box[0][1] + box[1][1] + 1):
                for i2 in range(box[0][2], box[0][2] + box[1][2] + 1):
                    self.assertTrue(not ((i0, i1, i2) in idx_set))

    def test_expand_idxs_float(self):
        for _ in range(4):
            r = random.randrange(2, 16)
            r2 = r**2
            ex = 2.0
            ex_r2 = (r + ex)**2
            idx_set = set()
            for idx in index.sphere_grid_index_iterator(r):
                idx_set.add(idx)
            ex_idx_set = index.expand_idxs_float(iter(idx_set), ex)
            for idx in ex_idx_set:
                self.assertTrue(idx[0]**2 + idx[1]**2 + idx[2]**2 <= ex_r2)
            self.assertTrue(len(ex_idx_set) >= len(idx_set))
            for idx in (ex_idx_set - idx_set):
                d2 = idx[0]**2 + idx[1]**2 + idx[2]**2
                self.assertTrue(r2 < d2)
                self.assertTrue(d2 <= ex_r2)
