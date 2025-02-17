# import random
import unittest
from src.main import multicluster


class Clustering(unittest.TestCase):

    def test_marge(self):
        cl0 = (1, 2, 3, 4, 5, 6)
        cl1 = (3, 4, 5, 6, 7, 8)
        cl2 = (13, 14, 15, 16, 17, 18)
        new_cls = multicluster.marge_clusters(
                ((cl0, 0), (cl1, 1), (cl2, 2)), 0.5)
        cl, i = next(new_cls)
        self.assertEqual(set(cl), set(cl0) | set(cl1))
        self.assertEqual(i, {0, 1})
        cl, i = next(new_cls)
        self.assertEqual(set(cl), set(cl2))
        self.assertEqual(i, {2, })
