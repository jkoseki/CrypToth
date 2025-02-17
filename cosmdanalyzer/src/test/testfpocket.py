"""fpocketの出力読み込みのユニットテスト"""
import unittest
from pathlib import Path
from src.fpocket import parser


class TestFpocket(unittest.TestCase):

    def test_fpocket(self):
        """fpocketの出力読み込みテスト"""
        with open(Path(__file__).parent / 'test_fpocket.pdb') as pdb, \
             open(Path(__file__).parent / 'test_fpocket_info.txt') as info:
            true_scores = ('0.413', '0.372', '0.335', '0.265')
            true_pdb_1 = ((22.161, 13.864, -10.205),
                          (18.871, 8.254, -7.981),
                          (22.110, 13.820, -10.123),)
            for i, (score, pos) in enumerate(parser.parse_fpocket(pdb, info)):
                self.assertEqual(score, true_scores[i])
                if i == 1:
                    for p, tp in zip(pos, true_pdb_1):
                        self.assertEqual(p[0], tp[0])
                        self.assertEqual(p[1], tp[1])
                        self.assertEqual(p[2], tp[2])
