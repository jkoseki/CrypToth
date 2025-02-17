import random
import unittest
from src import graph


class TestGraph(unittest.TestCase):

    def test_match(self):
        n_node = 8
        n_kind = 2
        node0 = tuple((i, random.randrange(0, n_kind)) for i in range(n_node))
        edge0 = create_edge(node0, 4)
        m = graph.match_all(node0[0], (lambda n: iter(edge0[n[0]])),
                            iter(node0), (lambda n: iter(edge0[n[0]])),
                            (lambda lhs, rhs: lhs[1] == rhs[1]))
        self.assertIsNotNone(m)
        for k, v in m.items():
            self.assertEqual(k[0], v[0])
        graph_conv = list(range(n_node))
        random.shuffle(graph_conv)
        node1 = list((graph_conv[n[0]], n[1]) for n in node0)
        node1.sort(key=(lambda n: n[0]))
        edge1 = dict()
        for k, v in edge0.items():
            buf = list()
            for n in v:
                buf.append(node1[graph_conv[n[0]]])
            edge1[graph_conv[k]] = buf
        m = graph.match_all(node0[0], (lambda n: iter(edge0[n[0]])),
                            iter(node1), (lambda n: iter(edge1[n[0]])),
                            (lambda lhs, rhs: lhs[1] == rhs[1]))
        self.assertIsNotNone(m)
        for k, v in m.items():
            self.assertTrue((graph_conv[k[0]] == v[0]) and (k[1] == v[1]))


def create_edge(nodes, max_edge: int) -> dict[int, tuple[int, int]]:
    buf_node = []
    for node in nodes:
        buf_node.append(node)
    n_node = len(nodes)
    random.shuffle(buf_node)
    in_nodes = [buf_node[0], ]
    edge = dict()
    for n in buf_node:
        edge[n[0]] = []
    for i in range(1, n_node):
        to = in_nodes[random.randrange(0, len(in_nodes))]
        edge[to[0]].append(buf_node[i])
        edge[buf_node[i][0]].append(to)
        if len(edge[to[0]]) >= max_edge:
            in_nodes.remove(to)
        in_nodes.append(buf_node[i])
    return edge
