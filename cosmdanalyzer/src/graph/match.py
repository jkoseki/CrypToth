"""グラフ構造のマッチング"""
from collections.abc import Callable, Iterable
import itertools
from typing import TypeVar


_N0 = TypeVar("_N0")
_N1 = TypeVar("_N1")


def match_all(node0: _N0,
              edge0: Callable[[_N0], Iterable[_N0]],
              nodes1: Iterable[_N1],
              edge1: Callable[[_N1], Iterable[_N1]],
              is_same_node: Callable[[_N0, _N1], bool]
              ) -> dict[_N0, _N1] | None:
    """2つのグラフが一致するか判定する.
    一致する場合はグラフ0からグラフ1のノード対応表を返す.

    Args:
        nodes0: グラフ0起点ノード
        edge0: グラフ0のノードから隣接ノードを返す関数
        nodes1: グラフ1のノード集合
        edge1: グラフ1のノードから隣接ノードを返す関数
        is_same_node: ノードが一致する場合はTrueを返す関数
    Returns:
        2つのグラフが一致する場合はグラフ0からグラフ1のノード対応表
        一致しない場合はNone
    """
    for node1 in nodes1:
        if is_same_node(node0, node1):
            same_node = {node0: node1}
            ret = _match_inner(node0, edge0, node1, edge1, is_same_node,
                               same_node, {node1, })
            if ret is not None:
                return same_node
    return None


def _match_inner(node0: _N0,
                 edge0: Callable[[_N0], Iterable[_N0]],
                 node1: _N1,
                 edge1: Callable[[_N1], Iterable[_N1]],
                 is_same_node: Callable[[_N0, _N1], bool],
                 same_node: dict[_N0, _N1],
                 searched1: set[_N1]) -> tuple[set[_N0], set[_N1]] | None:
    """一部のノード対応がわかっている2つのグラフが一致するか判定する.
    一致する場合はグラフ0からグラフ1のノード対応表を返す.

    Args:
        nodes0: グラフ0の起点ノード
        edge0: グラフ0のノードから隣接ノードを返す関数
        nodes1: node0に対応するグラフ1のノード
        edge1: グラフ1のノードから隣接ノードを返す関数
        is_same_node: ノードが一致する場合はTrueを返す関数
        same_node: グラフ0からグラフ1の同一ノードの対応表
        searched1: グラフ1の探索しないノード
    Returns:
        2つのグラフが一致する場合は
        探索した(グラフ0のノード集合, グラフ1のノード集合)
        一致しない場合はNone
    """
    neg0 = tuple(edge0(node0))
    neg1_p = tuple(edge1(node1))
    for neg1 in itertools.permutations(neg1_p):
        child_same0: set[_N0] = set()
        child_same1: set[_N1] = set()
        for n0, n1 in zip(neg0, neg1):
            if n0 in same_node:
                if n1 == same_node[n0]:
                    continue
                else:
                    break
            if n1 in searched1:
                break
            if is_same_node(n0, n1):
                child_same0.add(n0)
                child_same1.add(n1)
                same_node[n0] = n1
                searched1.add(n1)
                ret = _match_inner(n0, edge0, n1, edge1,
                                   is_same_node, same_node, searched1)
                if ret is None:
                    break
                child_same0.update(ret[0])
                child_same1.update(ret[1])
            else:
                break
        else:
            return (child_same0, child_same1)
        for n0 in child_same0:
            same_node.pop(n0)
        for n1 in child_same1:
            searched1.remove(n1)
    else:
        return None
