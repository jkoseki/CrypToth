"""VP木による近傍探索"""
from __future__ import annotations
from collections.abc import Callable, Iterable, Sequence
from sys import float_info
from typing import Generic, TypeVar


_V = TypeVar("_V")


class _DistanceAndValue(Generic[_V]):

    def __init__(self, distance: float, value: _V):
        self.d = distance
        self.v = value


class _SphereNode(Generic[_V]):

    def __init__(self, branch: _BranchNodeData[_V] | None):
        self.branch = branch


class _OtherNode(Generic[_V]):

    def __init__(self, branch: _BranchNodeData[_V] | None, value: _V):
        self.value = value
        self.branch = branch


class _BranchNodeData(Generic[_V]):

    def __init__(self, radius: float,
                 sphere_node: _SphereNode[_V], other_node: _OtherNode[_V]):
        self.radius = radius
        self.sphere_node = sphere_node
        self.other_node = other_node


class VpTree(Generic[_V]):
    """近傍探索用のVP木を表す"""

    def __init__(self, values: Iterable[_V],
                 distance_func: Callable[[_V, _V], float]):
        """VP木を作成する.

        Args:
            values: VP木の要素
            distance_func: 距離を定義する関数
        """
        values = iter(values)
        self._root: _OtherNode[_V] | None = None
        try:
            standard_v = next(values)
            self._distance_func = distance_func
            distance_array = []
            for v in values:
                distance_array.append(_DistanceAndValue[_V](0.0, v))
            self._root = _create_other_node(
                distance_func,
                standard_v,
                distance_array,
                0, len(distance_array),
            )
        except StopIteration:
            pass

    def nearest_neighbor(self, query: _V, thresthold: float | None = None
                         ) -> tuple[float, _V | None]:
        """最近傍点を検索する.

        Args:
            query: この要素の最近傍要素を探す
            thresthold: 指定距離未満の要素のみ探索する
        Returns:
            (最近傍要素の距離, 最近傍要素)
            thresthold未満の要素が見つからなかった場合は(thresthold, None)
        """
        if thresthold is None:
            thresthold = float_info.max
        if self._root is not None:
            update_func = _NearestNeighborUpdateFunc[_V]()
            min_d = _search_other_node(
                    self._distance_func, query,
                    self._root,
                    thresthold, update_func)
            return (min_d, update_func.get())
        return (thresthold, None)

    def neighbors(self, query: _V, thresthold: float
                  ) -> Sequence[tuple[float, _V]]:
        """queryから指定距離未満の要素を検索する.

        Args:
            query: この要素の最近傍要素を探す
            thresthold: 指定距離未満の要素のみ探索する
        Returns:
            条件を満たす要素の
            (queryからの距離, 要素)
        """
        update_func = _NeighborUpdateFunc[_V]()
        if self._root is not None:
            _search_other_node(
                    self._distance_func, query,
                    self._root,
                    thresthold, update_func)
        return update_func.get()

    def exists_neighbor(self, query: _V, thresthold: float
                        ) -> bool:
        """queryから指定距離未満の要素が存在する場合Trueを返す.

        Args:
            query: この要素の最近傍要素を探す
            thresthold: 指定距離未満の要素のみ探索する
        Returns:
            queryから指定距離未満の要素が存在する場合True
        """
        update_func = _ExistsUpdateFunc[_V]()
        if self._root is not None:
            _search_other_node(
                    self._distance_func, query,
                    self._root,
                    thresthold, update_func)
        return update_func.exists()


def _create_other_node(distance_func: Callable[[_V, _V], float],
                       standard_v: _V,
                       distance_array: list[_DistanceAndValue[_V]],
                       start: int, end: int,
                       ) -> _OtherNode[_V]:
    """VP木の球領域でないノードを作成する.

    Args:
        distance_func: 距離を定義する関数
        standard_v: このノードの基準要素
        distance_array: 距離と要素の配列(距離は使われない)
        start: distance_arrayの参照開始位置
        end: distance_arrayの参照終了位置
    Returns:
        球領域でないノード
    """
    if start == end:
        return _OtherNode[_V](None, standard_v)
    for i in range(start, end):
        dv = distance_array[i]
        dv.d = distance_func(standard_v, dv.v)
    distance_array[start:end] = sorted(
            distance_array[start:end], key=lambda v: v.d)
    mid_idx = int((end - start + 1) / 2) - 1 + start
    mid = distance_array[mid_idx]
    return _OtherNode[_V](
        _BranchNodeData[_V](
            mid.d,
            _create_sphere_node(distance_func, standard_v, distance_array,
                                start, mid_idx),
            _create_other_node(distance_func, mid.v, distance_array,
                               mid_idx + 1, end),
        ),
        standard_v,
    )


def _create_sphere_node(distance_func: Callable[[_V, _V], float],
                        standard_v: _V,
                        distance_array: list[_DistanceAndValue[_V]],
                        start: int, end: int) -> _SphereNode[_V]:
    """VP木の球領域のノードを作成する.

    Args:
        distance_func: 距離を定義する関数
        standard_v: このノードの基準要素
        distance_array: ソート済みの距離と要素の配列
        start: distance_arrayの参照開始位置
        end: distance_arrayの参照終了位置
    Returns:
        球領域のノード
    """
    if start == end:
        return _SphereNode[_V](None)
    mid_idx = int((end - start + 1) / 2) - 1 + start
    mid = distance_array[mid_idx]
    return _SphereNode[_V](
        _BranchNodeData[_V](
            mid.d,
            _create_sphere_node(distance_func, standard_v, distance_array,
                                start, mid_idx),
            _create_other_node(distance_func, mid.v, distance_array,
                               mid_idx + 1, end),
        ),
    )


def _search_other_node(
        distance_func: Callable[[_V, _V], float],
        query: _V, node: _OtherNode[_V],
        thresthold_d: float, update_func: Callable[[float, float, _V], float],
        ) -> float:
    """球領域でない枝ノードから最近傍点を検索する.

    Args:
        distance_func: 距離を定義する関数
        query: この要素の最近傍要素をノード以下から探す
        node: 球領域でない枝ノード
        thresthold_d: queryからのこの距離未満の要素のみ探索する
        update_func: 基準を満たす要素が見つかった場合に呼び出される関数
                     (しきい値, 距離, 要素)を受け取り新たなしきい値を返す.
    Returns:
        このノード以下のノードから計算された新たなthresthold_d
    """
    d = distance_func(query, node.value)
    if d < thresthold_d:
        thresthold_d = update_func(thresthold_d, d, node.value)
    if node.branch is None:
        return thresthold_d
    return _search_branch_node(
            distance_func, query, node.branch, d, thresthold_d, update_func)


def _search_sphere_node(
        distance_func: Callable[[_V, _V], float],
        query: _V, node: _SphereNode[_V],
        standard_d: float, thresthold_d: float,
        update_func: Callable[[float, float, _V], float],
        ) -> float:
    """球領域の枝ノードから最近傍点を検索する.

    Args:
        distance_func: 距離を定義する関数
        query: この要素の最近傍要素をノード以下から探す
        node: 球領域の枝ノード
        standard_d: nodeの基準点とqueryの距離
        thresthold_d: queryからのこの距離未満の要素のみ探索する
        update_func: 基準を満たす要素が見つかった場合に呼び出される関数
                     (しきい値, 距離, 要素)を受け取り新たなしきい値を返す.
    Returns:
        このノード以下のノードから計算された新たなthresthold_d
    """
    if node.branch is None:
        return thresthold_d
    return _search_branch_node(
            distance_func, query, node.branch, standard_d,
            thresthold_d, update_func)


def _search_branch_node(
        distance_func: Callable[[_V, _V], float], query: _V,
        branch: _BranchNodeData[_V],
        standard_d: float,
        thresthold_d: float,
        update_func: Callable[[float, float, _V], float],
        ) -> float:
    """枝ノードの子ノードから最近傍点を検索する.

    Args:
        distance_func: 距離を定義する関数
        query: この要素の最近傍要素をノード以下から探す
        branch: 枝ノードの情報
        standard_d: このノードの基準点とqueryの距離
        thresthold_d: queryからのこの距離未満の要素のみ探索する
        update_func: 基準を満たす要素が見つかった場合に呼び出される関数
                     (しきい値, 距離, 要素)を受け取り新たなしきい値を返す.
    Returns:
        このノード以下のノードから計算された新たなthresthold_d
    """
    if standard_d < branch.radius:
        thresthold_d = _search_sphere_node(
                distance_func, query, branch.sphere_node,
                standard_d, thresthold_d, update_func)
        if standard_d + thresthold_d >= branch.radius:
            thresthold_d = _search_other_node(
                distance_func, query, branch.other_node,
                thresthold_d, update_func)
    else:
        thresthold_d = _search_other_node(
            distance_func, query, branch.other_node, thresthold_d, update_func)
        if standard_d < thresthold_d + branch.radius:
            thresthold_d = _search_sphere_node(
                    distance_func, query, branch.sphere_node,
                    standard_d, thresthold_d, update_func)
    return thresthold_d


class _NearestNeighborUpdateFunc(Generic[_V]):
    """最近傍探索のupdate_func"""

    def __init__(self):
        self._v: _V | None = None

    def __call__(self, _: float, d: float, v: _V) -> float:
        self._v = v
        return d

    def get(self) -> _V | None:
        return self._v


class _NeighborUpdateFunc(Generic[_V]):
    """指定距離以内の要素を探索するupdate_func"""

    def __init__(self):
        self._ret: list[tuple[float, _V]] = list()

    def __call__(self, thresthold_d: float, d: float, v: _V) -> float:
        self._ret.append((d, v))
        return thresthold_d

    def get(self) -> list[tuple[float, _V]]:
        return self._ret


class _ExistsUpdateFunc(Generic[_V]):
    """指定距離以内の要素が存在するか調べるupdate_func"""

    def __init__(self):
        self._exists = False

    def __call__(self, thresthold_d: float, d: float, v: _V) -> float:
        self._exists = True
        return 0.0

    def exists(self) -> bool:
        return self._exists
