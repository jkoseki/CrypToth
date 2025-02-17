"""隣接点を同じクラスタと見なす"""
from collections import deque
from collections.abc import Callable, Iterable, Iterator
import itertools
from typing import Generic, TypeVar
from ..neighbors import vptree


_EL = TypeVar('_EL')


class _ElementInfo(Generic[_EL]):

    def __init__(self, element: _EL):
        self.element = element
        self.in_cluster = False


def neighbor(elements: Iterable[_EL],
             distance_func: Callable[[_EL, _EL], float],
             bandwidth: float,
             ) -> Iterator[Iterator[_EL]]:
    """隣接点を同じクラスタと見なす"""
    elements_info = tuple(map(lambda el: _ElementInfo(el), elements))
    tree = vptree.VpTree(elements_info, _InfoDistance(distance_func))
    for el in filter((lambda el: not el.in_cluster), elements_info):
        neg_d_els = tree.neighbors(el, bandwidth)
        el.in_cluster = True
        neg_els = map(_check_in_cluster,
                      filter(lambda el: not el.in_cluster,
                             map(lambda dv: dv[1], neg_d_els)))
        yield itertools.chain(
            iter((el.element, )),
            _neg_fill(deque(neg_els), tree, bandwidth)
            )


def _neg_fill(targets: deque[_ElementInfo[_EL]],
              neg_tree: vptree.VpTree[_ElementInfo[_EL]],
              bandwidth: float,
              ) -> Iterator[_EL]:
    while targets:
        cur_el = targets.popleft()
        yield cur_el.element
        neg_d_els = neg_tree.neighbors(cur_el, bandwidth)
        targets.extend(
            map(_check_in_cluster,
                filter(lambda el: not el.in_cluster,
                       map(lambda info: info[1], neg_d_els)
                       )
                )
            )


def _check_in_cluster(info: _ElementInfo[_EL]) -> _ElementInfo[_EL]:
    info.in_cluster = True
    return info


class _InfoDistance(Generic[_EL]):

    def __init__(self, distance_func: Callable[[_EL, _EL], float]):
        self._func = distance_func

    def __call__(self, el0: _ElementInfo[_EL], el1: _ElementInfo[_EL]
                 ) -> float:
        return self._func(el0.element, el1.element)
