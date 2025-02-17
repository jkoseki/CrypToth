"""汎用のイテレータ関係"""
from collections.abc import Callable, Iterable, Iterator
from typing import Any, TypeVar


_P = TypeVar("_P")
_R = TypeVar("_R")
_T = TypeVar("_T")


def deep2_map(func: Callable[[_P], _R], src: Iterable[Iterable[_P]]
              ) -> Iterator[Iterator[_R]]:
    """2重のIterableにmapを適用する.

    Args:
        func: 要素の変換関数
        src: 2重のIterableで要素を保持する
    Returns:
        2重のIterableにmapを適用したもの
    """
    for src_1 in src:
        yield map(func, src_1)


def len_iterator(itr: Iterable[Any]) -> int:
    """イテレータの要素総数を返す.

    Args:
        itr: イテレータ
    Returns:
        イテレータの要素総数
    """
    count = 0
    for _ in itr:
        count += 1
    return count
