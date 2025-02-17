from typing import Generic, TypeVar
from collections.abc import Callable, Hashable, Mapping

_K = TypeVar('_K', bound=Hashable)
_V = TypeVar('_V')


def dict_to_function(src: Mapping[_K, _V]) -> Callable[[_K], _V]:
    """辞書型からキーを入力として値を返す関数を生成する.

    Args:
        src: 生成元の辞書
    Returns:
        入力として値を返す関数
    """
    return lambda k: src[k]


class BufferdFunction(Generic[_K, _V]):
    """計算結果を辞書で保持する関数オブジェクト"""

    def __init__(self, func: Callable[[_K], _V]):
        self._func = func
        self._dict: dict[_K, _V] = dict()

    def __call__(self, key: _K) -> _V:
        if key not in self._dict:
            self._dict[key] = self._func(key)
        return self._dict[key]
