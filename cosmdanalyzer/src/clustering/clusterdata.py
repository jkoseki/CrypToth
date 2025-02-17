"""クラスタ情報を保持するデータ型を定義する. """
from collections.abc import Hashable, Iterable, Sequence
from typing import Generic, TypeVar

_EL = TypeVar('_EL', bound=Hashable)


class ClusterTable(Generic[_EL]):
    """要素の識別子からクラスタの識別子への変換表を
    アップデート可能な形式で保持する."""

    def __init__(self, elements: Iterable[_EL] | None = None):
        """初期化

        Args:
            elements: 要素の集合
        """
        self._element_links: dict[_EL, _EL] = dict()
        if elements is not None:
            for e in elements:
                self._element_links[e] = e

    def add_element(self, element: _EL) -> None:
        """要素を追加する

        Args:
            element: 追加する要素
        """
        self._element_links[element] = element

    def extend_elements(self, elements: Iterable[_EL]) -> None:
        """要素集合を追加する

        Args:
            elements: 追加する要素集合
        """
        for el in elements:
            self._element_links[el] = el

    def concat_cluster(self, element1_id: _EL, element2_id: _EL) -> bool:
        """異なるクラスタに含まれる要素を指定してクラスタを結合する.

        Args:
            element1_id: クラスタに含まれる要素の識別子
            element2_id: クラスタに含まれる要素の識別子
        Returns:
            クラスタの結合が行われた場合はTrue, それ以外はFalse
        """
        try:
            cluster1_id = self.get_cluster(element1_id)
            cluster2_id = self.get_cluster(element2_id)
            if cluster1_id == cluster2_id:
                return False
            self._element_links[cluster2_id] = cluster1_id
        except KeyError:
            return False
        return True

    def get_cluster(self, element_id: _EL) -> _EL:
        """要素の識別子から属するクラスタのに含まれる代表要素の識別子を返す.
        同じクラスタに属する要素は必ず同じ代表要素を返す.

        Args:
            elementId: 要素の識別子
        Returns:
            クラスタの代表要素の識別子
        """
        prev_ele = element_id
        while True:
            ele = self._element_links[prev_ele]
            if ele == prev_ele:
                return ele
            prev_ele = ele

    def create_element_to_cluster_dict(self) -> dict[_EL, int]:
        """要素の識別子からクラスタの識別子を返す辞書を生成する.

        Returns:
            要素の識別子からクラスタの識別子を返す辞書
        """
        element_to_index = self._give_element_to_index()
        cluster_dict = dict()
        for key in self._element_links.keys():
            cluster_dict[key] = element_to_index[self.get_cluster(key)]
        return cluster_dict

    def create_cluster_to_element_sequence(
            self) -> Sequence[Sequence[_EL]]:
        """クラスタ毎に所属する要素の配列を生成する.

        Returns:
            クラスタ毎に所属する要素の配列
        """
        element_to_index = self._give_element_to_index()
        clusters: tuple[list[_EL], ...] = tuple(
                list() for _ in range(len(element_to_index)))
        for key in self._element_links.keys():
            clusters[element_to_index[self.get_cluster(key)]].append(key)
        return clusters

    def _give_element_to_index(self) -> dict[_EL, int]:
        """クラスタの代表要素の識別子から
        クラスタのインデックスを返す辞書を生成する.

        Returns:
            クラスタの識別子からクラスタのインデックスを返す辞書
        """
        cluster_id = 0
        element_to_index: dict[_EL, int] = dict()
        for key, val in self._element_links.items():
            if key == val:
                element_to_index[key] = cluster_id
                cluster_id += 1
        return element_to_index
