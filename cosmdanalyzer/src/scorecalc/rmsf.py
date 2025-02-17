import math
import statistics
from collections.abc import Callable, Collection, Iterable, Iterator
from ..solidcalc.typehint import Vector3f
from ..solidcalc import vector3f


class AllPatchRmsfCalc:
    """パッチ集合のRMSFを計算する"""

    def __init__(self,
                 res_to_atoms: Callable[[int], Iterable[int]],
                 res_ids_set: Iterable[Collection[int]],
                 atom_to_weight: Callable[[int], float]):
        """

        Args:
            res_to_atoms: 残基IDを構成原子のID集合に変換する関数
            res_ids_set: 計算対象の残基ID集合
            atom_to_weight: 原子IDから質量を返す関数
        """
        self._patch_calc = tuple(
            PatchRmsfCalc(res_ids, res_to_atoms, atom_to_weight)
            for res_ids in res_ids_set)

    def add_frame(self, atom_to_pos: Callable[[int], Vector3f]) -> None:
        """1フレーム分の計算を行う.

        Args:
            atom_to_pos: 原子IDを原子座標に変換する関数
        """
        for calc in self._patch_calc:
            calc.add_frame(atom_to_pos)

    def get_result(self) -> Iterator[float]:
        """RMSFの計算結果を返す.

        Returns:
            パッチ毎の残基の平均RMSF, 残基がない場合は0.0
        """
        for calc in self._patch_calc:
            yield calc.get_result()


class PatchRmsfCalc:
    """残基集合の平均RMSFを計算する"""

    def __init__(self,
                 res_ids: Collection[int],
                 res_to_atoms: Callable[[int], Iterable[int]],
                 atom_to_weight: Callable[[int], float],
                 ):
        """
        Args:
            res_ids_set: 計算対象の残基ID集合
            res_to_atoms: 残基IDを構成原子のID集合に変換する関数
            atom_to_weight: 原子IDから質量を返す関数
        """
        self._res_calcs = tuple(
            ResidueRmsfCalc(tuple(res_to_atoms(res_id)), atom_to_weight)
            for res_id in res_ids)

    def add_frame(self, atom_to_pos: Callable[[int], Vector3f]) -> None:
        """1フレーム分の計算を原子座標から行う.

        Args:
            atom_to_pos: 原子indexから原子座標を返す関数
        """
        for calc in self._res_calcs:
            calc.add_frame(atom_to_pos)

    def get_result(self) -> float:
        """RMSFの計算結果を返す.

        Returns:
            残基の平均RMSF, 残基がない場合は0.0
        """
        if len(self._res_calcs) > 0:
            return statistics.mean(
                map(lambda calc: calc.get_result(), self._res_calcs))
        else:
            return 0.0


class ResidueRmsfCalc:
    """残基重心のRMSFを計算する"""

    def __init__(self, atom_ids: Collection[int],
                 atom_to_weight: Callable[[int], float]):
        self._rmsf_calc = PointRmsfCalc()
        self._atom_ids = atom_ids
        self._atom_to_weight = atom_to_weight

    def add_frame(self, atom_to_pos: Callable[[int], Vector3f]) -> None:
        cent = _calc_centroid((atom_to_pos(a), self._atom_to_weight(a))
                              for a in self._atom_ids)
        self._rmsf_calc.add_frame(cent)

    def get_result(self) -> float:
        return self._rmsf_calc.get_result()


class PointRmsfCalc:
    """1点のRMSFを計算する"""

    def __init__(self):
        """
        """
        self._sum = (0.0, 0.0, 0.0)
        self._sum2 = 0.0
        self._n_frames = 0

    def add_frame(self, pos: Vector3f) -> None:
        """1フレーム分の計算を行う.

        Args:
            pos: 追加する座標
        """
        self._n_frames += 1
        self._sum = vector3f.add(self._sum, pos)
        self._sum2 += vector3f.norm2(pos)

    def get_result(self) -> float:
        """RMSFの計算結果を返す.

        Returns:
            RMSF
        """
        return math.sqrt(
            (self._sum2 / self._n_frames)
            - sum((self._sum[i] / self._n_frames)**2 for i in range(3))
        )


def _mean_position(
        atom_ids: Iterable[int],
        atom_to_pos: Callable[[int], Vector3f]):
    count = 0
    v = (0.0, 0.0, 0.0)
    for pos in map(atom_to_pos, atom_ids):
        count += 1
        v = vector3f.add(v, pos)
    return vector3f.mul(v, 1.0 / count)


def _calc_centroid(points: Iterator[tuple[Vector3f, float]]) -> Vector3f:
    sum_w = 0.0
    sum_pos = (0.0, 0.0, 0.0)
    for pos, w in points:
        sum_w += w
        sum_pos = vector3f.add(sum_pos, vector3f.mul(pos, w))
    return vector3f.mul(sum_pos, 1.0 / sum_w)
