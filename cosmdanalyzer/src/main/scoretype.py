"""スコア計算の結果を保持する"""
from collections.abc import Callable, Collection, Iterable, Iterator, Sequence
import itertools
import math
import statistics
from .. import solidcalc
from ..scorecalc import convexity, compactness, protrusion, calcchargedensity
from ..solidcalc.typehint import Sphere, Vector3f
from . import fpocketscore


class ScoreSize:
    """複数のパッチの表面積のスコアを保持する"""

    def __init__(self,
                 patch_list: Sequence[Iterable[int]],
                 res_to_atoms: Callable[[int], Iterable[int]],
                 atom_to_vdw_radius: Callable[[int], float],
                 resolution: int,
                 calc_detail: bool = False):
        """

        Args:
            patch_list: パッチ毎に構成残基IDの集合を保持する
            res_to_atoms: 残基IDから構成原子IDの集合を返す関数
            atom_to_vdw_radius: 原子IDから原子半径を返す関数
            resolution: 原子表面を多面体で近似するときの頂点数
            calc_detail: 詳細情報を計算する場合はTrue
        """
        self._scores = tuple(MeanScore(calc_detail) for _ in patch_list)
        self._res_to_atoms = res_to_atoms
        self._atom_to_radius = atom_to_vdw_radius
        self._patch_list = patch_list
        self._resolution = resolution

    def add_frame(self, atom_to_pos: Callable[[int], Vector3f],
                  vdw_col_sphere: Callable[[int], Iterable[int]]) -> None:
        """1フレームのデータを追加する.

        Args:
            atom_to_pos: 原子IDから原子座標を返す関数
            vdw_col_sphere: 原子IDからファンデルワールス半径で見た場合に
                            衝突している原子ID集合を返す関数
        """
        for score, patch_res_ids in zip(self._scores, self._patch_list):
            score.add_score(sum(map(
                solidcalc.create_one_area_in_multi_spheres_func(
                    (lambda a: (atom_to_pos(a), self._atom_to_radius(a))),
                    vdw_col_sphere, self._resolution),
                _res_to_atom_iterator(patch_res_ids, self._res_to_atoms),
            )))

    def get_result(self) -> Iterator[float]:
        """スポット毎の平均スコアを返す

        Returns:
            スポット毎の平均スコア
        """
        for score in self._scores:
            yield score.get_result()

    def get_detail_result(self) -> Iterator[
            tuple[float, float, float, float, float, float, float]]:
        """スポット毎のスコアの詳細値を返す

        Returns:
            スポット毎のにススコアの以下の値
            (平均, 分散, 最小値, 第一四分位数, 中央値, 第三四分位数, 最大値)
        """
        for score in self._scores:
            yield score.get_detail_result()


class ScoreProtrusion:
    """複数のパッチのProtrusionのスコアを保持する"""

    def __init__(self,
                 patch_list: Sequence[Iterable[int]],
                 res_to_atoms: Callable[[int], Iterable[int]],
                 calc_detail: bool = False):
        """

        Args:
            patch_list: パッチ毎に構成残基IDの集合を保持する
            res_to_atoms: 残基IDから構成原子IDの集合を返す関数
            calc_detail: 詳細情報を計算する場合はTrue
        """
        self._scores = tuple(MeanScore(calc_detail) for _ in patch_list)
        self._res_to_atoms = res_to_atoms
        self._patch_list = patch_list

    def add_frame(self, atom_to_pos: Callable[[int], Vector3f],
                  distance_atom_in_sphere: Callable[[Sphere], Iterable[float]],
                  ) -> None:
        """1フレームのデータを追加する.

        Args:
            atom_to_pos: 原子IDから原子座標を返す関数
            distance_atom_in_sphere: 指定球内にある原子中心の
                                     球の中心からの距離を返す関数
        """
        for score, patch_res_ids in zip(self._scores, self._patch_list):
            p = protrusion.calc_patch_protrusion(
                    patch_res_ids,
                    atom_to_pos,
                    self._res_to_atoms,
                    distance_atom_in_sphere)
            if p is None:
                p = 0.0
            score.add_score(p)

    def get_result(self) -> Iterator[float]:
        """スポット毎の平均スコアを返す

        Returns:
            スポット毎の平均スコア
        """
        for score in self._scores:
            yield score.get_result()

    def get_detail_result(self) -> Iterator[
            tuple[float, float, float, float, float, float, float]]:
        """スポット毎のスコアの詳細値を返す

        Returns:
            スポット毎のにススコアの以下の値
            (平均, 分散, 最小値, 第一四分位数, 中央値, 第三四分位数, 最大値)
        """
        for score in self._scores:
            yield score.get_detail_result()


class ScoreConvexity:
    """複数のパッチのConvexityのスコアを保持する"""

    def __init__(self,
                 patch_list: Sequence[Iterable[int]],
                 res_to_atoms: Callable[[int], Iterable[int]],
                 res_to_ca: Callable[[int], int],
                 atom_to_res: Callable[[int], int],
                 atom_to_weight: Callable[[int], float],
                 neighbor_residue_distance: float,
                 calc_detail: bool = False):
        """

        Args:
            patch_list: パッチ毎に構成残基IDの集合を保持する
            res_to_atoms: 残基IDから構成原子IDの集合を返す関数
            res_to_ca: 残基IDからCα原子のIDを返す関数
            atom_to_res: 原子IDから所属する残基IDを返す関数
            atom_to_weight: 原子IDから原子量を返す関数
            neighbor_residue_distance: 残基を隣接していると見なす
                                       最近傍原子中心間の距離
            calc_detail: 詳細情報を計算する場合はTrue
        """
        self._scores = tuple(MeanScore(calc_detail) for _ in patch_list)
        self._res_to_atoms = res_to_atoms
        self._res_to_ca = res_to_ca
        self._patch_list = patch_list
        self._atom_to_residue = atom_to_res
        self._atom_to_weight = atom_to_weight
        self._neg_res_distance = neighbor_residue_distance

    def add_frame(self, atom_to_pos: Callable[[int], Vector3f],
                  atom_in_sphere: Callable[[Sphere], Iterable[int]],
                  is_exposed_atom: Callable[[int], bool],
                  ) -> None:
        """1フレームのデータを追加する.

        Args:
            atom_to_pos: 原子IDから原子座標を返す関数
            atom_in_sphere: 指定球内に原子中心ある原子のIDを返す関数
            is_exposed_atom: 原子が溶媒露出している場合はTrueを返す関数
        """
        for score, patch_res_ids in zip(self._scores, self._patch_list):
            score.add_score(convexity.calc_patch_convexity(
                patch_res_ids,
                self._res_to_atoms,
                self._res_to_ca,
                self._atom_to_residue,
                atom_to_pos,
                self._atom_to_weight,
                atom_in_sphere,
                is_exposed_atom,
                self._neg_res_distance))

    def get_result(self) -> Iterator[float]:
        """スポット毎の平均スコアを返す

        Returns:
            スポット毎の平均スコア
        """
        for score in self._scores:
            yield score.get_result()

    def get_detail_result(self) -> Iterator[
            tuple[float, float, float, float, float, float, float]]:
        """スポット毎のスコアの詳細値を返す

        Returns:
            スポット毎のにススコアの以下の値
            (平均, 分散, 最小値, 第一四分位数, 中央値, 第三四分位数, 最大値)
        """
        for score in self._scores:
            yield score.get_detail_result()


class ScoreCompactness:
    """複数のパッチのCompactnessのスコアを保持する"""

    def __init__(self,
                 patch_list: Sequence[Iterable[int]],
                 res_to_atoms: Callable[[int], Iterable[int]],
                 calc_detail: bool = False):
        """

        Args:
            patch_list: パッチ毎に構成残基IDの集合を保持する
            res_to_atoms: 残基IDから構成原子IDの集合を返す関数
            calc_detail: 詳細情報を計算する場合はTrue
        """
        self._scores = tuple(MeanScore(calc_detail) for _ in patch_list)
        self._res_to_atoms = res_to_atoms
        self._patch_list = patch_list

    def add_frame(self, atom_to_pos: Callable[[int], Vector3f],
                  is_exposed_atom: Callable[[int], bool],
                  ) -> None:
        """1フレームのデータを追加する.

        Args:
            atom_to_pos: 原子IDから原子座標を返す関数
            is_exposed_atom: 原子が溶媒露出している場合はTrueを返す関数
        """
        for score, patch_res_ids in zip(self._scores, self._patch_list):
            score.add_score(compactness.calc_patch_compactness(
                patch_res_ids,
                self._res_to_atoms,
                atom_to_pos,
                is_exposed_atom,
            ))

    def get_result(self) -> Iterator[float]:
        """スポット毎の平均スコアを返す

        Returns:
            スポット毎の平均スコア
        """
        for score in self._scores:
            yield score.get_result()

    def get_detail_result(self) -> Iterator[
            tuple[float, float, float, float, float, float, float]]:
        """スポット毎のスコアの詳細値を返す

        Returns:
            スポット毎のにススコアの以下の値
            (平均, 分散, 最小値, 第一四分位数, 中央値, 第三四分位数, 最大値)
        """
        for score in self._scores:
            yield score.get_detail_result()


class ScoreChargeDensity:
    """複数のパッチのcharge densityのスコアを保持する"""

    def __init__(self,
                 patch_list: Sequence[Iterable[int]],
                 res_to_atoms: Callable[[int], Iterable[int]],
                 atom_to_vdw_radius: Callable[[int], float],
                 solvent_radius: float,
                 atom_to_charge: Callable[[int], float],
                 resolution: int,
                 calc_detail: bool = False):
        """

        Args:
            patch_list: パッチ毎に構成残基IDの集合を保持する
            res_to_atoms: 残基IDから構成原子IDの集合を返す関数
            atom_to_vdw_radius: 原子IDから原子半径を返す関数
            solvent_radius: 溶媒半径
            atom_to_charge: 原子IDから電荷を返す関数
            resolution: 原子表面を多面体で近似するときの頂点数
            calc_detail: 詳細情報を計算する場合はTrue
        """
        self._scores = tuple(MeanScore(calc_detail) for _ in patch_list)
        self._res_to_atoms = res_to_atoms
        self._patch_list = patch_list
        self._atom_to_radius = atom_to_vdw_radius
        self._solvent_radius = solvent_radius
        self._atom_to_charge = atom_to_charge
        self._resolution = resolution

    def add_frame(self, atom_to_pos: Callable[[int], Vector3f],
                  is_exposed_atom: Callable[[int], bool],
                  as_col_sphere: Callable[[int], Iterable[int]],
                  ) -> None:
        """1フレームのデータを追加する.

        Args:
            atom_to_pos: 原子IDから原子座標を返す関数
            is_exposed_atom: 原子が溶媒露出している場合はTrueを返す関数
            as_col_sphere: 原子IDから溶媒接触面で見た場合に
                           衝突している原子ID集合を返す関数
        """
        for score, patch_res_ids in zip(self._scores, self._patch_list):
            patch_atom_ids = itertools.chain.from_iterable(
                map(self._res_to_atoms, patch_res_ids))
            score.add_score(calcchargedensity.calc_atoms_charge_density(
                filter(is_exposed_atom, patch_atom_ids),
                (lambda a: (atom_to_pos(a),
                            self._atom_to_radius(a)
                            + self._solvent_radius)),
                self._atom_to_charge,
                as_col_sphere,
                self._resolution,
            ))

    def get_result(self) -> Iterator[float]:
        """スポット毎の平均スコアを返す

        Returns:
            スポット毎の平均スコア
        """
        for score in self._scores:
            yield score.get_result()

    def get_detail_result(self) -> Iterator[
            tuple[float, float, float, float, float, float, float]]:
        """スポット毎のスコアの詳細値を返す

        Returns:
            スポット毎のにススコアの以下の値
            (平均, 分散, 最小値, 第一四分位数, 中央値, 第三四分位数, 最大値)
        """
        for score in self._scores:
            yield score.get_detail_result()


class ScoreFpocket:
    """複数のパッチのfpocketのスコアを保持する"""

    def __init__(self,
                 hotspot_list: Sequence[Collection[Vector3f]],
                 grid_size: float,
                 rate_threshold: float,
                 calc_detail: bool = False):
        """

        Args:
            hotspot_list: パッチ毎のホットスポット座標の集合の配列
            grid_size:
            rate_threshold: [0,1]hotspotに含まれるfpocketの座標点の割合が
                            rate未満の場合は無視する.
            calc_detail: 詳細情報を計算する場合はTrue
        """
        self._hotspot_list = hotspot_list
        self._grid_size = grid_size
        self._scores = tuple(MeanScore(calc_detail) for _ in hotspot_list)
        self._rate_threshold = rate_threshold

    def add_frame(self,
                  fpockets: Iterable[tuple[float, Iterable[Vector3f]]] | None
                  ) -> None:
        """1フレームのデータを追加する.

        Args:
            fpockets: ポケット毎のfpocketのスコアと座標集合,
                      fpocketの値が存在しない場合はNone
        """
        if fpockets is not None:
            for score, hotspot in zip(self._scores, self._hotspot_list):
                score.add_score(fpocketscore.fpocket_score(
                    hotspot, self._grid_size, fpockets, self._rate_threshold))
        else:
            for score in self._scores:
                score.add_score(0.0)

    def get_result(self) -> Iterator[float]:
        """スポット毎の平均スコアを返す

        Returns:
            スポット毎の平均スコア
        """
        for score in self._scores:
            yield score.get_result()

    def get_detail_result(self) -> Iterator[
            tuple[float, float, float, float, float, float, float]]:
        """スポット毎のスコアの詳細値を返す

        Returns:
            スポット毎のにススコアの以下の値
            (平均, 分散, 最小値, 第一四分位数, 中央値, 第三四分位数, 最大値)
        """
        for score in self._scores:
            yield score.get_detail_result()


class MeanScore:
    """複数のスコアの平均からなるスコアを管理する.
    Noneを無視する.
    """

    def __init__(self, calc_detail: bool):
        """

        Args:
            calc_detail: 詳細情報も計算する場合はTrue
        """
        self._detail_buf = self._detail_buf = list() if calc_detail else None
        self._num_scores = 0
        self._sum_scores = 0.0

    def add_score(self, new_score: float):
        """スコアを追加する

        Args:
            new_score: 追加するスコア
        """
        if new_score is not None:
            self._sum_scores += new_score
            self._num_scores += 1
        if self._detail_buf is not None:
            self._detail_buf.append(new_score)

    def get_result(self) -> float | None:
        """スコアの平均を返す.

        Returns:
            スコアの平均
        """
        if self._num_scores > 0:
            return self._sum_scores / self._num_scores
        return None

    def get_detail_result(self) -> tuple[
            float, float, float, float, float, float, float]:
        """スコアの傾向を表す詳細データを返す.
        コンストラクタでcalc_detailをTrueに指定していない場合はエラーを返す.

        Returns:
            (平均, 分散, 最小値, 第一四分位数, 中央値, 第三四分位数, 最大値)
            add_score()が呼び出されていない場合はNone
        """
        if self._num_scores == 0:
            return None
        sorted_score = sorted(self._detail_buf)
        return (self._sum_scores / self._num_scores,
                statistics.variance(sorted_score),
                sorted_score[0],
                sorted_percentile(sorted_score, 0.25),
                sorted_percentile(sorted_score, 0.5),
                sorted_percentile(sorted_score, 0.75),
                sorted_score[-1],
                )


def sorted_percentile(sorted_data: Sequence[float], per: float):
    """ソート済み配列のパーセンタイルを計算する.

    Args:
        sorted_data: ソート済み配列
        per: [0, 1]
    Returns:
        パーセンタイル
    """
    n_data = len(sorted_data)
    pos = (n_data - 1) * per
    first = int(pos)
    if first == n_data - 1:
        return sorted_data[-1]
    second_rate = pos - first
    return (sorted_data[first] * (1 - second_rate)
            + sorted_data[first + 1] * second_rate)


def _res_to_atom_iterator(
        res_ids: Iterable[int], res_to_atom: Callable[[int], Iterable[int]]
) -> Iterator[int]:
    """残基のインデックス集合から原子のインデックス集合に変換する."""
    for res_id in res_ids:
        for atom_id in res_to_atom(res_id):
            yield atom_id
