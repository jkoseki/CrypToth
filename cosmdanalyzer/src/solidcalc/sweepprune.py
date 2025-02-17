"""3次元のスイープアンドプルーンによる球同士の衝突判定"""
from collections.abc import Callable, Hashable, Iterable, Iterator
import operator
from typing import TypeVar
from .typehint import Sphere


_ID = TypeVar("_ID", bound=Hashable)


def create_strict_collided_dict(
        sphere_ids: Iterable[_ID],
        id_to_sphere: Callable[[_ID], Sphere],
        ) -> dict[_ID, list[_ID]]:
    """すべての球同士の厳密な衝突判定を行い結果を辞書で返す.

    Args:
        sphere_ids: 球の識別子集合
        id_to_sphere: 識別子から球情報への変換関数
    Return:
        {球の識別子, 衝突している球の識別子集合}
    """
    broad_col = sweep_and_prune(
            map(lambda i: (i, id_to_sphere(i)), sphere_ids))
    col_dict: dict[_ID, list[_ID]] = dict()
    for id0, id1 in strict_collision(broad_col, id_to_sphere):
        if id0 not in col_dict:
            col_dict[id0] = list()
        col_dict[id0].append(id1)
    return col_dict


def strict_collision(broad_col: Iterable[tuple[_ID, _ID]],
                     id_to_sphere: Callable[[_ID], Sphere],
                     ) -> Iterator[tuple[_ID, _ID]]:
    """球同士の厳密な衝突判定を候補ペアから行う.

    Args:
        broad_col: 衝突候補ペアを表す球の_ID対
        id_to_sphere: _IDから球情報への変換関数
    Return:
        厳密に衝突している球の_ID対を返すイテレータ
    """
    for o0, o1 in broad_col:
        p0, r0 = id_to_sphere(o0)
        p1, r1 = id_to_sphere(o1)
        if (((p0[0] - p1[0])**2 + (p0[1] - p1[1])**2 + (p0[2] - p1[2])**2)
                < (r0 + r1)**2):
            yield (o0, o1)


def sweep_and_prune(spheres: Iterable[tuple[_ID, Sphere]],
                    comp_lt: Callable[[_ID, _ID], bool] | None = None,
                    ) -> set[tuple[_ID, _ID]]:
    """球同士のおおまかな衝突判定をスイープアンドプルーンで行う.

    Args:
        spheres: 球の集合 (_ID, 球のサイズ)
        comp_lt: 1つめの要素が小さい時にTrueを返す関数
    Returns:
        衝突しているオブジェクトの_ID対の集合,
        comp_ltを設定した場合 _ID対は必ず_IDが小さいほうが前にくる.
        comp_ltがNoneの場合は順序を入れ替えた同じ組が計2回出力される.
    """
    ax = _create_sweep_and_prune_axis(spheres)
    col_obj = set(_sweep_and_prune_1d(map(lambda el: el[1:], ax[0]),
                                      comp_lt))
    for dim in range(1, 3):
        col_obj &= set(_sweep_and_prune_1d(
            map(lambda el: el[1:], ax[dim]), comp_lt))
    return col_obj


def _sweep_and_prune_1d(ax: Iterable[tuple[_ID, bool]],
                        comp_lt: Callable[[_ID, _ID], bool] | None,
                        ) -> Iterator[tuple[_ID, _ID]]:
    """1次元の整列済み軸情報を使用して
    1次元のスイープアンドプルーンによるおおまかな衝突判定を行う.

    Args:
        ax: 整列済み軸情報 (ObjectID, 開始点はTrue 終了点はFalse)
        comp_lt: 1つめの要素が小さい時にTrueを返す関数
    Returns:
        衝突しているオブジェクトの_ID対の集合,
        comp_ltを設定した場合 _ID対は必ず_IDが小さいほうが前にくる.
        comp_ltがNoneの場合は順序を入れ替えた同じ組が計2回出力される.
    """
    cur_opened: set[_ID] = set()
    for obj_id, is_open in ax:
        if is_open:
            for opened_obj_id in cur_opened:
                if comp_lt is not None:
                    yield ((obj_id, opened_obj_id)
                           if comp_lt(obj_id, opened_obj_id)
                           else (opened_obj_id, obj_id))
                else:
                    yield (obj_id, opened_obj_id)
                    yield (opened_obj_id, obj_id)
            cur_opened.add(obj_id)
        else:
            cur_opened.discard(obj_id)


def _create_sweep_and_prune_axis(
        spheres: Iterable[tuple[_ID, Sphere]],
        ) -> tuple[list[tuple[float, _ID, bool]],
                   list[tuple[float, _ID, bool]],
                   list[tuple[float, _ID, bool]]]:
    """各軸の整列済み軸情報を作成する.

    Args:
        spheres: 球の集合 (_ID, 球のサイズ)
    Returns:
        (座標, ObjectID, 開始点はTrue)の情報を各軸の座標でソートしたもの
    """
    ax: tuple[list[tuple[float, _ID, bool]],
              list[tuple[float, _ID, bool]],
              list[tuple[float, _ID, bool]]] = ([], [], [])
    for obj_id, (pos, r) in spheres:
        for dim in range(3):
            ax[dim].append((pos[dim] - r, obj_id, True))
            ax[dim].append((pos[dim] + r, obj_id, False))
    for dim in range(3):
        ax[dim].sort(key=operator.itemgetter(0))
    return ax
