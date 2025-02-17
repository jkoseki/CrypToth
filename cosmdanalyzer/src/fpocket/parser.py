"""fpocketの出力パーサー"""
import array
from collections.abc import Iterator, Sequence
from typing import Any, IO
from ..solidcalc.typehint import Vector3f


def parse_fpocket(src_pdb: IO[str], src_score: IO[str]
                  ) -> Iterator[tuple[float, Sequence[Vector3f]]]:
    """fpocketの出力からスコアとポケットのボクセル座標を取得する.

    Args:
        src_pdb: PDBの文字データを保持する入力ストリーム
        src_score: _info.txtファイルの文字列を保持する入力ストリーム
    Returns:
        各ポケット毎の(スコア, ボクセル座標)
    """
    scores = parse_score_from_info(src_score)
    pocket_ids = list(scores.keys())
    pocket_ids.sort()
    pos = parse_pocket_position_from_pdb(src_pdb)
    for pocket_id in pocket_ids:
        yield (scores[pocket_id], pos[pocket_id])


def parse_pocket_position_from_pdb(src: IO[str]
                                   ) -> dict[int, list[Vector3f]]:
    """PDBファイルから原子名APOLで記録された要素の座標を残基毎に返す.

    Args:
        src: PDBの文字データを保持する入力ストリーム
    Returns:
        key=残基番号, value=ボクセル座標
    """
    res_pockets: dict[int, list[Vector3f]] = dict()
    for line in src:
        if (line.startswith('HETATM')
                and len(line) >= 16 and line[12:16] == 'APOL'):
            ret = _read_pdb_atom_line(line)
            if ret is None:
                continue
            pocket_n, pos = ret
            if pocket_n in res_pockets:
                res_pockets[pocket_n].append(pos)
            else:
                res_pockets[pocket_n] = [pos, ]
    return res_pockets


def _read_pdb_atom_line(line: str) -> tuple[int, Vector3f] | None:
    """PDBの ATOM | HETATM 行から必要な情報のみを抽出する.

    Args:
        line: PDBの ATOM | HETATM 行の文字列
    Returns:
        (残基番号,3次元座標) 読み込みに失敗した場合はNone
    """
    # atom_id = int(line[6:11])
    res_id = int(line[22:26])
    try:
        x = float(line[30:38])
        y = float(line[38:46])
        z = float(line[46:54])
        return (res_id, (x, y, z))
    except (ValueError, IndexError):
        return None


def parse_score_from_info(src: IO[str]) -> dict[int, float]:
    """_info.txtファイルからスコア情報を読み出す.

    Args:
        src: _info.txtファイルの文字列を保持する入力ストリーム
    Returns:
        Pockert 1 から順にスコアを返す.
    """
    scores: dict[int, float] = dict()
    info = _parse_tab_indent(src)
    for pocket, v in info.items():
        if not pocket.startswith('Pocket'):
            continue
        try:
            pocket_number = int(pocket[6:])
        except ValueError:
            continue
        scores[pocket_number] = v['Score']
    return scores


def _parse_tab_indent(src: IO[str]) -> dict[str, Any]:
    """Tabインデントでデータ構造を表した文字データから構造データを読み出す.

    Args:
        src: 入力文字データ
    Returns:
        構造データ
    """
    dict_stack: list[dict[str, Any]] = [dict(), ]
    level = 0
    for line in src:
        if len(line.strip()) == 0:
            continue
        line_level = _count_tab_indent(line)
        for _ in range(line_level, level):
            dict_stack.pop()
        level = line_level
        name, value = _parse_line(line)
        if len(value) > 0:
            dict_stack[-1][name] = value
        else:
            new_dict: dict[str, Any] = dict()
            dict_stack[-1][name] = new_dict
            dict_stack.append(new_dict)
            level += 1
    return dict_stack[0]


def _parse_line(line: str) -> tuple[str, str]:
    """データ名:値の1行をパースする."""
    sp = line.strip().split(':', 2)
    return (sp[0].strip(), sp[1].strip())


def _count_tab_indent(line: str) -> int:
    """先頭のタブの数を数える."""
    i = 0
    for _ in range(len(line)):
        if line[i] != '\t':
            break
        i += 1
    return i
