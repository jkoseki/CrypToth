from collections.abc import Iterator
import csv
import os
from typing import IO


class LennardJonesCharge:
    """Lennard Jones Parametersの電荷テーブル"""

    def __init__(self, src_csv: str | bytes | os.PathLike):
        """Lennard Jones Parametersの電荷テーブルをCSVファイルから作成する.

        Args:
            src_csv: ヘッダにAtomTypeとeを含むCSVファイルのパス
        """
        with open(src_csv) as f:
            reader = csv.reader(f, skipinitialspace=True)
            header = next(reader, None)
            if header is not None:
                for i in range(len(header)):
                    header[i] = header[i].lower()
                type_idx = header.index('atomtype')
                charge_idx = header.index('e')
                self._type_to_charge = dict()
                for row in reader:
                    self._type_to_charge[row[type_idx]] = (
                            float(row[charge_idx]))

    def __call__(self, atom_type: str) -> float | None:
        """4文字原子タイプから電荷を返す.

        Args:
            atom_type: 4文字表記の原子タイプ
        Returns:
            電荷, 対応する値が存在しない場合はNone
        """
        if atom_type not in self._type_to_charge:
            for i in range(3, 0, -1):
                if not (atom_type[i].isdigit() or atom_type[i].isspace()):
                    break
            else:
                return None
            atom_type = atom_type[:i + 1] + (' ' * (3 - i))
        return self._type_to_charge.get(atom_type, None)


def load_glomacs_charge_from_file(
        src: str | bytes | os.PathLike,
        removes_h: bool
        ) -> Iterator[tuple[str,
                            list[tuple[int, float]],
                            dict[int, list[int]]]]:
    """Glomacsの電荷形式のファイルからグラフ構造を作成する.

    Args:
        src:
    Returns:
        残基名
        原子情報(原子番号, 電荷)
        結合情報(原子index,隣接原子index)
    """
    with open(src) as f:
        for v in load_charge(f, removes_h):
            yield v


_GRO_TO_PDB_NAME = {
    'HSD': 'HID',
    'HSE': 'HIE',
    'HSP': 'HIP',
}


def load_charge(src: IO[str], removes_h: bool,
                ) -> Iterator[tuple[str,
                                    list[tuple[int, float]],
                                    dict[int, list[int]]]]:
    """Glomacsの電荷形式の文字ストリームからグラフ構造を作成する.

    Args:
        src:
    Returns:
        残基名
        原子情報(原子番号, 電荷)
        結合情報(原子index,隣接原子index)
    """
    src_line = src.readline()
    while src_line:
        line = _clip_comment(src_line).rstrip()
        if line.startswith('['):
            name = line[1:line.find(']')].strip()
            atoms, bonds, src_line = _load_amino(src, removes_h)
            if bonds is not None:
                name = _GRO_TO_PDB_NAME.get(name, name)
                yield (name, atoms, bonds)
            continue
        src_line = src.readline()


def _load_amino(src: IO[str], removes_h: bool
                ) -> tuple[list[tuple[int, float]], dict[int, list[int]], str]:
    """

    Returns:
        原子情報(原子番号, 電荷)
        結合情報(原子index,隣接原子index)
        最後に読んだ行, 最後まで読んだ場合はNone
    """
    src_line = src.readline()
    atoms = None
    dict_atoms = dict()
    ret_atoms: list[tuple[int, float]] = list()
    bonds = None
    while src_line:
        line = _clip_comment(src_line)
        if line.startswith('['):
            break
        line = line.strip()
        if line.startswith('[ atoms ]'):
            atoms, src_line = _load_atoms(src, removes_h)
            for i, atom in enumerate(atoms):
                ret_atoms.append((atom[1], atom[2]))
                dict_atoms[atom[0]] = i
            continue
        if line.startswith('[ bonds ]'):
            bonds, src_line = _load_bonds(src, dict_atoms)
            continue
        src_line = src.readline()
    else:
        src_line = None
    return (ret_atoms, bonds, src_line)


def _load_atoms(src: IO[str], removes_h: bool
                ) -> tuple[list[tuple[str, int, float]], str]:
    """アミノ酸項目のatomsを読み込む.

    Args:
        src:
        removes_h: 水素を削除する場合はTrue
    Returns:
        (原子名, 原子番号, 電荷)の集合
        最後に読んだ行, 最後まで読んだ場合はNone
    """
    atoms = []
    for src_line in src:
        line = _clip_comment(src_line).strip()
        if len(line) == 0:
            continue
        if line.startswith('['):
            break
        sp = line.split()
        atomic_num = _symbol_to_atomic_num.get(_to_symbol(sp[1]), -1)
        if (atomic_num == -1) or (removes_h and (atomic_num == 1)):
            continue
        atoms.append((sp[0], atomic_num, float(sp[2])))
    else:
        src_line = None
    return (atoms, src_line)


def _load_bonds(src: IO[str], atoms: dict[str, int],
                ) -> tuple[dict[int, list[int]], str]:
    """アミノ酸項目のbondsを読み込む.

    Args:
        src:
        atoms: 原子名から原子indexへの変換表
    Returns:
        原子indexから隣接原子indexの対応表
        最後に読んだ行, 最後まで読んだ場合はNone
    """
    bonds = dict()
    for src_line in src:
        line = _clip_comment(src_line).strip()
        if len(line) == 0:
            continue
        if line.startswith('['):
            break
        sp = line.split()
        if (sp[0] not in atoms) or (sp[1] not in atoms):
            continue
        n0 = atoms[sp[0]]
        n1 = atoms[sp[1]]
        if n0 in bonds:
            bonds[n0].append(n1)
        else:
            bonds[n0] = [n1, ]
        if n1 in bonds:
            bonds[n1].append(n0)
        else:
            bonds[n1] = [n0, ]
    else:
        src_line = None
    return (bonds, src_line)


def _clip_comment(line: str) -> str:
    i = line.find(';')
    if i >= 0:
        return line[:i]
    return line


def _to_symbol(atom_type: str) -> str:
    return atom_type[0]


_symbol_to_atomic_num = {
        'H': 1,
        'He': 2,
        'Li': 3,
        'Be': 4,
        'B': 5,
        'C': 6,
        'N': 7,
        'O': 8,
        'F': 9,
        'Ne': 10,
        'Na': 11,
        'Mg': 12,
        'Al': 13,
        'Si': 14,
        'P': 15,
        'S': 16,
        'Cl': 17,
        'Ar': 18,
        'K': 19,
        'Ca': 20,
        'Sc': 21,
        'Ti': 22,
        'V': 23,
        'Cr': 24,
        'Mn': 25,
        'Fe': 26,
        'Co': 27,
        'Ni': 28,
        'Cu': 29,
        'Zn': 30,
        'Ga': 31,
        'Ge': 32,
        'As': 33,
        'Se': 34,
        'Br': 35,
        'Kr': 36,
        'Rb': 37,
        'Sr': 38,
        'Y': 39,
        'Zr': 40,
        'Nb': 41,
        }
