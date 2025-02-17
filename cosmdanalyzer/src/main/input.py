"""入力関係"""
from collections.abc import Iterable, Iterator
import itertools
import pathlib
import re
from operator import itemgetter
from os import PathLike
from typing import IO, NamedTuple


def trajectory_pdb_string_filter(
        pdb_lines: Iterable[str],
) -> tuple[str, int]:
    """PDBストリームからタンパク質のみのPDB文字列とプローブ重原子数を取得する

    Args:
        pdb_lines: PDB形式の文字列を行ごとに返すイテレータ
    Returns:
        (タンパク質のみのPDB, プローブ重原子数)
    """
    line_buf: list[str] = list()
    protain_max_id = _first_frame_protein_filter(pdb_lines, line_buf)
    n_probe_atoms = _first_frame_probe_counter(pdb_lines)
    line_buf.append('ENDMDL\n')
    for line in pdb_lines:
        if line.startswith('ATOM  ') or line.startswith('HETATM'):
            if int(line[6:11]) < protain_max_id:
                line_buf.append(line)
        if line.startswith('TER'):
            if int(line[6:11]) == protain_max_id:
                line_buf.append(line)
        elif line.startswith('MODEL') or line.startswith('ENDMDL'):
            line_buf.append(line)
    line_buf.append('END')
    return (''.join(line_buf), n_probe_atoms)


def trajectory_pdb_stream_filter(
        pdb_src: IO[str]
) -> tuple[str, int]:
    """PDBストリームからタンパク質のみのPDB文字列とプローブ重原子数を取得する

    Args:
        pdb_src: PDB形式の文字ストリーム
    Returns:
        (タンパク質のみのPDB, プローブ重原子数)
    """
    return trajectory_pdb_string_filter(_io_to_line_iterator(pdb_src))


def trajectory_pdb_files_filter(src: Iterable[str | bytes | PathLike]
                                ) -> tuple[str, int]:
    """PDBファイル集合からタンパク質のみのPDB文字列とプローブ重原子数を取得する

    Args:
        src: PDB形式のファイル集合, 同じ原子集合の座標のみ異なるデータをもつ.
    Returns:
        (タンパク質のみのPDB, プローブ重原子数)
    """
    return trajectory_pdb_string_filter(
        itertools.chain.from_iterable(
            map(_io_to_line_iterator, _path_to_io_iterator(src))
        )
    )


def _first_frame_protein_filter(
        pdb_lines: Iterable[str], out_buf: list[str]) -> int:
    """トラジェクトリPDBの先頭フレームのタンパク質を読み込む

    Args:
        pdb_lines: PDB形式の文字列を行ごとに返すイテレータ
        out_buf: タンパク質原子のPDB行の出力先
    Returns:
        タンパク質原子の最大ID + 1
    """
    for line in pdb_lines:
        if (line.startswith('ATOM  ') or line.startswith('HETATM')
                or line.startswith('MODEL')):
            out_buf.append(line)
        if line.startswith('TER'):
            out_buf.append(line)
            return int(line[6:11])


def _first_frame_probe_counter(pdb_lines: Iterable[str]) -> int:
    """
    Args:
        pdb_lines: PDB形式の文字列を行ごとに返すイテレータ
    Returns:
        プローブ重原子数
    """
    n_probe_atoms = 0
    line = next(pdb_lines)
    probe_name = line[17:20]
    if line[76:78] != ' H':
        n_probe_atoms += 1
    for line in pdb_lines:
        if line.startswith('ATOM  ') or line.startswith('HETATM'):
            if (line[17:20] == probe_name) and (line[76:78] != ' H'):
                n_probe_atoms += 1
        elif line.startswith('ENDMDL'):
            break
    return n_probe_atoms


def _io_to_line_iterator(src: IO[str]) -> Iterator[str]:
    for line in src:
        yield line


def _path_to_io_iterator(src: Iterable[str | bytes | PathLike]
                         ) -> Iterator[IO[str]]:
    for path in src:
        with open(path, mode='r') as f:
            yield f


class SystemInfo(NamedTuple):
    pdbs: Iterable[str | bytes | PathLike]
    dx: str | bytes | PathLike
    basename: str
    fpocket_pdb: str | bytes | PathLike
    fpocket_info: str | bytes | PathLike


def parse_src_dir(src_dir: str | bytes | PathLike
                  ) -> Iterator[SystemInfo]:
    """入力ディレクトリから使用するファイルのパスを抜き出す"""
    dx_pattern = re.compile(r'.*PMAP.*_nVH.dx')
    for dx_path in _search_rec_file(pathlib.Path(src_dir), dx_pattern):
        yield _parse_md_dir(dx_path.parent)


def _search_rec_file(src_dir: pathlib.Path, pattern: re.Pattern,
                     ) -> Iterator[pathlib.Path]:
    """指定された文字列で終わるファイル名を再帰的に検索する.
    ファイルが1つでも見つかったディレクトリは
    子ディレクトリ含めそれ以上探索しない.
    """
    for f in src_dir.iterdir():
        if pattern.match(f.name):
            yield f
            return
    for f in src_dir.iterdir():
        if f.is_dir():
            for p in _search_rec_file(f, pattern):
                yield p


def _get_system_name(src_path: pathlib.Path) -> str:
    """ファイルパスからシステム名を取得する"""
    def_name = src_path.name
    if src_path.is_file():
        src_path = src_path.parent
    while src_path.name.lower() == 'output':
        src_path = src_path.parent
        if not src_path.name:
            return def_name
    return src_path.name


def _parse_md_dir(src_dir: pathlib.Path) -> SystemInfo:
    dx_pattern = re.compile(r'.*PMAP.*_nVH.dx')
    system_buf: list[tuple[int, pathlib.Path]] = []
    dx_file: pathlib.Path = None
    cur_pdb: pathlib.Path = None
    fpocket_file: tuple[pathlib.Path, pathlib.Path] = (None, None)
    for p in src_dir.iterdir():
        if p.is_dir() and p.name.startswith('system'):
            system_buf.append((int(p.name[6:]), p))
        if p.is_dir() and p.name.startswith('fpocket'):
            fpocket_file = _parse_fpocket_dir(p)
        if p.is_file():
            if dx_pattern.match(p.name):
                dx_file = p
            elif p.name.endswith('_position_check2.pdb'):
                cur_pdb = p
    system_buf.sort(key=itemgetter(0))
    pdb_buf = tuple(filter((lambda p: p is not None),
                           itertools.chain(
                               (cur_pdb, ),
                               (_find_pdb(p) for _, p in system_buf))))
    return SystemInfo(
        pdbs=pdb_buf,
        dx=dx_file,
        basename=_get_system_name(src_dir),
        fpocket_pdb=fpocket_file[0],
        fpocket_info=fpocket_file[1],
    )


def _find_pdb(src_dir: pathlib.Path) -> pathlib.Path | None:
    for p in src_dir.iterdir():
        if p.name.endswith('_position_check2.pdb'):
            return p
    return None


def _parse_fpocket_dir(src_dir: pathlib.Path
                       ) -> tuple[pathlib.Path | None, pathlib.Path | None]:
    fpocket_pdb = None
    fpocket_info = None
    for p in src_dir.iterdir():
        if p.name.endswith('.pdb'):
            fpocket_pdb = p
        elif p.name.endswith('_info.txt'):
            fpocket_info = p
    return (fpocket_pdb, fpocket_info)
