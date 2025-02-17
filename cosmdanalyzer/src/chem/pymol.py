"""PyMol互換処理"""
from collections.abc import Iterable, Iterator
import math
import os
from rdkit import Chem
from typing import IO
from ..solidcalc.typehint import Vector3f


"""項目名,スコア"""
PatchScores = Iterable[tuple[str, float]]


def write_pymol_src(
        name: str,
        out_dir: str | bytes | os.PathLike,
        mol: Chem.rdchem.Mol,
        all_patch_info: Iterable[tuple[Iterable[Vector3f], Iterable[int],
                                       PatchScores]],
        ) -> None:
    """PyMol用の入力を出力する.

    Args:
        name: 名称
        mol: 分子オブジェクト
        all_patch_info: パッチの[hotspotの座標集合,原子のインデックス集合,
                                 [スコア名称, スコア]]
        out_dir: 出力ディレクトリ
    """
    spots_dir = os.path.join(out_dir, 'spots')
    os.makedirs(spots_dir, exist_ok=True)
    out_pdb_filename = name + '_out.pdb'
    out_pdb_path = os.path.join(out_dir, out_pdb_filename)
    out_pymol_script_path = os.path.join(out_dir, name + '.pml')
    out_info_path = os.path.join(out_dir, name + '_info.txt')

    with open(out_pymol_script_path, 'w') as out:
        _write_pymol_script(out, out_pdb_filename)
    # write pdb
    pdb_str = Chem.rdmolfiles.MolToPDBBlock(
            mol, mol.GetNumConformers() - 1, flavor=34)
    pdb_str_list = pdb_str.split('\n')
    with open(out_pdb_path, 'w') as pdb_out, \
         open(out_info_path, 'w') as info_out:
        for line in pdb_str_list:
            if line.startswith('ATOM  ') or line.startswith('HETATM'):
                pdb_out.write(line)
                pdb_out.write('\n')
        for i, patch_info in enumerate(all_patch_info, 1):
            part_pdb_path = os.path.join(spots_dir, 'spot{}.pdb'.format(i))
            with open(part_pdb_path, 'w') as part_pdb_out:
                _write_part_pdb(part_pdb_out, pdb_str_list, patch_info[1])
            for point_line in _create_point_pdb_lines(i, patch_info[0]):
                pdb_out.write(point_line)
                pdb_out.write('\n')
            _write_score_info(info_out, patch_info[2], i)


def write_score_info_file(
        out_file: str | bytes | os.PathLike,
        patches: Iterable[PatchScores]) -> None:
    """Score Infoファイルを出力する.

    Args:
        out_file: 出力ファイルのパス
        patches: パッチ毎のスコア名称とスコア
    """
    with open(out_file, 'w') as out:
        for i, patch in enumerate(patches, 1):
            _write_score_info(out, patch, i)


def _write_score_info(
        out: IO[str], patch: PatchScores, patch_number: int) -> None:
    out.write('Patch ')
    out.write(str(patch_number))
    out.write('\n')
    for name, score in patch:
        out.write('\t')
        out.write(name)
        out.write(' : ')
        out.write('\t')
        if (score is not None) and (not math.isnan(score)):
            out.write('{:.3f}'.format(score))
        else:
            out.write('n/a')
        out.write('\n')
    out.write('\n')


def _write_pymol_script(out: IO[str], pdb_filename: str):
    """PyMol用のスクリプトを出力する.

    Args:
        out: 出力先ストリーム
        pdb_filename: 全形PDBのファイル名
    """
    pymol_script = (
        'from pymol import cmd,stored\n'
        'load {}\n'
        '#select spots, resn STP\n'
        'stored.list=[]\n'
        'cmd.iterate("(resn STP)","stored.list.append(resi)")'
        '    #read info about residues STP\n'
        '#print stored.list\n'
        'lastSTP=stored.list[-1]	#get the index of the last residu\n'
        'hide lines, resn STP\n'
        '\n'
        '#show spheres, resn STP\n'
        'for my_index in range(1,int(lastSTP)+1): '
        'cmd.select("spot"+str(my_index), '
        '"resn STP and resi "+str(my_index))\n'
        'for my_index in range(2,int(lastSTP)+2): '
        'cmd.color(my_index,"spot"+str(my_index))\n'
        'for my_index in range(1,int(lastSTP)+1): '
        'cmd.show("spheres","spot"+str(my_index))\n'
        'for my_index in range(1,int(lastSTP)+1): '
        'cmd.set("sphere_scale","0.3","spot"+str(my_index))\n'
        'for my_index in range(1,int(lastSTP)+1): '
        'cmd.set("sphere_transparency","0.1","spot"+str(my_index))\n'
    )
    out.write(pymol_script.format(pdb_filename))


def _write_part_pdb(out: IO[str], pdb: Iterable[str], atom_idxs: Iterable[int]
                    ) -> None:
    """与えられたPDB文字列の一部を原子インデックスで指定して出力する.

    Args:
        out: 出力先ストリーム
        pdb: 1行毎のPDB文字列
        atom_idxs: 原子インデックスの集合
    """
    idx_set = set(atom_idxs)
    for line in _part_pdb(pdb, idx_set):
        out.write(line)
        out.write('\n')


def _part_pdb(pdb: Iterable[str], atom_idxs: set[int]) -> Iterator[str]:
    """与えられたPDB文字列の一部を原子インデックスで切り抜く.

    Args:
        pdb: 1行毎のPDB文字列
        atom_idxs: 原子インデックスの集合
    Returns:
        原子インデックスに対応する原子行のみのPDB
    """
    for line in pdb:
        if ((line.startswith('ATOM  ') or line.startswith('HETATM'))
                and (int(line[6:11]) in atom_idxs)):
            yield line


def _create_point_pdb_lines(residue_idx: int, positions: Iterable[Vector3f]
                            ) -> Iterator[str]:
    """PDB用の点描画用情報を生成する."""
    for i, pos in enumerate(positions):
        yield _create_point_pdb_line(i + 1, residue_idx, pos)


def _create_point_pdb_line(
        point_idx: int, residue_idx: int, pos: Vector3f) -> str:
    """PDB用の点描画用情報を1点分生成する."""
    template = "HETATM{:>5d} APOL STP C{:>4d}    {:8.3f}{:8.3f}{:8.3f}" \
               "  0.00  0.00          Ve  "
    return template.format(point_idx, residue_idx, *pos)
