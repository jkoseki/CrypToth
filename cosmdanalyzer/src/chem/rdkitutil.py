"""RDKit型のユーティリティー関数"""
from collections.abc import Iterable, Iterator
from os import PathLike
from rdkit import Chem
from .. import common
from ..solidcalc.typehint import Vector3f


class Mol:
    """RDKitの分子型から使用する情報のみのインターフェースに変換する"""

    def __init__(self, mol: Chem.rdchem.Mol):
        """

        Args:
            mol: RDKitの分子オブジェクト
        """
        self._mol = mol

    def get_rdkit_mol(self) -> Chem.rdchem.Mol:
        """RDKitのMolオブジェクトを返す.

        Returns:
            RDKitのMolオブジェクト
        """
        return self._mol

    def get_atom_idxs(self) -> Iterator[int]:
        """原子IDの集合を返す

        Returns:
            原子IDの集合
        """
        for atom in self._mol.GetAtoms():
            yield atom.GetIdx()

    def get_num_conformers(self) -> int:
        """コンフォーマー数を返す

        Returns:
            コンフォーマー数
        """
        return self._mol.GetNumConformers()

    def atom_to_atomic_number(self, atom_id: int) -> int:
        """原子IDに対応する原子番号を返す

        Returns:
            原子IDに対応する原子番号
        """
        return self._mol.GetAtomWithIdx(atom_id).GetAtomicNum()

    def atom_to_position(self, atom_idx: int, conformer_id: int) -> Vector3f:
        """原子IDに対応する原子座標を返す.

        Atoms:
            atom_idx: 原子ID
            conformer_id: コンフォーマーID
        Returns:
            中心座標
        """
        pos = self._mol.GetConformer(conformer_id).GetAtomPosition(atom_idx)
        return (pos.x, pos.y, pos.z)

    def atom_to_name(self, atom_idx: int) -> str:
        """原子IDに対応する原子名を返す.

        Atoms:
            atom_idx: 原子ID
        Returns:
            原子名
        """
        return self._mol.GetAtomWithIdx(atom_idx).GetMonomerInfo().GetName()

    def atom_to_vdw_radius(self, atom_idx: int) -> float:
        """原子IDに対応するファンデルワールス半径を返す.

        Atoms:
            atom_idx: 原子のID
        Returns:
            半径
        """
        return Chem.GetPeriodicTable().GetRvdw(
                self._mol.GetAtomWithIdx(atom_idx).GetAtomicNum())

    def atom_to_weight(self, atom_idx: int) -> float:
        """原子IDに対応する原子量を返す.

        Atoms:
            atom_idx: 原子ID
        Returns:
            原子量
        """
        return Chem.GetPeriodicTable().GetAtomicWeight(
                self._mol.GetAtomWithIdx(atom_idx).GetAtomicNum())

    def divide_to_residue(self, atom_idxs: Iterable[int]
                          ) -> dict[int, list[int]]:
        """原子ID集合を残基毎に分別する.

        Args:
            atom_idxs: 原子IDのイテレータ
        Returns:
            {残基番号: 残基に含まれる原子のインデックス集合}
        """
        res_atom_idxs = dict()
        for idx in atom_idxs:
            res_n = (self._mol.GetAtomWithIdx(idx).
                     GetPDBResidueInfo().GetResidueNumber())
            if res_n in res_atom_idxs:
                res_atom_idxs[res_n].append(idx)
            else:
                res_atom_idxs[res_n] = [idx,]
        return res_atom_idxs

    def atom_to_residue(self, atom_idx: int) -> int:
        """原子IDに対応する残基番号を返す.

        Rerutns:
            原子IDに対応する残基番号
        """
        return (self._mol.GetAtomWithIdx(atom_idx)
                .GetPDBResidueInfo().GetResidueNumber())

    def atom_to_residue_symbol(self, atom_idx: int):
        """原子IDに対応する3文字の残基名を返す.

        Args:
            atom_idx: 原子ID
        Rerutns:
            3文字の残基名
        """
        return (self._mol.GetAtomWithIdx(atom_idx)
                .GetPDBResidueInfo().GetResidueName())

    def get_neighbor_atoms(self, atom_idx: int) -> Iterator[int]:
        """指定原子と共有結合を持つ原子ID集合を返す.

        Args:
            atom_idx: 原子ID
        Returns:
            指定原子と共有結合を持つ原子ID集合
        """
        for atom in self._mol.GetAtomWithIdx(atom_idx).GetNeighbors():
            yield atom.GetIdx()


def create_mol_from_pdb_file(pdb_path: str | bytes | PathLike) -> Mol:
    """PDBファイルからMolオブジェクトを作成する.

    Args:
        pdb_path: pdbファイルのパス
    Return:
        Molオブジェクト
    """
    return Mol(Chem.MolFromPDBFile(common.path_to_str(pdb_path)))


def create_mol_from_pdb_str(pdb_str: str) -> Mol:
    return Mol(Chem.MolFromPDBBlock(pdb_str))
