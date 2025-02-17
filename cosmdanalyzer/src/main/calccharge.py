"""タンパク質原子の電荷対応表を作成する."""
from collections.abc import Callable, Iterable, Iterator, Sequence
from os import PathLike
from .. import chem
from .. import common
from .. import graph


def calc_atoms_charge_from_lenard_jones_file(
        atom_ids: Iterable[int],
        atom_to_name: Callable[[int], str],
        charge_file_path: str | bytes | PathLike,
        ) -> dict[int, float]:
    """電荷を原子名マッチングで付加する.

    Args:
        atom_ids: 原子IDの集合
        atom_to_name: 原子IDから原子名を返す関数
        charge_file_path: csv形式のファイル
    Returns:
        原子IDから電荷の対応表
    """
    id_to_charge: dict[int, float] = dict()
    charge_table = chem.LennardJonesCharge(charge_file_path)
    for atom_id in atom_ids:
        name = atom_to_name(atom_id)
        charge = charge_table(name)
        if charge is not None:
            id_to_charge[atom_id] = charge
    return id_to_charge


def calc_atoms_charge_from_rtp_file(
        res_ids: Iterable[int],
        res_to_atom_ids: Callable[[int], Iterable[int]],
        atom_to_atomic_number: Callable[[int], int],
        atom_to_neighbors: Callable[[int], Iterable[int]],
        atom_to_residue_symbol: Callable[[int], str],
        charge_file_path: str | bytes | PathLike,
        ) -> dict[int, float]:
    """

    Args:
        res_ids: 残基IDの集合
        res_to_atom_ids: 残基IDから構成原子IDを返す関数
        atom_to_atomic_number: 原子IDから原子番号を返す
        atom_to_neighbors: 原子IDから隣接原子ID集合を返す関数
        atom_to_residue_symbol: 原子IDから残基の3文字シンボルを返す関数
        charge_file_path: Glomacs電荷ファイルのパス
    Returns:
        原子IDから電荷の対応表
    """
    res_charge = dict()
    for symbol, ch, link in chem.load_glomacs_charge_from_file(
            charge_file_path, True):
        res_charge[symbol] = (ch, common.dict_to_function(link))
    return calc_atoms_charge(
            res_ids, res_to_atom_ids, atom_to_atomic_number,
            atom_to_neighbors, atom_to_residue_symbol,
            (lambda s: res_charge.get(s, (tuple(), lambda _: iter(tuple())))),
            )


def calc_atoms_charge(
        res_ids: Iterable[int],
        res_to_atom_ids: Callable[[int], Iterable[int]],
        atom_to_atomic_number: Callable[[int], int],
        atom_to_neighbors: Callable[[int], Iterable[int]],
        atom_to_residue_symbol: Callable[[int], str],
        res_atoms_charge: Callable[[str],
                                   tuple[Sequence[tuple[int, float]],
                                         Callable[[int], Iterable[int]]]]
        ) -> dict[int, float]:
    """

    Args:
        res_ids: 残基IDの集合
        res_to_atom_ids: 残基IDから構成原子IDを返す関数
        atom_to_atomic_number: 原子IDから原子番号を返す
        atom_to_neighbors: 原子IDから隣接原子ID集合を返す関数
        atom_to_residue_symbol: 原子IDから残基の3文字シンボルを返す関数
        res_atoms_charge: 3文字の残基名から(原子番号,電荷)の原子集合と
                          原子インデックスから隣接原子のインデックスを
                          返す関数を返す
    """
    atom_to_charge: dict[int, float] = dict()
    for res_id in res_ids:
        res_atoms = tuple(res_to_atom_ids(res_id))
        symbol = atom_to_residue_symbol(res_atoms[0])
        ret = _calc_residue_atoms_charge(
                res_atoms,
                atom_to_atomic_number,
                atom_to_neighbors,
                *res_atoms_charge(symbol))
        if ret is not None:
            for atom_id, charge in ret:
                atom_to_charge[atom_id] = charge
    return atom_to_charge


def _calc_residue_atoms_charge(
        res_atom_ids: Sequence[int],
        atom_to_atomic_number: Callable[[int], int],
        atom_to_neighbors: Callable[[int], Iterable[int]],
        template_atom_charge: Sequence[tuple[int, float]],
        template_atom_to_neighbors: Callable[[int], Iterable[int]]
        ) -> Iterator[tuple[int, float]] | None:
    """

    Args:
        res_atom_ids: 1つの残基の原子ID集合
        atom_to_atomic_number: 原子IDから原子番号を返す
        atom_to_neighbors: 原子IDから隣接原子ID集合を返す関数
        template_atom_charge: 基準残基の(原子番号,電荷)で表される原子集合
        template_atom_to_bonds: 基準残基の原子IDから
                                隣接原子ID集合を返す関数
    Returns:
        基準残基に一致する場合は(原子ID, 電荷)の集合
        一致しない場合はNone
    """
    res_atom_first = res_atom_ids[0]
    res_atom_set = set(iter(res_atom_ids))
    res_atom_to_neighbors = (lambda a: filter(lambda n: n in res_atom_set,
                                              atom_to_neighbors(a)))
    atom_to_cmp = (lambda a: (atom_to_atomic_number(a),
                              common.len_iterator(res_atom_to_neighbors(a))))
    template_to_cmp = (lambda a: (template_atom_charge[a][0],
                                  common.len_iterator(
                                      template_atom_to_neighbors(a))))
    m = graph.match_all(res_atom_first, res_atom_to_neighbors,
                        range(len(template_atom_charge)),
                        template_atom_to_neighbors,
                        (lambda r, t: atom_to_cmp(r) == template_to_cmp(t)))
    if m is not None:
        for atom_id in res_atom_ids:
            yield (atom_id, template_atom_charge[m[atom_id]][1])
    else:
        return None
