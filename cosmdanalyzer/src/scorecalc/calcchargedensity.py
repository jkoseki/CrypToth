"""charge densityの計算"""
from collections.abc import Callable, Hashable, Iterable
from .. import solidcalc
from ..solidcalc.typehint import Sphere


def calc_atoms_charge_density(
        exposed_atom_ids: Iterable[int],
        atom_to_as_sphere: Callable[[int], Sphere],
        atom_to_charge: Callable[[int], float],
        collided_atom_getter: Callable[[int], Iterable[Hashable]],
        resolution: int,
) -> float:
    """原子集合のcharge densityを計算する.

    Args:
        exposed_atom_ids: 溶媒露出原子のID集合
        atom_to_as_sphere: 原子の溶媒露出平面に相当する球を返す関数
        atom_to_charge: 原子の電荷を返す関数
        collided_atom_getter: 入力原子と溶媒露出平面が衝突している
                              他の原子IDの集合を返す関数
        resolution: 原子表面を多面体で近似するときの頂点数
    Returns:
        原子集合のcharge density
    """
    atoms = _SumChargeIterator(exposed_atom_ids, atom_to_charge)
    as_area = sum(map(
        solidcalc.create_one_area_in_multi_spheres_func(
            atom_to_as_sphere, collided_atom_getter, resolution),
        atoms
    ))
    charge = atoms.get_sum_charge()
    if as_area > 0:
        return charge / as_area
    else:
        return 0.0


class _SumChargeIterator:
    """イテレータ走査時に電荷の合計を計算する.
    電荷情報が無い原子は無視する.
    """

    def __init__(self, atom_ids: Iterable[int],
                 atom_to_charge: Callable[[int], float]):
        self._src = atom_ids
        self._atom_to_charge = atom_to_charge
        self._sum_charge = 0.0

    def __iter__(self):
        return self

    def __next__(self):
        while True:
            atom_id = next(self._src)
            try:
                charge = self._atom_to_charge(atom_id)
                self._sum_charge += charge
                return atom_id
            except KeyError:
                pass

    def get_sum_charge(self) -> float:
        return self._sum_charge
