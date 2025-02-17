import math
import operator
from rdkit import Chem
from gridData import Grid
from .. import chem
from .. import solidcalc
from .. import visualization
from .. import fpocket
# from . import testscikitlearn
from . import testclustering


def test():
    grid = solidcalc.CollisionGrid(1.0, (1.0, 1.0, 1.0, 10.0, 10.0, 10.0))
    obj1 = 'obj1'
    obj1_box = (2.0, 2.0, 2.0, 1.0, 1.0, 1.0)
    grid.add_object(obj1, obj1_box)
    obj2 = 'obj2'
    obj2_box = (7.0, 8.0, 9.0, 3.0, 3.0, 3.0)
    grid.add_object(obj2, obj2_box)
    obj3 = 'obj3'
    obj3_box = (0.0, 1.0, 0.0, 2.0, 2.0, 2.0)
    grid.add_object(obj3, obj3_box)
    for obj in grid.iterate_collided_objects(obj3):
        print(obj)


def plot_sphere():
    p = solidcalc.iterate_normalized_sphere_plot(256)
    visualization.plot_3d_points(p)


# def test_scikit_learn():
#     testscikitlearn.test_scikit_learn()


def test_single_linkage():
    testclustering.test_single_linkage()


def test_load_data():
    src_file = '../data/PMAP_TEST_PROJECT_nV.dx'
    g = Grid(src_file)
    visualization.plot_3d_points(itr_grid_0(g.grid))
    visualization.plot_3d_points(itr_grid_plus(g.grid))


def itr_grid_0(grid):
    for x in range(grid.shape[0]):
        for y in range(grid.shape[1]):
            for z in range(grid.shape[2]):
                if grid[x, y, z] == 0.0:
                    yield (x, y, z)


def itr_grid_plus(grid):
    for x in range(grid.shape[0]):
        for y in range(grid.shape[1]):
            for z in range(grid.shape[2]):
                if grid[x, y, z] > 0.0:
                    yield (x, y, z)


def test_load_fpoint():
    with open('../fspot_out/TEST_PROJECT_position_check_info.txt') as f:
        val = fpocket.parse_info(f)
    print(val)
