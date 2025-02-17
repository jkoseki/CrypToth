"""メインルーチン"""
import os
import sys
import tomli
from .. import arguments
from . import calcmain
from . import spot
from . import input


def load_setting(setting_file_path: str | bytes | os.PathLike) -> dict:
    with open(setting_file_path, 'rb') as f:
        toml_dict = tomli.load(f)
    return toml_dict


def main(root_dir: str | bytes | os.PathLike) -> None:
    """アプリケーションエントリーポイント

    Args:
        root_dir: このプロジェクトのルートディレクトリのパス
    """
    args = arguments.create_parser().parse_args()
    system_infos = input.parse_src_dir(args.src_dir)
    setting_path = args.setting
    if setting_path is None:
        setting_path = os.path.join(root_dir, 'data/setting.toml')
    setting = load_setting(setting_path)
    weight = setting['score']['weight']
    weight_array = (weight['gfe'],
                    weight['size'],
                    weight['protrusion'],
                    weight['convexity'],
                    weight['compactness'],
                    weight['hydrophobicity'],
                    weight['charge_density'],
                    weight['flexibility'],
                    weight['fpocket'],
                    )
    algo = (setting['clustering']['algorithm']).lower()
    if algo == 'single_linkage':
        clustering_input = spot.SingleLinkageInput(
                setting['clustering']['single_linkage']['threshold'],
                )
    elif algo == 'dbscan':
        clustering_input = spot.DbscanInput(
                setting['clustering']['dbscan']['epsilon'],
                setting['clustering']['dbscan']['min_pts'],
                )
    elif algo == 'mean_shift':
        clustering_input = spot.MeanShiftInput(
                setting['clustering']['mean_shift']['bandwidth'],
                )
    else:
        print('Unknown clustering algorithm {}'.format(algo), file=sys.stderr)
    calcmain.calc_main(system_infos,
                       args.out_dir,
                       setting['clustering']['occupancy'],
                       clustering_input,
                       setting['clustering']['extend'],
                       setting['score']['fpocket_threshold'],
                       os.path.join(root_dir, 'data/hydrophobicity.csv'),
                       os.path.join(root_dir, 'data/aminoacids.rtp'),
                       setting['score']['temperature'],
                       setting['score']['solvent_radius'],
                       weight_array,
                       setting['clustering']['spot_marge_rate'],
                       setting['score']['resolution'],
                       args.output_detail,
                       args.verbose,
                       )
