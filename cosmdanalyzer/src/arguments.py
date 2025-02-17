import argparse
import pathlib


def create_parser():
    """コマンドラインオプション設定"""
    parser = argparse.ArgumentParser(description='cosmdanalyzer')
    parser.add_argument('out_dir', help='出力ディレクトリ',
                        type=pathlib.Path)
    parser.add_argument('src_dir', help='入力ディレクトリ',
                        type=pathlib.Path)
    # parser.add_argument('src_pdb', help='トラジェクトリPDB',
    #                     type=pathlib.Path)
    # parser.add_argument('src_open_dx', help='存在確率',
    #                     nargs='+', type=pathlib.Path)
    parser.add_argument('-s', '--setting', help='設定ファイルのパス',
                        type=pathlib.Path)
    # parser.add_argument('--fpocket_info',
    #                     help='fpocketの出力*_info.txtファイルのパス',
    #                     type=pathlib.Path)
    # parser.add_argument('--fpocket_pdb', help='fpocketの出力PDBファイルのパス',
    #                     type=pathlib.Path)
    parser.add_argument('--output_detail',
                        help='スコアの傾向を表す詳細情報を出力する',
                        action='store_true')
    parser.add_argument('-v', '--verbose',
                        help='標準出力に詳細な処理情報を表示する',
                        action='store_true')
    return parser
