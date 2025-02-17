"""出力"""
from collections.abc import Iterable, Sequence
import itertools
import math
from typing import IO


def output_detail_score_csv(out: IO[str],
                            data: Iterable[Iterable[Iterable[float]]],
                            score_names: Sequence[str],
                            ) -> None:
    """スコア詳細情報をCSV形式でストリームに出力する

    Args:
        out: 出力先
        data: 0次元 スポット, 1次元 スコアの種類, 2次元 各分布の値
              分布の値は (spot,score_type,mean,variance,min,
                          0.25,0.50,0.75,max) の順
    """
    detail_types = ('spot', 'score_type', 'mean', 'variance', '0.0',
                    '0.25', '0.50', '0.75', '1.0')
    detail_itr = iter(detail_types)
    out.write(next(detail_itr))
    for d in detail_itr:
        out.write(',')
        out.write(d)
    out.write('\n')
    for i, spot_data in enumerate(data, 1):
        for name, score_data in zip(score_names, spot_data):
            out.write(str(i))
            out.write(',')
            out.write(name)
            if score_data is None:
                score_data = iter(tuple())
            score_data = itertools.islice(
                    itertools.chain(score_data, itertools.repeat(None)),
                    len(detail_types) - 2)
            for el in score_data:
                if (el is not None) and (not math.isnan(el)):
                    out.write(',{:.7g}'.format(el))
                else:
                    out.write(',n/a')
            out.write('\n')
