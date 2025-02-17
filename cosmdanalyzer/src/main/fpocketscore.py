from collections.abc import Collection, Iterable
from .. import solidcalc
from ..solidcalc.typehint import Vector3f


def fpocket_score(hotspot: Collection[Vector3f],
                  hotspot_size: float,
                  fpockets: Iterable[tuple[float, Iterable[Vector3f]]],
                  rate_threshold: float) -> float:
    """hotspotに対して最も重なりの大きいfpocketのスコアを求める.

    Args:
        hotspot: 1つのホットスポットの座標集合
        hotspot_size: ホットスポットの座標1点の幅
        fpockets: fpocketのスコアと座標集合
        rate_threshold: [0,1]hotspotに含まれるfpocketの座標点の割合が
                        rate未満の場合は無視する.
    Returns:
        hotspotに対して最も重なりの大きいfpocketのスコア
        rate_thresholdより重なりの大きいスコアがない場合は0.0
    """
    max_rate = 0
    max_score = 0.0
    for score, pos in fpockets:
        rate = solidcalc.calc_point_set_box_overwraped(
            pos,
            hotspot,
            hotspot_size)
        if (max_rate < rate) and rate:
            max_rate = rate
            max_score = score
    return max_score
