"""デバッグ用の可視化"""
from collections.abc import Collection, Iterable
from typing import Any
import plotly.graph_objects as go


def plot_3d_points(points: Iterable[tuple[float, float, float]]) -> None:
    """3次元離散プロット"""
    x = []
    y = []
    z = []
    for p in points:
        x.append(p[0])
        y.append(p[1])
        z.append(p[2])
    fig = go.Figure(data=[go.Scatter3d(x=x, y=y, z=z, mode='markers')])
    # fig = go.Figure(data=[go.Scatter3d(x=x, y=y, z=z, mode='lines')])
    lim = _same_aspect_range(x, y, z)
    fig.update_layout(scene={
        'xaxis': {'range': [lim[0], lim[1]]},
        'yaxis': {'range': [lim[2], lim[3]]},
        'zaxis': {'range': [lim[4], lim[5]]},
    })
    fig.show()


def plot_identified_3d_points(
        points: Iterable[tuple[tuple[float, float, float], Any]]) -> None:
    """識別子で色分けした3次元プロット"""
    x = []
    y = []
    z = []
    color = []
    for p, i in points:
        x.append(p[2])
        y.append(p[1])
        z.append(p[0])
        color.append(i)
    fig = go.Figure(data=[go.Scatter3d(x=x, y=y, z=z, mode='markers',
                                       marker=dict(color=color))])
    # fig.update_layout(scene = {'aspectratio': {'x': 1, 'y': 1, 'z': 1}})
    lim = _same_aspect_range(x, y, z)
    fig.update_layout(scene={
        'xaxis': {'range': [lim[0], lim[1]]},
        'yaxis': {'range': [lim[2], lim[3]]},
        'zaxis': {'range': [lim[4], lim[5]]},
    })
    fig.show()


def _same_aspect_range(x: Collection[float], y: Collection[float],
                       z: Collection[float]) \
                              -> tuple[float, float, float, float, float]:
    min_x = min(x)
    max_x = max(x)
    w_x = max_x - min_x
    min_y = min(y)
    max_y = max(y)
    w_y = max_y - min_y
    min_z = min(z)
    max_z = max(z)
    w_z = max_z - min_z
    half_w = 1.05 * max(w_x, w_y, w_z) / 2
    start_x = (min_x + max_x) / 2 - half_w
    end_x = start_x + half_w * 2
    start_y = (min_y + max_y) / 2 - half_w
    end_y = start_y + half_w * 2
    start_z = (min_z + max_z) / 2 - half_w
    end_z = start_z + half_w * 2
    return (start_x, end_x, start_y, end_y, start_z, end_z)
