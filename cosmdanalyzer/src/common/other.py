"""その他"""
import os


def path_to_str(path: str | bytes | os.PathLike) -> str:
    """ファイルパスを表す型をstr型に統一する.

    Args:
        path: ファイルパス
    Returns:
        ファイルパスを表す文字列
    """
    if isinstance(path, os.PathLike):
        return os.fsdecode(path)
    elif isinstance(path, bytes):
        return path.decode('utf-8')
    return path
