"""残基の情報"""
import csv
import os


class ResidueHydrophobicity:
    """残基の疎水性テーブル"""

    def __init__(self, src_csv: str | bytes | os.PathLike):
        """残基の疎水性テーブルをCSVファイルから作成する.

        Args:
            src_csv: ヘッダにsymbol3とhydrophobicityを含むCSVファイルのパス
        """
        with open(src_csv) as f:
            reader = csv.reader(f)
            header = next(reader, None)
            if header is not None:
                for i in range(len(header)):
                    header[i] = header[i].lower()
                symbol3_idx = header.index('symbol3')
                hydro_idx = header.index('hydrophobicity')
                self._symbol3_to_hydro = dict()
                for row in reader:
                    symbol3 = row[symbol3_idx].lower()
                    self._symbol3_to_hydro[symbol3] = float(row[hydro_idx])
                    if symbol3 == 'his':
                        for variant in ('hid', 'hie', 'hip'):
                            self._symbol3_to_hydro[variant] = (
                                    float(row[hydro_idx]))

    def __call__(self, symbol3: str) -> float:
        """残基の疎水性を返す.

        Args:
            symbol3: 3文字表記の残基シンボル
        Returns:
            疎水性
        """
        return self._symbol3_to_hydro[symbol3.lower()]
