# 構成

![モジュールの構成。公開 API の入口、解析の選択、自前の推定、読み書き、実行、描画](../diagrams/modules.svg)

Python から使うときの入口は `syutyu/__init__.py` だけです。`analysis/dispatch` は RDF・MSD・VACF を、LAMMPS の出力を読む方法（`lammps`）と軌跡から自分で計算する方法（`custom`）のどちらかを明示して選びます。自前の推定は `kernels.py` の短い NumPy の関数を使います。`cli.py` は公開 API を通さず、各モジュールを直接呼びます。`report/` は計算結果（`Result`）を受け取るだけで、ほかのモジュールを import しません。

図の元は `docs/diagrams/modules.dot` です。Graphviz で SVG を作ります。
