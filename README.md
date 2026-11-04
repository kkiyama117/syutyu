# syutyu

LAMMPS で分子動力学（MD）を実行し、その結果を Python で解析するライブラリです。動径分布関数（RDF）、配位数（CN）、平均二乗変位（MSD）、速度自己相関関数（VACF）、拡散係数を、LAMMPS 内蔵の計算結果から読む方法（`lammps`）と、軌跡から自分で計算する方法（`custom`）のどちらかを明示して求めます。

- [src/syutyu/](<src/syutyu/>)：数値実装・入出力・CLI。公開 API は `syutyu` パッケージから import します。
- [tests/](<tests/>)：数値・入出力・CLI のテスト。
- [docs/](<docs/index.md>)：[使い方](<docs/guide/usage.md>)・[計算方式](<docs/guide/analysis-methods.md>)・[API](<docs/reference/api.md>)。
- [tools/](<tools/>)：環境の検査、API リファレンスの生成、公開範囲の検査。
- `local/`：検証記録・閲覧用 LAMMPS ソース・キャッシュ（Git 管理外）。

```bash
uv sync --locked
uv run --locked syutyu --help
make check
```

LAMMPS 本体は `uv sync` では入りません。MD を実行するときは別に用意してください。
