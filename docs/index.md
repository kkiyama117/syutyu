# syutyu ドキュメント

LAMMPS で MD を実行し、RDF・MSD・VACF などを Python で解析するライブラリ `syutyu` の文書です。

- **使い方**：Python API と CLI、LAMMPS 内蔵計算と自前計算の違い。
- **API リファレンス**：公開クラスと関数の引数・戻り値・例外・単位。docstring から生成しています。
- **開発**：モジュールの構成と、閲覧用 LAMMPS ソースの版。

```{toctree}
:caption: 使い方
:maxdepth: 2

guide/usage
guide/analysis-methods
```

```{toctree}
:caption: リファレンス
:maxdepth: 2

reference/api
```

```{toctree}
:caption: 開発
:maxdepth: 2

development/structure
development/lammps-sources
```
