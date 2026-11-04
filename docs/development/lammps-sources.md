# 閲覧用ソース

## LAMMPS

- 場所：`local/lammps-src`（C++本体はその中の`src/`）
- 版：30 Sep 2026
- コミット：`8de817dd79bfe4525d5d39246a212d833e6dee07`
- 入手元：[LAMMPS公式リポジトリの該当コミット](https://github.com/lammps/lammps/tree/8de817dd79bfe4525d5d39246a212d833e6dee07)
- 再取得用：[同じコミットのソースアーカイブ](https://github.com/lammps/lammps/archive/8de817dd79bfe4525d5d39246a212d833e6dee07.tar.gz)
- ライセンス：[同梱LICENSE](<../../local/lammps-src/LICENSE>)を参照。

2026-10-08に、既存コンテナのビルドで取得済みだった同一コミットを`git archive`で書き出した。
Git履歴（`.git/`）を含まない閲覧用のソース一式で、コンパイルやコンパイラの追加はしていない。
実行用コンテナのLAMMPSと同じ版であり、ソースを編集しても既存の実行ファイルは変わらない。

外部ソース本体はこのリポジトリのGit管理外。この案内に版と入手先を残す。

## 読み始める場所

まず次のファイルから読むと、ソースの構成とポテンシャル・解析の実装が追いやすい。

| 内容 | この版のファイル |
| --- | --- |
| LJ＋長距離Coulombの例 | [pair_lj_cut_coul_long.cpp](<../../local/lammps-src/src/KSPACE/pair_lj_cut_coul_long.cpp>) |
| Nosé–Hooverの温度・圧力制御 | [fix_nh.cpp](<../../local/lammps-src/src/fix_nh.cpp>) |
| 平均二乗変位（MSD） | [compute_msd.cpp](<../../local/lammps-src/src/compute_msd.cpp>) |
| 動径分布関数（RDF） | [compute_rdf.cpp](<../../local/lammps-src/src/compute_rdf.cpp>) |
| バージョンの確認 | [version.h](<../../local/lammps-src/src/version.h>) |

対応する`.h`ファイルも同じ場所にある。まずクラス定義と、`.cpp`内の計算部分を見比べると追いやすい。
