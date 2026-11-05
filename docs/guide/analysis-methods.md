# LAMMPS内蔵計算と自前計算

[使い方](<usage.md>)。MDの時間発展はいずれもLAMMPSを使用し、解析の実装を選ぶ。

| 物理量 | method=lammps | method=custom |
| --- | --- | --- |
| RDF | compute rdfの規格化済み出力 | 最小像・i<j・ヒストグラム・体積規格化 |
| MSD | compute msd、単一原点 | 座標から単一／複数原点を選択 |
| VACF | compute vacf、単一原点 | 同一原子IDの速度内積を単一／複数原点で平均 |

```python
from syutyu import rdf, msd, vacf

custom_rdf = rdf('trajectory.lammpstrj', method='custom',
                 pairs=[('all','all')], rmax=12, bins=240, normalization='finite')
native_rdf = rdf('Ar_nvt.gr', method='lammps', pairs=[('Ar','Ar')])
custom_msd = msd('trajectory.lammpstrj', method='custom', dt_ps=.001,
                 start_step=1000, origin_mode='multiple')
native_msd = msd('Ar_nvt.msd', method='lammps', dt_ps=.001, origin_step=0)
custom_vacf = vacf('trajectory.lammpstrj', method='custom', dt_ps=.001,
                   velocity_unit='A/fs', origin_mode='multiple', max_lag_ps=1)
native_vacf = vacf('Ar_nvt.vacf', method='lammps', dt_ps=.001,
                   velocity_unit='A/fs', origin_step=0)
```

方式は明示指定で、ファイルがなければ別方式へ勝手に切り替えない。nativeには後から原子選択や平均方法を加え直せない。computeを作ったstepをorigin_stepとして指定し、dumpや表の先頭stepから推測しない。native VACFの規格化には原点での値が必要である。

## RDFの規格化

finiteでは単成分の一様・無相関な有限系で期待値が1となる：`g=H/[(N(N-1)/2)ΔV/V]`。populationでは密度N/Vを使い、同種の有限N極限が1−1/Nになる。新しい公開rdfはfiniteが既定。教材カーネルの [pair_rdf](<../../src/syutyu/kernels.py>) と内部の自前RDFには規格化を必ず明示する。呼ぶ入口によって知らないうちに規約が変わる状態を作らない。

nativeのブロック平均は元の入力に従い、Pythonでは保存されたフレームを使う。出力間隔や区間を揃えなければ、同じ物理量でも数値の一致は必要ない。CNは直接計数した平均近傍数で、平均密度から再構成しない。最終ビンを含めr<rmaxである。

## 時間原点・単位

customのmultipleは各lagで可能な全原点と原子を平均する。singleは最初の選択フレームを固定し、nativeとの照合用に使える。異原子相関は加えない。長いlagほど原子×原点数が減り、重なった原点は独立標本ではない。

1 step当たりのdt_psと、保存間隔×dt_psを混同しない。A/fs入力はA/psへ1000倍してから速度積を求め、VACFはÅ²/ps²、規格化後は無次元。Cv(0)=0なら規格化が定義できないので拒否する。VACFではドリフトを勝手に差し引かず、等質量系での平均速度RMSを診断記録する。これは一般の異質量系の質量重み付きCOMとは別のものである。

`com_masses` を与えると、全原子の質量で重み付けした系全体の重心の変位を各原子の変位から引いてから MSD を計算する。種ごとの重心を別々に引くことはしない。

## 検証

native/customは、同じ軌跡の同じ標本・単一原点で照合した。既存の検証でg/MSD/VACFは丸め誤差程度で一致し、複数原点との違いも確認した。数値一致は平衡化や物性値の収束の証明ではない。
