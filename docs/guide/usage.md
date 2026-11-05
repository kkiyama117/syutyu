# ライブラリとCLI

[計算方式](<analysis-methods.md>)／[API](<../reference/api.md>)。以下はこのリポジトリのルートから実行する。LAMMPS は `uv sync` では入らないので、MD を実行するときは別に用意する。

```bash
uv sync --locked
uv run --locked syutyu --help
```

## 公開API

```python
from syutyu import (
    Result, Trajectory, read_trajectory, read_thermo, save_result, load_result,
    Site, write_charge_data,
    rdf, msd, vacf, fit_diffusion, summarize_thermo,
    run_lammps, save_figure, export_vmd,
)
```

内部パッケージの`__init__`から同じ関数を再輸出しない。教材用カーネル（`syutyu.kernels`）は、そのモジュールを明示してimportする。

## 汎用の実行

```python
run = run_lammps('inputs/system.in', 'runs/my-system',
                 assets=['system.data', 'forcefield/params.inc'],
                 executable='lmp', ranks=1, timeout=180)
```

補助ファイルは主入力の親フォルダ基準で明示する。新しいrunフォルダへコピーし、入力ハッシュ・コマンド・成否・時間・警告・出力をrun.jsonへ残す。入力は信頼済みコードであることが前提で、LAMMPS自身のshellやファイルアクセスを隔離するsandboxではない。POSIXではtimeoutや割込み時にローカルのプロセス群を停止する。

```bash
uv run --locked syutyu run inputs/system.in --asset system.data \
  --output-root runs/my-system --ranks 1
```

## LAMMPS data の書き出し

```python
write_charge_data(
    'inputs/cscl.data', (4.115, 4.115, 4.115),
    [Site(1, 1.0, (0.0, 0.0, 0.0)), Site(2, -1.0, (0.5, 0.5, 0.5))],
    {1: 132.90543, 2: 35.4527}, title='CsCl unit cell')
```

電荷は素電荷 e、長さは Å、質量は g/mol。複製は LAMMPS の `replicate` で行う。

## 自前解析とnative出力

```python
trajectory = read_trajectory('trajectory.lammpstrj')
result = rdf(trajectory, method='custom', pairs=[('type:1', 'type:2')],
             rmax=10, bins=200, normalization='finite')
save_result(result, 'results/my-rdf')
save_figure(load_result('results/my-rdf'), 'results/my-rdf.png')
```

選択は`all`、`element:Ar`、`type:3`、`mol:12`。dumpにあるラベルだけを使い、元素から原子タイプや力場部位を推測しない。分子内除外は`intermolecular=True`で、非ゼロの同一mol IDを除く。mol=0同士は同一分子と扱わない。

RDFは各フレームの体積で規格化する。既定のfiniteは分母`N_A*N_B-N_overlap`、populationは`N_A*N_B`。後者の同種理想気体の有限N極限は1−1/Nである。CNは中心原子1個当たりの直接カウントであり、可変体積では平均gに平均密度を掛け直さない。rmaxは対象全フレームの最短辺の半分以下。結合近傍を除くLAMMPSのspecial_bonds処理までは再現しない。

```bash
uv run --locked syutyu rdf trajectory.lammpstrj --method custom \
  --pair type:1 type:2 --rmax 10 --bins 200 --normalization finite --output results/rdf
uv run --locked syutyu rdf native.gr --method lammps --output results/native-rdf
uv run --locked syutyu msd trajectory.lammpstrj --method custom \
  --select all --dt-ps .001 --origins multiple --output results/msd
uv run --locked syutyu vacf trajectory.lammpstrj --method custom \
  --dt-ps .001 --velocity-unit A/fs --max-lag-ps 1 --output results/vacf
```

`dt_ps`は1 LAMMPS step当たりのpsであり、dump間隔ではない。原子ID、セル、座標形式、速度の単位を明記する。直交・3方向周期セルのみを扱い、反応でIDやtypeが変わる系、傾斜セル等は受け付けない。距離はÅ、MSDはÅ²を前提とし、units ljの換算は行わない。

MSDはxu/yu/zu・image flagsを優先する。固定セルのwrapped座標だけなら、真の隣接フレーム変位の各成分がL/2未満という仮定でunwrapする。可変セルには明示unwrapped情報とlaboratory=Trueが必要で、その値はセル伸縮を含む実験室座標のMSDである。分子重心MSDは計算しない。系全体の重心ドリフトは、APIで `com_masses` を与えたときだけ引く（CLIにはない）。

## native表の段階的な読込

```python
from syutyu.io.native import read_rdf, read_msd, read_vacf

result = read_msd({'Na': 'Na.msd', 'Cl': 'Cl.msd'}, dt_ps=.001, origin_step=0)
```

MSD原点は必須。先頭の保存stepから推測しない。既存NPT dumpのunwrapが不確かなら、LAMMPSのcomputeで既に求めたMSDを読む。表の最後の列を全MSDとして扱い、必要ならcolumnを明示する。

```bash
uv run --locked syutyu native-msd --series Na Na.msd --series Cl Cl.msd \
  --dt-ps .001 --origin-step 0 --output results/native-msd
uv run --locked syutyu thermo log.lammps --block -1 --dt-ps .001 --output results/thermo
```

thermoはrun区間を分け、reset前後を混ぜない。平均と標準偏差は記述統計であり、相関した記録の平均値の標準誤差ではない。

## 描画と持ち出し

```bash
uv run --locked syutyu plot results/vacf --view vacf-normalized --output results/vacf.png
uv run --locked syutyu export-vmd trajectory.lammpstrj --output exports/vmd.zip
```

描画はRDF/CN/MSD/VACF/規格化VACF/温度。通常の結果・PNGの上書きにはoverwriteを明示し、配布ZIPは既存名へ上書きしない。

保存形式はschema3のResult（CSV tables + JSON summary）。方法、単位、原点、選択、限界を記録する。読込とMSDは軌跡をメモリに載せるため、大規模系のストリーミング処理ではない。原子MSDであり、粘性・熱伝導率・分子重心・傾斜セル等を対応済みとは扱わない。
