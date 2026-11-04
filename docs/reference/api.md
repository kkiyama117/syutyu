# APIリファレンス

このファイルは docstring から生成する。直接編集しない。

```bash
uv run --locked tools/docs/api_reference.py --check
```

- 対象は `syutyu` の公開クラスと関数。`cli.py` の内部、`syutyu.io.csvtable`、`_` で始まる名前は対象外。
- `from syutyu import ...` で使える18名には「公開API」と記す。それ以外は submodule を明示して import する。
- 単位・計算方法・実行例は [利用ガイド](../guide/usage.md)、[解析方法](../guide/analysis-methods.md)。

## 一覧

| 名前 | 種別 | 概要 |
| --- | --- | --- |
| [`Result`](#result) | クラス | 解析結果の表と要約をまとめた値。 |
| [`save_result`](#save_result) | 関数 | Result を schema3 の CSV と summary.json として保存する。 |
| [`load_result`](#load_result) | 関数 | schema3 の保存ディレクトリを読み、Result を返す。 |
| [`Trajectory`](#trajectory) | クラス | dump 1つ分のフレーム列と、原子ごとの固定ラベル。 |
| [`read_trajectory`](#read_trajectory) | 関数 | LAMMPS の custom dump を読み、Trajectory を返す。 |
| [`frame_indices`](#frame_indices) | 関数 | step の閉区間に入るフレーム位置を、間引きながら返す。 |
| [`read_thermo`](#read_thermo) | 関数 | 1つの thermo ブロックを選び、`time_ps` を加えて返す。 |
| [`read_thermo_blocks`](#read_thermo_blocks) | 関数 | `Step` で始まり `Loop time` で終わる表を、ブロックごとに分けて読む。 |
| [`read_rdf`](#read_rdf) | 関数 | .gr 形式の RDF を読み、選んだブロックの平均を等重みで返す。 |
| [`read_msd`](#read_msd) | 関数 | compute msd の fix ave/time 表を、系列名ごとに読む。 |
| [`read_vacf`](#read_vacf) | 関数 | compute vacf の表を読み、Å²/ps² へ換算して返す。 |
| [`Site`](#site) | クラス | 直方晶の単位胞にある1原子。 |
| [`write_charge_data`](#write_charge_data) | 関数 | 分率座標の単位胞から LAMMPS の atom_style charge 用 data ファイルを書く。 |
| [`rdf`](#rdf) | 関数 | 軌跡または LAMMPS の RDF 表から g(r) と積算配位数を返す。 |
| [`msd`](#msd) | 関数 | 自前計算の原子 MSD か、LAMMPS compute msd の出力を選んで返す。 |
| [`vacf`](#vacf) | 関数 | 自前計算の VACF か、LAMMPS compute vacf の出力を選んで返す。 |
| [`DiffusionFit`](#diffusionfit) | クラス | 1つの時間窓で MSD を直線近似した結果。 |
| [`fit_diffusion`](#fit_diffusion) | 関数 | 指定した時間窓で MSD を直線近似し、拡散係数を返す。 |
| [`summarize_thermo`](#summarize_thermo) | 関数 | thermo の区間について、平均・標準偏差・点数を返す。 |
| [`radial_distribution`](#radial_distribution) | 関数 | 中心と隣接の組ごとに、フレーム平均の RDF と積算配位数を返す。 |
| [`mean_squared_displacement`](#mean_squared_displacement) | 関数 | 選択ごとに原子 MSD を計算する。拡散フィットや重心除去は行わない。 |
| [`velocity_autocorrelation`](#velocity_autocorrelation) | 関数 | 選択ごとに速度自己相関を計算する。複数時間原点が Colab の定義。 |
| [`run_lammps`](#run_lammps) | 関数 | 入力と明示した補助ファイルを新しいディレクトリへ写し、実行してそのパスを返す。 |
| [`save_figure`](#save_figure) | 関数 | Result を PNG に描いて保存する。 |
| [`available_views`](#available_views) | 関数 | Result から描ける表示の一覧を返す。 |
| [`draw_view`](#draw_view) | 関数 | 1つの軸に表示を描き、図の下へ出す注記を返す。 |
| [`decorate_figure`](#decorate_figure) | 関数 | 図の題名・条件・注意書きを summary から書く。 |
| [`make_figure`](#make_figure) | 関数 | 描画済みの Figure を作って返す。保存はしない。 |
| [`export_vmd`](#export_vmd) | 関数 | 指定フレームと VMD スクリプト・manifest・説明を1つの ZIP にまとめる。 |
| [`iter_dump_blocks`](#iter_dump_blocks) | 関数 | 1フレーム分の生の dump を、ファイル全体を読み込まずに順に返す。 |
| [`minimum_image`](#minimum_image) | 関数 | 相対ベクトルを各方向 [-L/2, L/2] に戻す。 |
| [`pair_rdf`](#pair_rdf) | 関数 | 全フレームの i<j を数え、g(r) と積算配位数を返す。 |
| [`unwrap`](#unwrap) | 関数 | 隣接フレームの最小像変位を積算して非折り返し座標を作る。 |
| [`msd_multi_origin`](#msd_multi_origin) | 関数 | lag ごとに可能な全開始時刻と原子を平均した MSD を返す。 |

## syutyu.io.store

全解析が共有する現行の結果保存形式。過去形式の読み替えは行わない。

### Result

```python
class Result
```

解析結果の表と要約をまとめた値。

公開API: `from syutyu import Result`

**属性**

- `tables` — 表の名前（`rdf`、`msd`、`vacf`、`thermo`）から列名→配列への辞書。
- `summary` — 方法・単位・選択・限界を記録する辞書。保存時に `schema_version` と `tables` を加える。

### save_result

```python
save_result(result, output, *, overwrite=False)
```

Result を schema3 の CSV と summary.json として保存する。

公開API: `from syutyu import save_result`

**引数**

- `result` — 保存する `Result`。表は空でなく、`rdf`・`msd`・`vacf`・`thermo` のいずれかを含むこと。
- `output` — 保存先ディレクトリ。なければ作る。
- `overwrite` — 既存の非空ディレクトリへ書くことを明示する。省略時は拒否する。

**戻り値**

保存先の `Path`。

**例外**

- `ValueError` — 表が空、または未知の表名を含む。summary に JSON へ書けない値（NaN など）がある。
- `FileExistsError` — `output` が非空で `overwrite` を指定していない。

**注意**

書き込み前に既存ファイルを消さない。`overwrite=True` でも同じ名前の CSV は上書きされる。
`summary.json` は `schema_version: 3` と表の一覧を持ち、古い形式は読めない。

### load_result

```python
load_result(directory)
```

schema3 の保存ディレクトリを読み、Result を返す。

公開API: `from syutyu import load_result`

**引数**

- `directory` — `save_result` が作ったディレクトリ。

**戻り値**

`Result`。`tables` の各列は numpy 配列。

**例外**

- `ValueError` — `schema_version` が 3 でない、表の一覧が空、または未知の表名を含む。
- `FileNotFoundError` — `summary.json` がない。

**注意**

過去形式を読み替える処理はない。CSV は数値前提で、文字列の列は読めない。

## syutyu.io.dump

LAMMPS の custom dump を読み、種別・元素・分子で原子を明示的に選ぶ。

### Trajectory

```python
class Trajectory
```

dump 1つ分のフレーム列と、原子ごとの固定ラベル。

公開API: `from syutyu import Trajectory`

**属性**

- `steps` — 各フレームの LAMMPS step（整数）。
- `ids` — 原子 ID。全フレームで同じ昇順にそろえてある。
- `positions` — 各フレームの座標。Å、周期セル内の絶対座標。
- `boxes` — 各フレームの直方体セルの辺の長さ。Å。
- `origins` — 各フレームのセル原点。Å。
- `labels` — dump にあった `type`・`element`・`mol` の原子ごとの値（文字列）。
- `unwrapped` — 明示的な非折り返し座標。dump に xu/yu/zu か image flag があるときだけ入り、なければ None。
- `unwrap_source` — `unwrapped` の作り方（`xu/yu/zu` か `image flags`）。なければ None。
- `source_path` — 読み込んだファイルの絶対パス。
- `velocities` — 各フレームの速度。ID 順、保存時の単位のまま。vx/vy/vz がそろわなければ None。

#### select

```python
select(expression)
```

選択式に合う原子の真偽マスクを返す。

**引数**

- `expression` — `all`、`type:3`、`element:Na`、`mol:12` のいずれか。ラベルは dump にあった列だけを使う。

**戻り値**

原子数と同じ長さの真偽配列。

**例外**

- `ValueError` — 形式が `キー:値` でない、dump にないキーを指す、または1個も一致しない。

**注意**

元素から原子タイプや力場部位を推測しない。`mol:0` は未設定を表すことがあり、選んでも分子として意味があるとは限らない。

#### fixed_box

```python
fixed_box
```

全フレームでセルの辺と原点が変わっていないかを返す。

**戻り値**

変化が許容差 `1e-10` 未満なら True。可変体積の判定に使う。

### read_trajectory

```python
read_trajectory(path)
```

LAMMPS の custom dump を読み、Trajectory を返す。

公開API: `from syutyu import read_trajectory`

**引数**

- `path` — custom dump のファイル。`ITEM: TIMESTEP` で始まるフレームを並べたもの。

**戻り値**

`Trajectory`。座標は Å の絶対座標、セルは辺の長さ。速度は vx/vy/vz があれば ID 順で入る。

**例外**

- `ValueError` — 直方体・3方向周期（`pp pp pp`）以外のセル、ID やラベルの変化、重複 ID、step の減少、非有限値、vx/vy/vz の一部だけの出力。

**注意**

xu/yu/zu か image flag があれば非折り返し座標も作り、`unwrap_source` に残す。`q` などの追加列は読み捨てる。
単位は入力のまま扱い、Å 以外への換算はしない。`reset_timestep` をまたぐ軌跡は分割して読む。

### frame_indices

```python
frame_indices(trajectory, *, start_step=None, stop_step=None, stride=1)
```

step の閉区間に入るフレーム位置を、間引きながら返す。

`from syutyu.io.dump import frame_indices`

**引数**

- `trajectory` — `Trajectory`。`start_step`・`stop_step` の省略時は最初と最後の step を使う。
- `stride` — 選んだフレームを何個おきに残すか（1以上）。

**戻り値**

フレーム位置の整数配列。

**例外**

- `ValueError` — `stride` が正の整数でない、`start_step > stop_step`、該当フレームがない。

**注意**

間引きは区間で絞った後に行う。step が等間隔でない軌跡は、MSD や VACF 側の等間隔検査で止まる。

## syutyu.io.thermo

LAMMPS の thermo 出力を、ブロックごとに分けて読む。単位系は仮定しない。

### read_thermo

```python
read_thermo(path, *, block=-1, dt_ps)
```

1つの thermo ブロックを選び、`time_ps` を加えて返す。

公開API: `from syutyu import read_thermo`

**引数**

- `path` — LAMMPS の log または screen 出力。
- `block` — 選ぶブロックの番号。負なら末尾から数える。省略時は最後のブロック。
- `dt_ps` — LAMMPS の 1 step 当たりの ps。`timestep` コマンドの値そのものではない。

**戻り値**

列名→numpy 配列の辞書。`time_ps = Step * dt_ps` を加える。

**例外**

- `ValueError` — `dt_ps` が正でない、ブロックがない、選択ブロックが空・非有限・step が増加しない、`time_ps` 列と衝突する。
- `IndexError` — ブロック番号が範囲外。

**注意**

ブロックは結合しない。reset の前後を1つの時系列として扱わない。

### read_thermo_blocks

```python
read_thermo_blocks(path)
```

`Step` で始まり `Loop time` で終わる表を、ブロックごとに分けて読む。

公開API: `from syutyu import read_thermo_blocks`

**引数**

- `path` — LAMMPS の log または screen 出力。

**戻り値**

ブロックごとの辞書（列名→numpy 配列）のリスト。列名と step はブロック間で違ってよい。

**例外**

- `ValueError` — 列名の重複、列数や数値の不一致、`Loop time` のない未完了ブロック。

**注意**

WARNING などの表でない行は読み飛ばす。値は保存された単位のままで、換算しない。

## syutyu.io.native

LAMMPS の compute 出力（rdf・msd・vacf）を、再計算せずに読む。

### read_rdf

```python
read_rdf(path, *, pairs=None, start_step=None, stop_step=None)
```

.gr 形式の RDF を読み、選んだブロックの平均を等重みで返す。

`from syutyu.io.native import read_rdf`

**引数**

- `path` — `.gr` ファイル。`step 行数` の見出しの後に `index r g_0 CN_0 [g_1 CN_1 ...]` が続く。
- `pairs` — 表示用の (中心, 隣接) ラベル。省略時は `center_<i>` と `neighbor_<i>`。
- `start_step`、`stop_step` — 平均に使うブロックの step 範囲（両端を含む）。

**戻り値**

`Result`。`tables['rdf']` は `r_center_A`、`r_upper_A`、`g_<i>`、`CN_<i>`。

**例外**

- `ValueError` — 見出し・行数・列数・step の単調増加・ビン中心の等間隔のいずれかが形式に合わない。

**注意**

距離は Å とみなす。規格化と原子対の除外は LAMMPS の入力側の設定をそのまま残し、g を規格化し直したり、平均密度から CN を作り直したりしない。
`pairs` は表示名であり、元素やタイプの対応を推測するものではない。

### read_msd

```python
read_msd(files, *, dt_ps, origin_step, column=-1)
```

compute msd の fix ave/time 表を、系列名ごとに読む。

`from syutyu.io.native import read_msd`

**引数**

- `files` — 系列名→ファイルの辞書。全系列で step の刻みが同じであること。
- `dt_ps` — LAMMPS の 1 step 当たりの ps。
- `origin_step` — compute の原点にした step。省略できない。
- `column` — MSD の列。step の列は選べない。

**戻り値**

`Result`。`tables['msd']` は `lag_ps` と `MSD_<i>_A2`。

**例外**

- `ValueError` — 列数や行数が足りない、非有限値、step が原点より前か増加しない、系列間で刻みが違う、step 列を選んだ。

**注意**

選んだ列は Å² の全 MSD（成分や log ではない）を前提とする。座標の扱いと平均は LAMMPS の入力に従い、拡散係数や平衡性は推定しない。

### read_vacf

```python
read_vacf(path, *, dt_ps, origin_step, velocity_unit, max_lag_ps=None)
```

compute vacf の表を読み、Å²/ps² へ換算して返す。

`from syutyu.io.native import read_vacf`

**引数**

- `path` — step と vx*vx0、vy*vy0、vz*vz0、合計の5列を持つ表。
- `dt_ps` — LAMMPS の 1 step 当たりの ps。
- `origin_step` — compute の原点にした step。先頭の保存 step と一致すること。
- `velocity_unit` — 入力速度の単位。`A/fs` か `A/ps`。
- `max_lag_ps` — 残す lag の上限（両端を含む）。省略時は全 lag。

**戻り値**

`Result`。`tables['vacf']` は `lag_ps`、成分 `VACF_0_x_A2_ps2` など、`VACF_0_A2_ps2`、`VACF_0_normalized`。

**例外**

- `ValueError` — 単位が2つ以外、5列でない、step が原点から始まらないか増加しない、成分の和が合計と合わない、lag 上限が範囲外、Cv(0) が正でない。

**注意**

単一原点の値であり、複数原点の平均ではない。速度の平均や重心ドリフトの除去はしない。

## syutyu.io.lammps_data

LAMMPS の atom_style charge 用 data ファイルを、分率座標から書く。

### Site

```python
class Site
```

直方晶の単位胞にある1原子。

公開API: `from syutyu import Site`

**属性**

- `type` — LAMMPS の原子タイプ。1 以上。
- `charge` — 電荷。e。
- `fractional` — 分率座標。各成分は [0, 1)。セル長を掛けると Å の座標になる。

### write_charge_data

```python
write_charge_data(path, lengths, sites, masses, *, title='syutyu charge data', overwrite=False)
```

分率座標の単位胞から LAMMPS の atom_style charge 用 data ファイルを書く。

公開API: `from syutyu import write_charge_data`

**引数**

- `path` — 書き出すファイル。親ディレクトリは存在すること。
- `lengths` — 直方晶セルの (a, b, c)。Å。三つの有限な正の値。
- `sites` — 空でない `Site` の列。原子 ID はこの順に 1 から振る。
- `masses` — 原子タイプから質量（g/mol）への辞書。`sites` が使うタイプをすべて含み、使わないタイプを含まないこと。
- `title` — 1行目の見出し。改行を含まない空でない1行。
- `overwrite` — 既存ファイルを上書きすることを明示する。省略時は拒否する。

**戻り値**

書き出したファイルの `Path`。

**例外**

- `ValueError` — 長さ・原子タイプ・分率座標・電荷・質量・タイプの対応・`title` のいずれかが規則に合わない。
- `FileExistsError` — `path` が存在し `overwrite` を指定していない。
- `FileNotFoundError` — `path` の親ディレクトリがない。

**注意**

セルは原点 0 の直方晶で、角度は指定しない。電荷は e、長さは Å、質量は g/mol。

## syutyu.analysis.dispatch

LAMMPS と自前計算を明示的に選ぶ。別の推定量へ黙って切り替えない。

### rdf

```python
rdf(source, *, method, normalization='finite', **options)
```

軌跡または LAMMPS の RDF 表から g(r) と積算配位数を返す。

公開API: `from syutyu import rdf`

**引数**

- `source` — `method='custom'` なら `Trajectory` か dump のパス。`method='lammps'` なら `.gr` のパス。
- `method` — `custom`（自前計算）か `lammps`（compute rdf の出力を読む）。必須。
- `normalization` — `finite` か `population`。`lammps` では `finite` のみ。
- `options` — `custom` では `pairs`・`rmax`・`bins`・`start_step`・`stop_step`・`stride`・`intermolecular`。`lammps` では `pairs`・`start_step`・`stop_step`。

**戻り値**

`Result`。`tables['rdf']` に `r_center_A`、`r_upper_A`、`g_<i>`、`CN_<i>`。

**例外**

- `ValueError` — `method` が2つ以外、`lammps` で `population` を指定、その方法で使えないオプションを渡した。
- `TypeError` — `custom` で `pairs` や `rmax` を欠いた。

**注意**

この入口の既定は有限N規格化（Colab と LAMMPS の流儀）。下位の `radial_distribution` は `normalization` を必須にしており、既定を継承しない。

### msd

```python
msd(source, *, method, dt_ps, origin_mode=None, origin_step=None, **options)
```

自前計算の原子 MSD か、LAMMPS compute msd の出力を選んで返す。

公開API: `from syutyu import msd`

**引数**

- `source` — `custom` なら `Trajectory` か dump のパス。`lammps` ならファイル名か系列名→ファイル名の辞書。
- `method` — `custom` か `lammps`。必須。
- `dt_ps` — LAMMPS の 1 step 当たりの ps。
- `origin_mode` — `multiple`（全開始時刻を平均）か `single`（最初のフレームを原点）。`lammps` は単一原点のみ。
- `origin_step` — `lammps` では必須。`custom` では指定できない。
- `options` — `custom` では `selections`（既定 `['all']`）、`start_step`、`stop_step`、`stride`、`laboratory`。`lammps` では `column`。

**戻り値**

`Result`。`tables['msd']` に `lag_ps`、`MSD_<i>_A2`、`origin_atom_pairs_<i>`。

**例外**

- `ValueError` — `method` が2つ以外、`lammps` で `origin_step` がないか `multiple` を要求、`custom` で `origin_step` を渡した、使えないオプションを渡した。

**注意**

原子の MSD であり、分子重心の MSD ではない。重なった時間原点は独立でなく、拡散係数は推定しない。

### vacf

```python
vacf(source, *, method, dt_ps, velocity_unit, origin_mode=None, origin_step=None, **options)
```

自前計算の VACF か、LAMMPS compute vacf の出力を選んで返す。

公開API: `from syutyu import vacf`

**引数**

- `source` — `custom` なら速度を含む `Trajectory` か dump のパス。`lammps` なら fix ave/time の表。
- `method` — `custom` か `lammps`。必須。
- `dt_ps` — LAMMPS の 1 step 当たりの ps。
- `velocity_unit` — `A/fs` か `A/ps`。`A/fs` は二乗する前に `A/ps` へ換算する。
- `origin_mode` — `multiple`（Colab の定義）か `single`（フレーム0のみ）。`lammps` は単一原点のみ。
- `origin_step` — `lammps` では必須。`custom` では指定できない。
- `options` — `custom` では `selections`（既定 `['all']`）、`start_step`、`stride`、`max_lag_ps`。`lammps` では `max_lag_ps`。

**戻り値**

`Result`。`tables['vacf']` に `lag_ps`、成分 `VACF_<i>_x_A2_ps2` など、`VACF_<i>_A2_ps2`、`VACF_<i>_normalized`、`origin_atom_pairs_<i>`。

**例外**

- `ValueError` — 単位・原点モード・方法の指定が不正、`lammps` で `origin_step` がない、使えないオプションを渡した。

**注意**

単一原点と複数原点は別の推定量であり、同じ値として比較しない。速度の平均や重心ドリフトの除去はしない。

## syutyu.analysis.diffusion

明示した時間窓での3次元アインシュタイン拡散の見積もり。

### DiffusionFit

```python
class DiffusionFit
```

1つの時間窓で MSD を直線近似した結果。

公開API: `from syutyu import DiffusionFit`

**属性**

- `window_ps` — 近似に使った経過時間の範囲（両端を含む）。ps。
- `samples` — 窓に入った点の数。
- `slope_A2_ps` — MSD の傾き。Å²/ps。
- `intercept_A2` — 切片。Å²。
- `D_m2_s` — `slope / 6 × 1e-8` で求めた拡散係数。m²/s。
- `r_squared` — 窓内での決定係数。線形性の目安であり、平衡や収束の証拠ではない。
- `limitation` — この見積もりで主張しないことの既定の文。

#### predict

```python
predict(time_ps)
```

近似直線の値を返す。

**引数**

- `time_ps` — 経過時間。ps。スカラーか配列。

**戻り値**

`slope_A2_ps * time_ps + intercept_A2`。Å²。

### fit_diffusion

```python
fit_diffusion(time_ps, msd_A2, *, window_ps)
```

指定した時間窓で MSD を直線近似し、拡散係数を返す。

公開API: `from syutyu import fit_diffusion`

**引数**

- `time_ps` — MSD 原点からの経過時間。ps。狭義単調増加、非負。
- `msd_A2` — 同じ長さの MSD。Å²。非負。
- `window_ps` — 近似に使う範囲 `(下限, 上限)`。ps。データの範囲内で順序どおりであること。

**戻り値**

`DiffusionFit`。`D_m2_s` は `slope / 6 × 1e-8`。

**例外**

- `ValueError` — 長さ・形・有限性・単調性が条件を満たさない、窓が範囲外か逆順、窓内の点が3未満、傾きが負。

**注意**

窓は物理的な拡散領域を自分で選ぶ。最小二乗は自動で窓を選ばず、平衡や収束、独立な誤差を主張しない。

## syutyu.analysis.stats

thermo の記述統計。平衡に達したとは自動で主張しない。

### summarize_thermo

```python
summarize_thermo(table, *, start_ps=None, stop_ps=None)
```

thermo の区間について、平均・標準偏差・点数を返す。

公開API: `from syutyu import summarize_thermo`

**引数**

- `table` — `read_thermo` が返す辞書。`time_ps` 列が必要。
- `start_ps`、`stop_ps` — 集計する時間範囲（両端を含む）。省略時は全体。

**戻り値**

列名→`{'mean', 'std', 'samples'}` の辞書。`time_ps` と `Step` は含めない。

**例外**

- `ValueError` — `time_ps` が非有限か減少、区間が不正、該当サンプルがない、列の長さや値が不正。

**注意**

`std` は母集団標準偏差であり、平均の標準誤差ではない。記録されたサンプルは相関する。

## syutyu.analysis.rdf

種別とサイトを分けた RDF。力場に依存しない。

### radial_distribution

```python
radial_distribution(trajectory, pairs, *, rmax, bins=200, start_step=None, stop_step=None, stride=1, intermolecular=False, normalization)
```

中心と隣接の組ごとに、フレーム平均の RDF と積算配位数を返す。

`from syutyu.analysis.rdf import radial_distribution`

**引数**

- `trajectory` — `Trajectory`。
- `pairs` — `(中心, 隣接)` の選択式の並び。例 `('type:1', 'type:2')`。
- `rmax` — g(r) の上限。Å。対象全フレームの最短辺の半分以下。
- `bins` — ビン数。既定 200。
- `start_step`、`stop_step`、`stride` — 使うフレームの範囲と間引き。
- `intermolecular` — 同じ分子（mol が同じで 0 でない）の原子対を除く。
- `normalization` — `finite`（分母 `N_A*N_B-N_overlap`）か `population`（`N_A*N_B`）。必須。

**戻り値**

`Result`。`tables['rdf']` に `r_center_A`、`r_upper_A`、`g_<i>`、`CN_<i>`、`ordered_counts_<i>`。

**例外**

- `ValueError` — `normalization` が2つ以外、`rmax` が範囲外、`pairs` が空、`intermolecular` なのに mol 列がない。

**注意**

各フレームの体積で規格化し、CN は中心原子1個あたりの直接カウント。可変体積では平均 g に平均密度を掛け直さない。LAMMPS の `special_bonds` の除外は再現しない。

## syutyu.analysis.msd

力場に依存しない原子 MSD。

### mean_squared_displacement

```python
mean_squared_displacement(trajectory, selections, *, dt_ps, start_step=None, stop_step=None, stride=1, laboratory=False, origin_mode='multiple')
```

選択ごとに原子 MSD を計算する。拡散フィットや重心除去は行わない。

`from syutyu.analysis.msd import mean_squared_displacement`

**引数**

- `trajectory` — `Trajectory`。
- `selections` — 選択式の並び。各選択が1つの系列になる。
- `dt_ps` — LAMMPS の 1 step 当たりの ps。
- `start_step`、`stop_step`、`stride` — 使うフレームの範囲と間引き。等間隔であること。
- `laboratory` — 可変セルで、セル伸縮を含む実験室座標の MSD を明示的に求める。
- `origin_mode` — `multiple`（全開始時刻）か `single`（最初の選択フレーム）。

**戻り値**

`Result`。`tables['msd']` に `lag_ps`、`MSD_<i>_A2`、`origin_atom_pairs_<i>`。

**例外**

- `ValueError` — `dt_ps` が正でない、原点モードが不正、フレームが2未満か等間隔でない、選択が空、可変セルなのに非折り返し座標と `laboratory=True` がそろわない。

**注意**

原子の MSD であり、分子重心の MSD ではない。重なった時間原点は独立でなく、平均二乗変位から平衡を主張しない。

## syutyu.analysis.vacf

dump の速度から自前で VACF を求める。compute vacf の代替にはしない。

### velocity_autocorrelation

```python
velocity_autocorrelation(trajectory, selections, *, dt_ps, velocity_unit, start_step=None, stride=1, origin_mode='multiple', max_lag_ps=None)
```

選択ごとに速度自己相関を計算する。複数時間原点が Colab の定義。

`from syutyu.analysis.vacf import velocity_autocorrelation`

**引数**

- `trajectory` — 速度を含む `Trajectory`。
- `selections` — 選択式の並び。各選択が1つの系列になる。
- `dt_ps` — LAMMPS の 1 step 当たりの ps。dump の間隔ではない。
- `velocity_unit` — `A/fs` か `A/ps`。`A/fs` は二乗する前に換算する。
- `start_step`、`stride` — 使うフレームの範囲と間引き。等間隔であること。
- `origin_mode` — `multiple`（全開始時刻の平均）か `single`（フレーム0のみ）。
- `max_lag_ps` — 残す lag の上限（両端を含む）。省略時は全 lag。

**戻り値**

`Result`。`tables['vacf']` に `lag_ps`、成分 `VACF_<i>_x_A2_ps2` など、`VACF_<i>_A2_ps2`、`VACF_<i>_normalized`、`origin_atom_pairs_<i>`。

**例外**

- `ValueError` — `dt_ps` や単位・原点モードが不正、選択が空、フレームが2未満か等間隔でない、速度がないか非有限、lag 上限が範囲外、Cv(0) が0。

**注意**

速度は保存されたまま使い、平均も重心ドリフトも引かない。拡散係数・粘性・スペクトルの推定はしない。

## syutyu.run

明示された LAMMPS 入力だけを実行する。入力スクリプトは信頼済みコードとして扱う。

名前を挙げた補助ファイルだけをコピーする。LAMMPS の shell コマンドや include、ファイルアクセスを隔離する sandbox ではない。物理単位や系固有の設定は推測しない。

### run_lammps

```python
run_lammps(input_file, output_root, *, assets=(), ranks=2, timeout=180, executable='lmp')
```

入力と明示した補助ファイルを新しいディレクトリへ写し、実行してそのパスを返す。

公開API: `from syutyu import run_lammps`

**引数**

- `input_file` — LAMMPS の入力。親フォルダを基準に補助ファイルを解決する。
- `output_root` — run フォルダを作る場所。既存でもよい。
- `assets` — 入力の親フォルダからの相対パス（その配下の絶対パスも可）。サブフォルダは保たれる。
- `ranks` — MPI プロセス数。1 なら `lmp` を直接、2 以上なら `mpirun -np` で起動する。
- `timeout` — 1回の実行の上限秒。`lmp -help` は最大30秒。
- `executable` — LAMMPS の実行ファイル名。PATH から解決する。

**戻り値**

作った run フォルダの `Path`。`run.json` に入力ハッシュ・コマンド・成否・時間・警告・出力を残す。

**例外**

- `ValueError` — `ranks` や `timeout` が不正、`..` による移動、入力の外を指すパス、`run.json` などの予約名、重複した入力。
- `FileNotFoundError` — 実行ファイルや `mpirun` が見つからない。
- `RuntimeError` — `-help` の失敗、非ゼロ終了、出力の `ERROR` 行、時間切れ、その他の実行失敗。

**注意**

失敗した run もフォルダと `run.json` を残す。入力と過去の run は書き換えない。POSIX では時間切れや中断時にプロセスグループを止める。

## syutyu.report.figures

Result の表を共通に描く。pyplot を import しない headless 描画。

### save_figure

```python
save_figure(result, output, view=None, compare=True, *, overwrite=False)
```

Result を PNG に描いて保存する。

公開API: `from syutyu import save_figure`

**引数**

- `result` — `Result`。
- `output` — 保存先のパス。親フォルダは作る。
- `view` — `available_views` のキー。省略時は最初の表示。
- `compare` — 参照線や過去フレームとの比較を重ねる。
- `overwrite` — 既存ファイルへの上書きを明示する。

**戻り値**

保存先の `Path`。解像度は 160 dpi。

**例外**

- `FileExistsError` — 既存ファイルがあり `overwrite` を指定していない。
- `ValueError` — 描ける表がない、または `view` が使えない。

**注意**

GUI を開かない。

### available_views

```python
available_views(result)
```

Result から描ける表示の一覧を返す。

`from syutyu.report.figures import available_views`

**引数**

- `result` — `Result`。

**戻り値**

キー→表示名の辞書。`rdf`、`cn`、`msd`、`vacf`、`vacf-normalized`、`temperature` のうち描けるもの。

**例外**

- `ValueError` — 描ける表がない。

**注意**

`temperature` は `thermo` 表に `Temp` 列があるときだけ加わる。

### draw_view

```python
draw_view(ax, result, view=None, compare=True)
```

1つの軸に表示を描き、図の下へ出す注記を返す。

`from syutyu.report.figures import draw_view`

**引数**

- `ax` — 描画先の Matplotlib の軸。内容は消してから描く。
- `result` — `Result`。
- `view` — 表示のキー。省略時は最初の表示。
- `compare` — 参照線や過去フレームとの比較を重ねる。

**戻り値**

図に添える注記の文字列。規格化や原点の前提を短く書く。

**例外**

- `ValueError` — `view` が `available_views` にない。

### decorate_figure

```python
decorate_figure(fig, result)
```

図の題名・条件・注意書きを summary から書く。

`from syutyu.report.figures import decorate_figure`

**引数**

- `fig` — Matplotlib の Figure。
- `result` — `Result`。`summary` の `title`、`atoms`、`selected_steps`、`method`、`coordinate_frame` などを読む。

**注意**

可変セルの実験室座標 MSD では、補正済みの拡散係数ではないと明示する。

### make_figure

```python
make_figure(result, view=None, compare=True)
```

描画済みの Figure を作って返す。保存はしない。

`from syutyu.report.figures import make_figure`

**引数**

- `result` — `Result`。
- `view` — 表示のキー。省略時は最初の表示。
- `compare` — 参照線や過去フレームとの比較を重ねる。

**戻り値**

Matplotlib の `Figure`。10×6.5 インチ。

**例外**

- `ValueError` — 描ける表がない、または `view` が使えない。

## syutyu.report.bundle

VMD 用の持ち出し ZIP を作る。network・認証・GUI を使わない。

### export_vmd

```python
export_vmd(trajectory, output, *, stride=1)
```

指定フレームと VMD スクリプト・manifest・説明を1つの ZIP にまとめる。

公開API: `from syutyu import export_vmd`

**引数**

- `trajectory` — custom dump。
- `output` — 作る ZIP。既存のパスは拒否する。
- `stride` — 何フレームおきに保存するか（1以上）。1 なら全バイトをそのまま残す。

**戻り値**

manifest の辞書。元ファイルと保存フレームの SHA-256、枚数を含む。

**例外**

- `FileExistsError` — `output` が既にある。
- `FileNotFoundError` — 軌跡がない。
- `ValueError` — `stride` が正の整数でない、保存フレームが0。

**注意**

間引きは表示用で、動力学の解析に使うと元の軌跡と違う結果になる。完成した ZIP だけを公開し、途中のファイルは残さない。

### iter_dump_blocks

```python
iter_dump_blocks(path)
```

1フレーム分の生の dump を、ファイル全体を読み込まずに順に返す。

`from syutyu.report.bundle import iter_dump_blocks`

**引数**

- `path` — custom dump。

**戻り値**

フレームごとのバイト列。原子列は問わずそのまま保つ。

**例外**

- `ValueError` — 見出しが9行でそろわない、原子数が正でない、セル境界が2列か3列でない、原子行が足りない。

**注意**

持ち出しのための検査用で、数値の解析には使わない。一般三斜晶やベクトル形式のセルは扱わない。

## syutyu.kernels

解析と教材用の実装が共有する、NumPy だけの数値カーネル。

### minimum_image

```python
minimum_image(delta, box)
```

相対ベクトルを各方向 [-L/2, L/2] に戻す。

`from syutyu.kernels import minimum_image`

**引数**

- `delta` — 相対ベクトル。Å。末尾の軸が x, y, z。
- `box` — 直方体セルの辺の長さ。Å。

**戻り値**

最小像へ移した相対ベクトル。Å。

**注意**

直方体セルのみ。三斜晶の tilt は扱わない。

### pair_rdf

```python
pair_rdf(positions, box, edges, *, normalization)
```

全フレームの i<j を数え、g(r) と積算配位数を返す。

`from syutyu.kernels import pair_rdf`

**引数**

- `positions` — フレーム×原子×3 の座標。Å。
- `box` — 直方体セルの辺の長さ。Å。
- `edges` — ビンの境界。0 から始まり、最短辺の半分以下で終わること。
- `normalization` — `population`（密度 `N/V`）か `finite`（`(N-1)/V`、Colab の有限N補正）。必須。

**戻り値**

`(g, cn, counts)`。`counts` はビンごとの原子対数、`cn` はビン上端までの平均隣接原子数（自分自身を含まない）。

**例外**

- `ValueError` — `normalization` が2つ以外、フレームがないか原子が2個未満、ビンの範囲が不正。

**注意**

有限Nの理想気体で `population` は `1-1/N`、`finite` は 1 になる。規格化を既定値で継承させないため、キーワード必須にしている。

### unwrap

```python
unwrap(positions, box)
```

隣接フレームの最小像変位を積算して非折り返し座標を作る。

`from syutyu.kernels import unwrap`

**引数**

- `positions` — フレーム×原子×3 の折り返し座標。Å。
- `box` — 直方体セルの辺の長さ。Å。固定セルを前提とする。

**戻り値**

先頭フレームを基準に積算した座標。Å。

**注意**

各成分の真の変位が L/2 未満であることが前提。image flag のない dump では、この前提そのものを確かめられない。

### msd_multi_origin

```python
msd_multi_origin(unwrapped)
```

lag ごとに可能な全開始時刻と原子を平均した MSD を返す。

`from syutyu.kernels import msd_multi_origin`

**引数**

- `unwrapped` — フレーム×原子×3 の非折り返し座標。Å。

**戻り値**

`(msd, samples)`。`msd[lag]` は Å²、`samples[lag]` は使った (開始時刻, 原子) の組の数。

**注意**

重なった標本は独立でない。lag 0 は 0 として返す。
