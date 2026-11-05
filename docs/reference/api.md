# APIリファレンス

各項目の本文は docstring から Sphinx の autodoc と napoleon が生成する。docstring は Google スタイルで書く。このページには載せる名前と順番だけを書く。

- 対象は `syutyu` の公開クラスと関数。`cli.py` の内部、`_` で始まる名前は対象外。
- `syutyu.` の直下に書いた名前は `from syutyu import ...` で使える公開API（18名）。`syutyu.io.dump.frame_indices` のように submodule を含めて書いた名前は、その submodule から import する。
- 単位・計算方法・実行例は [利用ガイド](../guide/usage.md)、[解析方法](../guide/analysis-methods.md)。

## 一覧

```{eval-rst}
.. autosummary::
   :nosignatures:

   syutyu.Result
   syutyu.save_result
   syutyu.load_result
   syutyu.Trajectory
   syutyu.read_trajectory
   syutyu.io.dump.frame_indices
   syutyu.read_thermo
   syutyu.read_thermo_blocks
   syutyu.io.native.read_rdf
   syutyu.io.native.read_msd
   syutyu.io.native.read_vacf
   syutyu.Site
   syutyu.write_charge_data
   syutyu.rdf
   syutyu.msd
   syutyu.vacf
   syutyu.DiffusionFit
   syutyu.fit_diffusion
   syutyu.summarize_thermo
   syutyu.analysis.rdf.radial_distribution
   syutyu.analysis.msd.mean_squared_displacement
   syutyu.analysis.vacf.velocity_autocorrelation
   syutyu.run_lammps
   syutyu.save_figure
   syutyu.report.figures.available_views
   syutyu.report.figures.draw_view
   syutyu.report.figures.decorate_figure
   syutyu.report.figures.make_figure
   syutyu.export_vmd
   syutyu.report.bundle.iter_dump_blocks
   syutyu.kernels.minimum_image
   syutyu.kernels.pair_rdf
   syutyu.kernels.unwrap
   syutyu.kernels.msd_multi_origin
```

## syutyu.io.store

```{eval-rst}
.. automodule:: syutyu.io.store

.. autoclass:: syutyu.Result

.. autofunction:: syutyu.save_result

.. autofunction:: syutyu.load_result
```

## syutyu.io.dump

```{eval-rst}
.. automodule:: syutyu.io.dump

.. autoclass:: syutyu.Trajectory
   :members: select, fixed_box

.. autofunction:: syutyu.read_trajectory

.. autofunction:: syutyu.io.dump.frame_indices
```

## syutyu.io.thermo

```{eval-rst}
.. automodule:: syutyu.io.thermo

.. autofunction:: syutyu.read_thermo

.. autofunction:: syutyu.read_thermo_blocks
```

## syutyu.io.native

```{eval-rst}
.. automodule:: syutyu.io.native

.. autofunction:: syutyu.io.native.read_rdf

.. autofunction:: syutyu.io.native.read_msd

.. autofunction:: syutyu.io.native.read_vacf
```

## syutyu.io.lammps_data

```{eval-rst}
.. automodule:: syutyu.io.lammps_data

.. autoclass:: syutyu.Site

.. autofunction:: syutyu.write_charge_data
```

## syutyu.analysis.dispatch

```{eval-rst}
.. automodule:: syutyu.analysis.dispatch

.. autofunction:: syutyu.rdf

.. autofunction:: syutyu.msd

.. autofunction:: syutyu.vacf
```

## syutyu.analysis.diffusion

```{eval-rst}
.. automodule:: syutyu.analysis.diffusion

.. autoclass:: syutyu.DiffusionFit

.. autofunction:: syutyu.fit_diffusion
```

## syutyu.analysis.stats

```{eval-rst}
.. automodule:: syutyu.analysis.stats

.. autofunction:: syutyu.summarize_thermo
```

## syutyu.analysis.rdf

```{eval-rst}
.. automodule:: syutyu.analysis.rdf

.. autofunction:: syutyu.analysis.rdf.radial_distribution
```

## syutyu.analysis.msd

```{eval-rst}
.. automodule:: syutyu.analysis.msd

.. autofunction:: syutyu.analysis.msd.mean_squared_displacement
```

## syutyu.analysis.vacf

```{eval-rst}
.. automodule:: syutyu.analysis.vacf

.. autofunction:: syutyu.analysis.vacf.velocity_autocorrelation
```

## syutyu.run

```{eval-rst}
.. automodule:: syutyu.run

.. autofunction:: syutyu.run_lammps
```

## syutyu.report.figures

```{eval-rst}
.. automodule:: syutyu.report.figures

.. autofunction:: syutyu.save_figure

.. autofunction:: syutyu.report.figures.available_views

.. autofunction:: syutyu.report.figures.draw_view

.. autofunction:: syutyu.report.figures.decorate_figure

.. autofunction:: syutyu.report.figures.make_figure
```

## syutyu.report.bundle

```{eval-rst}
.. automodule:: syutyu.report.bundle

.. autofunction:: syutyu.export_vmd

.. autofunction:: syutyu.report.bundle.iter_dump_blocks
```

## syutyu.kernels

```{eval-rst}
.. automodule:: syutyu.kernels

.. autofunction:: syutyu.kernels.minimum_image

.. autofunction:: syutyu.kernels.pair_rdf

.. autofunction:: syutyu.kernels.unwrap

.. autofunction:: syutyu.kernels.msd_multi_origin
```
