# Correcting time-variation in ILL RunVsEnergy matrices

The ILL programs are ready-to-run drivers for applying CCM corrections to ROOT
matrices stored as run/time on the X axis and energy on the Y axis. The expected
ROOT layout is:

```text
RunVsEnergy_START_END/RunVsEnergy_detNUMBER
```

For example, detector 7 in runs 242083-243791 is read as:

```text
RunVsEnergy_242083_243791/RunVsEnergy_det7
```

Full option lists can be obtained with:

```bash
detectTimeEvo_ILL --help
solveTimeEvo_ILL --help
```

## Automated time evolution detection (not tested)

The `detectTimeEvo_ILL` executable checks whether one detector, or all detectors,
show an energy shift above a configurable threshold. It automatically chooses a
reference window around the first third of the matrix, calculates CCM shifts for
one ROI, and reports detectors whose shifts exceed the threshold.

Minimal example:

```bash
detectTimeEvo_ILL \
  --rootfile RunVsEnergy_all.root \
  --start_run 242083 \
  --end_run 243791 \
  --detector 7 \
  --ROIsource 66Ga
```

Useful options:

- `--alldet` checks detectors 0-40 and skips missing matrices.
- `--shift_threshold <value>` changes the reporting threshold, in energy units.
- `--draw` opens ROOT canvases for detectors above threshold.
- `--ROI <desired> <low> <high> <shift_low> <shift_high>` defines the ROI by hand.
- `--ROIsource <source>` defines a built-in ROI for sources such as `60Co`, `66Ga`,
  `133Ba`, and `226Ra`.

## Solving the time evolution

The `solveTimeEvo_ILL` executable calculates correction parameters for one
detector and writes a `TimeEvoCC.conf`-style output file. It expects the same ROOT
layout described above.

Minimal fast example:

```bash
solveTimeEvo_ILL \
  --rootfile RunVsEnergy_all.root \
  --start_run 242083 \
  --end_run 243791 \
  --detector 7 \
  --ROI_file solutions/ILL/ROI_det7.txt \
  --super_settings
```

Without `--super_settings`, the code runs a small grid search over hardcoded CCM
settings and chooses the setting with the best FWFM for a user-provided test
peak. In that mode, `--fit_peak <center> <low> <high>` is required. Use a peak
that is not already one of the correction ROIs, otherwise the optimization can
overfit the feature being corrected. (--fit_peak not tested)

Example with optimization:

```bash
solveTimeEvo_ILL \
  --rootfile RunVsEnergy_all.root \
  --start_run 242083 \
  --end_run 243791 \
  --detector 7 \
  --ROI_file solutions/ILL/ROI_det7.txt \
  --fit_peak 4295.187 4220 4360
```

## ROI input

ROIs can be supplied one at a time on the command line:

```bash
--ROI <desired> <low> <high> <shift_low> <shift_high>
```

or collected in a text file:

```text
ref_time 247800 247804
ROI 19600 19300 19800 -400 400
ROI 3845 3800 3900 -150 150
```

The ROI values are:

- `desired`: nominal energy of the peak or feature after correction.
- `low`, `high`: energy window used to build the reference and sample vectors.
- `shift_low`, `shift_high`: allowed displacement window around the ROI.

The solver adjusts the nominal ROI peak positions with a Theuerkauf peak fit
before running CCM.

## Correction functions

By default, `solveTimeEvo_ILL` fits a gain-only correction:

```text
[0]*x
```

Two alternate correction forms are available:

- `--gain_with_offset` uses `[0] + [1]*x` and requires at least two ROIs.
- `--gain_with_offset_and_sqrt` uses `[0] + [1]*sqrt(x) + [2]*x` and requires at
  least three ROIs.

These flags are mutually exclusive. The output table columns change to match the
selected correction function.

## Chained ranges

Use `--chain_ranges <start> <end> [...]` to apply the same optimized settings and
the same reference vectors to additional run ranges. The primary range is the one
given by `--start_run` and `--end_run`; chained ranges are then processed from
the same ROOT file using directories named `RunVsEnergy_START_END`.

Example:

```bash
solveTimeEvo_ILL \
  --rootfile RunVsEnergy_all.root \
  --start_run 242083 \
  --end_run 243791 \
  --detector 7 \
  --ROI_file solutions/ILL/ROI_det7.txt \
  --super_settings \
  --chain_ranges 243792 245000 245001 246000
```

## Outputs

Output files are written next to the input ROOT file, with names based on the run
range and detector:

```text
RunVsEnergy_START_END_detNUMBER_TimeEvoCC.conf
RunVsEnergy_START_END_detNUMBER_correctedTimeEvo.root
RunVsEnergy_START_END_detNUMBER_CCMconf.txt
```

The `TimeEvoCC.conf` file contains the correction parameters for each run/time
bin. The diagnostic ROOT file contains the original matrix, corrected matrix,
projection, ROI shift graphs, and interpolation graphs. The `CCMconf.txt` file is
written only when the optimizer mode is used and records the tested settings and
their costs.

## Simple examples

The `simple_ILL.cpp` file is a hardcoded example used for quick local tests and
plotting. It reads a specific matrix from `RunVsEnergy_all.root`, builds a CCM
object directly in C++, writes shift and fit tables, and saves diagnostic ROOT
output. For production use, prefer `detectTimeEvo_ILL` and `solveTimeEvo_ILL`.

## Helper scripts

Some simple ROOT scripts were added to `./timeevo_tools` to quickly visualize
and check the results of the time evolution correction.
