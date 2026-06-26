# Correcting time-variation in IMP energy-time matrices

The IMP programs are ready-to-run drivers for applying CCM corrections to
energy-time `TH2D` matrices in a flat ROOT file. The expected `HistAll.root`
layout contains histograms named:

```text
hetg0 ... hetg11
hetl0 ... hetl25
```

Full option lists can be obtained with:

```bash
detectTimeEvo_IMP --help
solveTimeEvo_IMP --help
```

## Automated time evolution detection

`detectTimeEvo_IMP` checks whether one IMP matrix, one group, or all IMP
matrices show an energy shift above a configurable threshold.

Minimal examples:

```bash
detectTimeEvo_IMP \
  --rootfile HistAll.root \
  --group l \
  --detector 7 \
  --ROIsource 66Ga

detectTimeEvo_IMP \
  --rootfile HistAll.root \
  --alldet \
  --ROI 19600 19300 19800 -400 400
```

Useful options:

- `--group <g|l>` selects `hetgNUMBER` or `hetlNUMBER`.
- `--alldet` checks `hetg0-11` and `hetl0-25`; with `--group`, it checks only
  that group.
- `--shift_threshold <value>` changes the reporting threshold, in energy units.
- `--draw` opens ROOT canvases for matrices above threshold.

## Solving the time evolution

`solveTimeEvo_IMP` calculates correction parameters for one IMP matrix and
writes a `TimeEvoCC.conf`-style output file.

Minimal fast example:

```bash
solveTimeEvo_IMP \
  --rootfile HistAll.root \
  --group l \
  --detector 7 \
  --ROI_file solutions/IMP/ROI_det7.txt \
  --super_settings
```

You can also select a histogram by exact name:

```bash
solveTimeEvo_IMP \
  --rootfile HistAll.root \
  --histogram hetg3 \
  --ROI 19600 19300 19800 -400 400 \
  --ref_time 243291 243295 \
  --super_settings
```

Without `--super_settings`, the code runs a small grid search over hardcoded CCM
settings and chooses the setting with the best FWFM for a user-provided test
peak. In that mode, `--fit_peak <center> <low> <high>` is required.

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

The solver adjusts nominal ROI peak positions with a Theuerkauf peak fit before
running CCM.

## Outputs

Output files are written next to the input ROOT file, with names based on the
selected histogram:

```text
HistAll_hetl7_TimeEvoCC.conf
HistAll_hetl7_correctedTimeEvo.root
HistAll_hetl7_CCMconf.txt
```

The diagnostic ROOT file contains the original matrix, corrected matrix,
projection, ROI shift graphs, and interpolation graphs.
