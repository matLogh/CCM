# Drift Correction Time-Evolution Tools

This directory contains two ROOT macros for combining and inspecting detector
time-evolution histograms written by `solveTimeEvo_ILL`:

```text
combine_offset_gain_timeevo.C
view_timeevo_det.C
```

`combine_offset_gain_timeevo.C` combines per-run-range ROOT files into one
output file per detector. `view_timeevo_det.C` opens those combined files and
draws detector comparisons in an interactive canvas.

## Combine Macro

`combine_offset_gain_timeevo.C` reads the per-run-range ROOT files produced by
the offset/gain correction workflows, rebins the energy axis by a factor of 2,
and writes wide run-vs-energy `TH2D` histograms covering the full run interval.

### Output Location

Combined files are written to:

```text
output/combined_offset_and_gain_timeevo/
```

For detector `det`, the output file is:

```text
RunVsEnergy_all_ranges_det<det>_correctedTimeEvo.root
```

Existing output files with the same name are overwritten.

### Input Directories

The default inputs are:

```text
../offset_and_gain/
../gain_only/
```

There is also an optional test input:

```text
../test/
```

All input directories are expected to use the same file naming pattern:

```text
RunVsEnergy_<runrange>_det<det>_correctedTimeEvo.root
```

The hard-coded run ranges are:

```text
242083_243791
243796_245219
245222_245290
246313_247251
247253_248623
```

### Expected Input Objects

From `../offset_and_gain/`, the macro reads:

```text
RunVsEnergy_det<det>
RunVsEnergy_det<det>_fixed
shift_ROI_0
shift_ROI_1
```

From `../gain_only/`, the macro reads:

```text
RunVsEnergy_det<det>_fixed
shift_ROI_0
```

When `include_test` is enabled, from `../test/` it reads:

```text
RunVsEnergy_det<det>_fixed
```

Missing gain-only or test files are skipped. A detector output is only skipped
entirely if none of the main `../offset_and_gain/` files can be used.

### Output Objects

Each output file contains:

```text
RunVsEnergy_det<det>_all_ranges
RunVsEnergy_det<det>_fixed_all_ranges
RunVsEnergy_det<det>_gain_only_all_ranges
shift_ROI_0_det<det>_all_ranges
shift_ROI_1_det<det>_all_ranges
ROI_0_gain_only
```

If `include_test` is enabled, the file also contains:

```text
RunVsEnergy_det<det>_test_all_ranges
```

### Binning

The combined histograms use one x bin per run number over:

```text
242083 to 248623
```

The input energy axis is assumed to have:

```text
32768 bins
```

The output energy axis is rebinned by:

```text
2
```

so the output histograms have `16384` y bins spanning `0` to `32768`.

The default display range is set to:

```text
3000 to 4000
```

This affects the axis range stored on the histogram, not the filled bin
contents.

### Running The Combine Macro

Load the macro in ROOT:

```cpp
.L combine_offset_gain_timeevo.C
```

Combine all detectors `0` through `31`:

```cpp
combine_offset_gain_timeevo();
```

Combine one detector:

```cpp
combine_offset_gain_timeevo(7, 7);
```

Combine a detector range:

```cpp
combine_offset_gain_timeevo(0, 15);
```

Include the optional `../test/` comparison histogram:

```cpp
combine_offset_gain_timeevo(0, 31, true);
```

Compile with ACLiC before running:

```cpp
.L combine_offset_gain_timeevo.C+
combine_offset_gain_timeevo(7, 7, true);
```

### Console Output

For each detector, the macro prints:

```text
output file path
input files being read
missing file/object warnings
number of main files used
number of gain-only files used
number of test files used, if enabled
combined histogram binning
combined graph point counts
```

Use the file counts to check that each detector used the expected run ranges.

## Viewer Macro

`view_timeevo_det.C` is an interactive ROOT viewer for the combined
time-evolution histograms produced by `combine_offset_gain_timeevo.C`.

The viewer opens one detector at a time and draws two stacked `TH2D` pads in a
single canvas. Both pads use log-z scale, crosshairs, and small left/right
margins.

### Viewer Input

The viewer reads combined files from:

```text
output/combined_offset_and_gain_timeevo/
```

For detector `det`, it opens:

```text
RunVsEnergy_all_ranges_det<det>_correctedTimeEvo.root
```

Expected histograms:

```text
RunVsEnergy_det<det>_all_ranges
RunVsEnergy_det<det>_fixed_all_ranges
RunVsEnergy_det<det>_gain_only_all_ranges
```

### Loading A Detector

Load the default comparison, raw on top and offset+gain fixed on bottom:

```cpp
LoadDet(7);
```

Load the gain-only comparison:

```cpp
LoadDet(7, true);
```

In gain-only mode, the pads are:

```text
top:    RunVsEnergy_det<det>_fixed_all_ranges
bottom: RunVsEnergy_det<det>_gain_only_all_ranges
```

This makes it easy to compare offset+gain correction against gain-only
correction.

### Canvas Size

Default canvas size:

```text
3840 x 864
```

This is intended for the wide two-monitor setup.

Compact canvas size:

```text
1512 x 820
```

Use compact mode with the third argument:

```cpp
LoadDet(7, false, true); // raw vs offset+gain fixed
LoadDet(7, true, true);  // offset+gain fixed vs gain-only
```

Use an exact custom width with:

```cpp
LoadDetWidth(7, false, 1400); // raw vs offset+gain fixed
LoadDetWidth(7, true, 1400);  // offset+gain fixed vs gain-only
```

`LoadDetWidth` keeps the standard height of `864`.

### Axis Control

Zoom both pads to the same energy range:

```cpp
ScaleY(3000, 4000);
```

Zoom both pads to the same run-number range:

```cpp
ScaleX(242083, 243791);
```

Reset both pads:

```cpp
ResetY();
```

If bad-run markers are currently drawn, they are redrawn after `ScaleY` and
`ScaleX`.

### Bad-Run Overlay

Overlay bad run numbers on the bottom pad as red `x` markers. The input file is
a text file containing run numbers, with comment/header lines allowed.

Use an explicit bad-run file and one y value:

```cpp
OverlayBadRuns("../bad_runs/det7_badruns.txt", 3600.0);
```

Use multiple y values:

```cpp
OverlayBadRuns("../bad_runs/det7_badruns.txt", 3200.0, 3600.0, 3900.0);
```

List-style syntax also works:

```cpp
OverlayBadRuns("../bad_runs/det7_badruns.txt", {3200.0, 3600.0, 3900.0});
```

There are shortcuts for the detector-7 bad-run file:

```cpp
OverlayBadRuns(3600.0);
OverlayBadRuns(3200.0, 3600.0, 3900.0);
OverlayBadRuns({3200.0, 3600.0, 3900.0});
```

Calling `LoadDet` clears any previous bad-run overlay. Call `OverlayBadRuns`
again after loading a new detector or changing between normal and gain-only
views.

### Typical Sessions

Raw vs offset+gain fixed, compact laptop view:

```cpp
.L view_timeevo_det.C
LoadDet(7, false, true);
ScaleY(3000, 4000);
OverlayBadRuns(3200.0, 3600.0, 3900.0);
```

Offset+gain fixed vs gain-only, compact laptop view:

```cpp
.L view_timeevo_det.C
LoadDet(7, true, true);
ScaleY(3000, 4000);
OverlayBadRuns("../bad_runs/det7_badruns.txt", 3600.0);
```

### Viewer Notes

`LoadDet` deletes and recreates only the canvas owned by this viewer. Other ROOT
windows, such as a `TBrowser`, should remain open unless they use the same ROOT
object name as the viewer canvas.

The viewer canvas is named:

```text
c_det<det>
```
