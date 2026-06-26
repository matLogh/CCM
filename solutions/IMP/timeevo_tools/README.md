# IMP Time-Evolution Helper Macros

These ROOT macros are helpers for inspecting outputs written by
`solveTimeEvo_IMP`.

## Viewing One Corrected Matrix

After running `solveTimeEvo_IMP`, the diagnostic output is written next to the
input file:

```text
HistAll_<histogram>_correctedTimeEvo.root
```

For example, `--group l --detector 7` writes:

```text
HistAll_hetl7_correctedTimeEvo.root
```

Open the viewer from the directory containing that output:

```cpp
.L solutions/IMP/timeevo_tools/view_timeevo_det.C
LoadDet(7, "l")
LoadDet(3, "g")
```

The macro draws the raw histogram (`hetl7`, for example) and the corrected
histogram (`hetl7_fixed`).

## Combining Outputs

IMP `HistAll.root` files use flat histograms (`hetg*`, `hetl*`) rather than ILL
run-range directories. The copied ILL combiner is therefore replaced with a small
stub that explains this and does not try to combine run ranges.
