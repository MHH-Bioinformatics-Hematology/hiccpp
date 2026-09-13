# Provenance

## Implementations this library follows

- **Reading** follows hicstraw 1.3.1, `pybind11_python/src/straw.cpp` of
  [aidenlab/straw](https://github.com/aidenlab/straw) at commit
  `82fba9cf7f323a432bad4ac9a33d1ee65de00a47` (MIT licence, Aiden Lab). The
  query semantics, including the float arithmetic of normalized and
  observed/expected values and the block selection of version 8 and version 9
  files, are re-implemented from that source; `src/reader.cpp` names the
  hicstraw function each part corresponds to.
- **The file format** is the Juicer `.hic` format as written by Juicer tools
  ([aidenlab/Juicebox](https://github.com/aidenlab/Juicebox), MIT licence):
  version 8 by release 1.22.01 (`juicebox.tools.utils.original.Preprocessor`
  with `VERSION = 8`) and version 9 by release 2.20.00 (`VERSION = 9`).
- hic2cool 1.0.1 ([4dn-dcic/hic2cool](https://github.com/4dn-dcic/hic2cool),
  MIT licence) served as a second reading of the version 8 layout.

## Test files in tests/data

`tests/data/generate_test_files.py` writes them from the 5 kb records of
HiCExplorer's `hicexplorer/test/test_data/hicHyperoptDetectLoopsHiCCUPS/SRR1791297_30.hic`,
read with hicstraw 1.3.1, on OpenJDK 11.0.1 (Zulu 11.2+3):

| File | Written by | Command |
|---|---|---|
| `SRR1791297_30.juicer_tools_1.22.01.v8.hic` | Juicer tools 1.22.01 | `pre -j 1 -r 1000000,250000,50000,10000 contacts.txt out.hic sacCer3.chrom.sizes` |
| `SRR1791297_30.juicer_tools_1.22.01.frag.v8.hic` | Juicer tools 1.22.01 | `pre -j 1 -f sites.txt -r 500000,50000,100f,20f contacts.txt out.hic sacCer3.chrom.sizes` |
| `SRR1791297_30.juicer_tools_2.20.00.v9.hic` | Juicer tools 2.20.00 | `pre -j 1 -k VC,VC_SQRT,KR,SCALE,INTER_SCALE,GW_SCALE -r 1000000,250000,50000,10000 contacts.txt out.hic sacCer3.chrom.sizes` |

Juicer tools 1.22.01 is the `share/juicer-1.6-0/juicer_tools.jar` of the
conda package `HCC::juicer-1.6-cpu_py310hc74bb38_0` (its manifest reports
1.22.01; the original download location no longer exists). Juicer tools
2.20.00 is `juicer_tools.2.20.00.jar` of the Juicebox v2.20.00 GitHub release.
Both jars target Java 8 (class file version 52).

`sites.txt` places a synthetic restriction site every 2 kb; the fragment maps
exist only to give the reader FRAG resolutions and FRAG vectors to read.
