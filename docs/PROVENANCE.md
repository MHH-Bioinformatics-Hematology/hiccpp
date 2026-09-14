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
- **Writing** follows Juicer tools 1.22.01 (Juicebox tag `1.22`) for version 8
  and 2.20.00 (tag `v2.20.00`) for version 9:
  `juicebox.tools.utils.original.Preprocessor`, `MatrixPP`,
  `MatrixZoomDataPP`, `BlockPP` and `ExpectedValueCalculation`, and
  `juicebox.tools.utils.norm.NormalizationCalculations`, `ZeroScale`,
  `final2.FinalScale`, `NormVectorUpdater` and `NormalizationVectorUpdater`,
  re-implemented in `src/writer.cpp` and `src/norms.cpp`, which name the Java
  method each part corresponds to.
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

### Versions 6 and 7

Reading versions 6 and 7 follows the same hicstraw 1.3.1 source, which reads
versions 6 to 9 and refuses anything below 6 (`readHeader`). Versions 6 and 7
share version 8's header, footer, zoom metadata, expected values and
normalization vectors; version 6 blocks are plain records of int32 binX,
int32 binY and float32 counts (`readBlock`, `version < 7`, as in Juicebox's
`DatasetReaderV2` and hic2cool's `read_block`).

No Juicer tools release that writes version 6 or 7 could be obtained on
2026-09-14, so hicfilecpp does not write them:

- Juicebox's `Preprocessor.java` has `VERSION = 8` in every commit of
  aidenlab/Juicebox, from the file's first commit there (2015-06-08) on, and
  in the tags 1.9.9, 1.22, v1.5.3 and v1.13.01; the version 7 writer predates
  that repository.
- The oldest jars still downloadable from the Aiden lab
  (`hicfiles.tc4ga.com/public/juicer/`: juicer_tools 1.6.2, 1.7.5 and 1.7.6)
  write version 8 (checked with `pre` on a small contact list); older names
  return 403. The GitHub releases start at 2.04.06, and conda offers juicer
  1.6 (Juicer tools 1.22.01) and juicertools 2.20.00.

The version 6 and 7 test files are therefore cut out of a deposited version 7
file by `tests/data/extract_legacy_subset.py`:

| File | Source | Command |
|---|---|---|
| `GM12878_combined_30.chr21_chr22.v7.hic` | `GSE63525_GM12878_insitu_primary+replicate_combined_30.hic` (GEO GSE63525, 39,901,821,731 bytes, version 7) | `extract_legacy_subset.py SOURCE tests/data` (chromosomes 21 and 22, 2.5 Mb, 1 Mb, 500 kb and 250 kb) |
| `GM12878_combined_30.chr21_chr22.v6.hic` | the same | the same run |

The version 7 file keeps the source's genome, attributes, zoom metadata,
compressed blocks, expected values and normalization vectors byte for byte
and rewrites only positions, sizes and chromosome indices (All, 21 and 22
become 0, 1 and 2; fragment resolutions are dropped). The version 6 file is
the same with every block re-encoded in the version 6 record layout and
deflated with zlib's default level; no real version 6 deposit was at hand,
so version 6 is validated on real pixels in the version 6 block layout
rather than on a file an old Juicer release wrote. Checked when the files
were made: hicstraw returns the same records from the source, the version 7
file and the version 6 file for every pair, resolution, normalization (NONE,
VC, VC_SQRT, KR, GW_KR, GW_VC, INTER_KR, INTER_VC) and matrix type (192
queries, 856,380 records), and hic2cool converts both files at 250 kb into
the same pixels and vectors.
