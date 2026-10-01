# hicstraw to C++ API mapping

hiccpp reproduces the reading API of hicstraw 1.3.1 and the writing of
Juicer tools pre and addNorm (1.22.01 for version 8, 2.20.00 for version 9).
Names, parameters and semantics follow those programs; the spelling follows
C++.

Conventions:

- `import hicstraw` corresponds to `#include <hiccpp/hiccpp.hpp>` and
  the namespace `hiccpp`.
- A numpy array of floats is a `std::vector<double>`; the float32 matrix of
  `getRecordsAsMatrix` is a `FloatMatrix` (row major, `rows`, `cols`,
  `values`).
- Where hicstraw prints a message and returns nothing, hiccpp returns the
  same nothing and keeps the message (`MatrixZoomData::message()`). Where
  hicstraw's behaviour is undefined, hiccpp throws `HicError`
  (see DEVIATIONS.md).

## Reading

| Python (hicstraw 1.3.1) | C++ (hiccpp) |
|---|---|
| `hicstraw.HiCFile(path)` | `HiCFile(const std::string& path)` |
| `f.getGenomeID()` | `std::string HiCFile::getGenomeID()` |
| `f.getResolutions()` | `std::vector<int32_t> HiCFile::getResolutions()` (BP resolutions) |
| `f.getChromosomes()` | `std::vector<Chromosome> HiCFile::getChromosomes()` |
| `chromosome.name`, `.index`, `.length` | `Chromosome::name`, `::index`, `::length` |
| `f.getMatrixZoomData(chr1, chr2, matrix_type, norm, unit, resolution)` | `MatrixZoomData HiCFile::getMatrixZoomData(chr1, chr2, matrixType, norm, unit, resolution)` |
| `mzd.getRecords(gx0, gx1, gy0, gy1)` | `std::vector<ContactRecord> MatrixZoomData::getRecords(gx0, gx1, gy0, gy1)` |
| `record.binX`, `.binY`, `.counts` | `ContactRecord::binX`, `::binY`, `::counts` |
| `mzd.getRecordsAsMatrix(gx0, gx1, gy0, gy1)` | `FloatMatrix MatrixZoomData::getRecordsAsMatrix(gx0, gx1, gy0, gy1)` |
| `mzd.getNormVector(index)` | `std::vector<double> MatrixZoomData::getNormVector(index)` |
| `mzd.getExpectedValues()` | `std::vector<double> MatrixZoomData::getExpectedValues()` |
| (C++ straw) `getNumberOfTotalRecords()` | `int64_t MatrixZoomData::getNumberOfTotalRecords()` |

`matrix_type` is `"observed"`, `"oe"` or `"expected"`; `norm` is `"NONE"` or
any normalization the file holds (`"VC"`, `"VC_SQRT"`, `"KR"`, `"SCALE"`,
`"GW_SCALE"`, ...); `unit` is `"BP"` or `"FRAG"`.

## Beyond hicstraw

| C++ | Returns |
|---|---|
| `HiCFile::version()` | 6 to 9 |
| `HiCFile::attributes()` | the header's key/value attributes, in file order |
| `HiCFile::getFragResolutions()` | FRAG resolutions |
| `HiCFile::fragmentSiteCounts()` | restriction sites per chromosome |
| `HiCFile::fragmentSites()` | the site positions themselves, per chromosome |
| `HiCFile::masterIndexPosition()`, `normVectorIndexHeader()` | footer positions |
| `HiCFile::hasMatrix(c1, c2)`, `matrixZoomHeaders(c1, c2)` | the matrices and their zoom headers |
| `HiCFile::expectedValuesKeys()`, `readExpectedValues(key)` | expected-value entries with their per-chromosome factors |
| `HiCFile::normVectorIndex()`, `getNormalizationTypes()`, `readNormVector(norm, chr, unit, resolution)` | normalization vectors |
| `HiCFile::hasNormalizedExpectedSection()` | false when the footer ends after the raw expected values |
| `MatrixZoomData::found()`, `message()` | whether the query found its data |
| `MatrixZoomData::zoomHeader()`, `blockIndex()`, `readBlock(entry)` | the zoom header and blocks, records in bins |
| `MatrixZoomData::forEachBlock(visit, threads)` | every block in ascending number, decoded on several threads |

Errors are `hiccpp::HicError`, derived from `std::runtime_error`.

## Writing

| Juicer tools | C++ (hiccpp) |
|---|---|
| `pre [options] <infile> <outfile> <genomeID>` | `writeHicFile(outfile, WriteOptions{.genomeId = genomeID, .chromosomes = ...}, source)` |
| `<infile>` contacts, "short with score" | a `PixelSource`: pixels of one chromosome pair at `sourceResolution`, each counting as one contact at the start of its bins |
| the chromosome sizes file `<genomeID>` | `WriteOptions::chromosomes`, in file order; "All" is added in front |
| (version) 1.22.01 writes 8, 2.20.00 writes 9 | `WriteOptions::version` |
| `-r <resolutions>` | `WriteOptions::resolutions` (base pair) and `fragResolutions` (the "f" suffixed entries) |
| `-f <sites file>` | `WriteOptions::fragmentSites`, one site list per chromosome |
| `-k <normalizations>` | `WriteOptions::normalizations` (VC, VC_SQRT, KR, SCALE) |
| `-n` | `WriteOptions::normalizations = {}` |
| `-j <threads>` | `WriteOptions::threads` (block compression; no effect on the output) |
| `addNorm <file>` with `-k` | `addNorm(file, normalizations, threads)` |
| `addNorm <file> <vector file>` (custom vectors) | `WriteOptions::providedNormalizations` and `normVector` at write time: vectors stored as given, under any label |

Not available: `-d`, `-m`, `-q`, `-c`, `-t`, `-s`, `-g`, `-z`, `-a`, position
randomization, `--v9-depth-base` (always 2), genome-wide and inter-chromosomal
normalizations (`-w`, GW_*, INTER_*) computed by the library (provided vectors
may carry these labels), custom expected value files, and `addNorm` on a file
with fragment resolutions.

`PixelSource::pixels(resolution, chr1, chr2, consume)` hands the pixels of one
chromosome pair to `consume` in batches of `Pixel{bin1, bin2, count}`.
`WriteOptions::sourceProvidesEveryResolution` makes the writer ask the source
for every resolution instead of binning `sourceResolution`.

### Fragment resolutions

`WriteOptions::fragResolutions` and `fragmentSites` make the writer add the
fragment maps of Juicer tools `pre -f`, for version 8 and version 9 alike. The
pixels come from `PixelSource::fragPixels`, whose bins count restriction
fragments: the writer asks for `sourceFragResolution` and bins it up to each
coarser fragment resolution, the way it bins `sourceResolution` for the base
pair ones, and `sourceProvidesEveryResolution` covers both. `fragPixels` has a
default implementation that hands over nothing, so a source written before this
option keeps compiling and then yields a file whose matrices hold base pair
zooms only.

| Option | Meaning |
|---|---|
| `fragResolutions` | fragment resolutions, written after the base pair ones, coarsest first |
| `fragmentSites` | the sites file: one ascending list of positions per entry of `chromosomes` |
| `sourceFragResolution` | the fragment resolution `fragPixels` is asked for; every entry of `fragResolutions` is a multiple of it |
| `fragNormVector` | the `normVector` of the fragment resolutions, for `providedNormalizations` |

A chromosome's fragment bin count at resolution N is its site count divided by
N plus one, and its block layout follows from that count the way a base pair
zoom's follows from the chromosome length. Expected values and normalization
vectors are written per fragment resolution under the unit "FRAG". The
whole-genome "All" matrix stays a base pair matrix, as it is in Juicer's files,
and "All" carries no sites.
