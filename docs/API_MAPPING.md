# hicstraw to C++ API mapping

hicfilecpp reproduces the reading API of hicstraw 1.3.1. Names, parameters
and semantics follow hicstraw; the spelling follows C++.

Conventions:

- `import hicstraw` corresponds to `#include <hicfilecpp/hicfilecpp.hpp>` and
  the namespace `hicfilecpp`.
- A numpy array of floats is a `std::vector<double>`; the float32 matrix of
  `getRecordsAsMatrix` is a `FloatMatrix` (row major, `rows`, `cols`,
  `values`).
- Where hicstraw prints a message and returns nothing, hicfilecpp returns the
  same nothing and keeps the message (`MatrixZoomData::message()`). Where
  hicstraw's behaviour is undefined, hicfilecpp throws `HicError`
  (see DEVIATIONS.md).

## Reading

| Python (hicstraw 1.3.1) | C++ (hicfilecpp) |
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
| `HiCFile::version()` | 8 or 9 |
| `HiCFile::attributes()` | the header's key/value attributes, in file order |
| `HiCFile::getFragResolutions()` | FRAG resolutions |
| `HiCFile::fragmentSiteCounts()` | restriction sites per chromosome |
| `HiCFile::masterIndexPosition()`, `normVectorIndexHeader()` | footer positions |
| `HiCFile::hasMatrix(c1, c2)`, `matrixZoomHeaders(c1, c2)` | the matrices and their zoom headers |
| `HiCFile::expectedValuesKeys()`, `readExpectedValues(key)` | expected-value entries with their per-chromosome factors |
| `HiCFile::normVectorIndex()`, `getNormalizationTypes()`, `readNormVector(norm, chr, unit, resolution)` | normalization vectors |
| `MatrixZoomData::found()`, `message()` | whether the query found its data |
| `MatrixZoomData::zoomHeader()`, `blockIndex()`, `readBlock(entry)` | the zoom header and blocks, records in bins |
| `MatrixZoomData::forEachBlock(visit, threads)` | every block in ascending number, decoded on several threads |

Errors are `hicfilecpp::HicError`, derived from `std::runtime_error`.
