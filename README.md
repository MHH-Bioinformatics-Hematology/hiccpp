# hicfilecpp

A C++20 library that reads the Juicer `.hic` format, versions 6 to 9, and
writes versions 8 and 9, without Python or Java.

- Reading reproduces hicstraw 1.3.1: `HiCFile`, `getMatrixZoomData`,
  `getRecords`, `getRecordsAsMatrix`, `getNormVector` and `getExpectedValues`
  return what hicstraw returns, record for record. Versions 6 and 7 are the
  older deposits (for example GEO GSE63525, Rao et al. 2014).
- Writing versions 6 and 7 is refused: no Juicer tools release that writes
  them can be obtained to validate against (docs/PROVENANCE.md).
- Writing reproduces Juicer tools `pre` and `addNorm`: version 8 files follow
  release 1.22.01, version 9 files release 2.20.00. The pixels, expected
  values and normalization vectors (VC, VC_SQRT, KR, SCALE) are Juicer's; the
  bytes are not. Blocks are compressed on several threads and the output does
  not depend on the number of threads.

```cpp
#include <hicfilecpp/hicfilecpp.hpp>

hicfilecpp::HiCFile hic("matrix.hic");
auto mzd = hic.getMatrixZoomData("chr1", "chr1", "observed", "KR", "BP", 10000);
auto records = mzd.getRecords(0, 5000000, 0, 5000000);      // binX, binY, counts
auto dense = mzd.getRecordsAsMatrix(0, 1000000, 0, 1000000);

hicfilecpp::WriteOptions options;
options.version = 9;
options.chromosomes = {{"chr1", 248956422}, {"chr2", 242193529}};
options.sourceResolution = 5000;                             // the pixels' resolution
options.resolutions = {5000, 10000, 100000, 1000000};
options.threads = 8;
hicfilecpp::writeHicFile("out.hic", options, source);        // source: a hicfilecpp::PixelSource
hicfilecpp::addNorm("other.hic", {"VC", "KR", "SCALE"});
```

Beyond hicstraw, `HiCFile` exposes the version, the attributes, FRAG
resolutions, the expected-value and normalization vector indexes, and
`MatrixZoomData::forEachBlock` decodes a matrix block by block on several
threads for bulk reading.

`docs/API_MAPPING.md` maps every hicstraw call and every Juicer tools option
to its C++ form; `docs/DEVIATIONS.md` lists the reproduced behaviour and the
deliberate deviations; `docs/PROVENANCE.md` names the implementations this
library follows and how the test files were made.

## Build and install

Dependencies: a C++20 compiler, CMake 3.21 or newer, and zlib. doctest is
fetched for the tests.

```sh
cmake -S . -B build -DCMAKE_PREFIX_PATH=/path/to/zlib/prefix
cmake --build build -j
ctest --test-dir build --output-on-failure
cmake --install build --prefix /opt/hicfilecpp
```

Downstream projects use either

```cmake
find_package(hicfilecpp 0.2 REQUIRED)
target_link_libraries(app PRIVATE hicfilecpp::hicfilecpp)
```

or `FetchContent_Declare(hicfilecpp ...)` with `FetchContent_MakeAvailable`.

## Equivalence harness

`harness/run.py` runs the reading cases through hicstraw and hicfilecpp and
compares the results: records and matrices bit for bit, vectors to three
significant digits or better, with CPU time and peak RSS gated against
hicstraw. The writing cases run Juicer tools pre on the same contacts, require
the source pixels back through hicstraw, hic2cool and Juicer tools dump,
compare pixels, expected values and normalization vectors with Juicer's, check
that the output is identical on one thread, on several and when repeated, and
gate CPU time against Juicer. hicstraw runs in its own Python environment:

```sh
python harness/run.py --driver build/harness/hicfilecpp-harness \
    --oracle-python /path/to/env/bin/python \
    --hicx-data /path/to/HiCExplorer/hicexplorer/test/test_data \
    --java /path/to/java --juicer8 juicer_tools_1.22.01.jar \
    --juicer9 juicer_tools.2.20.00.jar \
    --hic2cool-python /path/to/env-with-hic2cool/bin/python --out report
python harness/mutate.py --driver build/harness/hicfilecpp-harness \
    --oracle-python /path/to/env/bin/python --hicx-data ... \
    --straw-source /path/to/straw/pybind11_python --scratch scratch
```

`mutate.py` proves the reading cases can fail: it builds hicstraw with one
deliberate fault at a time and requires the targeted cases to fail.
