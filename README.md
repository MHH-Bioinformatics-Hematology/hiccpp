# hicfilecpp

A C++20 library that reads and writes the Juicer `.hic` format, versions 8
and 9, without Python or Java. Its reading API reproduces hicstraw 1.3.1:
`HiCFile`, `getMatrixZoomData`, `getRecords`, `getRecordsAsMatrix`,
`getNormVector` and `getExpectedValues` return what hicstraw returns, record
for record.

Status: reading. Writing (Juicer tools `pre` and `addNorm`) follows.

```cpp
#include <hicfilecpp/hicfilecpp.hpp>

hicfilecpp::HiCFile hic("matrix.hic");
auto mzd = hic.getMatrixZoomData("chr1", "chr1", "observed", "KR", "BP", 10000);
auto records = mzd.getRecords(0, 5000000, 0, 5000000);      // binX, binY, counts
auto dense = mzd.getRecordsAsMatrix(0, 1000000, 0, 1000000);
```

Beyond hicstraw, `HiCFile` exposes the version, the attributes, FRAG
resolutions, the expected-value and normalization vector indexes, and
`MatrixZoomData::forEachBlock` decodes a matrix block by block on several
threads for bulk reading.

`docs/API_MAPPING.md` maps every hicstraw call to its C++ form;
`docs/DEVIATIONS.md` lists reproduced hicstraw behaviour and the deliberate
deviations; `docs/PROVENANCE.md` names the implementations this library
follows and how the test files were made.

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
find_package(hicfilecpp 0.1 REQUIRED)
target_link_libraries(app PRIVATE hicfilecpp::hicfilecpp)
```

or `FetchContent_Declare(hicfilecpp ...)` with `FetchContent_MakeAvailable`.

## Equivalence harness

`harness/run.py` runs every case through hicstraw and through hicfilecpp and
compares the results: records and matrices bit for bit, vectors to three
significant digits or better. It also gates CPU time and peak RSS against
hicstraw. hicstraw runs in its own Python environment:

```sh
python harness/run.py --driver build/harness/hicfilecpp-harness \
    --oracle-python /path/to/env/bin/python \
    --hicx-data /path/to/HiCExplorer/hicexplorer/test/test_data --out report
python harness/mutate.py --driver build/harness/hicfilecpp-harness \
    --oracle-python /path/to/env/bin/python --hicx-data ... \
    --straw-source /path/to/straw/pybind11_python --scratch scratch
```

`mutate.py` proves the cases can fail: it builds hicstraw with one deliberate
fault at a time and requires the targeted cases to fail.
