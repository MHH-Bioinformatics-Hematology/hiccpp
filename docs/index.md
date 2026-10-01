# hiccpp

A C++20 library that reads the Juicer `.hic` format, versions 6 to 9, and writes
versions 8 and 9, without Python or Java.

Reading reproduces hicstraw 1.3.1 record for record. Writing reproduces Juicer
tools `pre` and `addNorm`: version 8 files follow release 1.22.01, version 9
files release 2.20.00, including the expected values, the VC, VC_SQRT, KR and
SCALE normalization vectors, and the restriction fragment resolutions of
`pre -f`. Blocks are compressed on several threads, and the output does not
depend on the number of threads.

```cpp
#include <iostream>
#include <vector>
#include <hiccpp/hiccpp.hpp>

int main() {
    const hiccpp::HiCFile hic("matrix.hic");
    hiccpp::MatrixZoomData mzd =
        hic.getMatrixZoomData("chr1", "chr1", "observed", "KR", "BP", 10000);
    const std::vector<hiccpp::ContactRecord> records = mzd.getRecords(0, 5000000, 0, 5000000);
    std::cout << records.size() << " records\n";
    return 0;
}
```

The programs under [`examples/`](https://github.com/MHH-Bioinformatics-Hematology/hiccpp/tree/main/examples)
are complete and build with the library: `read_hic` reads a file and
`write_hic` writes one. The pages that follow show them in full.

What the library does not do is listed in [the deviations page](DEVIATIONS.md): writing versions 6
and 7 is refused, since no obtainable Juicer tools release writes them to
validate against.
