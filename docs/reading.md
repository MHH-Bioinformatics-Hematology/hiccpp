# Reading .hic files

`HiCFile` opens a file of version 6 to 9 and answers the queries hicstraw
answers, with the same names, arguments and values, including hicstraw's float
arithmetic. Beyond hicstraw it exposes the file's metadata, its vectors and its
blocks for bulk reading.

```cpp
#include <hiccpp/hiccpp.hpp>

hiccpp::HiCFile hic("matrix.hic");
for (const auto& chromosome : hic.getChromosomes()) {
    // name, index, length
}

auto mzd = hic.getMatrixZoomData("chr1", "chr2", "observed", "NONE", "BP", 10000);
auto records = mzd.getRecords(0, 5000000, 0, 5000000);   // binX, binY, counts
auto dense = mzd.getRecordsAsMatrix(0, 1000000, 0, 1000000);
```

A query that finds nothing reports it rather than throwing, as hicstraw does:
`found()` is false and `message()` carries the text hicstraw prints.

## Reading a whole matrix

`forEachBlock` decodes every block of a zoom level on several threads and hands
each to the visitor in block order on the calling thread, which keeps memory to
roughly one block per thread.

```cpp
mzd.forEachBlock([](const hiccpp::BlockIndexEntry& entry,
                    std::vector<hiccpp::ContactRecord>& records) {
    // records hold bins and raw counts, in block order
}, 8);
```

## Fragment resolutions

`getFragResolutions`, `fragmentSiteCounts` and `fragmentSites` report the
restriction fragment maps of a file that has them; pass the unit `"FRAG"` to
`getMatrixZoomData` to read those matrices.
