# Reading .hic files

`HiCFile` opens a file of version 6 to 9 and answers the queries hicstraw
answers, with the same names, arguments and values, including hicstraw's float
arithmetic. Beyond hicstraw it exposes the file's metadata, its vectors and its
blocks for bulk reading.

```cpp
hiccpp::HiCFile hic("matrix.hic");
hic.version();
hic.getGenomeID();
hic.getResolutions();
hic.getChromosomes();                                    // name, index, length
hic.getNormalizationTypes();

auto mzd = hic.getMatrixZoomData("chr1", "chr2", "observed", "NONE", "BP", 10000);
auto records = mzd.getRecords(0, 5000000, 0, 5000000);   // binX, binY, counts
auto dense = mzd.getRecordsAsMatrix(0, 1000000, 0, 1000000);
```

`read_hic` on [the examples page](examples.md#reading-a-file) is a complete
program that prints this metadata, reads the records of one chromosome pair,
walks the same matrix block by block, and reads a normalization vector.

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
