# Reading .hic files

`HiCFile` opens a file of version 6 to 9 and answers the queries hicstraw
answers, with the same names, arguments and values, including hicstraw's float
arithmetic. Beyond hicstraw it exposes the file's metadata, its vectors and its
blocks for bulk reading.

`version`, `getGenomeID`, `getResolutions`, `getChromosomes` and
`getNormalizationTypes` report what the file holds; a `Chromosome` carries its
name, its index and its length. `getMatrixZoomData` opens one matrix, from which
`getRecords` returns the records as bins and counts and `getRecordsAsMatrix`
returns a dense block.

```cpp title="examples/read_hic.cpp"
--8<-- "examples/read_hic.cpp"
```

```
$ read_hic GM12878_combined_30.chr21_chr22.v6.hic 250000
version 6, genome hg19
resolutions: 2500000 1000000 500000 250000
2 chromosomes, first 21 of 48129895 bp
10319 records at 250000 bp, first 9250000 9250000 278
blocks hold 10319 records, 2.75158e+07 contacts
GW_KR vector of 193 values
```


A query that finds nothing reports it rather than throwing, as hicstraw does:
`found()` is false and `message()` carries the text hicstraw prints.

## Reading a whole matrix

`forEachBlock` decodes every block of a zoom level on several threads and hands
each to the visitor in block order on the calling thread, which keeps memory to
roughly one block per thread.

The visitor receives the index entry of the block and its records, which hold
bins and raw counts. `read_hic.cpp` above sums a whole zoom level that way.

## Fragment resolutions

`getFragResolutions`, `fragmentSiteCounts` and `fragmentSites` report the
restriction fragment maps of a file that has them; pass the unit `"FRAG"` to
`getMatrixZoomData` to read those matrices.
