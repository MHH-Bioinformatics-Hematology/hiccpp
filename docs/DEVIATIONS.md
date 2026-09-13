# Reproduced hicstraw behaviour and deliberate deviations

hicfilecpp reproduces what hicstraw 1.3.1 returns, including behaviour that
looks surprising. Each reproduced item is exercised by the harness cases named
in brackets, which run the same queries through hicstraw and hicfilecpp.

## hicstraw behaviour reproduced as is

1. `getRecords` reports `binX` and `binY` as genomic start positions (bin
   times resolution), not bins; the region bounds are inclusive on both ends.
   [R.SRR.sub.BP.10000.observed]
2. An intra-chromosomal query also returns the records of the mirrored region:
   a record matches when (x, y) or (y, x) lies inside the region.
   [R.J9.sub.BP.10000.observed, R.SRR.sub.BP.10000.observed]
3. The chromosome pair is ordered by chromosome index: `getMatrixZoomData`
   with the higher index first reads the same matrix, and the first region
   argument still applies to the lower-index chromosome.
   [R.SRR.sub.BP.10000.observed]
4. Normalized, observed/expected and expected values whose result is NaN or
   infinite are dropped, so normalized queries return fewer records than
   observed ones. [R.*.observed with norms, R.*.oe]
5. Intra-chromosomal "oe" divides by the expected value at the binned
   distance, clamped to the last entry; the expected vector is divided by the
   chromosome's normalization factor. Inter-chromosomal "oe" divides by
   `(sum / numBins1) / numBins2`, computed in float32 from the zoom header's
   float sum, with `numBins = length / resolution` (no rounding up).
   [R.*.oe, V.*]
6. A chromosome pair the file lacks, and an intra-chromosomal "oe" or
   "expected" query without expected values, return no records. hicstraw
   prints a message; hicfilecpp keeps it in `message()` and reports
   `found() == false`. [R.SRR.BP.2500000.oe, V.SRR.BP.2500000]
7. `getRecordsAsMatrix` fills the mirrored cell of intra-chromosomal records
   and returns the 1 by 1 matrix `[[0]]` when there are no records.
   [M.SRR.matrices, M.J9.matrices, M.J8F.matrices]
8. `getNormVector` returns the stored vector with its trailing entry past the
   last bin, and an empty vector for an index that is neither chromosome of
   the pair. [V.*]
9. For FRAG queries, positions are fragment bins times the resolution, while
   `numBins` still divides the chromosome's length in base pairs.
   [R.J8F.FRAG.*, M.J8F.matrices]
10. A repeated chromosome name resolves to its last occurrence.

## Deliberate deviations

1. Versions other than 8 and 9 raise `HicError("Version N is not supported:
   hicfilecpp reads .hic versions 8 and 9")`. hicstraw also reads versions 6
   and 7. (Unit test "versions other than 8 and 9 are refused".)
2. An unknown chromosome name raises `HicError("NAME not found in the
   file.")`. hicstraw silently reads the default chromosome, index 0 with
   length 0.
3. A normalization without vectors for both chromosomes raises `HicError`.
   hicstraw prints a message and then reads a vector from offset 0 of the
   file, which is undefined.
4. A resolution the matrix lacks raises `HicError("Error finding block data:
   ...")`. hicstraw prints the first part and continues with uninitialized
   block dimensions.
5. A file that cannot be opened, or is not a .hic file, raises `HicError`.
   hicstraw terminates the Python process or returns an empty object.
6. A normalization vector shorter than a record's bin yields NaN for that
   record, which is then dropped. hicstraw reads past the vector's end.
7. Genomic positions are computed in 64 bits. hicstraw multiplies bin and
   resolution in 32 bits, which overflows beyond 2^31 base pairs.
8. Blocks inflate to any size. hicstraw reserves ten times the compressed
   size and truncates larger blocks.
9. Remote (http) files are not supported.
