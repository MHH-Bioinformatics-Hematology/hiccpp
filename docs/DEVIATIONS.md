# Reproduced behaviour and deliberate deviations

hicfilecpp reproduces what hicstraw 1.3.1 returns and what Juicer tools pre and
addNorm write, including behaviour that looks surprising. Each reproduced item
is exercised by the harness cases named in brackets.

## Reading: hicstraw behaviour reproduced as is

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

## Reading: deliberate deviations

1. Versions below 6 and above 9 raise `HicError("Version N is not supported:
   hicfilecpp reads .hic versions 6 to 9")`. hicstraw refuses versions below
   6 with a message and reads a version above 9 as if it were 9. (Unit test
   "versions below 6 and above 9 are refused".)
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

## Writing: Juicer tools behaviour reproduced as is

Every write case runs Juicer tools pre on the same contacts and compares the
files through hicstraw [W.*]. The fragment items 8 to 11 are covered by the
unit test [F] instead, which writes a file from the FRAG pixels of the Juicer
written fragment reference and compares every FRAG record, expected value and
normalization vector with it.

1. A zoom level's header holds the sum of its contacts, counting
   off-diagonal intra-chromosomal contacts twice, and 0 for the occupied cell
   count and both percentiles: Juicer writes the header before the blocks
   that would fill them in. hicstraw's inter-chromosomal "oe" divides by that
   sum.
2. Normalization vectors hold as many entries as the matrix's grid axis has
   bins, which Juicer takes from the zoom level's block layout: block bin count
   times block column count. With a single block column that is
   `length / binSize + 2`. [W.GSM2644945_chr1.*, where chr1 needs three block
   columns at 100 kb]
3. VC is the row sum, VC_SQRT its square root, KR the Knight-Ruiz balancing
   with Juicer's row tossing (percentile thresholds of 1, 2, 3, 4 and 10
   percent of the non-zero row sums, at most six attempts, the fifth always
   discarded) and SCALE Juicer's FinalScale; every vector is scaled so the
   normalized matrix keeps the observed sum. A KR that never converges is
   written as all-NaN, as Juicer writes it.
4. Version 8 files use the double arithmetic of 1.22.01 and commons-math 2
   percentiles; version 9 files use the float vectors of 2.20.00 and the
   commons-math3 LEGACY percentile.
5. Expected values use Juicer's smoothing window of at least 400 contacts and
   per-chromosome scale factors, from the raw contacts for "NONE" and from the
   normalized records for each normalization.
6. The whole-genome "All" matrix bins every pixel at the kilobase position of
   its bin start, with a bin size of the genome length in kilobases divided by
   500.
7. Version 9 numbers intra-chromosomal blocks by their log2 distance from
   the diagonal and chooses 16 or 32 bit row and column indexes per block.
8. A fragment zoom takes the chromosome's site count where a base pair zoom
    takes its length, and its blocks are laid out without the version 9 column
    widening of fine resolutions: MatrixPP calls `getNumColumnsFromNumBins`
    with a cutoff of 0 there, so a fragment zoom always has
    `siteCount / binSize / 1000 + 1` block columns. The version 9 numbering of
    intra-chromosomal blocks by their distance from the diagonal does apply to
    fragment zooms. [F]
9. The fragment count a chromosome contributes to the expected values is its
    site count plus one in the raw section and its site count in the
    normalized one: Preprocessor builds its `fragmentCountMap` from
    `sites.length + 1`, while the normalization step reads the file back and
    DatasetReaderV2 puts `sites.length` in the map, so the two sections of one
    file disagree by one fragment. hicfilecpp reproduces both. [F]
10. A normalization vector of a fragment zoom holds `siteCount / binSize + 1`
    entries, the bin count of Juicer's HiCFragmentAxis, and not the
    `blockBinCount` times `blockColumnCount` of a base pair zoom, which is one
    more with a single block column. [F]
11. Juicer tools 2.20.00 writes no KR vectors unless `-k` asks for them, so a
    version 9 file written with its defaults holds VC, VC_SQRT and SCALE only,
    at the fragment resolutions as at the base pair ones. [F]

## Writing: deliberate deviations

1. The bytes differ from Juicer's: matrix headers follow their blocks, version
   8 blocks choose between the list-of-rows and dense layouts the way 1.22.01
   does but compression may differ, and the header carries the "software"
   attribute ("hicfilecpp <version>" unless set) and the caller's attributes,
   not Juicer's `hicFileScalingFactor`, `nviIndex` and `nviLength`.
2. Genome-wide and inter-chromosomal normalizations and the filters and
   statistics options of pre are not written.
3. Juicer skips a chromosome's normalization when its Java heap looks too
   small (records times 1000 at least the maximum heap); hicfilecpp always
   computes it. Juicer 2.20.00 spills blocks to temporary files and may then
   sum repeated non-integer counts in another order; such counts agree to
   floating point rounding.
4. Normalization vectors are kept in memory until the footer is written,
   eight bytes per bin per normalization.
5. addNorm refuses files with FRAG resolutions, where Juicer would also
   normalize the fragment maps; writeHicFile does write them. Like Juicer,
   addNorm replaces the normalizations a file already holds.
6. Non-finite counts, pixels outside their chromosome, a chromosome named
   "All", repeated chromosome names and lengths above 2147483647 raise
   `HicError`. So do fragment resolutions without a site list per chromosome,
   unsorted sites, and a fragment pixel past a chromosome's site count. Juicer
   checks none of these.
7. The fragment write cases are not in the equivalence harness: its writer
   cases derive one contact list from a base pair source, and a fragment map
   needs a sites file and the fragment columns of the contact format, which
   the case schema has no room for. [F] is the unit test "fragment resolutions
   reproduce the FRAG records, vectors and expected values of Juicer tools"
   instead, cross-checked outside the build by reading both files with
   hicstraw 1.3.1 and with Juicer tools dump.
