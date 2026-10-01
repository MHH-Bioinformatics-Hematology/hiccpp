# Writing .hic files

`writeHicFile` writes a version 8 or 9 file from a `PixelSource`, as Juicer
tools `pre` does, and computes the normalization vectors `addNorm` computes.
The source hands over the pixels of one chromosome pair at a time, so a whole
matrix never has to be in memory.

`write_hic` on [the examples page](examples.md#writing-a-file) is a complete
program that implements a `PixelSource` over a matrix held in memory and writes
it at two resolutions with the normalization vectors. A source that reads its
pixels from elsewhere, a database or another file, hands them over in as many
batches as it likes; the writer never holds more than one chromosome pair.

Versions 6 and 7 are refused with their own message: no Juicer tools release
that writes them can be obtained to validate against.

## Normalization vectors

By default the writer computes VC, VC_SQRT, KR and SCALE as Juicer tools does.
A caller that already has vectors, for example from another `.hic` file kept in
a cool file, supplies them instead through `providedNormalizations` and
`normVector`. Those are stored as given, under any label, and the normalized
expected values are computed from them.

`addNorm` computes and replaces the vectors of an existing file with base pair
resolutions, as Juicer tools `addNorm` does.

## Fragment resolutions

`fragResolutions` and `fragmentSites` add the restriction fragment matrices
`pre -f` writes, with the pixels coming from `PixelSource::fragPixels` at
`sourceFragResolution`.
