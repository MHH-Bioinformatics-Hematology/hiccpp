# Writing .hic files

`writeHicFile` writes a version 8 or 9 file from a `PixelSource`, as Juicer
tools `pre` does, and computes the normalization vectors `addNorm` computes.
The source hands over the pixels of one chromosome pair at a time, so a whole
matrix never has to be in memory.

```cpp
#include <hiccpp/hiccpp.hpp>

class MySource : public hiccpp::PixelSource {
  public:
    void pixels(std::int32_t resolution, std::int32_t chr1, std::int32_t chr2,
                const std::function<void(const hiccpp::Pixel*, std::size_t)>& consume) override {
        // hand over the pixels of this pair at this resolution, in any number of batches
    }
};

hiccpp::WriteOptions options;
options.version = 9;
options.genomeId = "hg38";
options.chromosomes = {{"chr1", 248956422}, {"chr2", 242193529}};
options.sourceResolution = 5000;
options.resolutions = {5000, 10000, 100000, 1000000};
options.threads = 8;

MySource source;
hiccpp::writeHicFile("out.hic", options, source);
```

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
