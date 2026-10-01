#ifndef HICCPP_WRITER_HPP
#define HICCPP_WRITER_HPP

// Writing .hic files, versions 8 and 9. writeHicFile corresponds to Juicer
// tools "pre" (followed by its normalization step) and addNorm to "addNorm":
// version 8 files follow the layout and arithmetic of Juicer tools 1.22.01,
// version 9 files those of 2.20.00. Expected values and normalization vectors
// agree with Juicer's to within floating point rounding; blocks may be
// encoded differently, so the bytes are not Juicer's. The output depends only
// on the input and the options, never on the number of threads.

#include <cstdint>
#include <functional>
#include <string>
#include <utility>
#include <vector>

#include "hiccpp/errors.hpp"

namespace hiccpp {

// A pixel handed to the writer: bins of the requested resolution, local to
// their chromosomes, and the count.
struct Pixel {
    int32_t bin1 = 0;
    int32_t bin2 = 0;
    float count = 0.0f;
};

class PixelSource {
public:
    virtual ~PixelSource() = default;

    // Hands every pixel of the chromosome pair (chr1, chr2) at `resolution`
    // to `consume`, in any number of batches. chr1 <= chr2 are 0-based
    // indexes into WriteOptions::chromosomes; bin1 lies on chr1 and bin2 on
    // chr2. An intra-chromosomal pixel is given once, in either triangle;
    // repeated pixels are summed. The writer asks for one pair and resolution
    // at a time, in pair order, and may ask for a pair once per resolution.
    virtual void pixels(int32_t resolution, int32_t chr1, int32_t chr2,
                        const std::function<void(const Pixel* pixels, size_t count)>& consume) = 0;

    // The same for fragment binned pixels, whose bins count restriction
    // fragments instead of base pairs: bin1 lies in 0 to the site count of
    // chr1 divided by `resolution`, as Juicer tools pre -f bins the fragment
    // numbers of its contacts. Only called when WriteOptions::fragResolutions
    // is not empty; the default hands over nothing, so a source written for
    // base pair resolutions keeps compiling.
    virtual void fragPixels(int32_t resolution, int32_t chr1, int32_t chr2,
                            const std::function<void(const Pixel* pixels, size_t count)>& consume) {
        (void)resolution;
        (void)chr1;
        (void)chr2;
        (void)consume;
    }
};

struct WriteOptions {
    // 8 or 9. Versions 6 and 7 are refused with their own message: no Juicer
    // tools release that writes them can be obtained to validate against
    // (docs/PROVENANCE.md).
    int32_t version = 9;
    std::string genomeId;
    // Name and length in base pairs, in file order. The whole-genome
    // pseudo-chromosome "All" is added in front, as Juicer tools does.
    std::vector<std::pair<std::string, int64_t>> chromosomes;
    // Base pair resolutions, written from the coarsest to the finest.
    std::vector<int32_t> resolutions;
    // Fragment resolutions, as Juicer tools pre takes them from the "f"
    // suffixed entries of -r (100f, 20f) when -f names a sites file. They are
    // written after the base pair resolutions, from the coarsest to the
    // finest, and need `fragmentSites`. The whole-genome "All" matrix stays a
    // base pair matrix, as it is in Juicer's files.
    std::vector<int32_t> fragResolutions;
    // The restriction sites of each chromosome, in ascending order: one entry
    // per entry of `chromosomes`, the lines of the sites file of Juicer tools
    // pre -f. Needed when `fragResolutions` is not empty, and empty for a
    // chromosome the sites file leaves out. A chromosome's fragment bin count
    // at resolution N is its site count divided by N plus one, the arithmetic
    // of Juicer's FragmentCalculation.getNumberFragments.
    std::vector<std::vector<int32_t>> fragmentSites;
    // The resolution of the pixels the source hands over. Each pixel counts
    // as one contact at the start of its bins, so every resolution, which
    // must be a multiple of this one, holds what Juicer tools pre makes of
    // the same contacts.
    int32_t sourceResolution = 0;
    // The fragment resolution of the pixels `PixelSource::fragPixels` hands
    // over, which every entry of `fragResolutions` must be a multiple of.
    int32_t sourceFragResolution = 0;
    // Ask the source for pixels at every resolution instead of binning the
    // source resolution. The whole-genome matrix then uses the finest base
    // pair one. This covers the fragment resolutions as well, so
    // `sourceFragResolution` is then unused.
    bool sourceProvidesEveryResolution = false;
    // Computed as Juicer tools addNorm does; any of VC, VC_SQRT, KR, SCALE.
    std::vector<std::string> normalizations{"VC", "VC_SQRT", "KR", "SCALE"};
    // Normalization vectors the caller supplies instead of having them
    // computed, such as the vectors of another .hic file kept in a cool file.
    // A name is any normalization type label (KR, GW_KR, a custom one) that
    // `normalizations` does not also list. For every chromosome with
    // intra-chromosomal pixels at a resolution, `normVector` is asked for the
    // vector of each name, chrIndex being the 0-based index into
    // `chromosomes`; an empty vector writes none for that chromosome. The
    // values are stored as given, without Juicer's scaling to the matrix sum,
    // and the normalized expected values are computed from them. A vector
    // shorter than the matrix grid is padded as Juicer pads its own: with 0
    // for VC types, NaN otherwise; a longer one is cut.
    std::vector<std::string> providedNormalizations;
    std::function<std::vector<double>(const std::string& name, int32_t chrIndex, int32_t resolution)> normVector;
    // The same for the fragment resolutions, asked for each name of
    // `providedNormalizations` at each entry of `fragResolutions`. Left unset,
    // no provided vector is written for the fragment resolutions.
    std::function<std::vector<double>(const std::string& name, int32_t chrIndex, int32_t resolution)> fragNormVector;
    // Header attributes written after "software", in order.
    std::vector<std::pair<std::string, std::string>> attributes;
    // The "software" attribute; empty means "hiccpp <version>".
    std::string software;
    // Threads for block encoding and compression.
    int threads = 1;
    // zlib compression level; -1 is zlib's default, which Juicer tools uses.
    int compressionLevel = -1;
};

// Writes a .hic file from the source's pixels.
void writeHicFile(const std::string& path, const WriteOptions& options, PixelSource& source);

// Computes the normalization vectors and normalized expected values of an
// existing version 8 or 9 file with base pair resolutions only, replacing any
// it holds, as Juicer tools addNorm does.
void addNorm(const std::string& path, const std::vector<std::string>& normalizations, int threads = 1);

}  // namespace hiccpp

#endif  // HICCPP_WRITER_HPP
