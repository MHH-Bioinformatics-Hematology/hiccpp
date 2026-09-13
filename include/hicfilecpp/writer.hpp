#ifndef HICFILECPP_WRITER_HPP
#define HICFILECPP_WRITER_HPP

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

#include "hicfilecpp/errors.hpp"

namespace hicfilecpp {

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
};

struct WriteOptions {
    // 8 or 9.
    int32_t version = 9;
    std::string genomeId;
    // Name and length in base pairs, in file order. The whole-genome
    // pseudo-chromosome "All" is added in front, as Juicer tools does.
    std::vector<std::pair<std::string, int64_t>> chromosomes;
    // Base pair resolutions, written from the coarsest to the finest.
    std::vector<int32_t> resolutions;
    // The resolution of the pixels the source hands over. Each pixel counts
    // as one contact at the start of its bins, so every resolution, which
    // must be a multiple of this one, holds what Juicer tools pre makes of
    // the same contacts.
    int32_t sourceResolution = 0;
    // Ask the source for pixels at every resolution instead of binning the
    // source resolution. The whole-genome matrix then uses the finest one.
    bool sourceProvidesEveryResolution = false;
    // Computed as Juicer tools addNorm does; any of VC, VC_SQRT, KR, SCALE.
    std::vector<std::string> normalizations{"VC", "VC_SQRT", "KR", "SCALE"};
    // Header attributes written after "software", in order.
    std::vector<std::pair<std::string, std::string>> attributes;
    // The "software" attribute; empty means "hicfilecpp <version>".
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

}  // namespace hicfilecpp

#endif  // HICFILECPP_WRITER_HPP
