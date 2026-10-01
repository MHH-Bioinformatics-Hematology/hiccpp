// Writing a .hic file: a PixelSource over a matrix held in memory, written at
// several resolutions with Juicer's normalization vectors, then read back.
//
//     write_hic out.hic 9

#include <cstdint>
#include <functional>
#include <iostream>
#include <map>
#include <string>
#include <utility>
#include <vector>

#include <hiccpp/hiccpp.hpp>

namespace {

constexpr std::int32_t kResolution = 10000;   // the resolution of the pixels below

// Two small chromosomes and the pixels of each pair, as bins of kResolution.
// A real program would read these from its own data structure instead.
struct Matrix {
    std::vector<std::pair<std::string, std::int64_t>> chromosomes{{"chr1", 400000}, {"chr2", 300000}};
    // (chromosome pair) -> pixels, the upper triangle for a pair with itself
    std::map<std::pair<std::int32_t, std::int32_t>, std::vector<hiccpp::Pixel>> pixels{
        {{0, 0}, {{0, 0, 120.0F}, {0, 1, 45.0F}, {1, 1, 90.0F}, {1, 2, 30.0F}, {2, 2, 75.0F}}},
        {{0, 1}, {{0, 0, 12.0F}, {1, 2, 8.0F}, {3, 1, 5.0F}}},
        {{1, 1}, {{0, 0, 60.0F}, {0, 1, 25.0F}, {1, 1, 80.0F}, {2, 2, 40.0F}}},
    };
};

// The writer asks the source for the pixels of one chromosome pair at a time.
class MatrixSource : public hiccpp::PixelSource {
  public:
    explicit MatrixSource(const Matrix& matrix) : matrix_(matrix) {}

    void pixels(std::int32_t resolution, std::int32_t chr1, std::int32_t chr2,
                const std::function<void(const hiccpp::Pixel*, std::size_t)>& consume) override {
        if (resolution != kResolution) {   // the writer only asks for sourceResolution
            return;
        }
        const auto found = matrix_.pixels.find({chr1, chr2});
        if (found == matrix_.pixels.end()) {
            return;                        // no contacts for this pair
        }
        consume(found->second.data(), found->second.size());
    }

  private:
    const Matrix& matrix_;
};

}  // namespace

int main(int argc, char** argv) {
    const std::string path = argc > 1 ? argv[1] : "out.hic";
    const std::int32_t version = argc > 2 ? std::stoi(argv[2]) : 9;

    const Matrix matrix;

    hiccpp::WriteOptions options;
    options.version = version;             // 8 as Juicer tools 1.22 writes, 9 as Juicer 2 does
    options.genomeId = "toy";
    options.chromosomes = matrix.chromosomes;
    options.sourceResolution = kResolution;
    options.resolutions = {kResolution, 50000};   // coarser levels are binned from the source
    options.normalizations = {"VC", "VC_SQRT", "SCALE"};
    options.threads = 4;                   // block compression; the output does not depend on it

    MatrixSource source(matrix);
    hiccpp::writeHicFile(path, options, source);

    const hiccpp::HiCFile written(path);
    std::cout << "wrote " << path << ": version " << written.version() << ", resolutions";
    for (const std::int32_t r : written.getResolutions()) {
        std::cout << ' ' << r;
    }
    std::cout << '\n';
    hiccpp::MatrixZoomData mzd =
        written.getMatrixZoomData("chr1", "chr1", "observed", "NONE", "BP", kResolution);
    std::cout << "chr1 holds " << mzd.getRecords(0, 400000, 0, 400000).size() << " records at "
              << kResolution << " bp\n";
    for (const std::string& norm : written.getNormalizationTypes()) {
        std::cout << "vector " << norm << '\n';
    }
    return 0;
}
