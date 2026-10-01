// A downstream program: opens a .hic file and prints its resolutions and the
// number of records of the first chromosome at the coarsest resolution.
#include <hiccpp/hiccpp.hpp>

#include <cstdio>

int main(int argc, char** argv) {
    if (argc != 2) {
        std::fprintf(stderr, "usage: consumer FILE.hic\n");
        return 2;
    }
    const hiccpp::HiCFile file(argv[1]);
    const auto chromosomes = file.getChromosomes();
    const int32_t resolution = file.getResolutions().front();
    const auto records = file.getMatrixZoomData(chromosomes[1].name, chromosomes[1].name, "observed", "NONE", "BP",
                                                resolution)
                             .getRecords(0, chromosomes[1].length, 0, chromosomes[1].length);
    std::printf("hiccpp %s: version %d, %zu records at %d bp\n", hiccpp::kVersion, file.version(),
                records.size(), resolution);
    return records.empty() ? 1 : 0;
}
