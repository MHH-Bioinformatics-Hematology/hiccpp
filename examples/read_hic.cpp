// Reading a .hic file of version 6 to 9: metadata, records of one chromosome
// pair, and every block of a zoom level.
//
//     read_hic matrix.hic 10000

#include <algorithm>
#include <cstdint>
#include <iostream>
#include <string>
#include <vector>

#include <hiccpp/hiccpp.hpp>

int main(int argc, char** argv) {
    if (argc != 3) {
        std::cerr << "usage: read_hic <file.hic> <resolution>\n";
        return 2;
    }
    const std::string path = argv[1];
    const std::int32_t resolution = std::stoi(argv[2]);

    const hiccpp::HiCFile hic(path);
    std::cout << "version " << hic.version() << ", genome " << hic.getGenomeID() << '\n';
    std::cout << "resolutions:";
    for (const std::int32_t r : hic.getResolutions()) {
        std::cout << ' ' << r;
    }
    std::cout << '\n';

    // A resolution the file does not hold has no matrix to read, so it is
    // checked before asking for one.
    const std::vector<std::int32_t> available = hic.getResolutions();
    if (std::find(available.begin(), available.end(), resolution) == available.end()) {
        std::cerr << "the file holds no " << resolution << " bp resolution\n";
        return 1;
    }

    // The chromosomes, without the "All" pseudo-chromosome the file begins with.
    std::vector<hiccpp::Chromosome> chromosomes;
    for (const hiccpp::Chromosome& c : hic.getChromosomes()) {
        if (c.index > 0) {
            chromosomes.push_back(c);
        }
    }
    const hiccpp::Chromosome& first = chromosomes.front();
    std::cout << chromosomes.size() << " chromosomes, first " << first.name << " of " << first.length
              << " bp\n";

    // The raw records of one intra-chromosomal matrix. The last two arguments
    // are genomic positions, as in hicstraw.
    hiccpp::MatrixZoomData mzd =
        hic.getMatrixZoomData(first.name, first.name, "observed", "NONE", "BP", resolution);
    if (!mzd.found()) {
        std::cerr << mzd.message() << '\n';
        return 1;
    }
    const std::vector<hiccpp::ContactRecord> records =
        mzd.getRecords(0, first.length, 0, first.length);
    std::cout << records.size() << " records at " << resolution << " bp";
    if (!records.empty()) {
        std::cout << ", first " << records.front().binX << ' ' << records.front().binY << ' '
                  << records.front().counts;
    }
    std::cout << '\n';

    // The same matrix block by block, decoded on four threads, which keeps
    // memory to about one block per thread.
    std::int64_t pixels = 0;
    double total = 0;
    mzd.forEachBlock(
        [&](const hiccpp::BlockIndexEntry&, std::vector<hiccpp::ContactRecord>& block) {
            pixels += static_cast<std::int64_t>(block.size());
            for (const hiccpp::ContactRecord& record : block) {
                total += record.counts;
            }
        },
        4);
    std::cout << "blocks hold " << pixels << " records, " << total << " contacts\n";

    // A normalization vector, when the file carries one.
    for (const std::string& norm : hic.getNormalizationTypes()) {
        if (const auto vector = hic.readNormVector(norm, first.index, "BP", resolution)) {
            std::cout << norm << " vector of " << vector->size() << " values\n";
            break;
        }
    }
    return 0;
}
