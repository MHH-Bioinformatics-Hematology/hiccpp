// The shortest program that opens a .hic file and reads records out of it.
//
//     quickstart matrix.hic chr1 10000

#include <cstdint>
#include <iostream>
#include <string>
#include <vector>

#include <hiccpp/hiccpp.hpp>

int main(int argc, char** argv) {
    if (argc != 4) {
        std::cerr << "usage: quickstart <file.hic> <chromosome> <resolution>\n";
        return 2;
    }
    const std::string chromosome = argv[2];
    const std::int32_t resolution = std::stoi(argv[3]);

    const hiccpp::HiCFile hic(argv[1]);

    // The length of the chromosome, to query the whole of it below.
    std::int64_t length = 0;
    for (const hiccpp::Chromosome& c : hic.getChromosomes()) {
        if (c.name == chromosome) {
            length = c.length;
        }
    }
    if (length == 0) {
        std::cerr << "no chromosome " << chromosome << " in this file\n";
        return 1;
    }

    hiccpp::MatrixZoomData mzd =
        hic.getMatrixZoomData(chromosome, chromosome, "observed", "NONE", "BP", resolution);
    // A query that finds nothing reports it rather than throwing, as hicstraw
    // does: found() is false and message() carries the text hicstraw prints.
    if (!mzd.found()) {
        std::cerr << mzd.message() << '\n';
        return 1;
    }

    // The two ranges are genomic positions, not bins, as in hicstraw.
    const std::vector<hiccpp::ContactRecord> records = mzd.getRecords(0, length, 0, length);
    double total = 0;
    for (const hiccpp::ContactRecord& record : records) {
        total += record.counts;
    }
    std::cout << records.size() << " records on " << chromosome << " at " << resolution
              << " bp, " << total << " contacts\n";
    return 0;
}
