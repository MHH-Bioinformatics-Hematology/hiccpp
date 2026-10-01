#include <doctest/doctest.h>

#include <algorithm>
#include <cstdio>
#include <filesystem>
#include <fstream>
#include <string>
#include <tuple>
#include <vector>

#include <hiccpp/hiccpp.hpp>

namespace {

const std::string kData = HICCPP_TEST_DATA;
const std::string kV8 = kData + "/SRR1791297_30.juicer_tools_1.22.01.v8.hic";
const std::string kFrag = kData + "/SRR1791297_30.juicer_tools_1.22.01.frag.v8.hic";
const std::string kV9 = kData + "/SRR1791297_30.juicer_tools_2.20.00.v9.hic";

std::string scratch(const std::string& name) {
    std::filesystem::create_directories(HICCPP_TEST_SCRATCH);
    return std::string(HICCPP_TEST_SCRATCH) + "/" + name;
}

std::vector<std::tuple<int32_t, int32_t, float>> sorted(const std::vector<hiccpp::ContactRecord>& records) {
    std::vector<std::tuple<int32_t, int32_t, float>> out;
    for (const auto& r : records) {
        out.emplace_back(r.binX, r.binY, r.counts);
    }
    std::sort(out.begin(), out.end());
    return out;
}

}  // namespace

TEST_CASE("headers of version 8 and version 9 files") {
    const hiccpp::HiCFile v8(kV8);
    CHECK(v8.version() == 8);
    CHECK(v8.getGenomeID() == "sacCer3.chrom.sizes");
    CHECK(v8.getResolutions() == std::vector<int32_t>{1000000, 250000, 50000, 10000});
    CHECK(v8.getChromosomes().size() == 17);
    CHECK(v8.getChromosomes()[0].name == "All");
    CHECK(v8.normVectorIndexHeader() == std::pair<int64_t, int64_t>{0, 0});

    const hiccpp::HiCFile v9(kV9);
    CHECK(v9.version() == 9);
    CHECK(v9.normVectorIndexHeader().first > 0);
    CHECK(v9.normVectorIndexHeader().second > 0);

    const hiccpp::HiCFile frag(kFrag);
    CHECK(frag.getFragResolutions() == std::vector<int32_t>{100, 20});
    CHECK(frag.fragmentSiteCounts().size() == 17);
}

TEST_CASE("records of a query equal the decoded blocks of the matrix") {
    for (const auto& path : {kV8, kV9}) {
        const hiccpp::HiCFile file(path);
        const auto chromosomes = file.getChromosomes();
        for (const auto& [a, b] : {std::pair{1, 1}, std::pair{2, 5}}) {
            const auto mzd = file.getMatrixZoomData(chromosomes[a].name, chromosomes[b].name, "observed", "NONE",
                                                    "BP", 10000);
            REQUIRE(mzd.found());
            const auto records =
                mzd.getRecords(0, chromosomes[a].length, 0, chromosomes[b].length);
            std::vector<hiccpp::ContactRecord> blocks;
            for (const auto& entry : mzd.blockIndex()) {
                for (auto r : mzd.readBlock(entry)) {
                    r.binX *= 10000;
                    r.binY *= 10000;
                    blocks.push_back(r);
                }
            }
            CHECK(!records.empty());
            CHECK(sorted(records) == sorted(blocks));
            CHECK(mzd.getNumberOfTotalRecords() == static_cast<int64_t>(blocks.size()));
        }
    }
}

TEST_CASE("forEachBlock gives the same blocks in the same order on any number of threads") {
    const hiccpp::HiCFile file(kV9);
    const auto mzd = file.getMatrixZoomData("NC_001136.10", "NC_001136.10", "observed", "NONE", "BP", 10000);
    auto collect = [&](int threads) {
        std::vector<std::tuple<int32_t, int32_t, int32_t, float>> out;
        mzd.forEachBlock(
            [&](const hiccpp::BlockIndexEntry& entry, std::vector<hiccpp::ContactRecord>& records) {
                for (const auto& r : records) {
                    out.emplace_back(entry.number, r.binX, r.binY, r.counts);
                }
            },
            threads);
        return out;
    };
    const auto one = collect(1);
    CHECK(!one.empty());
    CHECK(one == collect(4));
    CHECK(one == collect(32));
}

TEST_CASE("versions below 6 and above 9 are refused") {
    std::ifstream in(kV8, std::ios::binary);
    std::vector<char> bytes((std::istreambuf_iterator<char>(in)), std::istreambuf_iterator<char>());
    for (const int version : {-1, 0, 5, 10, 11}) {
        std::vector<char> copy = bytes;
        std::copy_n(reinterpret_cast<const char*>(&version), 4, copy.begin() + 4);
        const std::string path = scratch("version" + std::to_string(version) + ".hic");
        std::ofstream(path, std::ios::binary).write(copy.data(), static_cast<std::streamsize>(copy.size()));
        CHECK_THROWS_WITH_AS(hiccpp::HiCFile{path},
                             ("Version " + std::to_string(version) +
                              " is not supported: hiccpp reads .hic versions 6 to 9")
                                 .c_str(),
                             hiccpp::HicError);
    }
}

TEST_CASE("a version 8 file relabelled as version 7 reads the same records") {
    // Versions 7 and 8 share every layout hicstraw reads, so only the header
    // version differs. (Version 6 changes the block records, which the harness
    // checks on tests/data/GM12878_combined_30.chr21_chr22.v6.hic.)
    std::ifstream in(kV8, std::ios::binary);
    std::vector<char> bytes((std::istreambuf_iterator<char>(in)), std::istreambuf_iterator<char>());
    const int version = 7;
    std::copy_n(reinterpret_cast<const char*>(&version), 4, bytes.begin() + 4);
    const std::string path = scratch("relabelled_v7.hic");
    std::ofstream(path, std::ios::binary).write(bytes.data(), static_cast<std::streamsize>(bytes.size()));
    const hiccpp::HiCFile v7(path);
    const hiccpp::HiCFile v8(kV8);
    CHECK(v7.version() == 7);
    const auto records = [](const hiccpp::HiCFile& file) {
        const auto mzd = file.getMatrixZoomData("NC_001133.9", "NC_001133.9", "observed", "KR", "BP", 10000);
        std::vector<std::tuple<int32_t, int32_t, float>> out;
        for (const auto& r : mzd.getRecords(0, 230218, 0, 230218)) {
            out.emplace_back(r.binX, r.binY, r.counts);
        }
        return out;
    };
    CHECK(!records(v8).empty());
    CHECK(records(v7) == records(v8));
}

TEST_CASE("invalid files, chromosomes, zoom levels and vectors raise HicError") {
    const std::string text = scratch("not_hic.txt");
    std::ofstream(text) << "this is not a hic file\n";
    CHECK_THROWS_AS(hiccpp::HiCFile{text}, hiccpp::HicError);
    CHECK_THROWS_AS(hiccpp::HiCFile{scratch("does_not_exist.hic")}, hiccpp::HicError);

    const hiccpp::HiCFile file(kV8);
    CHECK_THROWS_AS(file.getMatrixZoomData("chrNone", "chrNone", "observed", "NONE", "BP", 10000),
                    hiccpp::HicError);
    CHECK_THROWS_AS(file.getMatrixZoomData("NC_001133.9", "NC_001133.9", "observed", "NONE", "BP", 12345),
                    hiccpp::HicError);
    CHECK_THROWS_AS(file.getMatrixZoomData("NC_001133.9", "NC_001133.9", "observed", "NO_SUCH_NORM", "BP", 10000),
                    hiccpp::HicError);
}

TEST_CASE("metadata accessors") {
    const hiccpp::HiCFile file(kV8);
    const auto types = file.getNormalizationTypes();
    CHECK(std::find(types.begin(), types.end(), "KR") != types.end());
    CHECK(std::find(types.begin(), types.end(), "SCALE") != types.end());
    const auto vector = file.readNormVector("VC", 1, "BP", 10000);
    REQUIRE(vector.has_value());
    // Juicer's vectors carry one entry past the last bin.
    CHECK(vector->size() >= static_cast<size_t>(file.getChromosomes()[1].length / 10000 + 1));
    const auto expected = file.readExpectedValues({"NONE", "BP", 10000});
    REQUIRE(expected.has_value());
    CHECK(!expected->values.empty());
    CHECK(!expected->normalizationFactors.empty());
    CHECK(file.hasMatrix(1, 2));
    CHECK(!file.matrixZoomHeaders(1, 1).empty());
}
