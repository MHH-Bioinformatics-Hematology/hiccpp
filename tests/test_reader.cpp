#include <doctest/doctest.h>

#include <algorithm>
#include <cstdio>
#include <filesystem>
#include <fstream>
#include <string>
#include <tuple>
#include <vector>

#include <hicfilecpp/hicfilecpp.hpp>

namespace {

const std::string kData = HICFILECPP_TEST_DATA;
const std::string kV8 = kData + "/SRR1791297_30.juicer_tools_1.22.01.v8.hic";
const std::string kFrag = kData + "/SRR1791297_30.juicer_tools_1.22.01.frag.v8.hic";
const std::string kV9 = kData + "/SRR1791297_30.juicer_tools_2.20.00.v9.hic";

std::string scratch(const std::string& name) {
    std::filesystem::create_directories(HICFILECPP_TEST_SCRATCH);
    return std::string(HICFILECPP_TEST_SCRATCH) + "/" + name;
}

std::vector<std::tuple<int32_t, int32_t, float>> sorted(const std::vector<hicfilecpp::ContactRecord>& records) {
    std::vector<std::tuple<int32_t, int32_t, float>> out;
    for (const auto& r : records) {
        out.emplace_back(r.binX, r.binY, r.counts);
    }
    std::sort(out.begin(), out.end());
    return out;
}

}  // namespace

TEST_CASE("headers of version 8 and version 9 files") {
    const hicfilecpp::HiCFile v8(kV8);
    CHECK(v8.version() == 8);
    CHECK(v8.getGenomeID() == "sacCer3.chrom.sizes");
    CHECK(v8.getResolutions() == std::vector<int32_t>{1000000, 250000, 50000, 10000});
    CHECK(v8.getChromosomes().size() == 17);
    CHECK(v8.getChromosomes()[0].name == "All");
    CHECK(v8.normVectorIndexHeader() == std::pair<int64_t, int64_t>{0, 0});

    const hicfilecpp::HiCFile v9(kV9);
    CHECK(v9.version() == 9);
    CHECK(v9.normVectorIndexHeader().first > 0);
    CHECK(v9.normVectorIndexHeader().second > 0);

    const hicfilecpp::HiCFile frag(kFrag);
    CHECK(frag.getFragResolutions() == std::vector<int32_t>{100, 20});
    CHECK(frag.fragmentSiteCounts().size() == 17);
}

TEST_CASE("records of a query equal the decoded blocks of the matrix") {
    for (const auto& path : {kV8, kV9}) {
        const hicfilecpp::HiCFile file(path);
        const auto chromosomes = file.getChromosomes();
        for (const auto& [a, b] : {std::pair{1, 1}, std::pair{2, 5}}) {
            const auto mzd = file.getMatrixZoomData(chromosomes[a].name, chromosomes[b].name, "observed", "NONE",
                                                    "BP", 10000);
            REQUIRE(mzd.found());
            const auto records =
                mzd.getRecords(0, chromosomes[a].length, 0, chromosomes[b].length);
            std::vector<hicfilecpp::ContactRecord> blocks;
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
    const hicfilecpp::HiCFile file(kV9);
    const auto mzd = file.getMatrixZoomData("NC_001136.10", "NC_001136.10", "observed", "NONE", "BP", 10000);
    auto collect = [&](int threads) {
        std::vector<std::tuple<int32_t, int32_t, int32_t, float>> out;
        mzd.forEachBlock(
            [&](const hicfilecpp::BlockIndexEntry& entry, std::vector<hicfilecpp::ContactRecord>& records) {
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

TEST_CASE("versions other than 8 and 9 are refused") {
    std::ifstream in(kV8, std::ios::binary);
    std::vector<char> bytes((std::istreambuf_iterator<char>(in)), std::istreambuf_iterator<char>());
    for (const int version : {6, 7, 10}) {
        std::vector<char> copy = bytes;
        std::copy_n(reinterpret_cast<const char*>(&version), 4, copy.begin() + 4);
        const std::string path = scratch("version" + std::to_string(version) + ".hic");
        std::ofstream(path, std::ios::binary).write(copy.data(), static_cast<std::streamsize>(copy.size()));
        CHECK_THROWS_WITH_AS(hicfilecpp::HiCFile{path},
                             ("Version " + std::to_string(version) +
                              " is not supported: hicfilecpp reads .hic versions 8 and 9")
                                 .c_str(),
                             hicfilecpp::HicError);
    }
}

TEST_CASE("invalid files, chromosomes, zoom levels and vectors raise HicError") {
    const std::string text = scratch("not_hic.txt");
    std::ofstream(text) << "this is not a hic file\n";
    CHECK_THROWS_AS(hicfilecpp::HiCFile{text}, hicfilecpp::HicError);
    CHECK_THROWS_AS(hicfilecpp::HiCFile{scratch("does_not_exist.hic")}, hicfilecpp::HicError);

    const hicfilecpp::HiCFile file(kV8);
    CHECK_THROWS_AS(file.getMatrixZoomData("chrNone", "chrNone", "observed", "NONE", "BP", 10000),
                    hicfilecpp::HicError);
    CHECK_THROWS_AS(file.getMatrixZoomData("NC_001133.9", "NC_001133.9", "observed", "NONE", "BP", 12345),
                    hicfilecpp::HicError);
    CHECK_THROWS_AS(file.getMatrixZoomData("NC_001133.9", "NC_001133.9", "observed", "NO_SUCH_NORM", "BP", 10000),
                    hicfilecpp::HicError);
}

TEST_CASE("metadata accessors") {
    const hicfilecpp::HiCFile file(kV8);
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
