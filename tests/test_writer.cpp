#include <doctest/doctest.h>

#include <algorithm>
#include <cmath>
#include <filesystem>
#include <fstream>
#include <iterator>
#include <map>
#include <string>
#include <tuple>
#include <vector>

#include <hicfilecpp/hicfilecpp.hpp>

namespace {

const std::string kData = HICFILECPP_TEST_DATA;
const std::string kJ8 = kData + "/SRR1791297_30.juicer_tools_1.22.01.v8.hic";
const std::string kJ9 = kData + "/SRR1791297_30.juicer_tools_2.20.00.v9.hic";

std::string scratch(const std::string& name) {
    std::filesystem::create_directories(HICFILECPP_TEST_SCRATCH);
    return std::string(HICFILECPP_TEST_SCRATCH) + "/" + name;
}

std::vector<char> bytesOf(const std::string& path) {
    std::ifstream in(path, std::ios::binary);
    return {std::istreambuf_iterator<char>(in), std::istreambuf_iterator<char>()};
}

using Key = std::tuple<int32_t, int32_t, float>;

std::vector<Key> recordsOf(const hicfilecpp::HiCFile& file, const std::string& a, const std::string& b,
                           int32_t resolution, const std::string& norm = "NONE") {
    const auto chromosomes = file.getChromosomes();
    int64_t la = 0;
    int64_t lb = 0;
    for (const auto& c : chromosomes) {
        la = c.name == a ? c.length : la;
        lb = c.name == b ? c.length : lb;
    }
    std::vector<Key> out;
    for (const auto& r :
         file.getMatrixZoomData(a, b, "observed", norm, "BP", resolution).getRecords(0, la, 0, lb)) {
        out.emplace_back(r.binX, r.binY, r.counts);
    }
    std::sort(out.begin(), out.end());
    return out;
}

// The 10 kb pixels of a Juicer-written test file.
class HicSource : public hicfilecpp::PixelSource {
public:
    HicSource(const std::string& path, int32_t resolution) : file_(path), resolution_(resolution) {
        for (const auto& c : file_.getChromosomes()) {
            if (c.index > 0) {
                chromosomes.emplace_back(c.name, c.length);
            }
        }
    }
    void pixels(int32_t resolution, int32_t chr1, int32_t chr2,
                const std::function<void(const hicfilecpp::Pixel*, size_t)>& consume) override {
        REQUIRE(resolution == resolution_);
        if (!file_.hasMatrix(chr1 + 1, chr2 + 1)) {
            return;
        }
        const auto& a = chromosomes[static_cast<size_t>(chr1)];
        const auto& b = chromosomes[static_cast<size_t>(chr2)];
        std::vector<hicfilecpp::Pixel> batch;
        for (const auto& r : file_.getMatrixZoomData(a.first, b.first, "observed", "NONE", "BP", resolution)
                                 .getRecords(0, a.second, 0, b.second)) {
            batch.push_back(hicfilecpp::Pixel{r.binX / resolution, r.binY / resolution, r.counts});
        }
        consume(batch.data(), batch.size());
    }
    std::vector<std::pair<std::string, int64_t>> chromosomes;

private:
    hicfilecpp::HiCFile file_;
    int32_t resolution_;
};

hicfilecpp::WriteOptions optionsFor(const HicSource& source, int32_t version) {
    hicfilecpp::WriteOptions options;
    options.version = version;
    options.genomeId = "sacCer3.chrom.sizes";
    options.chromosomes = source.chromosomes;
    options.resolutions = {10000, 50000};
    options.sourceResolution = 10000;
    return options;
}

double worstRelative(const std::vector<double>& a, const std::vector<double>& b) {
    REQUIRE(a.size() == b.size());
    double worst = 0;
    for (size_t i = 0; i < a.size(); ++i) {
        if (std::isnan(a[i]) || std::isnan(b[i])) {
            CHECK(std::isnan(a[i]) == std::isnan(b[i]));
            continue;
        }
        if (b[i] == 0) {
            CHECK(a[i] == 0);
            continue;
        }
        worst = std::max(worst, std::abs(a[i] - b[i]) / std::abs(b[i]));
    }
    return worst;
}

}  // namespace

TEST_CASE("writing versions 6 and 7 is refused with the reason, other versions as before") {
    HicSource source(kJ8, 10000);
    for (const int32_t version : {6, 7}) {
        const std::string path = scratch("refused.v" + std::to_string(version) + ".hic");
        CHECK_THROWS_WITH_AS(hicfilecpp::writeHicFile(path, optionsFor(source, version), source),
                             ("writing .hic version " + std::to_string(version) +
                              " is not supported: no Juicer tools release that writes it can be "
                              "obtained to validate against; hicfilecpp reads versions 6 to 9 and "
                              "writes 8 and 9")
                                 .c_str(),
                             hicfilecpp::HicError);
    }
    for (const int32_t version : {5, 10}) {
        const std::string path = scratch("refused.v" + std::to_string(version) + ".hic");
        CHECK_THROWS_WITH_AS(hicfilecpp::writeHicFile(path, optionsFor(source, version), source),
                             ("version must be 8 or 9, got " + std::to_string(version)).c_str(),
                             hicfilecpp::HicError);
    }
}

TEST_CASE("written files return the source pixels and Juicer's coarser pixels") {
    for (const int32_t version : {8, 9}) {
        HicSource source(kJ8, 10000);
        const std::string path = scratch("roundtrip.v" + std::to_string(version) + ".hic");
        hicfilecpp::writeHicFile(path, optionsFor(source, version), source);
        const hicfilecpp::HiCFile written(path);
        const hicfilecpp::HiCFile juicer(kJ8);
        CHECK(written.version() == version);
        CHECK(written.getResolutions() == std::vector<int32_t>{50000, 10000});
        CHECK(written.getChromosomes().size() == juicer.getChromosomes().size());
        for (const auto& [a, b] : {std::pair{"NC_001133.9", "NC_001133.9"}, std::pair{"NC_001136.10", "NC_001136.10"},
                                   std::pair{"NC_001134.8", "NC_001145.3"}}) {
            for (const int32_t resolution : {10000, 50000}) {
                const auto mine = recordsOf(written, a, b, resolution);
                CHECK(!mine.empty());
                CHECK(mine == recordsOf(juicer, a, b, resolution));
            }
        }
        const auto wholeMine = recordsOf(written, "All", "All", 24);
        CHECK(!wholeMine.empty());
    }
}

TEST_CASE("normalization vectors and expected values agree with Juicer tools") {
    for (const auto& [version, reference] : {std::pair{8, kJ8}, std::pair{9, kJ9}}) {
        HicSource source(kJ8, 10000);
        const std::string path = scratch("norms.v" + std::to_string(version) + ".hic");
        hicfilecpp::writeHicFile(path, optionsFor(source, version), source);
        const hicfilecpp::HiCFile written(path);
        const hicfilecpp::HiCFile juicer(reference);
        for (const std::string norm : {"VC", "VC_SQRT", "KR", "SCALE"}) {
            for (const int32_t resolution : {10000, 50000}) {
                for (int32_t chr = 1; chr <= 16; ++chr) {
                    const auto mine = written.readNormVector(norm, chr, "BP", resolution);
                    const auto theirs = juicer.readNormVector(norm, chr, "BP", resolution);
                    REQUIRE(mine.has_value() == theirs.has_value());
                    if (mine) {
                        CAPTURE(norm);
                        CAPTURE(resolution);
                        CAPTURE(chr);
                        CHECK(worstRelative(*mine, *theirs) <= 1e-3);
                    }
                }
                const auto mine = written.readExpectedValues({norm, "BP", resolution});
                const auto theirs = juicer.readExpectedValues({norm, "BP", resolution});
                REQUIRE(mine.has_value());
                REQUIRE(theirs.has_value());
                CHECK(worstRelative(mine->values, theirs->values) <= 1e-3);
            }
        }
        for (const int32_t resolution : {10000, 50000}) {
            const auto mine = written.readExpectedValues({"NONE", "BP", resolution});
            const auto theirs = juicer.readExpectedValues({"NONE", "BP", resolution});
            REQUIRE(mine.has_value());
            CHECK(worstRelative(mine->values, theirs->values) <= 1e-6);
        }
    }
}

TEST_CASE("the output does not depend on the number of threads or on the run") {
    for (const int32_t version : {8, 9}) {
        HicSource source(kJ8, 10000);
        auto options = optionsFor(source, version);
        const std::string one = scratch("threads1.hic");
        const std::string many = scratch("threads8.hic");
        const std::string again = scratch("threads8.again.hic");
        options.threads = 1;
        hicfilecpp::writeHicFile(one, options, source);
        options.threads = 8;
        hicfilecpp::writeHicFile(many, options, source);
        hicfilecpp::writeHicFile(again, options, source);
        CHECK(bytesOf(one) == bytesOf(many));
        CHECK(bytesOf(many) == bytesOf(again));
    }
}

TEST_CASE("addNorm on a file written without normalizations equals writing with them") {
    for (const int32_t version : {8, 9}) {
        HicSource source(kJ8, 10000);
        auto options = optionsFor(source, version);
        const std::string with = scratch("with_norms.hic");
        const std::string without = scratch("added_norms.hic");
        hicfilecpp::writeHicFile(with, options, source);
        options.normalizations.clear();
        hicfilecpp::writeHicFile(without, options, source);
        CHECK(hicfilecpp::HiCFile(without).normVectorIndex().empty());
        hicfilecpp::addNorm(without, {"VC", "VC_SQRT", "KR", "SCALE"}, 4);
        CHECK(bytesOf(with) == bytesOf(without));
    }
}

TEST_CASE("invalid writer input raises HicError") {
    HicSource source(kJ8, 10000);
    auto options = optionsFor(source, 9);
    const std::string path = scratch("invalid.hic");
    auto bad = options;
    bad.version = 7;
    CHECK_THROWS_AS(hicfilecpp::writeHicFile(path, bad, source), hicfilecpp::HicError);
    bad = options;
    bad.resolutions = {15000};
    CHECK_THROWS_AS(hicfilecpp::writeHicFile(path, bad, source), hicfilecpp::HicError);
    bad = options;
    bad.normalizations = {"GW_KR"};
    CHECK_THROWS_AS(hicfilecpp::writeHicFile(path, bad, source), hicfilecpp::HicError);
    bad = options;
    bad.chromosomes[0].second = 1000;
    CHECK_THROWS_AS(hicfilecpp::writeHicFile(path, bad, source), hicfilecpp::HicError);
    CHECK_THROWS_AS(hicfilecpp::addNorm(kData + "/SRR1791297_30.juicer_tools_1.22.01.frag.v8.hic", {"VC"}),
                    hicfilecpp::HicError);
}
