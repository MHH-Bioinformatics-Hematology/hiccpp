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
const std::string kJ8Frag = kData + "/SRR1791297_30.juicer_tools_1.22.01.frag.v8.hic";
const std::string kJ9Frag = kData + "/SRR1791297_30.juicer_tools_2.20.00.frag.v9.hic";

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

// The pixels of a Juicer-written file with fragment maps: base pair pixels at
// the finest base pair resolution and fragment pixels at the finest fragment
// one, with the file's own restriction sites.
class FragSource : public hicfilecpp::PixelSource {
public:
    FragSource(const std::string& path, int32_t resolution, int32_t fragResolution)
        : file_(path), resolution_(resolution), fragResolution_(fragResolution) {
        for (const auto& c : file_.getChromosomes()) {
            if (c.index > 0) {
                chromosomes.emplace_back(c.name, c.length);
            }
        }
        sites.assign(file_.fragmentSites().begin() + 1, file_.fragmentSites().end());
    }

    // A fragment bin region spans resolution fragments, and a fragment number
    // may be the site count itself, so the query reaches one bin past the last.
    int64_t fragExtent(int32_t chr, int32_t resolution) const {
        return static_cast<int64_t>(sites[static_cast<size_t>(chr)].size()) * resolution + resolution;
    }

    void emit(int32_t resolution, int32_t chr1, int32_t chr2, const std::string& unit,
              const std::function<void(const hicfilecpp::Pixel*, size_t)>& consume) {
        if (!file_.hasMatrix(chr1 + 1, chr2 + 1)) {
            return;
        }
        const auto& a = chromosomes[static_cast<size_t>(chr1)];
        const auto& b = chromosomes[static_cast<size_t>(chr2)];
        const int64_t ea = unit == "FRAG" ? fragExtent(chr1, resolution) : a.second;
        const int64_t eb = unit == "FRAG" ? fragExtent(chr2, resolution) : b.second;
        std::vector<hicfilecpp::Pixel> batch;
        for (const auto& r :
             file_.getMatrixZoomData(a.first, b.first, "observed", "NONE", unit, resolution).getRecords(0, ea, 0, eb)) {
            batch.push_back(hicfilecpp::Pixel{r.binX / resolution, r.binY / resolution, r.counts});
        }
        consume(batch.data(), batch.size());
    }

    void pixels(int32_t resolution, int32_t chr1, int32_t chr2,
                const std::function<void(const hicfilecpp::Pixel*, size_t)>& consume) override {
        REQUIRE(resolution == resolution_);
        emit(resolution, chr1, chr2, "BP", consume);
    }
    void fragPixels(int32_t resolution, int32_t chr1, int32_t chr2,
                    const std::function<void(const hicfilecpp::Pixel*, size_t)>& consume) override {
        REQUIRE(resolution == fragResolution_);
        emit(resolution, chr1, chr2, "FRAG", consume);
    }

    std::vector<std::pair<std::string, int64_t>> chromosomes;
    std::vector<std::vector<int32_t>> sites;

private:
    hicfilecpp::HiCFile file_;
    int32_t resolution_;
    int32_t fragResolution_;
};

// The resolutions of the FRAG test files: pre -r 500000,50000,100f,20f.
hicfilecpp::WriteOptions fragOptionsFor(const FragSource& source, int32_t version) {
    hicfilecpp::WriteOptions options;
    options.version = version;
    options.genomeId = "sacCer3.chrom.sizes";
    options.chromosomes = source.chromosomes;
    options.resolutions = {500000, 50000};
    options.fragResolutions = {100, 20};
    options.fragmentSites = source.sites;
    options.sourceResolution = 50000;
    options.sourceFragResolution = 20;
    // Juicer tools 2.20.00 writes no KR unless -k asks for it, so the version 9
    // reference holds none.
    if (version == 9) {
        options.normalizations = {"VC", "VC_SQRT", "SCALE"};
    }
    return options;
}

std::vector<Key> fragRecordsOf(const hicfilecpp::HiCFile& file, const FragSource& source, int32_t chr1, int32_t chr2,
                               int32_t resolution) {
    std::vector<Key> out;
    for (const auto& r : file
                             .getMatrixZoomData(source.chromosomes[static_cast<size_t>(chr1)].first,
                                                source.chromosomes[static_cast<size_t>(chr2)].first, "observed",
                                                "NONE", "FRAG", resolution)
                             .getRecords(0, source.fragExtent(chr1, resolution), 0,
                                         source.fragExtent(chr2, resolution))) {
        out.emplace_back(r.binX, r.binY, r.counts);
    }
    std::sort(out.begin(), out.end());
    return out;
}

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

TEST_CASE("provided normalization vectors are stored as given with their expected values") {
    for (const auto& [version, reference] : {std::pair{8, kJ8}, std::pair{9, kJ9}}) {
        HicSource source(kJ8, 10000);
        auto options = optionsFor(source, version);
        const hicfilecpp::HiCFile juicer(reference);
        options.normalizations = {"VC"};
        options.providedNormalizations = {"KR", "GW_KR"};
        int32_t asked = 0;
        options.normVector = [&](const std::string& name, int32_t chrIndex, int32_t resolution) {
            ++asked;
            // GW_KR is a label Juicer tools also uses; its values are KR's.
            auto vector = juicer.readNormVector(name == "GW_KR" ? "KR" : name, chrIndex + 1, "BP", resolution);
            REQUIRE(vector.has_value());
            if (chrIndex == 3 && name == "GW_KR") {
                return std::vector<double>{};  // no vector for this chromosome
            }
            vector->resize(vector->size() - 1);  // the padding is restored
            return *vector;
        };
        const std::string path = scratch("provided.v" + std::to_string(version) + ".hic");
        hicfilecpp::writeHicFile(path, options, source);
        CHECK(asked == 2 * 2 * 16);
        const hicfilecpp::HiCFile written(path);
        CHECK(written.hasNormalizedExpectedSection());
        for (const int32_t resolution : {10000, 50000}) {
            for (int32_t chr = 1; chr <= 16; ++chr) {
                const auto theirs = juicer.readNormVector("KR", chr, "BP", resolution);
                const auto kr = written.readNormVector("KR", chr, "BP", resolution);
                REQUIRE(kr.has_value());
                REQUIRE(kr->size() == theirs->size());
                for (size_t i = 0; i < kr->size(); ++i) {
                    CHECK((std::isnan((*kr)[i]) ? std::isnan((*theirs)[i]) : (*kr)[i] == (*theirs)[i]));
                }
                CHECK(written.readNormVector("GW_KR", chr, "BP", resolution).has_value() == (chr != 4));
                CHECK(written.readNormVector("VC", chr, "BP", resolution).has_value());
            }
            const auto mine = written.readExpectedValues({"KR", "BP", resolution});
            const auto theirs = juicer.readExpectedValues({"KR", "BP", resolution});
            REQUIRE(mine.has_value());
            CHECK(worstRelative(mine->values, theirs->values) <= 1e-3);
            CHECK(written.readExpectedValues({"GW_KR", "BP", resolution}).has_value());
        }
        const auto types = written.getNormalizationTypes();
        CHECK(std::find(types.begin(), types.end(), "GW_KR") != types.end());

        auto bad = options;
        bad.providedNormalizations = {"VC"};
        CHECK_THROWS_AS(hicfilecpp::writeHicFile(scratch("bad.hic"), bad, source), hicfilecpp::HicError);
        bad = options;
        bad.normVector = nullptr;
        CHECK_THROWS_AS(hicfilecpp::writeHicFile(scratch("bad.hic"), bad, source), hicfilecpp::HicError);
        bad = options;
        bad.providedNormalizations = {"K R"};
        CHECK_THROWS_AS(hicfilecpp::writeHicFile(scratch("bad.hic"), bad, source), hicfilecpp::HicError);
    }
}

TEST_CASE("a footer without a normalized expected-value section is reported") {
    HicSource source(kJ8, 10000);
    auto options = optionsFor(source, 8);
    options.normalizations.clear();
    const std::string path = scratch("no_norm_section.hic");
    hicfilecpp::writeHicFile(path, options, source);
    CHECK(hicfilecpp::HiCFile(path).hasNormalizedExpectedSection());
    // Version 8 without normalizations ends with two zero counts: the
    // normalized expected values and the vector index.
    std::filesystem::resize_file(path, std::filesystem::file_size(path) - 8);
    const hicfilecpp::HiCFile truncated(path);
    CHECK_FALSE(truncated.hasNormalizedExpectedSection());
    CHECK(truncated.getNormalizationTypes().empty());
    CHECK(hicfilecpp::HiCFile(kJ9).hasNormalizedExpectedSection());
}

TEST_CASE("fragment resolutions reproduce the FRAG records, vectors and expected values of Juicer tools") {
    for (const auto& [version, reference] : {std::pair{8, kJ8Frag}, std::pair{9, kJ9Frag}}) {
        FragSource source(reference, 50000, 20);
        const std::string path = scratch("frag.v" + std::to_string(version) + ".hic");
        hicfilecpp::writeHicFile(path, fragOptionsFor(source, version), source);
        const hicfilecpp::HiCFile written(path);
        const hicfilecpp::HiCFile juicer(reference);
        CAPTURE(version);
        CHECK(written.getResolutions() == juicer.getResolutions());
        CHECK(written.getFragResolutions() == juicer.getFragResolutions());
        CHECK(written.fragmentSiteCounts() == juicer.fragmentSiteCounts());
        CHECK(written.fragmentSites() == juicer.fragmentSites());

        const auto nChromosomes = static_cast<int32_t>(source.chromosomes.size());
        // Every zoom header: the unit, the zoom index, the sum and the block
        // layout Juicer's MatrixPP arithmetic gives the fragment zooms.
        for (int32_t chr = 0; chr < nChromosomes; ++chr) {
            const auto mine = written.matrixZoomHeaders(chr + 1, chr + 1);
            const auto theirs = juicer.matrixZoomHeaders(chr + 1, chr + 1);
            REQUIRE(mine.size() == theirs.size());
            for (size_t k = 0; k < mine.size(); ++k) {
                CAPTURE(chr);
                CAPTURE(mine[k].unit);
                CAPTURE(mine[k].binSize);
                CHECK(mine[k].unit == theirs[k].unit);
                CHECK(mine[k].binSize == theirs[k].binSize);
                CHECK(mine[k].zoomIndex == theirs[k].zoomIndex);
                CHECK(mine[k].blockBinCount == theirs[k].blockBinCount);
                CHECK(mine[k].blockColumnCount == theirs[k].blockColumnCount);
                CHECK(mine[k].sumCounts == theirs[k].sumCounts);
            }
        }

        // Every FRAG record of every pair, and the base pair records too.
        int64_t records = 0;
        for (int32_t c1 = 0; c1 < nChromosomes; ++c1) {
            for (int32_t c2 = c1; c2 < nChromosomes; ++c2) {
                for (const int32_t resolution : {100, 20}) {
                    const auto theirs = fragRecordsOf(juicer, source, c1, c2, resolution);
                    CAPTURE(c1);
                    CAPTURE(c2);
                    CAPTURE(resolution);
                    CHECK(fragRecordsOf(written, source, c1, c2, resolution) == theirs);
                    records += static_cast<int64_t>(theirs.size());
                }
            }
        }
        CHECK(records > 40000);
        for (const int32_t resolution : {500000, 50000}) {
            const auto& name = source.chromosomes[1].first;
            CHECK(recordsOf(written, name, name, resolution) == recordsOf(juicer, name, name, resolution));
        }

        // The FRAG expected values, their per-chromosome factors and the FRAG
        // normalization vectors, for every entry the reference holds.
        int32_t expectedKeys = 0;
        for (const auto& key : juicer.expectedValuesKeys()) {
            if (key.unit != "FRAG") {
                continue;
            }
            const auto mine = written.readExpectedValues(key);
            const auto theirs = juicer.readExpectedValues(key);
            REQUIRE(theirs.has_value());
            REQUIRE(mine.has_value());
            CAPTURE(key.normalization);
            CAPTURE(key.binSize);
            CHECK(worstRelative(mine->values, theirs->values) <= 1e-3);
            REQUIRE(mine->normalizationFactors.size() == theirs->normalizationFactors.size());
            for (size_t k = 0; k < mine->normalizationFactors.size(); ++k) {
                CHECK(mine->normalizationFactors[k].first == theirs->normalizationFactors[k].first);
                CHECK(mine->normalizationFactors[k].second ==
                      doctest::Approx(theirs->normalizationFactors[k].second).epsilon(1e-3));
            }
            ++expectedKeys;
        }
        // Two raw entries and one per normalization per fragment resolution.
        CHECK(expectedKeys == (version == 8 ? 10 : 8));
        int32_t vectors = 0;
        for (const auto& entry : juicer.normVectorIndex()) {
            if (entry.unit != "FRAG") {
                continue;
            }
            const auto mine = written.readNormVector(entry.normalization, entry.chrIndex, "FRAG", entry.resolution);
            const auto theirs = juicer.readNormVector(entry.normalization, entry.chrIndex, "FRAG", entry.resolution);
            REQUIRE(theirs.has_value());
            REQUIRE(mine.has_value());
            CAPTURE(entry.normalization);
            CAPTURE(entry.chrIndex);
            CAPTURE(entry.resolution);
            REQUIRE(mine->size() == theirs->size());
            CHECK(worstRelative(*mine, *theirs) <= 1e-3);
            ++vectors;
        }
        CHECK(vectors == (version == 8 ? 128 : 96));
    }
}

TEST_CASE("a fragment file does not depend on the number of threads") {
    for (const int32_t version : {8, 9}) {
        FragSource source(kJ8Frag, 50000, 20);
        auto options = fragOptionsFor(source, version);
        const std::string one = scratch("frag.threads1.hic");
        const std::string many = scratch("frag.threads8.hic");
        options.threads = 1;
        hicfilecpp::writeHicFile(one, options, source);
        options.threads = 8;
        hicfilecpp::writeHicFile(many, options, source);
        CHECK(bytesOf(one) == bytesOf(many));
    }
}

TEST_CASE("fragment options are checked and a base pair source needs none") {
    FragSource source(kJ8Frag, 50000, 20);
    const auto options = fragOptionsFor(source, 9);
    const std::string path = scratch("frag.invalid.hic");
    auto bad = options;
    bad.fragResolutions = {30};  // not a multiple of the source fragment resolution
    CHECK_THROWS_AS(hicfilecpp::writeHicFile(path, bad, source), hicfilecpp::HicError);
    bad = options;
    bad.sourceFragResolution = 0;
    CHECK_THROWS_AS(hicfilecpp::writeHicFile(path, bad, source), hicfilecpp::HicError);
    bad = options;
    bad.fragmentSites.pop_back();
    CHECK_THROWS_AS(hicfilecpp::writeHicFile(path, bad, source), hicfilecpp::HicError);
    bad = options;
    bad.fragResolutions.clear();
    CHECK_THROWS_AS(hicfilecpp::writeHicFile(path, bad, source), hicfilecpp::HicError);
    bad = options;
    bad.fragmentSites[0] = {1000, 500};  // out of order
    CHECK_THROWS_AS(hicfilecpp::writeHicFile(path, bad, source), hicfilecpp::HicError);
    bad = options;
    bad.fragResolutions = {100, 100};
    CHECK_THROWS_AS(hicfilecpp::writeHicFile(path, bad, source), hicfilecpp::HicError);

    // A source that does not override fragPixels hands over no fragment pixels,
    // so the matrices hold base pair zooms only.
    HicSource plain(kJ8, 10000);
    auto onlyBp = optionsFor(plain, 9);
    onlyBp.fragResolutions = {100};
    onlyBp.sourceFragResolution = 20;
    onlyBp.fragmentSites = source.sites;
    const std::string onlyBpPath = scratch("frag.none.hic");
    hicfilecpp::writeHicFile(onlyBpPath, onlyBp, plain);
    const hicfilecpp::HiCFile written(onlyBpPath);
    CHECK(written.getFragResolutions() == std::vector<int32_t>{100});
    CHECK(written.fragmentSiteCounts().size() == source.sites.size() + 1);
    for (const auto& header : written.matrixZoomHeaders(1, 1)) {
        CHECK(header.unit == "BP");
    }
}

TEST_CASE("provided normalization vectors are stored for the fragment resolutions too") {
    FragSource source(kJ8Frag, 50000, 20);
    auto options = fragOptionsFor(source, 8);
    const hicfilecpp::HiCFile juicer(kJ8Frag);
    options.normalizations = {"VC"};
    options.providedNormalizations = {"KR"};
    options.normVector = [&](const std::string& name, int32_t chrIndex, int32_t resolution) {
        auto vector = juicer.readNormVector(name, chrIndex + 1, "BP", resolution);
        return vector.value_or(std::vector<double>{});
    };
    int32_t askedFrag = 0;
    options.fragNormVector = [&](const std::string& name, int32_t chrIndex, int32_t resolution) {
        ++askedFrag;
        auto vector = juicer.readNormVector(name, chrIndex + 1, "FRAG", resolution);
        return vector.value_or(std::vector<double>{});
    };
    const std::string path = scratch("frag.provided.hic");
    hicfilecpp::writeHicFile(path, options, source);
    CHECK(askedFrag == 2 * 16);
    const hicfilecpp::HiCFile written(path);
    for (const int32_t resolution : {100, 20}) {
        for (int32_t chr = 1; chr <= 16; ++chr) {
            const auto theirs = juicer.readNormVector("KR", chr, "FRAG", resolution);
            const auto mine = written.readNormVector("KR", chr, "FRAG", resolution);
            REQUIRE(theirs.has_value());
            REQUIRE(mine.has_value());
            CAPTURE(resolution);
            CAPTURE(chr);
            REQUIRE(mine->size() == theirs->size());
            for (size_t i = 0; i < mine->size(); ++i) {
                CHECK((std::isnan((*mine)[i]) ? std::isnan((*theirs)[i]) : (*mine)[i] == (*theirs)[i]));
            }
        }
        // The expected values the provided vectors give match Juicer's KR ones,
        // which came from the same vectors.
        const auto mine = written.readExpectedValues({"KR", "FRAG", resolution});
        const auto theirs = juicer.readExpectedValues({"KR", "FRAG", resolution});
        REQUIRE(mine.has_value());
        REQUIRE(theirs.has_value());
        CHECK(worstRelative(mine->values, theirs->values) <= 1e-3);
    }
    // Without fragNormVector the fragment resolutions hold the computed
    // normalizations only.
    options.fragNormVector = nullptr;
    const std::string other = scratch("frag.provided.none.hic");
    hicfilecpp::writeHicFile(other, options, source);
    const hicfilecpp::HiCFile plain(other);
    CHECK_FALSE(plain.readNormVector("KR", 1, "FRAG", 20).has_value());
    CHECK(plain.readNormVector("VC", 1, "FRAG", 20).has_value());
    CHECK(plain.readNormVector("KR", 1, "BP", 50000).has_value());
}
