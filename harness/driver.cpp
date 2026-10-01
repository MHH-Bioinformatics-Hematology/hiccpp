// hiccpp-harness: runs one harness case through hiccpp.
//
//   hiccpp-harness CASE.json OUT_DIR
//
// Writes OUT_DIR/result.json ({"op_seconds": s, "result": {...}}) and the
// arrays the result refers to as .npy files, in the layout oracle.py writes
// for the same case.

#include <algorithm>
#include <chrono>
#include <cstdio>
#include <fstream>
#include <functional>
#include <map>
#include <sstream>
#include <stdexcept>
#include <string>
#include <utility>
#include <vector>

#include <hiccpp/hiccpp.hpp>

#include "json.hpp"

namespace {

using harness::Json;
using harness::quote;

void writeNpy(const std::string& path, const char* descr, const std::vector<int64_t>& shape, const void* data,
              size_t bytes) {
    std::string dims;
    for (size_t i = 0; i < shape.size(); ++i) {
        dims += std::to_string(shape[i]);
        dims += (shape.size() == 1 || i + 1 < shape.size()) ? "," : "";
        if (i + 1 < shape.size()) {
            dims += " ";
        }
    }
    std::string header = std::string("{'descr': '") + descr + "', 'fortran_order': False, 'shape': (" + dims + "), }";
    const size_t total = 10 + header.size() + 1;
    header += std::string((64 - total % 64) % 64, ' ');
    header += '\n';
    std::ofstream out(path, std::ios::binary);
    out.write("\x93NUMPY\x01\x00", 8);
    const auto length = static_cast<uint16_t>(header.size());
    out.write(reinterpret_cast<const char*>(&length), 2);
    out.write(header.data(), static_cast<std::streamsize>(header.size()));
    out.write(static_cast<const char*>(data), static_cast<std::streamsize>(bytes));
    if (!out) {
        throw std::runtime_error("cannot write " + path);
    }
}

std::string opHeader(const Json& c) {
    const hiccpp::HiCFile file(c["file"].str());
    std::ostringstream out;
    out << "{\"kind\": \"header\", \"genome\": " << quote(file.getGenomeID()) << ", \"resolutions\": [";
    const auto resolutions = file.getResolutions();
    for (size_t i = 0; i < resolutions.size(); ++i) {
        out << (i ? ", " : "") << resolutions[i];
    }
    out << "], \"chromosomes\": [";
    const auto chromosomes = file.getChromosomes();
    for (size_t i = 0; i < chromosomes.size(); ++i) {
        out << (i ? ", " : "") << "[" << quote(chromosomes[i].name) << ", " << chromosomes[i].index << ", "
            << chromosomes[i].length << "]";
    }
    out << "]}";
    return out.str();
}

std::string opRecords(const Json& c, const std::string& dir) {
    const hiccpp::HiCFile file(c["file"].str());
    std::vector<int32_t> binX;
    std::vector<int32_t> binY;
    std::vector<float> counts;
    std::string sizes;
    const Json& queries = c["queries"];
    for (size_t k = 0; k < queries.size(); ++k) {
        const Json& q = queries[k];
        const auto mzd = file.getMatrixZoomData(q[1].str(), q[2].str(), c["matrix_type"].str(), q[0].str(),
                                                c["unit"].str(), static_cast<int32_t>(c["resolution"].i64()));
        const auto records = mzd.getRecords(q[3].i64(), q[4].i64(), q[5].i64(), q[6].i64());
        for (const auto& record : records) {
            binX.push_back(record.binX);
            binY.push_back(record.binY);
            counts.push_back(record.counts);
        }
        sizes += (k ? ", " : "") + std::to_string(records.size());
    }
    const auto n = static_cast<int64_t>(counts.size());
    writeNpy(dir + "/binX.npy", "<i4", {n}, binX.data(), binX.size() * sizeof(int32_t));
    writeNpy(dir + "/binY.npy", "<i4", {n}, binY.data(), binY.size() * sizeof(int32_t));
    writeNpy(dir + "/counts.npy", "<f4", {n}, counts.data(), counts.size() * sizeof(float));
    return "{\"kind\": \"records\", \"n\": [" + sizes + "]}";
}

std::string opMatrices(const Json& c, const std::string& dir) {
    const hiccpp::HiCFile file(c["file"].str());
    std::string shapes;
    const Json& queries = c["queries"];
    for (size_t k = 0; k < queries.size(); ++k) {
        const Json& q = queries[k];
        const auto mzd = file.getMatrixZoomData(q[0].str(), q[1].str(), q[2].str(), q[3].str(), q[4].str(),
                                                static_cast<int32_t>(q[5].i64()));
        const auto matrix = mzd.getRecordsAsMatrix(q[6].i64(), q[7].i64(), q[8].i64(), q[9].i64());
        writeNpy(dir + "/matrix_" + std::to_string(k) + ".npy", "<f4", {matrix.rows, matrix.cols},
                 matrix.values.data(), matrix.values.size() * sizeof(float));
        shapes += (k ? ", " : "") + std::string("[") + std::to_string(matrix.rows) + ", " +
                  std::to_string(matrix.cols) + "]";
    }
    return "{\"kind\": \"matrices\", \"shapes\": [" + shapes + "]}";
}

std::string opVectors(const Json& c, const std::string& dir) {
    const hiccpp::HiCFile file(c["file"].str());
    std::string lengths;
    const Json& queries = c["queries"];
    for (size_t k = 0; k < queries.size(); ++k) {
        const Json& q = queries[k];
        const auto mzd = file.getMatrixZoomData(q[1].str(), q[1].str(), c["matrix_type"].str(), q[0].str(),
                                                c["unit"].str(), static_cast<int32_t>(c["resolution"].i64()));
        const auto expected = mzd.getExpectedValues();
        const auto norm = mzd.getNormVector(static_cast<int32_t>(q[2].i64()));
        writeNpy(dir + "/expected_" + std::to_string(k) + ".npy", "<f8", {static_cast<int64_t>(expected.size())},
                 expected.data(), expected.size() * sizeof(double));
        writeNpy(dir + "/norm_" + std::to_string(k) + ".npy", "<f8", {static_cast<int64_t>(norm.size())},
                 norm.data(), norm.size() * sizeof(double));
        lengths += (k ? ", " : "") + std::string("[") + std::to_string(expected.size()) + ", " +
                   std::to_string(norm.size()) + "]";
    }
    return "{\"kind\": \"vectors\", \"lengths\": [" + lengths + "]}";
}

std::string readText(const std::string& path) {
    std::ifstream in(path);
    std::stringstream buffer;
    buffer << in.rdbuf();
    return buffer.str();
}

// The contacts writer_io.py prep writes: packed (chr1, chr2, bin1, bin2,
// count) records grouped by chromosome pair, read one pair at a time.
class PackedSource : public hiccpp::PixelSource {
public:
    explicit PackedSource(const std::string& dir) : path_(dir + "/pixels.bin") {
        const Json pairs = Json::parse(readText(dir + "/pairs.json"));
        for (size_t i = 0; i < pairs.size(); ++i) {
            ranges_[{static_cast<int32_t>(pairs[i][0].i64()), static_cast<int32_t>(pairs[i][1].i64())}] = {
                pairs[i][2].i64(), pairs[i][3].i64()};
        }
        const Json chroms = Json::parse(readText(dir + "/chroms.json"));
        for (size_t i = 0; i < chroms.size(); ++i) {
            chromosomes.emplace_back(chroms[i][0].str(), chroms[i][1].i64());
        }
    }
    void pixels(int32_t, int32_t chr1, int32_t chr2,
                const std::function<void(const hiccpp::Pixel*, size_t)>& consume) override {
        const auto it = ranges_.find({chr1, chr2});
        if (it == ranges_.end()) {
            return;
        }
        struct Packed {
            int32_t chr1;
            int32_t chr2;
            int32_t bin1;
            int32_t bin2;
            float count;
        };
        std::ifstream in(path_, std::ios::binary);
        in.seekg(it->second.first * static_cast<int64_t>(sizeof(Packed)));
        std::vector<Packed> packed;
        std::vector<hiccpp::Pixel> batch;
        for (int64_t done = it->second.first; done < it->second.second;) {
            const int64_t n = std::min<int64_t>(it->second.second - done, 1 << 20);
            packed.resize(static_cast<size_t>(n));
            in.read(reinterpret_cast<char*>(packed.data()), n * static_cast<int64_t>(sizeof(Packed)));
            if (!in) {
                throw std::runtime_error("short read in " + path_);
            }
            batch.resize(packed.size());
            for (size_t k = 0; k < packed.size(); ++k) {
                batch[k] = hiccpp::Pixel{packed[k].bin1, packed[k].bin2, packed[k].count};
            }
            consume(batch.data(), batch.size());
            done += n;
        }
    }
    std::vector<std::pair<std::string, int64_t>> chromosomes;

private:
    std::string path_;
    std::map<std::pair<int32_t, int32_t>, std::pair<int64_t, int64_t>> ranges_;
};

// The observed pixels of a .hic file at one resolution, decoded block by block.
class HicSource : public hiccpp::PixelSource {
public:
    HicSource(const std::string& path, int threads) : file_(path), threads_(threads) {
        for (const auto& c : file_.getChromosomes()) {
            if (c.index > 0) {
                chromosomes.emplace_back(c.name, c.length);
            }
        }
    }
    void pixels(int32_t resolution, int32_t chr1, int32_t chr2,
                const std::function<void(const hiccpp::Pixel*, size_t)>& consume) override {
        if (!file_.hasMatrix(chr1 + 1, chr2 + 1)) {
            return;
        }
        const auto& a = chromosomes[static_cast<size_t>(chr1)].first;
        const auto& b = chromosomes[static_cast<size_t>(chr2)].first;
        const auto mzd = file_.getMatrixZoomData(a, b, "observed", "NONE", "BP", resolution);
        std::vector<hiccpp::Pixel> batch;
        mzd.forEachBlock(
            [&](const hiccpp::BlockIndexEntry&, std::vector<hiccpp::ContactRecord>& records) {
                batch.resize(records.size());
                for (size_t k = 0; k < records.size(); ++k) {
                    batch[k] = hiccpp::Pixel{records[k].binX, records[k].binY, records[k].counts};
                }
                consume(batch.data(), batch.size());
            },
            threads_);
    }
    std::vector<std::pair<std::string, int64_t>> chromosomes;

private:
    hiccpp::HiCFile file_;
    int threads_;
};

std::string opWrite(const Json& c) {
    hiccpp::WriteOptions options;
    options.version = static_cast<int32_t>(c["version"].i64());
    options.genomeId = c["genome"].str();
    for (size_t i = 0; i < c["resolutions"].size(); ++i) {
        options.resolutions.push_back(static_cast<int32_t>(c["resolutions"][i].i64()));
    }
    options.sourceResolution = static_cast<int32_t>(c["source_resolution"].i64());
    options.normalizations.clear();
    for (size_t i = 0; i < c["normalizations"].size(); ++i) {
        options.normalizations.push_back(c["normalizations"][i].str());
    }
    options.threads = static_cast<int>(c["threads"].i64());
    const bool viaAddNorm = c.has("via_addnorm") && c["via_addnorm"].truthy();
    const std::vector<std::string> norms = options.normalizations;
    if (viaAddNorm) {
        options.normalizations.clear();
    }
    const std::string& output = c["output"].str();
    if (c.has("source_hic")) {
        HicSource source(c["source_hic"].str(), options.threads);
        options.chromosomes = source.chromosomes;
        hiccpp::writeHicFile(output, options, source);
    } else {
        PackedSource source(c["inputs"].str());
        options.chromosomes = source.chromosomes;
        hiccpp::writeHicFile(output, options, source);
    }
    if (viaAddNorm) {
        hiccpp::addNorm(output, norms, options.threads);
    }
    return "{\"kind\": \"written\"}";
}

// Juicer tools pre input from a .hic file: one "short with score" line per
// pixel at the start of its bins, and the chromosome sizes.
std::string opContacts(const Json& c) {
    HicSource source(c["source_hic"].str(), static_cast<int>(c["threads"].i64()));
    const auto resolution = static_cast<int32_t>(c["resolution"].i64());
    {
        std::ofstream sizes(c["sizes"].str());
        for (const auto& [name, length] : source.chromosomes) {
            sizes << name << '\t' << length << '\n';
        }
    }
    std::FILE* out = std::fopen(c["output"].str().c_str(), "w");
    if (out == nullptr) {
        throw std::runtime_error("cannot write " + c["output"].str());
    }
    std::vector<char> buffer(size_t{1} << 22);
    std::setvbuf(out, buffer.data(), _IOFBF, buffer.size());
    int64_t lines = 0;
    const auto n = static_cast<int32_t>(source.chromosomes.size());
    for (int32_t a = 0; a < n; ++a) {
        for (int32_t b = a; b < n; ++b) {
            const char* name1 = source.chromosomes[static_cast<size_t>(a)].first.c_str();
            const char* name2 = source.chromosomes[static_cast<size_t>(b)].first.c_str();
            source.pixels(resolution, a, b, [&](const hiccpp::Pixel* pixels, size_t count) {
                for (size_t k = 0; k < count; ++k) {
                    std::fprintf(out, "0 %s %lld 0 0 %s %lld 1 %.9g\n", name1,
                                 static_cast<long long>(pixels[k].bin1) * resolution, name2,
                                 static_cast<long long>(pixels[k].bin2) * resolution,
                                 static_cast<double>(pixels[k].count));
                }
                lines += static_cast<int64_t>(count);
            });
        }
    }
    if (std::fclose(out) != 0) {
        throw std::runtime_error("cannot write " + c["output"].str());
    }
    return "{\"kind\": \"contacts\", \"lines\": " + std::to_string(lines) + "}";
}

}  // namespace

int main(int argc, char** argv) {
    if (argc != 3) {
        std::fprintf(stderr, "usage: hiccpp-harness CASE.json OUT_DIR\n");
        return 2;
    }
    const std::string text = readText(argv[1]);
    const std::string dir = argv[2];
    const auto start = std::chrono::steady_clock::now();
    std::string result;
    try {
        const Json c = Json::parse(text);
        const std::string& op = c["op"].str();
        if (op == "header") {
            result = opHeader(c);
        } else if (op == "records") {
            result = opRecords(c, dir);
        } else if (op == "matrices") {
            result = opMatrices(c, dir);
        } else if (op == "write") {
            result = opWrite(c);
        } else if (op == "contacts") {
            result = opContacts(c);
        } else if (op == "vectors") {
            result = opVectors(c, dir);
        } else {
            std::fprintf(stderr, "unknown op %s\n", op.c_str());
            return 2;
        }
    } catch (const hiccpp::HicError& error) {
        result = "{\"kind\": \"error\", \"message\": " + quote(error.what()) + "}";
    }
    const double seconds = std::chrono::duration<double>(std::chrono::steady_clock::now() - start).count();
    std::ofstream out(dir + "/result.json");
    out << "{\"op_seconds\": " << seconds << ", \"result\": " << result << "}\n";
    return out ? 0 : 1;
}
