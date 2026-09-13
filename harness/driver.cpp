// hicfilecpp-harness: runs one harness case through hicfilecpp.
//
//   hicfilecpp-harness CASE.json OUT_DIR
//
// Writes OUT_DIR/result.json ({"op_seconds": s, "result": {...}}) and the
// arrays the result refers to as .npy files, in the layout oracle.py writes
// for the same case.

#include <chrono>
#include <cstdio>
#include <fstream>
#include <sstream>
#include <string>
#include <vector>

#include <hicfilecpp/hicfilecpp.hpp>

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
    const hicfilecpp::HiCFile file(c["file"].str());
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
    const hicfilecpp::HiCFile file(c["file"].str());
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
    const hicfilecpp::HiCFile file(c["file"].str());
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
    const hicfilecpp::HiCFile file(c["file"].str());
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

}  // namespace

int main(int argc, char** argv) {
    if (argc != 3) {
        std::fprintf(stderr, "usage: hicfilecpp-harness CASE.json OUT_DIR\n");
        return 2;
    }
    std::string text;
    {
        std::ifstream in(argv[1]);
        std::stringstream buffer;
        buffer << in.rdbuf();
        text = buffer.str();
    }
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
        } else if (op == "vectors") {
            result = opVectors(c, dir);
        } else {
            std::fprintf(stderr, "unknown op %s\n", op.c_str());
            return 2;
        }
    } catch (const hicfilecpp::HicError& error) {
        result = "{\"kind\": \"error\", \"message\": " + quote(error.what()) + "}";
    }
    const double seconds = std::chrono::duration<double>(std::chrono::steady_clock::now() - start).count();
    std::ofstream out(dir + "/result.json");
    out << "{\"op_seconds\": " << seconds << ", \"result\": " << result << "}\n";
    return out ? 0 : 1;
}
