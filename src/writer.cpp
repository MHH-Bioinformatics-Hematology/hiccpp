// Writing .hic files: Juicer tools pre (juicebox.tools.utils.original.
// Preprocessor, MatrixPP, MatrixZoomDataPP) and addNorm
// (juicebox.tools.utils.norm.NormalizationVectorUpdater, NormVectorUpdater),
// following release 1.22.01 for version 8 and 2.20.00 for version 9.
//
// Layout of a written file: the header; for every chromosome pair with
// pixels, the blocks of every resolution followed by the pair's matrix header
// and block indexes; the whole-genome matrix the same way; the footer with the
// matrix index, the raw expected values, the normalized expected values, the
// normalization vector index and the vectors. Juicer places each matrix
// header before its blocks; readers locate both through the index, so the
// order does not matter to them.

#include "hicfilecpp/writer.hpp"

#include <algorithm>
#include <climits>
#include <functional>
#include <cmath>
#include <cstring>
#include <limits>
#include <map>
#include <optional>
#include <set>
#include <string>
#include <vector>

#include "codec.hpp"
#include "hicfilecpp/hicfile.hpp"
#include "hicfilecpp/version.hpp"
#include "norms.hpp"
#include "output.hpp"
#include "state.hpp"

namespace hicfilecpp {

namespace detail {

namespace {

const std::vector<std::string> kNormOrder{"VC", "VC_SQRT", "KR", "SCALE"};
constexpr int32_t kBlockSize = 1000;  // Preprocessor.BLOCK_SIZE
constexpr int32_t kMaxSqrt = 46340;   // (int) Math.sqrt(Integer.MAX_VALUE)

struct Layout {
    int32_t binSize = 0;
    int32_t blockBinCount = 0;
    int32_t blockColumnCount = 0;
    bool intra = false;
    int32_t version = 9;
};

// MatrixPP.getNumColumnsFromNumBins. The cutoff is 0 for a fragment zoom,
// which is how MatrixPP calls it there, so the version 9 column widening of
// fine base pair zooms never applies to fragment resolutions.
int32_t numColumns(int32_t version, int32_t nBins, int32_t binSize, int32_t cutoff) {
    int32_t nColumns = nBins / kBlockSize + 1;
    if (version == 8) {
        if (nColumns > std::sqrt(static_cast<double>(INT32_MAX))) {
            nColumns = kMaxSqrt - 1;
        }
        return nColumns;
    }
    if (binSize < cutoff) {
        const int64_t numerator = static_cast<int64_t>(nBins) * binSize;
        const int64_t denominator = static_cast<int64_t>(kBlockSize) * cutoff;
        nColumns = static_cast<int32_t>(numerator / denominator) + 1;
    }
    return std::min(nColumns, kMaxSqrt - 1);
}

// The MatrixPP and MatrixZoomDataPP constructors. For a fragment zoom
// `extent` is the chromosome's site count (MatrixPP takes
// FragmentCalculation.getNumberFragments where it takes the length for a base
// pair zoom) and the cutoff is 0.
Layout pairLayout(int32_t version, int64_t extent1, int64_t extent2, int32_t binSize, bool intra, bool frag) {
    Layout layout;
    layout.binSize = binSize;
    layout.intra = intra;
    layout.version = version;
    const int64_t extent = std::max(extent1, extent2);
    const auto nBins = static_cast<int32_t>(extent / binSize + 1);
    const int32_t cutoff = frag ? 0 : (intra ? 500 : 5000);
    layout.blockColumnCount = numColumns(version, nBins, binSize, cutoff);
    layout.blockBinCount = nBins / layout.blockColumnCount + 1;
    return layout;
}

// Preprocessor.getInitialGenomeWideMatrixPP
Layout wholeGenomeLayout(int32_t version, int32_t genomeLengthKb) {
    Layout layout;
    layout.binSize = std::max(genomeLengthKb / 500, 1);
    layout.intra = true;
    layout.version = version;
    const int32_t nBinsX = genomeLengthKb / layout.binSize + 1;
    layout.blockColumnCount = nBinsX / kBlockSize + 1;
    layout.blockBinCount = nBinsX / layout.blockColumnCount + 1;
    return layout;
}

// MatrixZoomDataPP.incrementCount; version 9 numbers intra-chromosomal blocks
// by their distance from the diagonal (LogDepth, base 2).
int32_t blockNumber(const Layout& layout, int32_t x, int32_t y) {
    if (layout.intra && layout.version > 8) {
        const double v = std::abs(x - y) / std::sqrt(2) / layout.blockBinCount;
        const auto depth = static_cast<int32_t>(std::log(1 + v) / std::log(2));
        const int32_t position = (x + y) / 2 / layout.blockBinCount;
        return depth * layout.blockColumnCount + position;
    }
    return layout.blockColumnCount * (y / layout.blockBinCount) + x / layout.blockBinCount;
}

struct Entry {
    int32_t block;
    int32_t x;
    int32_t y;
    float count;
};

// A stable LSD radix sort on the block number.
void sortByBlock(std::vector<Entry>& entries) {
    std::vector<Entry> scratch(entries.size());
    std::vector<size_t> counts(size_t{1} << 16);
    for (int pass = 0; pass < 2; ++pass) {
        const int shift = pass * 16;
        std::fill(counts.begin(), counts.end(), 0);
        for (const Entry& e : entries) {
            counts[(static_cast<uint32_t>(e.block) >> shift) & 0xFFFF]++;
        }
        size_t total = 0;
        for (auto& count : counts) {
            const size_t n = count;
            count = total;
            total += n;
        }
        for (const Entry& e : entries) {
            scratch[counts[(static_cast<uint32_t>(e.block) >> shift) & 0xFFFF]++] = e;
        }
        entries.swap(scratch);
    }
}

// MatrixZoomDataPP.writeBlock, before compression. Records are sorted by y,
// then x. Version 8 picks the smaller of the list-of-rows and dense
// representations as 1.22.01 does; version 9 always writes list of rows, as
// 2.20.00 does.
void encodeBlock(int32_t version, const std::vector<NormRecord>& records, ByteWriter& out) {
    out.bytes.clear();
    out.bytes.reserve(records.size() * 6 + 32);
    int32_t binXOffset = INT32_MAX;
    int32_t binYOffset = INT32_MAX;
    int32_t binXMax = 0;
    int32_t binYMax = 0;
    bool isInteger = true;
    float maxCounts = 0;
    size_t rows = 0;
    for (size_t i = 0; i < records.size(); ++i) {
        const NormRecord& r = records[i];
        binXOffset = std::min(binXOffset, r.x);
        binYOffset = std::min(binYOffset, r.y);
        binXMax = std::max(binXMax, r.x);
        binYMax = std::max(binYMax, r.y);
        isInteger = isInteger && std::floor(r.counts) == r.counts;
        maxCounts = std::max(r.counts, maxCounts);
        if (i == 0 || r.y != records[i - 1].y) {
            rows++;
        }
    }
    const bool useShort = isInteger && maxCounts < 32767;
    const int64_t valueSize = useShort ? 2 : 4;
    const int64_t lorSize = static_cast<int64_t>(rows) * 4 + static_cast<int64_t>(records.size()) * valueSize;
    auto putValue = [&](float counts) {
        if (useShort) {
            out.put<int16_t>(static_cast<int16_t>(counts));
        } else {
            out.put<float>(counts);
        }
    };
    out.put<int32_t>(static_cast<int32_t>(records.size()));
    out.put<int32_t>(binXOffset);
    out.put<int32_t>(binYOffset);
    out.put<int8_t>(useShort ? 0 : 1);

    bool rowShort = true;
    bool colShort = true;
    if (version > 8) {
        colShort = binXMax - binXOffset + 1 < 32767;
        rowShort = binYMax - binYOffset + 1 < 32767;
        out.put<int8_t>(colShort ? 0 : 1);
        out.put<int8_t>(rowShort ? 0 : 1);
    } else {
        const auto w = static_cast<int16_t>(binXMax - binXOffset + 1);
        const NormRecord& last = records.back();
        const int64_t nDensePts = static_cast<int64_t>(last.y - binYOffset) * w + (last.x - binXOffset) + 1;
        if (lorSize >= nDensePts * valueSize) {
            out.put<int8_t>(2);
            out.put<int32_t>(static_cast<int32_t>(nDensePts));
            out.put<int16_t>(w);
            int64_t lastIdx = 0;
            for (const NormRecord& r : records) {
                const int64_t idx = static_cast<int64_t>(r.y - binYOffset) * w + (r.x - binXOffset);
                for (int64_t i = lastIdx; i < idx; ++i) {
                    if (useShort) {
                        out.put<int16_t>(INT16_MIN);
                    } else {
                        out.put<float>(std::numeric_limits<float>::quiet_NaN());
                    }
                }
                putValue(r.counts);
                lastIdx = idx + 1;
            }
            return;
        }
    }
    out.put<int8_t>(1);
    auto putRowIndex = [&](int64_t value) {
        if (rowShort) {
            out.put<int16_t>(static_cast<int16_t>(value));
        } else {
            out.put<int32_t>(static_cast<int32_t>(value));
        }
    };
    auto putColIndex = [&](int64_t value) {
        if (colShort) {
            out.put<int16_t>(static_cast<int16_t>(value));
        } else {
            out.put<int32_t>(static_cast<int32_t>(value));
        }
    };
    putRowIndex(static_cast<int64_t>(rows));
    for (size_t begin = 0; begin < records.size();) {
        size_t end = begin + 1;
        while (end < records.size() && records[end].y == records[begin].y) {
            end++;
        }
        putRowIndex(records[begin].y - binYOffset);
        putColIndex(static_cast<int64_t>(end - begin));
        for (size_t i = begin; i < end; ++i) {
            putColIndex(records[i].x - binXOffset);
            putValue(records[i].counts);
        }
        begin = end;
    }
}

struct EncodedBlock {
    int32_t number = 0;
    std::vector<char> bytes;
    std::vector<NormRecord> records;
};

// Sorts the entries of one zoom level into blocks, merges repeated pixels in
// input order, encodes and compresses the blocks on the pool and writes them
// in ascending block number. With `keep`, also returns the merged records in
// that order, which is the order Juicer reads them back for normalization.
std::vector<BlockIndexEntry> writeBlocks(OutputFile& out, ThreadPool& pool, int32_t version, int level,
                                         std::vector<Entry>& entries, std::vector<NormRecord>* keep) {
    sortByBlock(entries);
    std::vector<std::pair<size_t, size_t>> ranges;
    for (size_t i = 0; i < entries.size();) {
        size_t j = i + 1;
        while (j < entries.size() && entries[j].block == entries[i].block) {
            j++;
        }
        ranges.emplace_back(i, j);
        i = j;
    }
    std::vector<BlockIndexEntry> index;
    index.reserve(ranges.size());
    const size_t window = static_cast<size_t>(pool.size()) * 32;
    std::vector<EncodedBlock> encoded;
    for (size_t begin = 0; begin < ranges.size(); begin += window) {
        const size_t end = std::min(ranges.size(), begin + window);
        encoded.assign(end - begin, EncodedBlock{});
        pool.parallelFor(end - begin, [&](size_t i) {
            const auto [first, last] = ranges[begin + i];
            std::stable_sort(entries.begin() + static_cast<std::ptrdiff_t>(first),
                             entries.begin() + static_cast<std::ptrdiff_t>(last),
                             [](const Entry& a, const Entry& b) { return a.y != b.y ? a.y < b.y : a.x < b.x; });
            std::vector<NormRecord> records;
            records.reserve(last - first);
            for (size_t k = first; k < last; ++k) {
                const Entry& e = entries[k];
                if (!records.empty() && records.back().x == e.x && records.back().y == e.y) {
                    records.back().counts += e.count;
                } else {
                    records.push_back(NormRecord{e.x, e.y, e.count});
                }
            }
            ByteWriter raw;
            encodeBlock(version, records, raw);
            EncodedBlock& block = encoded[i];
            block.number = entries[first].block;
            deflateZlib(raw.bytes.data(), raw.bytes.size(), level, block.bytes);
            if (keep != nullptr) {
                block.records = std::move(records);
            }
        });
        for (auto& block : encoded) {
            const int64_t position = out.position();
            out.write(block.bytes);
            index.push_back(BlockIndexEntry{block.number, position, static_cast<int32_t>(block.bytes.size())});
            if (keep != nullptr) {
                keep->insert(keep->end(), block.records.begin(), block.records.end());
            }
        }
    }
    return index;
}

struct ZoomWritten {
    std::string unit = "BP";
    int32_t zoomIndex = 0;
    double sum = 0;
    Layout layout;
    std::vector<BlockIndexEntry> index;
};

// Preprocessor.writeMatrix and writeZoomHeader. Juicer writes the header
// before the blocks, so the stored sum holds the contacts (off-diagonal
// intra-chromosomal ones twice) and the cell count and percentiles are 0.
MatrixEntry writeMatrixHeader(OutputFile& out, int32_t chr1, int32_t chr2, const std::vector<ZoomWritten>& zooms) {
    ByteWriter m;
    m.put<int32_t>(chr1);
    m.put<int32_t>(chr2);
    m.put<int32_t>(static_cast<int32_t>(zooms.size()));
    for (const auto& zoom : zooms) {
        m.cstr(zoom.unit);
        m.put<int32_t>(zoom.zoomIndex);
        m.put<float>(static_cast<float>(zoom.sum));
        m.put<float>(0.0f);
        m.put<float>(0.0f);
        m.put<float>(0.0f);
        m.put<int32_t>(zoom.layout.binSize);
        m.put<int32_t>(zoom.layout.blockBinCount);
        m.put<int32_t>(zoom.layout.blockColumnCount);
        m.put<int32_t>(static_cast<int32_t>(zoom.index.size()));
        for (const auto& entry : zoom.index) {
            m.put<int32_t>(entry.number);
            m.put<int64_t>(entry.position);
            m.put<int32_t>(entry.size);
        }
    }
    MatrixEntry entry;
    entry.position = out.position();
    entry.size = static_cast<int32_t>(m.size());
    out.write(m);
    return entry;
}

void putValues(ByteWriter& w, int32_t version, const std::vector<double>& values) {
    if (version > 8) {
        w.put<int64_t>(static_cast<int64_t>(values.size()));
        for (const double v : values) {
            w.put<float>(static_cast<float>(v));
        }
    } else {
        w.put<int32_t>(static_cast<int32_t>(values.size()));
        for (const double v : values) {
            w.put<double>(v);
        }
    }
}

void putExpected(ByteWriter& w, int32_t version, ExpectedValueCalculation& ev, bool withType) {
    ev.computeDensity();
    if (withType) {
        w.cstr(ev.type());
    }
    w.cstr(ev.unit());
    w.put<int32_t>(ev.gridSize());
    putValues(w, version, ev.densityAvg());
    w.put<int32_t>(static_cast<int32_t>(ev.chrScaleFactors().size()));
    for (const auto& [chr, factor] : ev.chrScaleFactors()) {
        w.put<int32_t>(chr);
        if (version > 8) {
            w.put<float>(static_cast<float>(factor));
        } else {
            w.put<double>(factor);
        }
    }
}

void validateNormalizations(const std::vector<std::string>& normalizations) {
    std::set<std::string> seen;
    for (const auto& norm : normalizations) {
        if (std::find(kNormOrder.begin(), kNormOrder.end(), norm) == kNormOrder.end()) {
            throw HicError("unsupported normalization " + norm + "; hicfilecpp computes VC, VC_SQRT, KR and SCALE");
        }
        if (!seen.insert(norm).second) {
            throw HicError("normalization " + norm + " is listed twice");
        }
    }
}

// NormalizationVectorUpdater: the vectors of every chromosome at every
// resolution and the normalized expected values they give.
class NormBuilder {
public:
    using Provider = std::function<std::vector<double>(const std::string&, int32_t, int32_t)>;
    // One zoom level of the file, in the order Juicer's
    // Dataset.getAllPossibleResolutions returns them: every base pair zoom,
    // then every fragment zoom. `lengths` are the chromosome lengths for a
    // base pair zoom and the fragment counts the file holds (the site counts,
    // which is what DatasetReaderV2 puts in the fragmentCountMap) for a
    // fragment one.
    struct Zoom {
        std::string unit = "BP";
        int32_t resolution = 0;
        std::vector<int64_t> lengths;
        Provider provider;
    };

    NormBuilder(int32_t version, std::vector<Zoom> zooms, const std::vector<std::string>& requested,
                std::vector<std::string> provided = {})
        : version_(version), zooms_(std::move(zooms)), names_(kNormOrder), provided_(std::move(provided)) {
        for (const auto& type : kNormOrder) {
            if (std::find(requested.begin(), requested.end(), type) != requested.end()) {
                types_.push_back(type);
            }
        }
        // Provided types are indexed after the computed ones, in the order given.
        for (const auto& name : provided_) {
            if (std::find(names_.begin(), names_.end(), name) == names_.end()) {
                names_.push_back(name);
            }
        }
    }

    bool active() const { return !types_.empty() || !provided_.empty(); }

    // vectorLength is the bin count of the matrix's grid axis. For a base pair
    // zoom Juicer takes it from the block layout (MatrixZoomData:
    // blockBinCount times blockColumnCount), not from the chromosome length;
    // for a fragment zoom the axis is a HiCFragmentAxis, whose bin count is
    // the chromosome's site count divided by the resolution plus one.
    void addChromosome(int32_t chrIndex, int32_t zoom, const std::vector<NormRecord>& records,
                       int64_t vectorLength) {
        if (records.empty() || !active()) {
            return;
        }
        if (version_ > 8) {
            compute<float>(chrIndex, zoom, records, vectorLength);
            provide<float>(chrIndex, zoom, records, vectorLength);
        } else {
            compute<double>(chrIndex, zoom, records, vectorLength);
            provide<double>(chrIndex, zoom, records, vectorLength);
        }
    }

    // Writes the normalized expected values, the vector index and the vectors
    // at the output's position; returns the index's position and size.
    std::pair<int64_t, int64_t> write(OutputFile& out) {
        ByteWriter e;
        int32_t count = 0;
        for (auto& [key, ev] : expected_) {
            count += ev.hasData() ? 1 : 0;
        }
        e.put<int32_t>(count);
        for (auto& [key, ev] : expected_) {
            if (ev.hasData()) {
                putExpected(e, version_, ev, true);
            }
        }
        out.write(e);
        std::stable_sort(vectors_.begin(), vectors_.end(), [](const StoredVector& a, const StoredVector& b) {
            return std::tie(a.zoomOrder, a.chr, a.typeOrder) < std::tie(b.zoomOrder, b.chr, b.typeOrder);
        });
        const int64_t indexPosition = out.position();
        int64_t indexSize = 4;
        for (const auto& v : vectors_) {
            indexSize += static_cast<int64_t>(names_[static_cast<size_t>(v.typeOrder)].size()) + 1 + 4 +
                         static_cast<int64_t>(zooms_[static_cast<size_t>(v.zoomOrder)].unit.size()) + 1 + 4 + 8 +
                         (version_ > 8 ? 8 : 4);
        }
        ByteWriter index;
        index.put<int32_t>(static_cast<int32_t>(vectors_.size()));
        int64_t position = indexPosition + indexSize;
        for (const auto& v : vectors_) {
            const auto n = static_cast<int64_t>(v.values.size());
            const int64_t size = version_ > 8 ? 8 + 4 * n : 4 + 8 * n;
            index.cstr(names_[static_cast<size_t>(v.typeOrder)]);
            index.put<int32_t>(v.chr);
            index.cstr(zooms_[static_cast<size_t>(v.zoomOrder)].unit);
            index.put<int32_t>(zooms_[static_cast<size_t>(v.zoomOrder)].resolution);
            index.put<int64_t>(position);
            if (version_ > 8) {
                index.put<int64_t>(size);
            } else {
                index.put<int32_t>(static_cast<int32_t>(size));
            }
            position += size;
        }
        out.write(index);
        for (const auto& v : vectors_) {
            ByteWriter w;
            putValues(w, version_, v.values);
            out.write(w);
        }
        return {indexPosition, indexSize};
    }

private:
    struct StoredVector {
        int32_t zoomOrder = 0;
        int32_t chr = 0;
        int32_t typeOrder = 0;
        std::vector<double> values;
    };

    int32_t typeOrder(const std::string& type) const {
        return static_cast<int32_t>(std::find(names_.begin(), names_.end(), type) - names_.begin());
    }
    bool wants(const std::string& type) const {
        return std::find(types_.begin(), types_.end(), type) != types_.end();
    }
    ExpectedValueCalculation& expected(int32_t zoom, const std::string& type) {
        const auto key = std::make_pair(zoom, typeOrder(type));
        auto it = expected_.find(key);
        if (it == expected_.end()) {
            const Zoom& z = zooms_[static_cast<size_t>(zoom)];
            it = expected_.emplace(key, ExpectedValueCalculation(z.lengths, z.resolution, type, z.unit)).first;
        }
        return it->second;
    }

    template <class T>
    void compute(int32_t chr, int32_t zoom, const std::vector<NormRecord>& records, int64_t size) {
        const bool vc = wants("VC");
        const bool vcSqrt = wants("VC_SQRT");
        if (vc || vcSqrt) {
            std::vector<T> values = computeVC<T>(records, size);
            std::vector<T> roots(values.size(), T(0));
            if (vcSqrt) {
                for (size_t i = 0; i < values.size(); ++i) {
                    roots[i] = static_cast<T>(std::sqrt(static_cast<double>(values[i])));
                }
            }
            if (vc) {
                update<T>(chr, zoom, "VC", std::move(values), records);
            }
            if (vcSqrt) {
                update<T>(chr, zoom, "VC_SQRT", std::move(roots), records);
            }
        }
        if (wants("KR")) {
            update<T>(chr, zoom, "KR", computeKR<T>(records, size, version_ > 8), records);
        }
        if (wants("SCALE")) {
            std::vector<T> scale = computeScale<T>(records, size);
            if (!scale.empty()) {
                update<T>(chr, zoom, "SCALE", std::move(scale), records);
            }
        }
    }

    // The caller's vectors, stored as given (writeHicFile's chromosome indexes
    // start at 1 for the whole-genome pseudo-chromosome; the caller's at 0).
    template <class T>
    void provide(int32_t chr, int32_t zoom, const std::vector<NormRecord>& records, int64_t size) {
        const Zoom& z = zooms_[static_cast<size_t>(zoom)];
        if (!z.provider) {
            return;
        }
        for (const auto& name : provided_) {
            const std::vector<double> given = z.provider(name, chr - 1, z.resolution);
            if (given.empty()) {
                continue;
            }
            const bool vcType = name == "VC" || name == "VC_SQRT" ||
                                (name.size() > 3 && name.compare(name.size() - 3, 3, "_VC") == 0) ||
                                (name.size() > 8 && name.compare(name.size() - 8, 8, "_VC_SQRT") == 0);
            const T pad = vcType ? T(0) : std::numeric_limits<T>::quiet_NaN();
            std::vector<T> vector(static_cast<size_t>(std::max<int64_t>(size, 0)), pad);
            const size_t n = std::min(vector.size(), given.size());
            for (size_t i = 0; i < n; ++i) {
                vector[i] = static_cast<T>(given[i]);
            }
            store<T>(chr, zoom, name, std::move(vector), records);
        }
    }

    // NormalizationVectorUpdater.updateExpectedValueCalculationForChr
    template <class T>
    void update(int32_t chr, int32_t zoom, const std::string& type, std::vector<T> vector,
                const std::vector<NormRecord>& records) {
        const double factor = sumFactor(records, vector);
        for (auto& value : vector) {
            value = static_cast<T>(value * factor);
        }
        store<T>(chr, zoom, type, std::move(vector), records);
    }

    template <class T>
    void store(int32_t chr, int32_t zoom, const std::string& type, std::vector<T> vector,
               const std::vector<NormRecord>& records) {
        addDistancesFromRecords(expected(zoom, type), chr, records, vector);
        StoredVector stored;
        stored.zoomOrder = zoom;
        stored.chr = chr;
        stored.typeOrder = typeOrder(type);
        stored.values.assign(vector.begin(), vector.end());
        vectors_.push_back(std::move(stored));
    }

    int32_t version_;
    std::vector<Zoom> zooms_;
    std::vector<std::string> types_;
    std::vector<std::string> names_;
    std::vector<std::string> provided_;
    std::map<std::pair<int32_t, int32_t>, ExpectedValueCalculation> expected_;
    std::vector<StoredVector> vectors_;
};

void validate(const WriteOptions& options) {
    if (options.version == 6 || options.version == 7) {
        throw HicError("writing .hic version " + std::to_string(options.version) +
                       " is not supported: no Juicer tools release that writes it can be obtained to "
                       "validate against; hicfilecpp reads versions 6 to 9 and writes 8 and 9");
    }
    if (options.version != 8 && options.version != 9) {
        throw HicError("version must be 8 or 9, got " + std::to_string(options.version));
    }
    if (options.chromosomes.empty()) {
        throw HicError("no chromosomes given");
    }
    std::set<std::string> names;
    for (const auto& [name, length] : options.chromosomes) {
        std::string lower = name;
        std::transform(lower.begin(), lower.end(), lower.begin(), [](unsigned char c) { return std::tolower(c); });
        if (name.empty() || lower == "all") {
            throw HicError("invalid chromosome name '" + name + "'");
        }
        if (!names.insert(name).second) {
            throw HicError("chromosome " + name + " is listed twice");
        }
        if (length <= 0 || length > INT32_MAX) {
            throw HicError("chromosome " + name + " has length " + std::to_string(length) +
                           "; lengths must lie in 1 to 2147483647");
        }
    }
    if (options.resolutions.empty()) {
        throw HicError("no resolutions given");
    }
    std::set<int32_t> resolutions;
    for (const int32_t resolution : options.resolutions) {
        if (resolution <= 0 || !resolutions.insert(resolution).second) {
            throw HicError("resolutions must be positive and distinct");
        }
    }
    std::set<int32_t> fragResolutions;
    for (const int32_t resolution : options.fragResolutions) {
        if (resolution <= 0 || !fragResolutions.insert(resolution).second) {
            throw HicError("fragResolutions must be positive and distinct");
        }
    }
    if (!options.fragResolutions.empty()) {
        if (options.fragmentSites.size() != options.chromosomes.size()) {
            throw HicError("fragmentSites holds " + std::to_string(options.fragmentSites.size()) +
                           " entries but there are " + std::to_string(options.chromosomes.size()) +
                           " chromosomes; fragResolutions needs one site list per chromosome");
        }
        for (size_t i = 0; i < options.fragmentSites.size(); ++i) {
            const auto& sites = options.fragmentSites[i];
            for (size_t k = 0; k < sites.size(); ++k) {
                if (sites[k] < 0 || (k > 0 && sites[k] < sites[k - 1])) {
                    throw HicError("the restriction sites of " + options.chromosomes[i].first +
                                   " must be non-negative and in ascending order");
                }
            }
        }
    } else if (!options.fragmentSites.empty()) {
        throw HicError("fragmentSites is given without fragResolutions");
    }
    if (!options.sourceProvidesEveryResolution) {
        if (options.sourceResolution <= 0) {
            throw HicError("sourceResolution must be positive");
        }
        for (const int32_t resolution : options.resolutions) {
            if (resolution % options.sourceResolution != 0) {
                throw HicError("resolution " + std::to_string(resolution) + " is not a multiple of the source resolution " +
                               std::to_string(options.sourceResolution));
            }
        }
        if (!options.fragResolutions.empty()) {
            if (options.sourceFragResolution <= 0) {
                throw HicError("sourceFragResolution must be positive");
            }
            for (const int32_t resolution : options.fragResolutions) {
                if (resolution % options.sourceFragResolution != 0) {
                    throw HicError("fragment resolution " + std::to_string(resolution) +
                                   " is not a multiple of the source fragment resolution " +
                                   std::to_string(options.sourceFragResolution));
                }
            }
        }
    }
    validateNormalizations(options.normalizations);
    std::set<std::string> provided;
    for (const auto& name : options.providedNormalizations) {
        if (name.empty() || name == "NONE" ||
            std::any_of(name.begin(), name.end(), [](unsigned char c) { return c <= ' ' || c >= 127; })) {
            throw HicError("invalid provided normalization name '" + name + "'");
        }
        if (!provided.insert(name).second) {
            throw HicError("provided normalization " + name + " is listed twice");
        }
        if (std::find(options.normalizations.begin(), options.normalizations.end(), name) !=
            options.normalizations.end()) {
            throw HicError("normalization " + name + " is both computed and provided");
        }
    }
    if (!options.providedNormalizations.empty() && !options.normVector) {
        throw HicError("providedNormalizations needs a normVector function");
    }
    if (options.threads < 1) {
        throw HicError("threads must be at least 1");
    }
    if (options.compressionLevel < -1 || options.compressionLevel > 9) {
        throw HicError("compressionLevel must lie in -1 to 9");
    }
}

}  // namespace

}  // namespace detail

void writeHicFile(const std::string& path, const WriteOptions& options, PixelSource& source) {
    using namespace detail;
    validate(options);
    const int32_t version = options.version;
    std::vector<int32_t> resolutions = options.resolutions;
    std::sort(resolutions.begin(), resolutions.end(), std::greater<>());
    // Preprocessor.setResolutions sorts the fragment resolutions of -r
    // descending as well, and writes them after the base pair ones.
    std::vector<int32_t> fragResolutions = options.fragResolutions;
    std::sort(fragResolutions.begin(), fragResolutions.end(), std::greater<>());
    const bool frag = !fragResolutions.empty();
    const auto nChromosomes = static_cast<int32_t>(options.chromosomes.size());

    std::vector<int64_t> lengths(static_cast<size_t>(nChromosomes) + 1);
    std::vector<int64_t> offsets(static_cast<size_t>(nChromosomes) + 1, 0);
    int64_t genomeLength = 0;
    for (int32_t i = 0; i < nChromosomes; ++i) {
        offsets[static_cast<size_t>(i) + 1] = genomeLength;
        lengths[static_cast<size_t>(i) + 1] = options.chromosomes[static_cast<size_t>(i)].second;
        genomeLength += options.chromosomes[static_cast<size_t>(i)].second;
    }
    const auto genomeLengthKb = static_cast<int32_t>(genomeLength / 1000);
    lengths[0] = genomeLengthKb;

    // The site count of each chromosome, "All" first with none, as the header
    // stores them. The matrix layout and the normalization vectors count
    // fragments with this (FragmentCalculation.getNumberFragments), while the
    // raw expected values count them with one more: Preprocessor builds its
    // fragmentCountMap from sites.length + 1, and DatasetReaderV2, which the
    // normalization step reads the file back with, from sites.length.
    std::vector<int64_t> siteCounts(static_cast<size_t>(nChromosomes) + 1, 0);
    std::vector<int64_t> siteCountsPlusOne(static_cast<size_t>(nChromosomes) + 1, 0);
    if (frag) {
        for (int32_t i = 0; i < nChromosomes; ++i) {
            const auto count = static_cast<int64_t>(options.fragmentSites[static_cast<size_t>(i)].size());
            siteCounts[static_cast<size_t>(i) + 1] = count;
            siteCountsPlusOne[static_cast<size_t>(i) + 1] = count + 1;
        }
    }

    OutputFile out(path, std::nullopt);
    ByteWriter header;
    header.cstr("HIC");
    header.put<int32_t>(version);
    const auto masterOffset = static_cast<int64_t>(header.size());
    header.put<int64_t>(0);
    header.cstr(options.genomeId);
    const auto nviOffset = static_cast<int64_t>(header.size());
    if (version > 8) {
        header.put<int64_t>(0);
        header.put<int64_t>(0);
    }
    header.put<int32_t>(static_cast<int32_t>(options.attributes.size()) + 1);
    header.cstr("software");
    header.cstr(options.software.empty() ? std::string("hicfilecpp ") + kVersion : options.software);
    for (const auto& [key, value] : options.attributes) {
        header.cstr(key);
        header.cstr(value);
    }
    header.put<int32_t>(nChromosomes + 1);
    header.cstr("All");
    if (version > 8) {
        header.put<int64_t>(genomeLengthKb);
    } else {
        header.put<int32_t>(genomeLengthKb);
    }
    for (const auto& [name, length] : options.chromosomes) {
        header.cstr(name);
        if (version > 8) {
            header.put<int64_t>(length);
        } else {
            header.put<int32_t>(static_cast<int32_t>(length));
        }
    }
    header.put<int32_t>(static_cast<int32_t>(resolutions.size()));
    for (const int32_t resolution : resolutions) {
        header.put<int32_t>(resolution);
    }
    // Preprocessor.writeHeader: the fragment resolutions, then the restriction
    // sites of every chromosome in header order, "All" included with none.
    header.put<int32_t>(static_cast<int32_t>(fragResolutions.size()));
    for (const int32_t resolution : fragResolutions) {
        header.put<int32_t>(resolution);
    }
    if (frag) {
        header.put<int32_t>(0);  // "All" is absent from the sites file
        for (const auto& sites : options.fragmentSites) {
            header.put<int32_t>(static_cast<int32_t>(sites.size()));
            for (const int32_t site : sites) {
                header.put<int32_t>(site);
            }
        }
    }
    out.write(header);

    ThreadPool pool(options.threads);
    // The zoom levels of the file, in the order Juicer holds them: every base
    // pair resolution, then every fragment resolution.
    struct ZoomSpec {
        std::string unit;
        int32_t resolution = 0;
        int32_t zoomIndex = 0;
        bool frag = false;
    };
    std::vector<ZoomSpec> zoomSpecs;
    for (size_t z = 0; z < resolutions.size(); ++z) {
        zoomSpecs.push_back(ZoomSpec{"BP", resolutions[z], static_cast<int32_t>(z), false});
    }
    for (size_t z = 0; z < fragResolutions.size(); ++z) {
        zoomSpecs.push_back(ZoomSpec{"FRAG", fragResolutions[z], static_cast<int32_t>(z), true});
    }
    std::vector<ExpectedValueCalculation> expectedNone;
    std::vector<NormBuilder::Zoom> normZooms;
    for (const auto& spec : zoomSpecs) {
        expectedNone.emplace_back(spec.frag ? siteCountsPlusOne : lengths, spec.resolution, "NONE", spec.unit);
        normZooms.push_back(NormBuilder::Zoom{spec.unit, spec.resolution, spec.frag ? siteCounts : lengths,
                                             spec.frag ? options.fragNormVector : options.normVector});
    }
    NormBuilder norms(version, normZooms, options.normalizations, options.providedNormalizations);

    const Layout wholeLayout = wholeGenomeLayout(version, genomeLengthKb);
    const int64_t wholeBins = genomeLengthKb / wholeLayout.binSize + 1;
    std::vector<float> wholeCells(static_cast<size_t>(wholeBins * wholeBins), 0.0f);
    std::vector<uint8_t> wholePresent(wholeCells.size(), 0);
    double wholeSum = 0;

    std::vector<std::pair<std::string, MatrixEntry>> matrices;
    std::vector<Entry> entries;
    std::vector<NormRecord> records;
    for (int32_t c1 = 0; c1 < nChromosomes; ++c1) {
        for (int32_t c2 = c1; c2 < nChromosomes; ++c2) {
            const bool intra = c1 == c2;
            const int64_t length1 = lengths[static_cast<size_t>(c1) + 1];
            const int64_t length2 = lengths[static_cast<size_t>(c2) + 1];
            std::vector<ZoomWritten> zooms;
            for (size_t z = 0; z < zoomSpecs.size(); ++z) {
                const ZoomSpec& spec = zoomSpecs[z];
                const int32_t resolution = spec.resolution;
                const int32_t sourceStep = spec.frag ? options.sourceFragResolution : options.sourceResolution;
                const int32_t request = options.sourceProvidesEveryResolution ? resolution : sourceStep;
                // A fragment zoom takes the site counts where a base pair one
                // takes the chromosome lengths (MatrixPP's fragment branch).
                const int64_t extent1 = spec.frag ? siteCounts[static_cast<size_t>(c1) + 1] : length1;
                const int64_t extent2 = spec.frag ? siteCounts[static_cast<size_t>(c2) + 1] : length2;
                const Layout layout = pairLayout(version, extent1, extent2, resolution, intra, spec.frag);
                // The whole-genome matrix is a base pair matrix, so it is
                // filled on the finest base pair pass.
                const bool wholeGenomePass = !spec.frag && z + 1 == resolutions.size();
                const bool firstOfUnit = spec.frag ? z == resolutions.size() : z == 0;
                ExpectedValueCalculation& ev = expectedNone[z];
                entries.clear();
                double sum = 0;
                const auto consume = [&](const Pixel* pixels, size_t count) {
                    for (size_t k = 0; k < count; ++k) {
                        const Pixel& p = pixels[k];
                        const int64_t pos1 = static_cast<int64_t>(p.bin1) * request;
                        const int64_t pos2 = static_cast<int64_t>(p.bin2) * request;
                        // A fragment number may equal the site count: the
                        // position past the last site lies on one fragment more.
                        const bool outside = spec.frag ? pos1 > extent1 || pos2 > extent2
                                                       : pos1 >= extent1 || pos2 >= extent2;
                        if (p.bin1 < 0 || p.bin2 < 0 || outside) {
                            throw HicError("pixel (" + std::to_string(p.bin1) + ", " + std::to_string(p.bin2) +
                                           ") at " + std::to_string(request) + " " + spec.unit +
                                           " lies outside chromosomes " +
                                           options.chromosomes[static_cast<size_t>(c1)].first + " and " +
                                           options.chromosomes[static_cast<size_t>(c2)].first);
                        }
                        if (!std::isfinite(p.count)) {
                            throw HicError("pixel count is not finite");
                        }
                        const float score = p.count;
                        auto x = static_cast<int32_t>(pos1 / resolution);
                        auto y = static_cast<int32_t>(pos2 / resolution);
                        sum += score;
                        if (intra) {
                            if (x > y) {
                                std::swap(x, y);
                            }
                            if (x != y) {
                                sum += score;
                            }
                            ev.addDistance(c1 + 1, x, y, score);
                        }
                        entries.push_back(Entry{blockNumber(layout, x, y), x, y, score});
                        if (wholeGenomePass) {
                            // Preprocessor.getGenomicPosition, in kilobases
                            const auto g1 = static_cast<int32_t>((offsets[static_cast<size_t>(c1) + 1] + pos1) / 1000);
                            const auto g2 = static_cast<int32_t>((offsets[static_cast<size_t>(c2) + 1] + pos2) / 1000);
                            wholeSum += score;
                            int32_t wx = g1 / wholeLayout.binSize;
                            int32_t wy = g2 / wholeLayout.binSize;
                            if (wx > wy) {
                                std::swap(wx, wy);
                            }
                            if (wx != wy) {
                                wholeSum += score;
                            }
                            const auto cell = static_cast<size_t>(wx * wholeBins + wy);
                            if (wholePresent[cell] != 0) {
                                wholeCells[cell] += score;
                            } else {
                                wholePresent[cell] = 1;
                                wholeCells[cell] = score;
                            }
                        }
                    }
                };
                if (spec.frag) {
                    source.fragPixels(request, c1, c2, consume);
                } else {
                    source.pixels(request, c1, c2, consume);
                }
                if (entries.empty()) {
                    if (firstOfUnit) {
                        // No pixels of this unit for the pair: skip its zooms,
                        // as Juicer writes no matrix for a pair without
                        // contacts.
                        z = spec.frag ? zoomSpecs.size() : resolutions.size() - 1;
                        continue;
                    }
                    throw HicError("the source gave pixels for " + options.chromosomes[static_cast<size_t>(c1)].first +
                                   " and " + options.chromosomes[static_cast<size_t>(c2)].first +
                                   " at some " + spec.unit + " resolutions only");
                }
                const bool keep = intra && norms.active();
                records.clear();
                ZoomWritten zoom;
                zoom.unit = spec.unit;
                zoom.zoomIndex = spec.zoomIndex;
                zoom.sum = sum;
                zoom.layout = layout;
                zoom.index = writeBlocks(out, pool, version, options.compressionLevel, entries,
                                         keep ? &records : nullptr);
                zooms.push_back(std::move(zoom));
                if (keep) {
                    // The grid axis a normalization vector spans: the block
                    // layout for a base pair zoom, the fragment axis, one bin
                    // per resolution fragments, for a fragment one.
                    const int64_t vectorLength =
                        spec.frag ? siteCounts[static_cast<size_t>(c1) + 1] / resolution + 1
                                  : static_cast<int64_t>(layout.blockBinCount) * layout.blockColumnCount;
                    norms.addChromosome(c1 + 1, static_cast<int32_t>(z), records, vectorLength);
                }
            }
            if (!zooms.empty()) {
                matrices.emplace_back(std::to_string(c1 + 1) + "_" + std::to_string(c2 + 1),
                                      writeMatrixHeader(out, c1 + 1, c2 + 1, zooms));
            }
        }
    }

    // The whole-genome matrix, Preprocessor.computeWholeGenomeMatrix.
    entries.clear();
    for (int64_t x = 0; x < wholeBins; ++x) {
        for (int64_t y = x; y < wholeBins; ++y) {
            const auto cell = static_cast<size_t>(x * wholeBins + y);
            if (wholePresent[cell] != 0) {
                entries.push_back(Entry{blockNumber(wholeLayout, static_cast<int32_t>(x), static_cast<int32_t>(y)),
                                        static_cast<int32_t>(x), static_cast<int32_t>(y), wholeCells[cell]});
            }
        }
    }
    if (!entries.empty()) {
        ZoomWritten zoom;
        zoom.sum = wholeSum;
        zoom.layout = wholeLayout;
        zoom.index = writeBlocks(out, pool, version, options.compressionLevel, entries, nullptr);
        matrices.insert(matrices.begin(), {"0_0", writeMatrixHeader(out, 0, 0, {zoom})});
    }
    std::vector<Entry>().swap(entries);

    // Preprocessor.writeFooter
    const int64_t masterPosition = out.position();
    ByteWriter footer;
    if (version > 8) {
        footer.put<int64_t>(0);
    } else {
        footer.put<int32_t>(0);
    }
    const size_t counted = footer.size();
    footer.put<int32_t>(static_cast<int32_t>(matrices.size()));
    for (const auto& [key, entry] : matrices) {
        footer.cstr(key);
        footer.put<int64_t>(entry.position);
        footer.put<int32_t>(entry.size);
    }
    footer.put<int32_t>(static_cast<int32_t>(expectedNone.size()));
    for (auto& ev : expectedNone) {
        putExpected(footer, version, ev, false);
    }
    const int64_t nBytes = static_cast<int64_t>(footer.size() - counted);
    if (version > 8) {
        std::memcpy(footer.bytes.data(), &nBytes, sizeof(int64_t));
    } else {
        const auto nBytes32 = static_cast<int32_t>(nBytes);
        std::memcpy(footer.bytes.data(), &nBytes32, sizeof(int32_t));
    }
    out.write(footer);
    const auto [indexPosition, indexSize] = norms.write(out);
    if (version > 8) {
        out.patchValue<int64_t>(nviOffset, indexPosition);
        out.patchValue<int64_t>(nviOffset + 8, indexSize);
    }
    out.patchValue<int64_t>(masterOffset, masterPosition);
    out.close();
}

void addNorm(const std::string& path, const std::vector<std::string>& normalizations, int threads) {
    using namespace detail;
    validateNormalizations(normalizations);
    if (threads < 1) {
        throw HicError("threads must be at least 1");
    }
    int64_t truncateAt = 0;
    int32_t version = 0;
    int64_t nviOffset = 0;
    std::optional<NormBuilder> builder;
    {
        const HiCFile file(path);
        const FileState& state = stateOf(file);
        if (!state.fragResolutions.empty()) {
            throw HicError("addNorm supports base pair resolutions only; " + path + " has FRAG resolutions");
        }
        version = state.version;
        truncateAt = state.normalizedSectionPosition;
        nviOffset = 4 + 4 + 8 + static_cast<int64_t>(state.genome.size()) + 1;
        std::vector<int64_t> lengths;
        for (const auto& chromosome : state.chromosomes) {
            lengths.push_back(chromosome.length);
        }
        std::vector<NormBuilder::Zoom> zooms;
        for (const int32_t resolution : state.bpResolutions) {
            zooms.push_back(NormBuilder::Zoom{"BP", resolution, lengths, {}});
        }
        builder.emplace(version, zooms, normalizations);
        for (size_t z = 0; z < state.bpResolutions.size(); ++z) {
            const int32_t resolution = state.bpResolutions[z];
            for (size_t chr = 1; chr < state.chromosomes.size(); ++chr) {
                const auto index = static_cast<int32_t>(chr);
                const auto headers = file.matrixZoomHeaders(index, index);
                const auto header = std::find_if(headers.begin(), headers.end(), [&](const ZoomHeader& h) {
                    return h.unit == "BP" && h.binSize == resolution;
                });
                if (header == headers.end()) {
                    continue;
                }
                const auto mzd = file.getMatrixZoomData(state.chromosomes[chr].name, state.chromosomes[chr].name,
                                                        "observed", "NONE", "BP", resolution);
                std::vector<NormRecord> records;
                mzd.forEachBlock(
                    [&](const BlockIndexEntry&, std::vector<ContactRecord>& block) {
                        for (const auto& r : block) {
                            records.push_back(NormRecord{r.binX, r.binY, r.counts});
                        }
                    },
                    threads);
                builder->addChromosome(index, static_cast<int32_t>(z), records,
                                       static_cast<int64_t>(header->blockBinCount) * header->blockColumnCount);
            }
        }
    }
    OutputFile out(path, truncateAt);
    const auto [indexPosition, indexSize] = builder->write(out);
    if (version > 8) {
        out.patchValue<int64_t>(nviOffset, indexPosition);
        out.patchValue<int64_t>(nviOffset + 8, indexSize);
    }
    out.close();
}

}  // namespace hicfilecpp
