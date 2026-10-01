#ifndef HICCPP_DETAIL_CODEC_HPP
#define HICCPP_DETAIL_CODEC_HPP

// zlib streams, the compression of .hic blocks.

#include <cstddef>
#include <vector>

namespace hiccpp::detail {

// Inflates one zlib stream. Like hicstraw, a stream that ends early yields
// the bytes decoded so far.
void inflateZlib(const char* data, size_t size, std::vector<char>& out);

// Deflates into a zlib stream at the given level; the output depends only on
// the input, the level and the zlib version.
void deflateZlib(const char* data, size_t size, int level, std::vector<char>& out);

}  // namespace hiccpp::detail

#endif  // HICCPP_DETAIL_CODEC_HPP
