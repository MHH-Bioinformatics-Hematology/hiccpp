#ifndef HICCPP_ERRORS_HPP
#define HICCPP_ERRORS_HPP

#include <stdexcept>
#include <string>

namespace hiccpp {

// Every failure hiccpp reports: an unreadable or malformed file, an
// unsupported version, an unknown chromosome, a missing zoom level or vector,
// or invalid writer input.
class HicError : public std::runtime_error {
public:
    using std::runtime_error::runtime_error;
};

}  // namespace hiccpp

#endif  // HICCPP_ERRORS_HPP
