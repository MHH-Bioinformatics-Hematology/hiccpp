#ifndef HICFILECPP_HARNESS_JSON_HPP
#define HICFILECPP_HARNESS_JSON_HPP

// A small JSON reader for the harness case files, and string quoting for the
// result documents the driver writes.

#include <cmath>
#include <cstdint>
#include <cstdio>
#include <stdexcept>
#include <string>
#include <utility>
#include <vector>

namespace harness {

struct Json {
    enum class Type { Null, Bool, Number, String, Array, Object };

    Type type = Type::Null;
    bool boolean = false;
    double number = 0.0;
    std::string string;
    std::vector<Json> array;
    std::vector<std::pair<std::string, Json>> object;

    bool has(const std::string& key) const {
        for (const auto& item : object) {
            if (item.first == key) {
                return true;
            }
        }
        return false;
    }
    const Json& operator[](const std::string& key) const {
        for (const auto& item : object) {
            if (item.first == key) {
                return item.second;
            }
        }
        throw std::runtime_error("case is missing key " + key);
    }
    const Json& operator[](size_t index) const { return array.at(index); }
    size_t size() const { return array.size(); }
    const std::string& str() const {
        if (type != Type::String) {
            throw std::runtime_error("expected a JSON string");
        }
        return string;
    }
    int64_t i64() const {
        if (type != Type::Number) {
            throw std::runtime_error("expected a JSON number");
        }
        return static_cast<int64_t>(std::llround(number));
    }
    double f64() const { return number; }
    bool truthy() const { return type == Type::Bool ? boolean : type != Type::Null; }

    static Json parse(const std::string& text) {
        size_t pos = 0;
        Json value = parseValue(text, pos);
        skipSpace(text, pos);
        if (pos != text.size()) {
            throw std::runtime_error("trailing characters in JSON");
        }
        return value;
    }

private:
    static void skipSpace(const std::string& t, size_t& p) {
        while (p < t.size() && (t[p] == ' ' || t[p] == '\n' || t[p] == '\r' || t[p] == '\t')) {
            ++p;
        }
    }
    static std::string parseString(const std::string& t, size_t& p) {
        if (t[p] != '"') {
            throw std::runtime_error("expected a string in JSON");
        }
        ++p;
        std::string out;
        while (p < t.size() && t[p] != '"') {
            char c = t[p++];
            if (c == '\\') {
                char e = t.at(p++);
                switch (e) {
                    case 'n': out += '\n'; break;
                    case 't': out += '\t'; break;
                    case 'r': out += '\r'; break;
                    case 'b': out += '\b'; break;
                    case 'f': out += '\f'; break;
                    case 'u': {
                        const unsigned code = static_cast<unsigned>(std::stoul(t.substr(p, 4), nullptr, 16));
                        p += 4;
                        if (code < 0x80) {
                            out += static_cast<char>(code);
                        } else if (code < 0x800) {
                            out += static_cast<char>(0xC0 | (code >> 6));
                            out += static_cast<char>(0x80 | (code & 0x3F));
                        } else {
                            out += static_cast<char>(0xE0 | (code >> 12));
                            out += static_cast<char>(0x80 | ((code >> 6) & 0x3F));
                            out += static_cast<char>(0x80 | (code & 0x3F));
                        }
                        break;
                    }
                    default: out += e; break;
                }
            } else {
                out += c;
            }
        }
        if (p >= t.size()) {
            throw std::runtime_error("unterminated JSON string");
        }
        ++p;
        return out;
    }
    static Json parseValue(const std::string& t, size_t& p) {
        skipSpace(t, p);
        if (p >= t.size()) {
            throw std::runtime_error("unexpected end of JSON");
        }
        Json v;
        const char c = t[p];
        if (c == '{') {
            v.type = Type::Object;
            ++p;
            skipSpace(t, p);
            if (t[p] == '}') {
                ++p;
                return v;
            }
            while (true) {
                skipSpace(t, p);
                std::string key = parseString(t, p);
                skipSpace(t, p);
                if (t.at(p++) != ':') {
                    throw std::runtime_error("expected ':' in JSON");
                }
                v.object.emplace_back(std::move(key), parseValue(t, p));
                skipSpace(t, p);
                const char next = t.at(p++);
                if (next == '}') {
                    return v;
                }
                if (next != ',') {
                    throw std::runtime_error("expected ',' in JSON object");
                }
            }
        }
        if (c == '[') {
            v.type = Type::Array;
            ++p;
            skipSpace(t, p);
            if (t[p] == ']') {
                ++p;
                return v;
            }
            while (true) {
                v.array.push_back(parseValue(t, p));
                skipSpace(t, p);
                const char next = t.at(p++);
                if (next == ']') {
                    return v;
                }
                if (next != ',') {
                    throw std::runtime_error("expected ',' in JSON array");
                }
            }
        }
        if (c == '"') {
            v.type = Type::String;
            v.string = parseString(t, p);
            return v;
        }
        if (t.compare(p, 4, "true") == 0) {
            v.type = Type::Bool;
            v.boolean = true;
            p += 4;
            return v;
        }
        if (t.compare(p, 5, "false") == 0) {
            v.type = Type::Bool;
            p += 5;
            return v;
        }
        if (t.compare(p, 4, "null") == 0) {
            p += 4;
            return v;
        }
        size_t used = 0;
        v.type = Type::Number;
        v.number = std::stod(t.substr(p), &used);
        p += used;
        return v;
    }
};

inline std::string quote(const std::string& text) {
    std::string out = "\"";
    for (const char c : text) {
        switch (c) {
            case '"': out += "\\\""; break;
            case '\\': out += "\\\\"; break;
            case '\n': out += "\\n"; break;
            case '\t': out += "\\t"; break;
            case '\r': out += "\\r"; break;
            default:
                if (static_cast<unsigned char>(c) < 0x20) {
                    char buffer[8];
                    std::snprintf(buffer, sizeof(buffer), "\\u%04x", static_cast<unsigned>(c));
                    out += buffer;
                } else {
                    out += c;
                }
        }
    }
    return out + "\"";
}

}  // namespace harness

#endif  // HICFILECPP_HARNESS_JSON_HPP
