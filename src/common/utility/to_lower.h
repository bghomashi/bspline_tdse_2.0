#ifndef __TO_LOWER_H__
#define __TO_LOWER_H__

#include <string>
#include <algorithm>

inline std::string ToLower(std::string data) {
    std::transform(data.begin(), data.end(), data.begin(),
        [](unsigned char c){ return std::tolower(c); });
    return data;
}


#endif