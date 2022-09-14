#pragma once

#include "cl_gpu/cl_kernel.hpp"

auto sum_src =  "\n\n__kernel void sum(__global complex *a, int len) {\n" 
                "    int i = get_global_id(0);\n" 
                "    int size = get_global_size(0);\n" 
                "    if (i < size) {\n" \
                "        a[i] += a[len - i - 1];\n" 
                "    }\n" 
                "}\n\n\n";