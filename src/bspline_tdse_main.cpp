#include <iostream>
#include <unordered_map>
#include <functional>
#include "input/validate_input.h"


int main(int argc, char **args) {
    if (!ValidateAndRunTDSE(argc, args, "input.json"))
        return -1;
    return 0;
}