#ifndef __VALIDATE_INPUT_H__
#define __VALIDATE_INPUT_H__

#include <string>

bool ValidateAndRunTDSE(int argc, char **args, const std::string& input_file);
bool ValidateRunTISE(int argc, char **args, const std::string& input_file);

#endif