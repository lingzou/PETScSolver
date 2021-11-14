#ifndef UTILS_H
#define UTILS_H

#include <vector>
#include <cmath>

static const double PI = asin(1) * 2;

// Some helper functions
void sysError(std::string);

enum TimeScheme
{
  BDF1    = 0,
  BDF2    = 1,
  CN      = 2,
  INVALID = 99
};

namespace UTILS
{
  // convert string to TimeScheme Enum
  TimeScheme StringToEnum(std::string str);

  // convert string to bool
  bool StringToBool(std::string str);

  std::string trim_file_name(std::string file_name_with_path);
  void linearReconstruction(double, double, std::vector<double> &, std::vector<double> &, std::vector<double> &);
  void linearReconstruction(double u_W, double u_E, double u_P, double &u_west, double &u_east);
}

#endif /*UTILS_H*/
