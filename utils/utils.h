#pragma once

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
  double bilinear(double x, double y, double x0, double x1, double y0, double y1, double z00, double z01, double z10, double z11);

  double BDF2Tran(double u, double u_o, double u_oo, double dt, double dt_o);
}
