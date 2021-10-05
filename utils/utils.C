#include <cmath>
#include <iostream>
#include "utils.h"

void sysError(std::string message)
{ std::cerr << message << std::endl; exit(1); }

TimeScheme
UTILS::StringToEnum(std::string str)
{
  if      (str.compare("BDF1") == 0)  return BDF1;
  else if (str.compare("BDF2") == 0)  return BDF2;
  else if (str.compare("CN")   == 0)  return CN;
  else    {sysError("ERROR: UNKNOWN TimeScheme: " + str); return INVALID;}
}

bool
UTILS::StringToBool(std::string str)
{
  std::string temp_str = str;

  // Keep the input str for later error message, this temp_str transformed to be all upper case
  std::transform(temp_str.begin(), temp_str.end(), temp_str.begin(), ::toupper);

  if(temp_str.compare("TRUE") == 0)         return true;
  else if(temp_str.compare("FALSE") == 0)   return false;
  else    { sysError("Input string '" + str + "' cannot be converted as a boolean value."); }

  return false;
}

std::string
UTILS::trim_file_name(std::string full_file_name)
{
  unsigned int pos_of_point = full_file_name.find_last_of(".");
  unsigned int pos_of_slash = full_file_name.find_last_of("\\/");
  unsigned int length = pos_of_point - pos_of_slash - 1;
  if(length < 1)
  {
    sysError("ERROR: The input file name, '" + full_file_name + "', cannot be properly trimmed.");
    return std::string("");
  }
  else
    return full_file_name.substr(pos_of_slash + 1, length);
}

void
UTILS::linearReconstruction(double l_ghost, double r_ghost,
  std::vector<double> &u, std::vector<double> &u_w, std::vector<double> &u_e)
{
  for(int i = 0; i < u.size(); i++)
  {
    double u_P = u[i];
    double u_W = (i == 0) ? l_ghost : u[i - 1];
    double u_E = (i == u.size() - 1) ? r_ghost : u[i + 1];

    // I think this is a more decent implementation of TVD limiter
    // Reference:
    // [1] Berger, M., Aftosmis, M.J., Murman, S.M., 2005. Analysis of slope limiters on irreg- ular grids.
    //     In: AIAA Paper 2005-0490, 43rd AIAA Aerospace Sciences Meeting, Jan. 10e13, Reno, NV, 2005.
    double wf = 0.0;
    if ((u_E - u_P)*(u_P - u_W) > 0.0)    // u_P is in between u_W and u_E, i.e., not a local extremum
      if (std::fabs(u_E - u_W) > 1.e-10)  // The difference is large enough, so reconstruction is meaningful
      {
        double f = (u_P - u_W) / (u_E - u_W);
        wf = 2.0 * f * (1.0 - f) / ((1.0 - f) * (1.0 - f) + f * f);   // Eqn. (10) of Ref. [1]
      }
    // Linear reconstructed values
    u_w[i] = u_P - 0.25 * wf * (u_E - u_W);   // Eqn. (5) of Ref. [1]
    u_e[i] = u_P + 0.25 * wf * (u_E - u_W);
  }
}

void
UTILS::linearReconstruction(double u_W, double u_E, double u_P, double &u_west, double &u_east)
{
  double wf = 0.0;
  if ((u_E - u_P)*(u_P - u_W) > 0.0)    // u_P is in between u_W and u_E, i.e., not a local extremum
    if (std::fabs(u_E - u_W) > 1.e-10)  // The difference is large enough, so reconstruction is meaningful
    {
      double f = (u_P - u_W) / (u_E - u_W);
      wf = 2.0 * f * (1.0 - f) / ((1.0 - f) * (1.0 - f) + f * f);   // Eqn. (10) of Ref. [1]
    }
  // Linear reconstructed values
  u_west = u_P - 0.25 * wf * (u_E - u_W);   // Eqn. (5) of Ref. [1]
  u_east = u_P + 0.25 * wf * (u_E - u_W);
}
