#include <cmath>
#include "utils.h"

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
