/*
 * kmaplib_common.cpp
 *
 * This file contains the common utility functions that are shared across different
 * parts of the kinetic modeling and optimization process. These include functions
 * for vector norms, setting and getting kinetic parameters, and convolution operations.
 */

#include "kmaplib.h"
// #include "mex.h"
#include <cmath>
#include <cstdlib>
#include <algorithm>
#include <climits>
#include <limits>


//------------------------------------------------------------------------------
// frame
//------------------------------------------------------------------------------
// Averages the time activity curve (TAC) over the frame duration.
// This function calculates the average activity in the specified time frames.
void frame(double *scant, double td, double *c_t, int num_frm, int num_c, double *ct)
{
   int      m, k;
   double   tmp, dt;
   int      num_time;
   int      j, jstart, jend;

   if (!scant || !c_t || !ct || num_frm <= 0 || num_c <= 0 ||
       !std::isfinite(td) || td <= 0.0 ||
       !std::isfinite(scant[2*num_frm-1]) ||
       scant[2*num_frm-1] / td < 1.0 ||
       scant[2*num_frm-1] / td >= static_cast<double>(INT_MAX)) {
      // Caller must provide a valid, sampled time grid. Never divide by zero
      // or index outside a buffer on a malformed acquisition.
      if (ct && num_frm > 0 && num_c > 0)
         std::fill_n(ct, num_frm * num_c, std::numeric_limits<double>::quiet_NaN());
      return;
   }
   num_time = static_cast<int>(scant[2*num_frm-1] / td);

   for (m = 0; m < num_frm; m++) {
      const double start = scant[m];
      const double end = scant[m + num_frm];
      if (!std::isfinite(start) || !std::isfinite(end) || start < 0.0 ||
          end <= start || end > scant[2*num_frm-1] ||
          end / td >= static_cast<double>(INT_MAX)) {
         for (k = 0; k < num_c; ++k)
            ct[m + k * num_frm] = std::numeric_limits<double>::quiet_NaN();
         continue;
      }
      jstart = static_cast<int>(start / td);
      jend = static_cast<int>(end / td);
      if (jstart < 0 || jend > num_time) {
         for (k = 0; k < num_c; ++k)
            ct[m + k * num_frm] = std::numeric_limits<double>::quiet_NaN();
         continue;
      }
      dt = static_cast<double>(jend - jstart);
      if (jend == jstart) {
         // Frame shorter than the integration grid (including its final
         // partial cell): interpolate at its
         // midpoint instead of computing 0/0. This is a numerical sampling
         // approximation; for quantitative short frames use a finer td.
         const double position = std::max(0.0, std::min(
             static_cast<double>(num_time - 1), (start + end)/(2.0*td) - 0.5));
         const int left = static_cast<int>(position);
         const int right = std::min(num_time - 1, left + 1);
         const double fraction = position - left;
         for (k = 0; k < num_c; ++k) {
            ct[m + k * num_frm] = (1.0 - fraction) * c_t[left + k*num_time]
                + fraction * c_t[right + k*num_time];
         }
         continue;
      }
      for (k = 0; k < num_c; k++) {
         tmp = 0.0;
         for (j = jstart; j < jend; j++) {
            tmp += c_t[j + k * num_time];
         }
         ct[m + k * num_frm] = tmp / dt;
      }
   }
}

//------------------------------------------------------------------------------
// kconv_exp
//------------------------------------------------------------------------------
// Computes the convolution of an exponential function with the input function.
// c(t) = k1 * exp(-k2 * t) ⊗ u(t) where ⊗ represents the convolution operator.
void kconv_exp(double k1, double k2, double *u, int num_time, double td, double *c)
{
   int      i;
   double   ek2, prev, tmp;

   k1 *= td / 60.0;  // Convert td from seconds to minutes
   k2 *= td / 60.0;
   ek2 = exp(-k2);
   prev = 0;

   if (ek2 == 1.0) {  // Special case when k2 is 0
      for (i = 0; i < num_time; i++) {
         prev += u[i];
         c[i] = k1 * prev;
      }
   } else {
      tmp = (1.0 - ek2) / k2;
      for (i = 0; i < num_time; i++) {
         prev = prev * ek2 + k1 * tmp * u[i];
         c[i] = prev;
      }
   }
}

//------------------------------------------------------------------------------
// setkin
//------------------------------------------------------------------------------
// Assigns sensitive parameters from x to x0 based on the sensitivity indicators.
void setkin(double *x, int num_x, int *xsens, double *x0)
{
   int j = 0;
   int k = 0;
   while (j < num_x) {
      if (xsens[k] == 0)
         ++k;
      else {
         x0[k] = x[j];
         ++j;
         ++k;
      }
   }
}

//------------------------------------------------------------------------------
// getkin
//------------------------------------------------------------------------------
// Retrieves sensitive parameters from x0 to x based on the sensitivity indicators.
void getkin(double *x0, int *xsens, int num_x, double *x)
{
   int j = 0;
   int k = 0;
   while (j < num_x) {
      if (xsens[k] == 0)
         ++k;
      else {
         x[j] = x0[k];
         ++j;
         ++k;
      }
   }
}

//------------------------------------------------------------------------------
// vecnorm2
//------------------------------------------------------------------------------
// Computes the L2 norm (Euclidean norm) of vector x.
double vecnorm2(double *x, int num)
{
   int      i;
   double   temp = 0.0;

   for (i = 0; i < num; i++)
      temp += x[i] * x[i];
   return sqrt(temp);
}

//------------------------------------------------------------------------------
// vecnormw
//------------------------------------------------------------------------------
// Computes the weighted norm of vector x.
double vecnormw(double *w, double *x, int num)
{
   int      i;
   double   temp = 0.0;

   for (i = 0; i < num; i++)
      temp += w[i] * x[i] * x[i];
   return temp;
}
