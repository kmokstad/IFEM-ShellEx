// $Id$
//==============================================================================
//!
//! \file WaveFunc.C
//!
//! \date Oct 31 2025
//!
//! \author Knut Morten Okstad / SINTEF
//!
//! \brief Wave function implementations.
//!
//==============================================================================

#include "WaveFunc.h"
#include <fstream>
#include <sstream>
#ifdef USE_OPENMP
#include <omp.h>
#endif


WaveSpectrum::WaveSpectrum (const char* file, double angle, double g, double t0)
  : alpha(angle), grav(g), time(t0)
{
  std::ifstream is(file);
  if (!is)
  {
    std::cerr <<"\n *** WaveSpectrum: Failed to open file "<< file << std::endl;
    return;
  }

  char temp[1024];
  while (is.good() && is.getline(temp,1024))
  {
    if (temp[0] == '#') continue;
    std::stringstream str(temp);
    int i;
    Component wc;
    str >> i >> wc.A >> wc.omega >> wc.eps;
    wave.emplace_back(wc);
  }

  alpha *= M_PI/180.0; // Convert to radians

#ifdef USE_OPENMP
  delta.resize(omp_get_max_threads());
#else
  delta.resize(1);
#endif
}


void WaveSpectrum::setParam (const std::string& name, double value)
{
#ifdef USE_OPENMP
  const size_t i = omp_get_thread_num();
#else
  const size_t i = 0;
#endif

  if (name == "ux")
    delta[i].x = value;
  else if (name == "uy")
    delta[i].y = value;
  else if (name == "uz")
    delta[i].z = value;
}


double WaveSpectrum::evaluate (const Vec3& X) const
{
  const Vec4* Xt = dynamic_cast<const Vec4*>(&X);
  const double t = Xt ? Xt->t : 0.0;
#ifdef USE_OPENMP
  const size_t i = omp_get_thread_num();
#else
  const size_t i = 0;
#endif

  double XoG = ((X.x+delta[i].x)*cos(alpha) +
                (X.y+delta[i].y)*sin(alpha)) / grav;

  double res = 0.0;
  for (const Component& w : wave)
    res += w.A*sin(w.omega*t - w.omega*XoG + w.eps);

  if (time > 0.0 && t < time)
    res *= t/time;

  return res;
}


HydroStaticPressure::HydroStaticPressure (RealFunc& h, double z0,
                                          double g, double d)
  : zeta(h), grav(g), rhow(d), w_line(z0)
{
#ifdef USE_OPENMP
  dZ.resize(omp_get_max_threads(),0.0);
#else
  dZ.resize(1,0.0);
#endif
}


void HydroStaticPressure::setParam (const std::string& name, double value)
{
  if (name == "uz")
  {
#ifdef USE_OPENMP
    const size_t i = omp_get_thread_num();
#else
    const size_t i = 0;
#endif
    dZ[i] = value;
  }
  zeta.setParam(name,value);
}


double HydroStaticPressure::evaluate (const Vec3& X) const
{
#ifdef USE_OPENMP
  const size_t i = omp_get_thread_num();
#else
  const size_t i = 0;
#endif
  const double depth = X.z + dZ[i] - (w_line + zeta(X));
  return depth < 0.0 ? -rhow*grav*depth : 0.0;
}
