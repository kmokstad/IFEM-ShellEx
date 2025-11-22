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
#include "Utilities.h"
#include "IFEM.h"
#include "tinyxml2.h"
#include <fstream>
#include <sstream>
#include <cstdlib>
#ifdef USE_OPENMP
#include <omp.h>
#endif


RealFunc* WaveSpectrum::parse (const tinyxml2::XMLElement* elem, double g)
{
  double Hs = 0.0, Tp = 1.0, Twm = 100.0, fmin = 0.1, fmax = 10.0;
  if (elem->FirstChild())
    IFEM::cout <<": file = "<< elem->FirstChild()->Value();
  else if (utl::getAttribute(elem,"Hs",Hs))
  {
    IFEM::cout <<": Hs = "<< Hs;
    if (utl::getAttribute(elem,"Tp",Tp))
      IFEM::cout <<" Tp = "<< Tp;
    if (utl::getAttribute(elem,"Twm",Twm))
      IFEM::cout <<" Twm = "<< Twm;
    if (utl::getAttribute(elem,"fmin",fmin))
      IFEM::cout <<" fmin = "<< fmin;
    if (utl::getAttribute(elem,"fmax",fmax))
      IFEM::cout <<" fmax = "<< fmax;
  }
  else
  {
    IFEM::cout <<" (invalid specification, ignored).\n";
    return nullptr;
  }

  double angle = 180.0, vs = 0.0, rampT = 0.0;
  if (utl::getAttribute(elem,"angle",angle))
    IFEM::cout <<", wave angle = "<< angle;
  if (utl::getAttribute(elem,"speed",vs))
    IFEM::cout <<", vessel speed = "<< vs;
  if (utl::getAttribute(elem,"ramp",rampT) && rampT > 0.0)
    IFEM::cout <<", ramp up time = "<< rampT;

  if (elem->FirstChild())
    return new WaveSpectrum(elem->FirstChild()->Value(),angle,vs,g,rampT);
  else
    return new PiersonMoskovitz(Hs,Tp,Twm,fmin,fmax,angle,vs,g,rampT);
}


WaveSpectrum::WaveSpectrum (const char* file,
                            double a, double v, double g, double t0)
  : alpha(a), vs(v), grav(g), time(t0)
{
  this->init();

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
    Component wc;
    str >> wc.id >> wc.Ampl >> wc.omega >> wc.eps;
    wave.emplace_back(wc);
  }
}


PiersonMoskovitz::PiersonMoskovitz (double Hs, double Tp, double Twm,
                                    double fmin, double fmax,
                                    double a, double v, double g, double t0)
  : WaveSpectrum(a,v,g,t0)
{
  this->init();

  const size_t Nspectr = (fmax - fmin) * Twm;
  const double twoPi   = 2.0 * M_PI;
  const double omega_p = twoPi / Tp;
  const double omega_4 = pow(omega_p,4.0);

  // Lambda function generating a Pierson-Moskowitz wave height component.
  auto S = [Hs,omega_4] (double omega)
  {
    return (0.3125 * Hs*Hs * omega_4 * pow(omega,-5.0) *
            exp(-1.25*pow(omega,-4.0)/omega_4));
  };

  wave.reserve(Nspectr);
  double omega = twoPi * fmin;
  double domega = twoPi / Twm;

  srand(0); // we want the same random series for each run
  for (size_t i = 0; i < Nspectr; i++, omega += domega)
  {
    double phi = twoPi * rand() / static_cast<double>(RAND_MAX);
    wave.emplace_back(sqrt(2.0 * S(omega) * domega), omega, phi);
  }
}


void WaveSpectrum::init ()
{
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
                (X.y+delta[i].y)*sin(alpha) + vs*t) / grav;

  double res = 0.0;
  for (const Component& w : wave)
    res += w.Ampl*sin(w.omega*t - w.omega*XoG + w.eps);

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


/*!
  The pressure is assumed in the direction of the outward-directed
  normal vector of the shell surface at the evaluation point.
  Therefore, it will always have a non-positive value.
*/

double HydroStaticPressure::evaluate (const Vec3& X) const
{
#ifdef USE_OPENMP
  const size_t i = omp_get_thread_num();
#else
  const size_t i = 0;
#endif
  const double depth = X.z + dZ[i] - (w_line + zeta(X));
  return depth < 0.0 ? rhow*grav*depth : 0.0;
}
