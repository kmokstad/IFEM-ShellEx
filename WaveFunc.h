// $Id$
//==============================================================================
//!
//! \file WaveFunc.h
//!
//! \date Oct 31 2025
//!
//! \author Knut Morten Okstad / SINTEF
//!
//! \brief Wave function implementations.
//!
//==============================================================================

#ifndef WAVE_FUNC_H
#define WAVE_FUNC_H

#include "Function.h"


/*!
  \brief A scalar-valued spatial function, wave spectrum (sum of sines).
*/

class WaveSpectrum : public RealFunc
{
  //! \brief A struct describing a single wave component of a spectrum.
  struct Component
  {
    double A = 0.0;     //!< Amplitude
    double omega = 0.0; //!< Angular frequency
    double eps = 0.0;   //!< phase shift
  };

  std::vector<Component> wave;  //!< Wave component container
  double                 alpha; //!< Wave direction angle (w.r.t. global X-axis)
  double                 grav;  //!< Gravitation constant
  double                 time;  //!< Ramp-up time
  std::vector<Vec3>      delta; //!< Spatial point offsets

public:
  //! \brief Constructor initializing the function parameters from a file.
  //! \param[in] file Name of file to read wave spectrum from
  //! \param[in] angle Wave direction angle (w.r.t. to positive global X-axis)
  //! \param[in] g Gravitation constant
  //! \param[in] t0 Ramp-up time
  WaveSpectrum(const char* file, double angle, double g, double t0 = 0.0);

  //! \brief Returns whether the function is identically zero or not.
  bool isZero() const override { return wave.empty(); }
  //! \brief Returns whether the function is time-independent or not.
  bool isConstant() const override { return false; }

  //! \brief Sets an additional parameter in the function.
  void setParam(const std::string& name, double value) override;

protected:
  //! \brief Evaluates the function.
  double evaluate(const Vec3& X) const override;
};


/*!
  \brief A scalar-valued spatial function, hydrostatic pressure
*/

class HydroStaticPressure : public RealFunc
{
  RealFunc&           zeta;   //!< Spatial function describing the water surface
  double              grav;   //!< Gravitation constant
  double              rhow;   //!< Water mass density
  double              w_line; //!< Initial waterline z-position
  std::vector<double> dZ;     //!< Vertical point offsets

public:
  //! \brief Constructor initializing the function parameters.
  //! \param[in] h Sea elevation function
  //! \param[in] z0 Initial vertical waterline position
  //! \param[in] g Gravition constant
  //! \param[in] d Mass density of sea water
  HydroStaticPressure(RealFunc& h, double z0 = 0.0,
                      double g = 9.81, double d = 1000.0);

  //! \brief Returns whether the function is identically zero or not.
  bool isZero() const override { return zeta.isZero() && w_line == 0.0; }
  //! \brief Returns whether the function is time-independent or not.
  bool isConstant() const override { return false; }

  //! \brief Sets an additional parameter in the function.
  void setParam(const std::string& name, double value) override;

protected:
  //! \brief Evaluates the function.
  double evaluate(const Vec3& X) const override;
};

#endif
