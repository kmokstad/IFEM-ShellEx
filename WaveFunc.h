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

namespace tinyxml2 { class XMLElement; }


/*!
  \brief A scalar-valued spatial function, wave spectrum (sum of sines).
*/

class WaveSpectrum : public RealFunc
{
protected:
  //! \brief A struct describing a single wave component of a spectrum.
  struct Component
  {
    int    id;    //!< Component id
    double Ampl;  //!< Amplitude
    double omega; //!< Angular frequency
    double eps;   //!< phase shift

    //! \brief Default constructor.
    Component(double h = 0.0, double w = 0.0, double e = 0.0)
      : id(0), Ampl(h), omega(w), eps(e) {}
  };

  std::vector<Component> wave;  //!< Wave component container
  double                 alpha; //!< Wave direction angle (w.r.t. global X-axis)
  double                 vs;    //!< Vessel speed (against the wave direction)
  double                 grav;  //!< Gravitation constant
  double                 time;  //!< Ramp-up time
  std::vector<Vec3>      delta; //!< Spatial point offsets

  //! \brief The default constructor is protected, used by sub-classes only.
  WaveSpectrum(double a = 0.0, double v = 0.0, double g = 9.81, double t0 = 0.0)
    : alpha(a), vs(v), grav(g), time(t0) {}

public:
  //! \brief Constructor initializing the function parameters from a file.
  //! \param[in] file Name of file to read wave spectrum from
  //! \param[in] a Wave direction angle (w.r.t. to positive global X-axis)
  //! \param[in] v Vessel speed (against the wave direction)
  //! \param[in] g Gravitation constant
  //! \param[in] t0 Ramp-up time
  WaveSpectrum(const char* file, double a, double v, double g, double t0 = 0.0);

  //! \brief Returns whether the function is identically zero or not.
  bool isZero() const override { return wave.empty(); }
  //! \brief Returns whether the function is time-independent or not.
  bool isConstant() const override { return false; }

  //! \brief Sets an additional parameter in the function.
  void setParam(const std::string& name, double value) override;

  //! \brief Static method creating a WaveSpectrum instance.
  //! \param[in] elem XML-element to parse function parameters from.
  //! \param[in] g Gravitation constant
  static RealFunc* parse(const tinyxml2::XMLElement* elem, double g);

protected:
  //! \brief Common initializer (used by constructors only).
  void init();

  //! \brief Evaluates the sea elevation function.
  double evaluate(const Vec3& X) const override;
};


/*!
  \brief A scalar-valued spatial function, Pierson-Moskovitz wave spectrum.
*/

class PiersonMoskovitz : public WaveSpectrum
{
public:
  //! \brief The constructor realizes the Pierson-Moskowitz wave spectrum.
  //! \param[in] Hs Significant wave height
  //! \param[in] Tp Spectral peak period
  //! \param[in] Twm Total length of wave maker time span
  //! \param[in] fmin Minimum frequency to consider
  //! \param[in] fmax Maximum frequency to consider
  //! \param[in] a Wave direction angle (w.r.t. to positive global X-axis)
  //! \param[in] v Vessel speed (against the wave direction)
  //! \param[in] g Gravitation constant
  //! \param[in] t0 Ramp-up time
  PiersonMoskovitz(double Hs, double Tp, double Twm,
                   double fmin, double fmax,
                   double a, double v, double g, double t0);
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
  //! \brief Evaluates the hydrostatic pressure function.
  double evaluate(const Vec3& X) const override;
};

#endif
