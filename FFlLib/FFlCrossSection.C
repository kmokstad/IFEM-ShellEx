// SPDX-FileCopyrightText: 2023 SAP SE
//
// SPDX-License-Identifier: Apache-2.0
//
// This file is part of FEDEM - https://openfedem.org
////////////////////////////////////////////////////////////////////////////////

#include "FFlLib/FFlIOAdaptors/FFlCrossSection.H"
#include "FFaLib/FFaDefinitions/FFaMsg.H"
#include "FFaLib/FFaAlgebra/FFaMath.H"
#include <functional>


/*!
  Formulas for the cross section parameters for the various section types
  can be found in the SIEMENS Element Library Reference documentation, e.g.,

  https://docs.plm.automation.siemens.com/data_services/resources/nxnastran/10/help/en_US/tdocExt/pdf/element.pdf

  See Section 8.1 Using Supplied Beam and Bar Libraries.
*/

FFlCrossSection::FFlCrossSection (const std::string& Type,
                                  const std::vector<double>& Dim)
{
  name = Type;

  // Lambda function for calculating x^2
  std::function<double(double)> pow2 = [](double x) -> double { return x*x; };
  // Lambda function for calculating x^3
  std::function<double(double)> pow3 = [](double x) -> double { return x*x*x; };

  if (Type == "ROD")
  {
    double R2 = pow2(Dim[0]);

    A   = M_PI*R2;
    Iyy = Izz = 0.25*A*R2;
    J   = Iyy + Izz;
    K1  = K2  = 0.9;
  }
  else if (Type == "TUBE")
  {
    double Ro2 = pow2(Dim[0]);
    double Ri2 = pow2(Dim[1]);

    A   =       M_PI*(Ro2 - Ri2);
    Iyy = Izz = M_PI*(pow2(Ro2) - pow2(Ri2))/4.0;
    J   = Iyy + Izz;
    K1  = K2  = 0.5;
  }
  else if (Type == "BAR")
  {
    double b = Dim[0];
    double h = Dim[1];

    A   = b*h;
    Iyy = pow2(b)*A/12.0;
    Izz = pow2(h)*A/12.0;
    if (h > b) std::swap(b,h);
    J   = pow2(h)*A*(1.0 - 0.63*(h/b)*(1.0-pow2(pow2(h/b))/12.0))/3.0;
    K1  = K2 = 5.0/6.0;
  }
  else if (Type == "BOX")
  {
    double b  = Dim[0];
    double h  = Dim[1];
    double t1 = Dim[2];
    double t2 = Dim[3];

    double bi = b - 2.0*t2;
    double hi = h - 2.0*t1;

    A   = b*h - bi*hi;
    Iyy = (h*pow3(b) - hi*pow3(bi))/12.0;
    Izz = (b*pow3(h) - bi*pow3(hi))/12.0;
    J   = 2.0*t2*t1*pow2((b-t2)*(h-t1))/(b*t2 + h*t1 - pow2(t2) - pow2(t1));
    K1  = 2.0*hi*t2/A;
    K2  = 2.0*bi*t1/A;
  }
  else if (Type == "I")
  {
    double a  = Dim[1];
    double b  = Dim[2];
    double tw = Dim[3];
    double ta = Dim[4];
    double tb = Dim[5];
    double hw = Dim[0] - (ta + tb);
    double hf = hw + 0.5*(ta + tb);

    A = a*ta + hw*tw + b*tb;

    // I-profile centroid location w.r.t lower flange centroid
    double ya = (0.5*(hw+ta)*hw*tw + hf*tb*b)/A;
    // Web centroid location w.r.t. I-profile centroid
    double yw = 0.5*(hw+ta) - ya;
    // Upper flange centroid location w.r.t. I-profile centroid
    double yb = hf - ya;

    Iyy = (ta*pow3(a) + hw*pow3(tw) + tb*pow3(b))/12.0;
    Izz = (a*pow3(ta) + tw*pow3(hw) + b*pow3(tb))/12.0
      +    pow2(ya)*a*ta + pow2(yw)*hw*tw + pow2(yb)*b*tb;
    J  = (a*pow3(ta) + hf*pow3(tw) + b*pow3(tb))/3.0;
    K1 = hw*tw/A;
    K2 = 5.0*(a*ta+b*tb)/(6.0*A);
    S1 = hf*tb*pow3(b)/(ta*pow3(a) + tb*pow3(b)) - ya;

#ifdef FFL_DEBUG
    std::cout <<"I-profile: a="<< a <<" b="<< b <<" h="<< hw+ta+tb
              <<" tw="<< tw <<" ta="<< ta <<" tb="<< tb
              <<"\n           Iyy="<< Iyy << " Izz="<< Izz
              <<"\n           A="<< A << " J="<< J << std::endl;
#endif
  }
  else if (Type == "T")
  {
    double bf = Dim[0];
    double tf = Dim[2];
    double tw = Dim[3];
    double hw = Dim[1] - tf;
    double hf = hw + 0.5*tf;

    A = bf*tf + hw*tw;

    // T-profile centroid location w.r.t. flange centroid
    double yf = 0.5*(hw+tf)*hw*tw/A;
    // Web centroid location w.r.t. T-profile centroid
    double yw = 0.5*(hw+tf) - yf;

    Iyy = (tf*pow3(bf) + hw*pow3(tw))/12.0;
    Izz = (bf*pow3(tf) + tw*pow3(hw))/12.0 + pow2(yw)*hw*tw + pow2(yf)*bf*tf;
    J   = (hf*pow3(tw) + bf*pow3(tf))/3.0;
    K1  = hw*tw/A;
    K2  = bf*tf/A;
    S1  = yf;

#ifdef FFL_DEBUG
    std::cout <<"T-profile: bf="<< bf <<" h="<< hw+tf
              <<" tf="<< tf <<" tw="<< tw
              <<"\n           Iyy="<< Iyy << " Izz="<< Izz
              <<"\n           A="<< A << " J="<< J << std::endl;
#endif
  }
  else if (Type == "L")
  {
    double b  = Dim[0];
    double h  = Dim[1];
    double t1 = Dim[2];
    double t2 = Dim[3];
    double h2 = h - t1;
    double b1 = b - t2;

    A = b*t1 + h2*t2;

    // L-profile centroid location w.r.t. outer (lower left) corner
    double yc = (b1*pow2(t1) + t2*pow2(h))*0.5/A;
    double zc = (t1*pow2(b) + h2*pow2(t2))*0.5/A;

    Iyy = (t1*pow3(b) + h2*pow3(t2))/3.0 - A*pow2(zc);
    Izz = (b1*pow3(t1) + t2*pow3(h))/3.0 - A*pow2(yc);
    Izy = (pow2(b*t1) + pow2(h*t2) - pow2(t1*t2))/4.0 - A*yc*zc;
    J   = ((b-0.5*t2)*pow3(t1) + (h-0.5*t1)*pow3(t2))/3.0;
    K1  = h2*t2/A;
    K2  = b1*t1/A;
    S1  = 0.5*t1 - yc;
    S2  = 0.5*t2 - zc;

#ifdef FFL_DEBUG
    std::cout <<"L-profile: b="<< b <<" h="<< h <<" t1="<< t1 <<" t2="<< t2
              <<"\n           yc="<< yc << " zc="<< zc
              <<"\n           Iyy="<< Iyy << " Izz="<< Izz << " Izy="<< Izy
              <<"\n           A="<< A << " J="<< J << std::endl;
#endif
  }
  else
    ListUI <<"\n *** FFlCrossSection: Type \"" << Type <<"\" is not supported."
           <<"\n            Replace it with a general cross section entry.\n";
}


double FFlCrossSection::findMainAxes ()
{
  if (Izy == 0.0) return 0.0;

  double fi = 0.5*atan2(-2.0*Izy,Iyy-Izz);
  double cf = cos(fi);
  double sf = sin(fi);
  double s2 = sin(2.0*fi);

  // Find the principal moments of inertia
  double Iy = Izz*sf*sf + Iyy*cf*cf - Izy*s2;
  double Iz = Izz*cf*cf + Iyy*sf*sf + Izy*s2;
#ifdef FFL_DEBUG
  std::cout <<"           phi="<< 180.0*fi/M_PI
            <<" I1="<< Iy <<" I2="<< Iz << std::endl;
#endif
  Iyy = Iy;
  Izz = Iz;
  Izy = 0.0;

  // Transform the shear centre offset
  s2 = cf*S2 - sf*S1;
  S1 = cf*S1 + sf*S2;
  S2 = s2;

  // Transform the shear stiffness factors
  s2 = cf*cf*K2 + sf*sf*K1;
  K1 = cf*cf*K1 + sf*sf*K2;
  K2 = s2;

  return 180.0*fi/M_PI;
}


std::ostream& operator<< (std::ostream& os, const FFlCrossSection& cs)
{
  os << cs.name <<":";
  if (fabs(cs.A) > 1.0e-16) os <<" A="<< cs.A;
  if (fabs(cs.Izz) > 1.0e-16) os <<" Izz="<< cs.Izz;
  if (fabs(cs.Iyy) > 1.0e-16) os <<" Iyy="<< cs.Iyy;
  if (fabs(cs.Izy) > 1.0e-16) os <<" Izy="<< cs.Izy;
  if (fabs(cs.J) > 1.0e-16) os <<" J="<< cs.J;
  if (fabs(cs.K1)+fabs(cs.K2) > 1.0e-16) os <<" K1="<< cs.K1 <<" K2="<< cs.K2;
  if (fabs(cs.S1)+fabs(cs.S2) > 1.0e-16) os <<" S1="<< cs.S1 <<" S2="<< cs.S2;
  if (fabs(cs.NSM) > 1.0e-16) os <<" NSM="<< cs.NSM;
  return os << std::endl;
}
