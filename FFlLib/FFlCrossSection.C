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


FFlCrossSection::FFlCrossSection (const std::string& Type,
                                  const std::vector<double>& Dim)
{
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
    J   = pow2(h)*A*(1.0 - 0.6*h/b)/3.0; // TODO: verify this
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
    double minThk = t1 < t2 ? t1 : t2;

    A   = b*h - bi*hi;
    Iyy = (h*pow3(b) - hi*pow3(bi))/12.0;
    Izz = (b*pow3(h) - bi*pow3(hi))/12.0;
    J   = 2.0*pow2(b*h)*minThk/(b+h); // TODO: verify this, probably assumes t1 << h and t2 << b
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

    A = a*ta + hw*tw + b*tb;

    double y1 = 0.5*ta;
    double y2 = y1 + 0.5*(ta+hw);
    double y3 = y2 + 0.5*(tb+hw);

    double ya = y1 - (a*ta*y1 + tw*hw*y2 + b*tb*y3)/A;
    double yw = ya + 0.5*(ta+hw);
    double yb = yw + 0.5*(tb+hw);

    Iyy = (ta*pow3(a) + hw*pow3(tw) + tb*pow3(b))/12.0;
    Izz = (a*pow3(ta) + tw*pow3(hw) + b*pow3(tb))/12.0
      +    pow2(ya)*a*ta + pow2(yw)*hw*tw + pow2(yb)*b*tb;
    J  = (a*pow3(ta) + hw*pow3(tw) + b*pow3(tb))/3.0;
    K1 = hw*tw/A;
    K2 = 5.0*(a*ta+b*tb)/(6.0*A);

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

    A = bf*tf + hw*tw;

    double yw = 0.5*hw - (0.5*tw*pow2(hw) + bf*(tf*hw + 0.5*pow2(tf)))/A;
    double yf = yw + 0.5*(hw+tf);

    Iyy = (tf*pow3(bf) + hw*pow3(tw))/12.0;
    Izz = (bf*pow3(tf) + tw*pow3(hw))/12.0 + pow2(yw)*hw*tw + pow2(yf)*bf*tf;
    J   = (hw*pow3(tw) + bf*pow3(tf))/3.0;
    K1  = hw*tw/A;
    K2  = bf*tf/A;

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


std::ostream& operator<< (std::ostream& os, const FFlCrossSection& cs)
{
  if (fabs(cs.A) > 1.0e-16) os <<" A="<< cs.A;
  if (fabs(cs.Izz) > 1.0e-16) os <<" Izz="<< cs.Izz;
  if (fabs(cs.Iyy) > 1.0e-16) os <<" Iyy="<< cs.Iyy;
  if (fabs(cs.Izy) > 1.0e-16) os <<" Izy="<< cs.Izy;
  if (fabs(cs.J) > 1.0e-16) os <<" J="<< cs.J;
  if (fabs(cs.K1)+fabs(cs.K2) > 1.0e-16) os <<" K1="<< cs.K1 <<" K2="<< cs.K2;
  if (fabs(cs.S1)+fabs(cs.S2) > 1.0e-16) os <<" S1="<< cs.S1 <<" S2="<< cs.S2;
  if (fabs(cs.NSM) > 1.0e-16) os <<" NSM="<< cs.NSM;
  return os;
}
