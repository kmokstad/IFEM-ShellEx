// $Id$
//==============================================================================
//!
//! \file AndesShell.h
//!
//! \brief Class representing a linear elastic shell with rotational DOFs.
//!
//! \author Knut Morten Okstad / SINTEF
//!
//! \date Feb 15 2024
//!
//==============================================================================

#ifndef _ANDES_SHELL_H
#define _ANDES_SHELL_H

#include "ElasticBase.h"
#include <set>

class ASMu2DNastran;
class ScalarFunc;
class RealFunc;


/*!
  \brief Class representing the integrand of an ANDES shell.
*/

class AndesShell : public ElasticBase
{
public:
  //! \brief Default constructor.
  //! \param[in] n Number of consequtive solution vectors to reside in core
  //! \param[in] modal If \e true, a modal dynamics simulation is performed
  explicit AndesShell(unsigned short int n = 1, bool modal = false);
  //! \brief The destructor writes out the ignored bad elements.
  virtual ~AndesShell();

  using ElasticBase::parseMatProp;
  //! \brief Parses material properties from an XML-element.
  virtual Material* parseMatProp(const tinyxml2::XMLElement* elem, bool);

  //! \brief Prints out the problem definition to the log stream.
  virtual void printLog() const;

  //! \brief Checks if the point \a Xc is inside the thickness loss area or not.
  bool isInLossArea(const Vec3& Xc) const;

  //! \brief Defines the solution mode before the element assembly is started.
  //! \param[in] mode The solution mode to use
  virtual void setMode(SIM::SolutionMode mode);

  //! \brief Initializes and toggles the use of left-hand-side matrix buffers.
  //! \param[in] nEl Number of elements in the model/toggle.
  //! If larger than 1, element matrix buffers are allocated to given size.
  //! If equal to 1, element matrices are recomputed.
  //! If equal to 0, reuse buffered element matrices.
  virtual void initLHSbuffers(size_t nEl);

  //! \brief Initialization of integrand with patch-specific data.
  virtual void initForPatch(const ASMbase* pch);

  //! \brief Defines the pressure field.
  //! \param[in] pf Spatial function describing the pressure field
  //! \param[in] code Pressure function code
  //! \param[in] sName Name of element set the pressure applies to
  //! \param[in] pch The patch containing those elements
  bool setPressure(RealFunc* pf, int code,
                   const std::string& sName = "", const ASMbase* pch = nullptr);

  using ElasticBase::initIntegration;
  //! \brief Initializes the integrand with the number of integration points.
  //! \param[in] nGp Total number of interior integration points
  virtual void initIntegration(size_t nGp, size_t);

  using ElasticBase::getLocalIntegral;
  //! \brief Returns a local integral container for the given element.
  //! \param[in] nen Number of nodes on element
  virtual LocalIntegral* getLocalIntegral(size_t nen, size_t, bool) const;

  //! \brief Returns which integrand to be used.
  virtual int getIntegrandType() const;

  using ElasticBase::initElement;
  //! \brief Initializes current element for numerical integration.
  //! \param[in] MNPC Matrix of nodal point correspondance for current element
  //! \param[in] fe Nodal and integration point data for current element
  //! \param[in] Xc Cartesian coordinates of the element center
  //! \param elmInt The local integral object for current element
  virtual bool initElement(const std::vector<int>& MNPC,
                           const FiniteElement& fe, const Vec3& Xc, size_t,
                           LocalIntegral& elmInt);

  using ElasticBase::evalInt;
  //! \brief Evaluates the integrand at an interior point.
  //! \param elmInt The local integral object to receive the contributions
  //! \param[in] fe Finite element data of current integration point
  //! \param[in] X Cartesian coordinates of current integration point
  virtual bool evalInt(LocalIntegral& elmInt, const FiniteElement& fe,
                       const Vec3& X) const;

  using ElasticBase::finalizeElement;
  //! \brief Finalizes the element matrices after the numerical integration.
  //! \param elmInt The local integral object to receive the contributions
  //! \param[in] fe Nodal and integration point data for current element
  //! \param[in] time Parameters for nonlinear and time-dependent simulations
  virtual bool finalizeElement(LocalIntegral& elmInt, const FiniteElement& fe,
                               const TimeDomain& time, size_t);

  //! \brief Returns whether there are any load values to write to VTF.
  virtual bool hasTractionValues() const;

  //! \brief Writes the surface pressure for a given time step to VTF-file.
  //! \param vtf The VTF-file object to receive the pressure vectors
  //! \param[in] iStep Load/time step identifier
  //! \param geoBlk Running geometry block counter
  //! \param nBlock Running result block counter
  virtual bool writeGlvT(VTF* vtf, int iStep, int& geoBlk, int& nBlock) const;

  //! \brief Returns the number of primary/secondary solution field components.
  //! \param[in] fld which field set to consider (1=primary, 2=secondary)
  virtual size_t getNoFields(int fld) const { return fld < 2 ? npv+1 : n2v; }
  //! \brief Returns the name of a secondary solution field component.
  //! \param[in] i Field component index
  //! \param[in] prefix Name prefix for all components
  virtual std::string getField2Name(size_t i, const char* prefix) const;

  //! \brief Computes some derived primary solution quantities.
  virtual void primaryScalarFields(Matrix& field);

  //! \brief Specifies the recovery of von Mises stresses only.
  void vonMisesOnly() { n2v = 2; }

protected:
  //! \brief Evaluates the secondary solution at a result point.
  //! \param[out] s The solution field values at current point
  //! \param[in] eV Element-level primary solution vectors
  //! \param[in] fe Finite element data at current point
  virtual bool evalSol2(Vector& s, const Vectors& eV,
                        const FiniteElement& fe, const Vec3&) const;

  //! \brief Returns whether element \a iel has pressure loads or not.
  bool havePressure(int iel = -1) const;
  //! \brief Evaluates the surface pressure function(s) at specified point.
  //! \param p Updated surface pressue value at evaluation point
  //! \param[in] X Cartesian coordinates of evaluation point
  //! \param[in] n Shell surface normal vector at evaluation point
  //! \param[in] iel External index of element containing the evaluation point
  void addPressure(Vec3& p, const Vec3& X, const Vec3& n, int iel) const;

private:
  double Thck0; //!< Initial (uniform) shell thickness
  double Thick; //!< Current shell thickness
  double Emod;  //!< Current Young's modules
  double Rny;   //!< Current Poisson's ratio
  double Rho;   //!< Current mass density
  bool ovrMat;  //!< If \e true, patch-level material properties are overridden

  double trInside;  //!< Thickness loss level inside given box domain
  double trOutside; //!< Thickness loss level outside given box domain
  Vec3 Xlow; //!< Lower bound of thickness loss box
  Vec3 Xupp; //!< Uppoer bound of thickness loss box
  ScalarFunc* thickLoss; //!< Thickness loss function

  const ASMu2DNastran* currentPatch; //!< Pointer to underlying FE model

  std::map<int,RealFunc*> presFld; //!< Pointers to pressure field functions

  mutable std::vector<Vec3Pair> presVal; //!< Pressure field point values

  std::set<int> degenerated;  //!< List of detected degenerated elements
  std::set<int> straightline; //!< List of detected straight line elements

  Matrices myKmats; //!< Element stiffness matrix buffer
  Matrices myMmats; //!< Element mass matrix buffer

  size_t n2v; //!< Number of secondary variables to output

  bool isModal; //!< Flag for modal dynamics simulation
};

#endif
