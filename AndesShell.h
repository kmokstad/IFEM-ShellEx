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
#include "BeamProperty.h"
#include <set>

class ASMu2DNastran;
class ASMuBeam;
class ElasticBeam;
class ScalarFunc;
class RealFunc;


/*!
  \brief Class representing the integrand of an ANDES shell.
*/

class AndesShell : public ElasticBase
{
public:
  //! \brief Default constructor.
  //! \param[in] ns Number of consecutive solution states to reside in core
  //! \param[in] modal If \e true, a modal dynamics simulation is performed
  //! \param[in] withBeams If \e true, the model also contains beam elements
  explicit AndesShell(unsigned short int ns = 1, bool modal = false,
                      bool withBeams = false);
  //! \brief The destructor writes out the ignored bad elements.
  virtual ~AndesShell();

  using ElasticBase::parseMatProp;
  //! \brief Parses material properties from an XML-element.
  virtual Material* parseMatProp(const tinyxml2::XMLElement* elem);

  //! \brief Prints out the problem definition to the log stream.
  virtual void printLog() const;

  //! \brief Checks if the point \a Xc is inside the thickness loss area or not.
  bool isInLossArea(const Vec3& Xc) const;

  //! \brief Initializes a time integration parameter for the integrand.
  //! \param[in] i Index of the integration parameter to define
  //! \param[in] prm The parameter value to assign
  virtual void setIntegrationPrm(unsigned short int i, double prm);

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
  //! \param[in] iEl Global element number (1-based)
  virtual LocalIntegral* getLocalIntegral(size_t nen, size_t iEl, bool) const;

  //! \brief Defines the global integral for calculating reaction forces only.
  virtual void setSecondaryInt(GlobalIntegral* gq);
  //! \brief Returns the system quantity to be integrated by \a *this.
  virtual GlobalIntegral& getGlobalInt(GlobalIntegral* gq) const;

  //! \brief Defines which FE quantities are needed by the integrand.
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

  //! \brief Returns \e true if no elements failed during assembly.
  bool allElementsOK() const { return failedElements.empty(); }

  //! \brief Assigns parameter values to the pressure functions.
  virtual void setParam(const std::string& name, const Vec3& value);

protected:
  //! \brief Evaluates the secondary solution at a result point.
  //! \param[out] s The solution field values at current point
  //! \param[in] eV Element-level primary solution vectors
  //! \param[in] fe Finite element data at current point
  virtual bool evalSol2(Vector& s, const Vectors& eV,
                        const FiniteElement& fe, const Vec3&) const;

  //! \brief Returns whether global element \a iEl has pressure loads or not.
  bool havePressure(int iEl = -1) const;
  //! \brief Evaluates the surface pressure function(s) at specified point.
  //! \param p Updated surface pressue value at evaluation point
  //! \param[in] X Cartesian coordinates of evaluation point
  //! \param[in] n Shell surface normal vector at evaluation point
  //! \param[in] iEl Global element number (1-based) containing evaluation point
  void addPressure(Vec3& p, const Vec3& X, const Vec3& n, int iEl) const;

private:
  double Thck0; //!< Initial (uniform) shell thickness
  double Thick; //!< Current shell thickness
  double Emod;  //!< Current Young's modules
  double GorNu; //!< Current Poisson's ratio or shear modulus (G)
  double Rho;   //!< Current mass density
  bool ovrMat;  //!< If \e true, patch-level material properties are overridden

  double trInside;  //!< Thickness loss level inside given box domain
  double trOutside; //!< Thickness loss level outside given box domain
  Vec3 Xlow; //!< Lower bound of thickness loss box
  Vec3 Xupp; //!< Uppoer bound of thickness loss box
  ScalarFunc* thickLoss; //!< Thickness loss function

  const ASMu2DNastran* currentPatch; //!< Pointer to underlying FE model
  const ASMuBeam*      beamPatch;    //!< Pointer to underlying beam model
  ElasticBeam*         beamProblem;  //!< Pointer to beam problem integrand
  BeamProperty         myBeamProps;  //!< Beam cross section properties
  GlobalIntegral*      myReacI;      //!< Reaction-forces-only integral

  std::map<int,RealFunc*> presFld; //!< Pointers to pressure field functions

  mutable std::vector<Vec3Pair> presVal; //!< Pressure field point values

  std::set<int> degenerated;    //!< List of detected degenerated elements
  std::set<int> straightline;   //!< List of detected straight line elements
  std::set<int> failedElements; //!< List of element with assembly failure

  Matrices myKmats; //!< Element stiffness matrix buffer
  Matrices myMmats; //!< Element mass matrix buffer

  size_t n2v; //!< Number of secondary variables to output

  bool isModal; //!< Flag for modal dynamics simulation
};

#endif
