// SPDX-FileCopyrightText: 2023 SAP SE
//
// SPDX-License-Identifier: Apache-2.0
//
// This file is part of FEDEM - https://openfedem.org
////////////////////////////////////////////////////////////////////////////////

#ifndef FFL_ELEMENT_BASE_H
#define FFL_ELEMENT_BASE_H

#include "FFlLib/FFlPartBase.H"
#include "FFlLib/FFlAttributeBase.H"
#include "FFlLib/FFlFEAttributeRefs.H"
#include "FFlLib/FFlFENodeRefs.H"
#ifdef FT_USE_VISUALS
#include "FFlLib/FFlVisualRefs.H"
#endif
#include "FFlLib/FFlTypeInfoSpec.H"
#include "FFaLib/FFaPatterns/FFaGenericFactory.H"

class FFlFEResultBase;
class FFlFEElmResult;
class FFaCheckSum;
class FFaTensor3;
class FaMat33;
class FaVec3;

#ifdef FF_NAMESPACE
namespace FF_NAMESPACE {
#endif

class FFlLinkHandler;


class FFlElementBase : public FFlPartBase,
		       public FFlFEAttributeRefs,
#ifdef FT_USE_VISUALS
		       public FFlFENodeRefs,
		       public FFlVisualRefs
#else
		       public FFlFENodeRefs
#endif
{
protected:
  using Elements = std::vector<FFlElementBase*>;

  FFlElementBase(int ID = -1);
  FFlElementBase(const FFlElementBase& obj);
  FFlElementBase& operator=(const FFlElementBase&) = delete;

public:
  virtual ~FFlElementBase();

  virtual void calculateChecksum(FFaCheckSum* cs, int cstype) const;
  virtual FFlElementBase* clone() const = 0;

  const std::string& getTypeName() const
  { return this->getTypeInfoSpec()->getTypeName(); }

  FFlTypeInfoSpec::Cathegory getCathegory() const
  { return this->getTypeInfoSpec()->getCathegory(); }

  // Reference to underlying finite element (for strain coat elements only):
  virtual void setFElement(FFlElementBase*) {}
  virtual void setFElement(const int) {}
  virtual bool resolveElmRef(const Elements&, bool = false) { return true; }
  virtual FFlElementBase* getFElement() const { return NULL; }

  //! \brief Checks that the element have positive volume, try to fix it if not.
  virtual int checkOrientation(bool = false) { return 1; }
  //! \brief Splits a parabolic element into linear elements.
  virtual bool split(Elements&, FFlLinkHandler*, int = 0) { return false; }

  // Structural properties:
  virtual double getMassDensity() const;
  virtual bool getVolumeAndInertia(double& volume, FaVec3& cog,
				   FFaTensor3& inertia) const;
  bool getMassProperties(double& mass, FaVec3& cog, FFaTensor3& inertia) const;

  // Geometrical mapping:
  virtual bool getLocalSystem(double*) const { return false; }
  virtual FaMat33 getGlobalizedElmCS() const;
  virtual FaVec3 interpolate(const double*, const std::vector<FaVec3>&) const;
  virtual FaVec3 mapping(double, double = 0.0, double = 0.0) const;
  virtual bool invertMapping(const FaVec3&, double*) const;

  // Calculation options - used for e.g. stress calculations and result print:
  bool doCalculations() const { return calculateResults; }
  void setUpForCalculations(bool trueFalse) { calculateResults = trueFalse; }

  // Result management:
  bool hasResults() const { return myResults != NULL; }
  FFlFEElmResult* getResults();
  void deleteResults();

#ifndef FT_USE_VISUALS
  bool isVisible() const { return true; }
#endif

private:
  FFlFEResultBase* myResults;
  bool             calculateResults;
};


using ElementFactory = FFaGenericFactory<FFlElementBase>;

#define FFL_FE_ELEMENT_FACTORY_INIT(clName) public: \
  static void create(int ID, FFlElementBase*& elm) { elm = new clName(ID); } \
  clName& operator=(const clName&) = delete; \
  virtual FFlElementBase* clone() const { return new clName(*this); }

#ifdef FF_NAMESPACE
} // namespace
#endif
#endif
