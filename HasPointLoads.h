// $Id$
//==============================================================================
//!
//! \file HasPointLoads.h
//!
//! \date Jun 12 2025
//!
//! \author Knut Morten Okstad / SINTEF
//!
//! \brief Class representing a nodal point load container.
//!
//==============================================================================

#ifndef _HAS_POINT_LOADS_H
#define _HAS_POINT_LOADS_H

#include <vector>
#include <map>

class ScalarFunc;
class AlgEqSystem;
namespace tinyxml2 { class XMLElement; }


/*!
  \brief Class representing a nodal point load container.
*/

class HasPointLoads
{
protected:
  //! \brief Struct defining a nodal point load.
  struct PointLoad
  {
    int         inod; //!< Node index
    int         ldof; //!< Local DOF number
    ScalarFunc* func; //!< The load magnitude as function of time
    //! \brief Default constructor.
    explicit PointLoad(int n = 0, int d = 0, ScalarFunc* f = nullptr)
      : inod(n), ldof(d), func(f) {}
  };

  //! \brief Empty default constructor.
  HasPointLoads() {}
  //! \brief The destructor deletes the load functions.
  virtual ~HasPointLoads();

  //! \brief Parses a point load definition from an XML-element
  bool parseLoad(const tinyxml2::XMLElement* elem);

  //! \brief Renumbers the global node numbers of the point loads.
  bool renumberLoadedNodes(const std::map<int,int>& nodeMap);

  //! \brief Assemble right-hand-side vector contributions from the point loads.
  bool assembleLoads(AlgEqSystem& eqSys, double t) const;

private:
  std::vector<PointLoad> myLoads; //!< Nodal point load container
};

#endif
