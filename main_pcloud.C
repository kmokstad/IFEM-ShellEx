// $Id$
//==============================================================================
//!
//! \file main_pcloud.C
//!
//! \date Mar 17 2026
//!
//! \author Knut Morten Okstad / SINTEF
//!
//! \brief Utility for extracting a point cloud from a FE model.
//!
//==============================================================================

#include "IFEMNastranReader.h"
#include "FFlLinkHandler.H"
#include "FFlGroup.H"
#include "FFlFEParts/FFlAllFEParts.H"

#include <algorithm>
#include <numeric>
#include <fstream>
#include <sstream>
#include <cstring>
#include <cstdlib>


/*!
  \brief Auxiliary class to handle initialization of heap-allocated singeltons.
*/

class FFlInit
{
public:
  //! \brief The constructor initializes the finite element library.
  FFlInit() { FFl::initAllElements(); }
  //! \brief The destructor releases the finite element library.
  ~FFlInit() { FFl::releaseAllElements(); }
};


/*!
  \brief Reads a Nastran bulk data file into a FFlLinkHandler object.
*/

FFlLinkHandler* readNastran (const char* fileName)
{
  std::ifstream is(fileName);
  std::stringstream sets;

  int lCount = IFEMNastranReader::parseSets(is,sets,fileName);
  if (lCount < 0) return nullptr; // invalid Nastran file

  FFlLinkHandler* fem = new FFlLinkHandler();
  if (IFEMNastranReader reader(*fem,lCount);
      reader.readFE(is,sets) && fem->resolve())
    std::cout <<"\nParsing Nastran bulk data file succeeded."<< std::endl;
  else
  {
    delete fem;
    std::cerr <<"\n *** Parsing/resolving FE data failed.\n"
              <<"     The FE model is probably not consistent and has not been"
              <<" resolved completely."<< std::endl;
    return nullptr;
  }

  int nnod = fem->getNodeCount(FFlLinkHandler::FFL_FEM);
  int nel  = fem->getElementCount(FFlTypeInfoSpec::SHELL_ELM);
  int nBel = fem->getElementCount(FFlTypeInfoSpec::BEAM_ELM);
  std::cout <<"\nTotal number of nodes:          "<< nnod
            <<"\nNumber of shell elements:       "<< nel;
  if (nBel > 0)
    std::cout <<"\nNumber of beam elements:        "<< nBel;
  std::cout <<"\nNumber of constraint elements:  "
            << fem->getElementCount(FFlTypeInfoSpec::CONSTRAINT_ELM)
            <<"\nNumber of other elements:       "
            << fem->getElementCount(FFlTypeInfoSpec::OTHER_ELM)
            << std::endl;
  if (size_t allN = fem->getNodeCount(FFlLinkHandler::FFL_ALL); allN > nnod)
    std::cout <<"\n  ** Warning: This model contains "<< allN-nnod
              <<" node(s) without any element connections (ignored)."
              <<"\n     Please check the FE data file.\n"<< std::endl;

  return fem;
}


/*!
  \brief Main program for the point cloud extraction utility.
*/

int main (int argc, char** argv)
{
  if (argc < 2)
  {
    std::cout <<"usage: "<< argv[0] <<" <bdf-file>";
    for (char d = 'x'; d <= 'z'; d++)
      std::cout <<" [-"<< d <<"min <"<< d <<"0>] [-"<< d <<"max <"<< d << "1>]";
    std::cout <<" [-group <ID>] [-out <filename>] [-distancesort]"<< std::endl;
    return 0;
  }

  FFlInit _initializer;
  FFlLinkHandler* fem = readNastran(argv[1]);
  if (!fem) return 1;

  FaVec3 max, min;
  fem->getExtents(max,min);
  std::cout <<"Model extension (diameter):     "
            << (max-min).length() << std::endl;

  // Lambda function parsing the point cloud bound.
  auto&& parseRange = [](const char* option, const char* value,
                         double& x0, double& x1, char d)
  {
    static char xmin[6] = "-xmin";
    static char xmax[6] = "-xmax";
    xmin[1] = d;
    xmax[1] = d;
    if (!strcmp(option,xmin))
    {
      if (double xval = atof(value); xval > x0)
        x0 = xval;
      return true;
    }
    else if (!strcmp(option,xmax))
    {
      if (double xval = atof(value); xval < x1)
        x1 = xval;
      return true;
    }
    return false;
  };

  bool distansort = !strncmp(argv[argc-1],"-dist",5);
  const char* out = nullptr;
  FFlGroup* group = nullptr;
  for (int i = 2; i+1 < argc; i += 2)
    if (!strcmp(argv[i],"-out"))
      out = argv[i+1];
    else if (!strncmp(argv[i],"-dist",5))
      distansort = true, --i;
    else if (!strcmp(argv[i],"-group"))
    {
      group = fem->getGroup(atoi(argv[i+1]));
      std::cout <<"Element groups:";
      for (GroupCIter it = fem->groupsBegin(); it != fem->groupsEnd(); ++it)
        std::cout << (group == it->second ? "\n*\t" : "\n\t")
                  << it->first <<"\t"<< it->second->getName();
      std::cout << std::endl;
    }
    else for (int d = 0; d < 3; d++)
      if (parseRange(argv[i],argv[i+1],min[d],max[d],'x'+d))
        break;

  std::cout <<"\nExtracting a point cloud for the domain ["<< min
            <<"] x ["<< max <<"]";
  if (group)
    std::cout <<"\nfrom the element group "
              << group->getID() <<": "<< group->getName();
  std::cout << std::endl;

  std::vector<FaVec3> cloudPts;
  FaVec3 X0;

  // Lambda function adding an element center to the point cloud.
  auto&& addElmCenter = [&cloudPts,&X0,&min,&max](const FFlElementBase* elm)
  {
    FaVec3 X;
    if (elm->getVolumeAndCoG(X) <= 0.0)
      return; // zero-olum element, ignore

    for (int i = 0; i < 3; i++)
      if (X[i] < min[i] || X[i] > max[i])
        return; // point is outisde the bounding box

    cloudPts.push_back(X);
    X0 += X;
  };

  if (group)
    for (const GroupElemRef& elm : *group)
      addElmCenter(elm.getReference());
  else // no group specified, use all shell elements instead
    for (ElementsCIter e = fem->elementsBegin(); e != fem->elementsEnd(); ++e)
      if ((*e)->getCathegory() == FFlTypeInfoSpec::SHELL_ELM)
        addElmCenter(*e);

  delete fem; // release the FE model
  if (cloudPts.empty())
    return 2; // no points, invalid bounding box perhaps

  X0 /= cloudPts.size();
  std::cout <<"\nFound "<< cloudPts.size() <<" points, centered at "<< X0
            << std::endl;

  // Compute the theta angle of a cylindrical coordinate system,
  // with the global X-axis as the cylinder axis and where
  // zero angle is for points along the positive global Z-axis.
  // Therefore, the angles are shifted by 90 degrees.
  std::vector<double> angles;
  angles.reserve(cloudPts.size());
  for (const FaVec3& X : cloudPts)
    if (double theta = (X-X0).getAsCylCoords(VX).y(); theta < 0.5*M_PI)
      angles.push_back(theta*180.0/M_PI + 270.0);
    else
      angles.push_back(theta*180.0/M_PI - 90.0);

  // Index-sort the points w.r.t. the angles
  std::vector<int> indices(cloudPts.size());
  std::iota(indices.begin(),indices.end(),0);
  std::sort(indices.begin(),indices.end(),[&angles](int a, int b)
            {
              return angles[a] < angles[b];
            });

  std::cout <<"End points: "<< cloudPts[indices.front()]
            <<" (theta=" << angles[indices.front()]
            <<") and "<< cloudPts[indices.back()]
            <<" (theta=" << angles[indices.back()] <<")"<< std::endl;

  if (distansort)
  {
    // Reorder the points using the smallest-distance-to-neighbor approach
    std::vector<int> newOrder;
    newOrder.reserve(indices.size());
    newOrder.push_back(indices.front());
    indices.erase(indices.begin());
    while (!indices.empty())
    {
      size_t indxmin = 0;
      double distmin = (max-min).length();
      const FaVec3& Xlast = cloudPts[newOrder.back()];
      for (size_t i = 0; i < indices.size(); i++)
        if (double dist = (cloudPts[indices[i]]-Xlast).length(); dist < distmin)
        {
          distmin = dist;
          indxmin = i;
        }
      newOrder.push_back(indices[indxmin]);
      indices.erase(indices.begin()+indxmin);
    }
    indices.swap(newOrder);
  }

  if (out)
  {
    // Print out the point cloud in sorted order
    std::ofstream fs(out);
    for (int idx : indices)
      fs << cloudPts[idx] <<"\n";
  }

  return 0;
}
