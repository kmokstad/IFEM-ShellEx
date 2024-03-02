// SPDX-FileCopyrightText: 2023 SAP SE
//
// SPDX-License-Identifier: Apache-2.0
//
// This file is part of FEDEM - https://openfedem.org
////////////////////////////////////////////////////////////////////////////////

#include <algorithm>
#include <iostream>

#include "FFlLib/FFlFEElementTopSpec.H"

void
FFlFEElementTopSpec::dump() const
{
  std::vector<FaceType>::const_iterator fit;
  std::vector<EdgeType>::const_iterator eit;

  std::cout <<"Topology";
  for (fit = myFaces.begin(); fit != myFaces.end(); fit++)
  {
    std::cout <<"\n\tFace";
    for (eit = fit->myEdges.begin(); eit != fit->myEdges.end(); eit++)
      std::cout <<" ["<< eit->first <<" "<< eit->second <<"]";
  }
  std::cout << std::endl;
}


void
FFlFEElementTopSpec::setTopology(int cFace, int* vFaces, bool expanded)
{
  int faces, i;
  for (faces = i = 0; faces < cFace; faces++, i++)
  {
    FaceType aFace;
    int faceStartIndex = i;
    // create new edges
    while (vFaces[++i] != -1)
      aFace.addEdge(vFaces[i-1],vFaces[i]);
    // create closing edge
    if (i-1 > faceStartIndex)
      aFace.addEdge(vFaces[i-1],vFaces[faceStartIndex]);
    else
      i--;

    // add face
    if (expanded)
      this->addExpandedFace(aFace);
    else
      this->addFace(aFace);
  }
}


void
FFlFEElementTopSpec::setExpandedTopology(int cFace, int* vFaces)
{
  this->setTopology(cFace,vFaces,true);
}


/*!
  Returns the local element node ids for each vertex in face.
*/

bool
FFlFEElementTopSpec::getFaceTopology(short int faceNumber,
					  bool isExpandedFace,
					  bool switchNormal,
					  short int idxOffset,
					  std::vector<int>& topology) const
{
  topology.clear();

  const std::vector<FaceType>& faces = isExpandedFace ? myExpandedFaces : myFaces;
  if (faceNumber < 0 || (size_t)faceNumber >= faces.size()) {
#ifdef FFL_DEBUG
    std::cerr <<"FFlFEElementTopSpec::getFaceTopology: faceNumber "<< faceNumber
	 <<" out of range."<< std::endl;
#endif


    return false;
  }

  std::vector<EdgeType> edges = faces[faceNumber].myEdges;

  if ( idxOffset > 0 && (size_t)idxOffset < edges.size() )
  {
    rotate(edges.begin(), edges.begin()+idxOffset, edges.end());
  }
  else if ( idxOffset != 0 )
  {
     std::cerr << "FFlFEElementTopSpec::getFaceTopology: edge offset "
               << idxOffset
	             << " out of range.\n\t\tIn element with node count: "
	             << myNodeCount
	             << " node DOFS: "
	             << myNodeDOFs
	             << std::endl;
  }
  if ( switchNormal )
  {
    reverse(edges.begin()+1, edges.end());
  }

  topology.resize(edges.size());
  for ( size_t i = 0; i < edges.size(); ++i )
  {
    topology[i] = switchNormal ? edges[i].second : edges[i].first;
  }

  return true;
}


/*!
  Get the number of vertices in the given local face.
*/

int FFlFEElementTopSpec::getFaceVertexCount(short int faceNumber,
					    bool isExpandedFace) const
{
  std::vector<FaceType>const & faces = isExpandedFace ? myExpandedFaces : myFaces;
  if (faceNumber < 0 || (size_t)faceNumber >= faces.size()) return 0;

  return faces[faceNumber].myEdges.size();
}



/*!
  Returns the local face number of a solid defined by two local node ids.
  The way we determine the unique face number depends on the element type,
  and matches the description of the Nastran Bulk Data File PLOAD4 entry.

  Tetrahedrons: node1 must be on the face and node2 must not be on the face.
  Pentahedrons: if node2 != 0, both node1 and node2 must be on a quad face.
                if node2 == 0, node1 must be on a triangular face.
  Hexahedrons: both node1 and node2 must be on a face.
*/

short int
FFlFEElementTopSpec::getFaceNum(int node1, int node2) const
{
  if (node1 <= 0 || node1 == node2) return 0;
  if (node1 > myNodeCount || node2 > myNodeCount) return 0;

  size_t i, j, k, nFaces = myFaces.size();
  switch (myFaces.size())
    {
  case 1: // Shell surface, just check that node1 is among the face nodes
    for (j = 0; j < myFaces.front().myEdges.size(); j++)
    	if (myFaces.front().myEdges[j].first == node1)
    		return 1;
    break;

    case 4: // Tetrahedron
      for (i = 0; i < nFaces; i++)
      {
	// node1 must be on the face, and node2 must not be on the face
	short int theFace = 0;
	for (j = 0; j < myFaces[i].myEdges.size() && theFace >= 0; j++)
	  if (myFaces[i].myEdges[j].first == node1 && theFace == 0)
	    theFace = i+1;
	  else if (myFaces[i].myEdges[j].first == node2)
	    theFace = -1;

	if (theFace > 0) return theFace;
      }
      break;

    case 5: // Pentahedron
      if (node2 > 0)
	nFaces = 3; // assuming that the 3 quadrilateral faces are ordered first
      else
      {
	// node1 must be on a triangular face
	for (i = 3; i < nFaces; i++)
	  for (j = 0; j < myFaces[i].myEdges.size(); j++)
	    if (myFaces[i].myEdges[j].first == node1)
	      return i+1;

	break;
      }

    case 6: // Hexahedron, or quadrilateral face of Pentahedron
      for (i = 0; i < nFaces; i++)
	// both node1 and node 2 must be on the face
	for (j = 0; j < myFaces[i].myEdges.size(); j++)
	  if (myFaces[i].myEdges[j].first == node1)
	    for (k = 0; k < myFaces[i].myEdges.size(); k++)
	      if (myFaces[i].myEdges[k].first == node2)
		return i+1;
      break;
    }

  return 0; // no match
}
