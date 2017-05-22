/**
 * \file VertexProximityRemoval.h
 *
 * \ingroup Clusterer
 * 
 * \brief Class def header for a class VertexProximityRemoval
 *
 * @author david caratelli
 */

/** \addtogroup Clusterer

    @{*/

#ifndef LARLITE_VERTEXPROXIMITYREMOVAL_H
#define LARLITE_VERTEXPROXIMITYREMOVAL_H

#include "TwoDimTools/Linearity.h"

#include "HitRemovalBase.h"

#include <map>

namespace larlite {
  /**
     \class VertexProximityRemoval
     User custom analysis class made by SHELL_USER_NAME
   */
  class VertexProximityRemoval : public HitRemovalBase {
  
  public:

    /// Default constructor
    VertexProximityRemoval();

    /// Default destructor
    virtual ~VertexProximityRemoval(){}

    /** IMPLEMENT in VertexProximityRemoval.cc!
        Initialization method to be called before the analysis event loop.
    */ 
    virtual bool initialize();

    /** IMPLEMENT in VertexProximityRemoval.cc! 
        Analyze a data event-by-event  
    */
    virtual bool analyze(storage_manager* storage);

    /// set vertex radius
    void setVtxRadius(double r) { _vtx_rad = r; }
    /// set maximum linearity
    void setMaxLin(double l) { _lin_max = l; }
    
    /// Set Producers
    void setClusterProducer(std::string s) { _clusterProducer = s; }
    void setVertexProducer (std::string s) { _vertexProducer  = s; }

  protected:

    /// ROI radius
    double _vtx_rad;

    /// max linearity
    double _lin_max;
    
    /// producers
    std::string _clusterProducer;
    std::string _vertexProducer;

  };
}
#endif

//**************************************************************************
// 
// For Analysis framework documentation, read Manual.pdf here:
//
// http://microboone-docdb.fnal.gov:8080/cgi-bin/ShowDocument?docid=3183
//
//**************************************************************************

/** @} */ // end of doxygen group 
