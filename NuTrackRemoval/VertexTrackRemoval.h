/**
 * \file VertexTrackRemoval.h
 *
 * \ingroup Clusterer
 * 
 * \brief Class def header for a class VertexTrackRemoval
 *
 * @author david caratelli
 */

/** \addtogroup Clusterer

    @{*/

#ifndef LARLITE_VERTEXTRACKREMOVAL_H
#define LARLITE_VERTEXTRACKREMOVAL_H

#include "TwoDimTools/Linearity.h"

#include "HitRemovalBase.h"

#include <map>

namespace larlite {
  /**
     \class VertexTrackRemoval
     User custom analysis class made by SHELL_USER_NAME
   */
  class VertexTrackRemoval : public HitRemovalBase {
  
  public:

    /// Default constructor
    VertexTrackRemoval();

    /// Default destructor
    virtual ~VertexTrackRemoval(){}

    /** IMPLEMENT in VertexTrackRemoval.cc!
        Initialization method to be called before the analysis event loop.
    */ 
    virtual bool initialize();

    /** IMPLEMENT in VertexTrackRemoval.cc! 
        Analyze a data event-by-event  
    */
    virtual bool analyze(storage_manager* storage);

    /// set maximum linearity allowed for this
    void setMaxLinearity(double l) { _max_lin_v.push_back( l ); }
    /// set minimum number of hits
    void setMinNHits(int n) { _min_n_hits_v.push_back( n ); }
    /// set maximum proton linearity
    void setMaxProtonLin(double l)  { _max_proton_lin = l; }
    /// set max proton distance
    void setMaxProtonDist(double d) { _max_proton_dist = d; }
    
    /// set vertex radius
    void setVtxRadius(double r) { _vtx_rad = r; }
    
    /// Set Producers
    void setClusterProducer(std::string s) { _clusterProducer = s; }
    void setVertexProducer (std::string s) { _vertexProducer  = s; }

  protected:

    /// maximum linearity for hits
    std::vector<double> _max_lin_v;
    // min number of hits for cluster to be considered for removal
    std::vector<int>    _min_n_hits_v;

    /// ROI radius
    double _vtx_rad;
    /// max proton distance
    double _max_proton_dist;
    /// max proton linearity
    double _max_proton_lin;
    
    /// producers
    std::string _clusterProducer;
    std::string _vertexProducer;

    int _nhits;
    double _lin;
    double _local_lin_truncated;
    double _local_lin_avg;
    
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
