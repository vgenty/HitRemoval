/**
 * \file RemoveDeltaRays.h
 *
 * \ingroup Clustering
 * 
 * \brief Class def header for a class RemoveDeltaRays
 *
 * @author david
 */

/** \addtogroup Clustering

    @{*/

#ifndef LARLITE_REMOVEDELTARAYS_H
#define LARLITE_REMOVEDELTARAYS_H

#include "HitRemovalBase.h"

#include "TwoDimTools/Linearity.h"

#include <map>

struct BBox {
  double wmin;
  double wmax;
  double tmin;
  double tmax;
};

namespace larlite {
  /**
     \class RemoveDeltaRays
     User custom analysis class made by SHELL_USER_NAME
   */
  class RemoveDeltaRays : public HitRemovalBase {
  
  public:

    /// Default constructor
    RemoveDeltaRays();

    /// Default destructor
    virtual ~RemoveDeltaRays(){}

    /** IMPLEMENT in RemoveDeltaRays.cc!
        Initialization method to be called before the analysis event loop.
    */ 
    virtual bool initialize();

    /** IMPLEMENT in RemoveDeltaRays.cc! 
        Analyze a data event-by-event  
    */
    virtual bool analyze(storage_manager* storage);

    /// set delta-ray distances
    void setDeltaRayDistMin(double d) { _d_delta_min = d; }
    void setDeltaRayDistMax(double d) { _d_delta_max = d; }
    void setMaxDeltaHits(int n) { _max_delta_hits = n; }
    void setROI(double r) { _roi = r; }

    /// Set Producers
    void setClusterProducer(std::string s) { _clusProducer = s; }
    void setVertexProducer(std::string s) { _vertexProducer = s; }

  protected:
    
    /// function that determines if cluster is a delta-rays
    bool DeltaRay(const std::vector<unsigned int>& muon,
		  const std::vector<unsigned int>& deltaray);

    // distance between hits
    double _distSq_(const larlite::hit& h1, const larlite::hit& h2);

    // keep track of event hits
    event_hit* _ev_hit;

    /// Producers
    std::string _clusProducer;
    std::string _vertexProducer;

    // map connecting each cluster index to the cluster bounding box
    // each entry is [wmin, wmax, tmin, tmax]
    std::map< size_t, BBox > _clus_bbox;
    
    // minimum & max dist to delta-ray
    double _d_delta_min, _d_delta_max;
    int _max_delta_hits;
    double _roi;
    
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
