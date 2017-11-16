/**
 * \file DRRemoval.h
 *
 * \ingroup Clustering
 * 
 * \brief Class def header for a class DRRemoval
 *
 * @author david caratelli [dcaratelli@nevis.columbia.edu]
 */

/** \addtogroup Clustering

    @{*/

#ifndef LARLITE_DRREMOVAL_H
#define LARLITE_DRREMOVAL_H

#include "HitRemovalBase.h"

#include "TwoDimTools/Linearity.h"
#include "GeoAlgo/GeoAlgo.h"

#include <map>

typedef std::vector< std::pair<double,double> > Traj2D;

namespace larlite {
  /**
     \class DRRemoval
     User custom analysis class made by SHELL_USER_NAME
   */
  class DRRemoval : public HitRemovalBase {
  
  public:

    /// Default constructor
    DRRemoval();

    /// Default destructor
    virtual ~DRRemoval(){}

    /** IMPLEMENT in DRRemoval.cc!
        Initialization method to be called before the analysis event loop.
    */ 
    virtual bool initialize();

    /** IMPLEMENT in DRRemoval.cc! 
        Analyze a data event-by-event  
    */
    virtual bool analyze(storage_manager* storage);

    /// set delta-ray distances
    void setDeltaRayDistMin(double d) { _dmin = d; }
    void setDeltaRayDistMax(double d) { _dmax = d; }
    void setMaxDeltaHits(int n) { _max_delta_hits = n; }
    void setROI(double r) { _roi = r; }

    /// Set Producers
    void setClusterProducer(std::string s) { _clusProducer   = s; }
    void setVertexProducer (std::string s) { _vertexProducer = s; }
    void setTrackProducer  (std::string s) { _trackProducer  = s; }

  protected:

    /// vector storing vectors of 2D trajectory points
    /// one vector per track, for each plane
    std::vector< std::vector< Traj2D > > _trk2D_v;
    
    /// function that determines if cluster is a delta-rays
    bool IsDeltaRay(const std::vector<unsigned int>& hitidx_v, const int& pl);

    /// is a point contained within certain bounds
    bool Contained(const std::pair<double,double>& pt2D,
		   const double& tmin, const double& tmax,
		   const double& wmin, const double& wmax);

    // keep track of event hits
    event_hit* _ev_hit;

    ::geoalgo::GeoAlgo _geoAlgo;

    /// Producers
    std::string _trackProducer;
    std::string _clusProducer;
    std::string _vertexProducer;

    // minimum & max dist to delta-ray
    double _dmin, _dmax;
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
