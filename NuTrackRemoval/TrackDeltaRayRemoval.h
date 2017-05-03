/**
 * \file TrackDeltaRayRemoval.h
 *
 * \ingroup Clustering
 * 
 * \brief Class def header for a class TrackDeltaRayRemoval
 *
 * @author david
 */

/** \addtogroup Clustering

    @{*/

#ifndef LARLITE_TRACKDELTARAYREMOVAL_H
#define LARLITE_TRACKDELTARAYREMOVAL_H

#include "Analysis/ana_base.h"
#include "DataFormat/hit.h"

#include "TwoDimTools/Linearity.h"

namespace larlite {
  /**
     \class TrackDeltaRayRemoval
     User custom analysis class made by SHELL_USER_NAME
   */
  class TrackDeltaRayRemoval : public ana_base{
  
  public:

    /// Default constructor
    TrackDeltaRayRemoval();

    /// Default destructor
    virtual ~TrackDeltaRayRemoval(){}

    /** IMPLEMENT in TrackDeltaRayRemoval.cc!
        Initialization method to be called before the analysis event loop.
    */ 
    virtual bool initialize();

    /** IMPLEMENT in TrackDeltaRayRemoval.cc! 
        Analyze a data event-by-event  
    */
    virtual bool analyze(storage_manager* storage);

    /** IMPLEMENT in TrackDeltaRayRemoval.cc! 
        Finalize method to be called after all events processed.
    */
    virtual bool finalize();

    /// Verbosity setter
    void setVerbose(bool on) { _verbose = on; }

    /// set delta-ray distances
    void setDeltaRayDistMin(double d) { _d_delta_min = d; }
    void setDeltaRayDistMax(double d) { _d_delta_max = d; }
    // set max # hits for a delta-ray
    void setMaxNHitsDelta(int n) { _nhitmax = n; }

    /// Set Producers
    void setClusterProducer(std::string s) { _clusProducer = s; }

  protected:

    /// function that determines if cluster is a delta-rays
    bool DeltaRay(const std::vector<unsigned int>& muon,
		  const std::vector<unsigned int>& deltaray);

    // distance between hits
    double _distSq_(const larlite::hit& h1, const larlite::hit& h2);

    // keep track of event hits
    event_hit* _ev_hit;

    /// verbosity flag
    bool _verbose;

    /// conversion factors for hits
    double _wire2cm, _time2cm;

    /// Producers
    std::string _clusProducer;

    // minimum & max dist to delta-ray
    double _d_delta_min, _d_delta_max;

    // maximum # of hits for a delta-ray
    int _nhitmax;
    
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
