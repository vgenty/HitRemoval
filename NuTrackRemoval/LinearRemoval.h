/**
 * \file LinearRemoval.h
 *
 * \ingroup NuTrackRemoval
 * 
 * \brief Class def header for a class LinearRemoval
 *
 * @author david
 */

/** \addtogroup NuTrackRemoval

    @{*/

#ifndef LARLITE_LINEARREMOVAL_H
#define LARLITE_LINEARREMOVAL_H

#include "Analysis/ana_base.h"

#include "HitRemovalBase.h"

#include "TwoDimTools/Linearity.h"

namespace larlite {
  /**
     \class LinearRemoval
     User custom analysis class made by SHELL_USER_NAME
   */
  class LinearRemoval : public ana_base, HitRemovalBase {
  
  public:

    /// Default constructor
    LinearRemoval(){ _name="LinearRemoval"; _fout=0;}

    /// Default destructor
    virtual ~LinearRemoval(){}

    virtual bool initialize();

    virtual bool analyze(storage_manager* storage);

    virtual bool finalize();

    /// Set Producers
    void setClusterProducer(std::string s) { _clusterProducer = s; }

    /// set maximum linearity allowed for this
    void setMaxLinearity(double l) { _max_lin_v.push_back( l ); }
    /// set minimum number of hits
    void setMinNHits(int n) { _min_n_hits_v.push_back( n ); }

    /// set maximum SSV
    void setMaxSSV(double s) { _ssv = s; }

  protected:

    /// producers
    std::string _clusterProducer;

    /// maximum linearity for hits
    std::vector<double> _max_lin_v;
    // min number of hits for cluster to be considered for removal
    std::vector<int>    _min_n_hits_v;
    // maximum ssv allowed
    double _ssv;
    
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
