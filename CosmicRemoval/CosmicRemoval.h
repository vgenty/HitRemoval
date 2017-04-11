/**
 * \file CosmicRemoval.h
 *
 * \ingroup Clusterer
 * 
 * \brief Class def header for a class CosmicRemoval
 *
 * @author david caratelli
 */

/** \addtogroup Clusterer

    @{*/

#ifndef LARLITE_COSMICREMOVAL_H
#define LARLITE_COSMICREMOVAL_H

#include "Analysis/ana_base.h"
#include "DataFormat/hit.h"

#include "TwoDimTools/Linearity.h"

#include <map>

namespace larlite {
  /**
     \class CosmicRemoval
     User custom analysis class made by SHELL_USER_NAME
   */
  class CosmicRemoval : public ana_base{
  
  public:

    /// Default constructor
    CosmicRemoval();

    /// Default destructor
    virtual ~CosmicRemoval(){}

    /** IMPLEMENT in CosmicRemoval.cc!
        Initialization method to be called before the analysis event loop.
    */ 
    virtual bool initialize();

    /** IMPLEMENT in CosmicRemoval.cc! 
        Analyze a data event-by-event  
    */
    virtual bool analyze(storage_manager* storage);

    /** IMPLEMENT in CosmicRemoval.cc! 
        Finalize method to be called after all events processed.
    */
    virtual bool finalize();

    /// set maximum linearity allowed for this
    void setMaxLinearity(double l) { _max_lin_v.push_back( l ); }
    /// Verbosity setter
    void setVerbose(bool on) { _verbose = on; }

    /// set minimum number of hits
    void setMinNHits(int n) { _min_n_hits_v.push_back( n ); }
    
    /// Set Producers
    void setClusterProducer(std::string s) { _clusterProducer = s; }
    void setVertexProducer (std::string s) { _vtxProducer     = s; }

  protected:

    /// vertex coordinates
    std::vector<double> _vtx_w_cm;
    std::vector<double> _vtx_t_cm;

    /// maximum linearity for hits
    std::vector<double> _max_lin_v;

    // min number of hits for cluster to be considered for removal
    std::vector<int>    _min_n_hits_v;
    
    /// verbosity flag
    bool _verbose;

    /// conversion factors for hits
    double _wire2cm, _time2cm;

    /// Producers
    std::string _clusterProducer;
    std::string _vtxProducer;
    std::string _out_hitProducer;

    TTree* _tree;
    int _nhits;
    double _lin;
    double _angle;
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
