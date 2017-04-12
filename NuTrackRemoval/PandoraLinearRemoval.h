/**
 * \file PandoraLinearRemoval.h
 *
 * \ingroup Clusterer
 * 
 * \brief Class def header for a class PandoraLinearRemoval
 *
 * @author david caratelli
 */

/** \addtogroup Clusterer

    @{*/

#ifndef LARLITE_PANDORALINEARREMOVAL_H
#define LARLITE_PANDORALINEARREMOVAL_H

#include "Analysis/ana_base.h"
#include "DataFormat/hit.h"

#include "TwoDimTools/Linearity.h"

#include <map>

namespace larlite {
  /**
     \class PandoraLinearRemoval
     User custom analysis class made by SHELL_USER_NAME
   */
  class PandoraLinearRemoval : public ana_base{
  
  public:

    /// Default constructor
    PandoraLinearRemoval();

    /// Default destructor
    virtual ~PandoraLinearRemoval(){}

    /** IMPLEMENT in PandoraLinearRemoval.cc!
        Initialization method to be called before the analysis event loop.
    */ 
    virtual bool initialize();

    /** IMPLEMENT in PandoraLinearRemoval.cc! 
        Analyze a data event-by-event  
    */
    virtual bool analyze(storage_manager* storage);

    /** IMPLEMENT in PandoraLinearRemoval.cc! 
        Finalize method to be called after all events processed.
    */
    virtual bool finalize();

    /// set maximum linearity allowed for this
    void setMaxLinearity(double l) { _max_lin_v.push_back( l ); }
    /// set minimum number of hits
    void setMinNHits(int n) { _min_n_hits_v.push_back( n ); }
    
    /// Verbosity setter
    void setVerbose(bool on) { _verbose = on; }
    void setDebug  (bool on) { _debug   = on; }
    

    /// set max distance pandora cluster can have to vertex
    /// in order to be considered for linear removal at this
    /// stage.
    void setMaxDVtx(double d) { _dvtx_max = d; }

    /// Set Producers
    void setClusterProducer(std::string s) { _clusterProducer = s; }
    void setVertexProducer (std::string s) { _vertexProducer  = s; }

  protected:

    /// vertex coordinates
    std::vector<double> _vtx_w_cm;
    std::vector<double> _vtx_t_cm;

    /// maximum linearity for hits
    std::vector<double> _max_lin_v;
    // min number of hits for cluster to be considered for removal
    std::vector<int>    _min_n_hits_v;

    /// max vtx-clus distance allowed
    double _dvtx_max;
    
    /// verbosity flag
    bool _verbose;
    /// debug flag
    bool _debug;

    /// conversion factors for hits
    double _wire2cm, _time2cm;

    /// producers
    std::string _clusterProducer;
    std::string _vertexProducer;

    TTree* _tree;
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
