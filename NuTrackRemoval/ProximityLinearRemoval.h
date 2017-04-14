/**
 * \file ProximityLinearRemoval.h
 *
 * \ingroup NuTrackRemoval
 * 
 * \brief Class def header for a class ProximityLinearRemoval
 *
 * @author david
 */

/** \addtogroup NuTrackRemoval

    @{*/

#ifndef LARLITE_PROXIMITYLINEARREMOVAL_H
#define LARLITE_PROXIMITYLINEARREMOVAL_H

#include "Analysis/ana_base.h"

#include "DataFormat/vertex.h"

#include "TwoDimTools/Linearity.h"

namespace larlite {
  /**
     \class ProximityLinearRemoval
     User custom analysis class made by SHELL_USER_NAME
   */
  class ProximityLinearRemoval : public ana_base{
  
  public:

    /// Default constructor
    ProximityLinearRemoval(){ _name="ProximityLinearRemoval"; _fout=0;}

    /// Default destructor
    virtual ~ProximityLinearRemoval(){}

    virtual bool initialize();

    virtual bool analyze(storage_manager* storage);

    virtual bool finalize();

    /// Set Producers
    void setClusterProducer(std::string s) { _clusterProducer = s; }
    void setVertexProducer (std::string s) { _vertexProducer  = s; }

    /// set maximum linearity allowed for this
    void setMaxLinearity(double l) { _max_lin_v.push_back( l ); }
    /// set minimum number of hits
    void setMinNHits(int n) { _min_n_hits_v.push_back( n ); }

  protected:

    // load vertex
    bool loadVertex(event_vertex* ev_vtx);

    /// producers
    std::string _clusterProducer;
    std::string _vertexProducer;

    /// vertex coordinates
    std::vector<double> _vtx_w_cm;
    std::vector<double> _vtx_t_cm;

    /// maximum linearity for hits
    std::vector<double> _max_lin_v;
    // min number of hits for cluster to be considered for removal
    std::vector<int>    _min_n_hits_v;
    
    /// conversion factors for hits
    double _wire2cm, _time2cm;
    
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
