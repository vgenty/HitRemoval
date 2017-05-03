/**
 * \file RemoveHitsNearVtx.h
 *
 * \ingroup Clusterer
 * 
 * \brief Class def header for a class RemoveHitsNearVtx
 *
 * @author david caratelli
 */

/** \addtogroup Clusterer

    @{*/

#ifndef LARLITE_REMOVEHITSNEARVERTEX_H
#define LARLITE_REMOVEHITSNEARVERTEX_H

#include "Analysis/ana_base.h"

#include "DataFormat/hit.h"
#include "DataFormat/vertex.h"

#include "TwoDimTools/Linearity.h"

#include <map>

namespace larlite {
  /**
     \class RemoveHitsNearVtx
     User custom analysis class made by SHELL_USER_NAME
   */
  class RemoveHitsNearVtx : public ana_base{
  
  public:

    /// Default constructor
    RemoveHitsNearVtx();

    /// Default destructor
    virtual ~RemoveHitsNearVtx(){}

    /** IMPLEMENT in RemoveHitsNearVtx.cc!
        Initialization method to be called before the analysis event loop.
    */ 
    virtual bool initialize();

    /** IMPLEMENT in RemoveHitsNearVtx.cc! 
        Analyze a data event-by-event  
    */
    virtual bool analyze(storage_manager* storage);

    /** IMPLEMENT in RemoveHitsNearVtx.cc! 
        Finalize method to be called after all events processed.
    */
    virtual bool finalize();

    /// Verbosity setter
    void setVerbose(bool on) { _verbose = on; }
    void setDebug  (bool on) { _debug   = on; }

    /// vertex buffer region
    void setVtxRadius(double r) { _vtx_rad = r; }
    
    /// Set Producers
    void setHitProducer    (std::string s) { _hitProducer     = s; }
    void setVertexProducer (std::string s) { _vertexProducer  = s; }

  protected:

    bool loadVertex(event_vertex* ev_vtx);

    /// vertex coordinates
    std::vector<double> _vtx_w_cm;
    std::vector<double> _vtx_t_cm;

    /// ROI radius
    double _vtx_rad;
    
    /// verbosity flag
    bool _verbose;
    /// debug flag
    bool _debug;

    /// conversion factors for hits
    double _wire2cm, _time2cm;

    /// producers
    std::string _hitProducer;
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
