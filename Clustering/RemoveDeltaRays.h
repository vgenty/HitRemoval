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

#include "Analysis/ana_base.h"
#include "DataFormat/hit.h"

#include "TwoDimTools/Linearity.h"

namespace larlite {
  /**
     \class RemoveDeltaRays
     User custom analysis class made by SHELL_USER_NAME
   */
  class RemoveDeltaRays : public ana_base{
  
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

    /** IMPLEMENT in RemoveDeltaRays.cc! 
        Finalize method to be called after all events processed.
    */
    virtual bool finalize();

    /// Verbosity setter
    void setVerbose(bool on) { _verbose = on; }

    /// set ROI size (cm away from vtx)
    void setROIRadius(double r) { _roi_radius = r; }

    /// Set Producers
    void setClusterProducer(std::string s) { _clusProducer = s; }
    void setVertexProducer (std::string s) { _vtxProducer  = s; }

  protected:

    /// vertex coordinates
    std::vector<double> _vtx_w_cm;
    std::vector<double> _vtx_t_cm;

    /// verbosity flag
    bool _verbose;

    /// conversion factors for hits
    double _wire2cm, _time2cm;

    /// Producers
    std::string _clusProducer;
    std::string _vtxProducer;

    /// ROI radius
    double _roi_radius;
    
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
