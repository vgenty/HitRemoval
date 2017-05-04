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
#include "DataFormat/vertex.h"

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

    void addSlopePt(double slope, double lin) {
      _pts_x_v.push_back(slope);
      _pts_y_v.push_back(lin);
    }
    
    /// Verbosity setter
    void setVerbose(bool on) { _verbose = on; }
    void setDebug  (bool on) { _debug   = on; }
    

    /// set max distance pandora cluster can have to vertex
    /// in order to be considered for linear removal at this
    /// stage.
    void setMaxDVtx(double d) { _dvtx_max = d; }
    /// set ROI bound. Clusters escaping this region will be removed
    void setROIRadius(double r) { _roi_rad = r; }
    /// set maximum vtx. distance for a possible proton
    void setProtonDMax(double d) { _proton_dmax = d; }
    /// maximum SSV for protons
    void setMaxSSV(double l) { _ssv_max = l; }
    /// set max slope
    void setSlopeMin(double s) { _slope_min = s; }
    /// set maximum local linearity for low-slope clusters
    void setLLT(double l ) { _llt_min = l; }
    
    /// Set Producers
    void setClusterProducer(std::string s) { _clusterProducer = s; }
    void setVertexProducer (std::string s) { _vertexProducer  = s; }

  protected:

    bool loadVertex(event_vertex *ev_vtx);

    bool lineCut(const double& x, const double&y);

    /// vertex coordinates
    std::vector<double> _vtx_w_cm;
    std::vector<double> _vtx_t_cm;

    // slope separation value
    double _slope_min;

    double _llt_min;

    // vector of x/y points on slope/intercept curve
    std::vector<double> _pts_x_v;
    std::vector<double> _pts_y_v;
    std::vector<double> _slope;
    std::vector<double> _intercept;
    
    /// max vtx-clus distance allowed
    double _dvtx_max;

    /// PROTON CUTS
    /// maximum proton length
    double _proton_dmax;
    /// maximum SSV
    double _ssv_max;

    /// ROI radius
    double _roi_rad;
    
    /// verbosity flag
    bool _verbose;
    /// debug flag
    bool _debug;

    /// conversion factors for hits
    double _wire2cm, _time2cm;

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
