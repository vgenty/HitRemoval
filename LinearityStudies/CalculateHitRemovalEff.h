/**
 * \file CalculateHitRemovalEff.h
 *
 * \ingroup LinearityStudies
 * 
 * \brief Class def header for a class CalculateHitRemovalEff
 *
 * @author david caratelli
 */

/** \addtogroup LinearityStudies

    @{*/

#ifndef LARLITE_CALCULATEHITREMOVALEFF_H
#define LARLITE_CALCULATEHITREMOVALEFF_H

#include "Analysis/ana_base.h"
#include "TTree.h"

namespace larlite {
  /**
     \class CalculateHitRemovalEff
     User custom analysis class made by SHELL_USER_NAME
   */
  class CalculateHitRemovalEff : public ana_base{
  
  public:

    /// Default constructor
    CalculateHitRemovalEff();

    /// Default destructor
    virtual ~CalculateHitRemovalEff(){}

    /** IMPLEMENT in CalculateHitRemovalEff.cc!
        Initialization method to be called before the analysis event loop.
    */ 
    virtual bool initialize();

    /** IMPLEMENT in CalculateHitRemovalEff.cc! 
        Analyze a data event-by-event  
    */
    virtual bool analyze(storage_manager* storage);

    /** IMPLEMENT in CalculateHitRemovalEff.cc! 
        Finalize method to be called after all events processed.
    */
    virtual bool finalize();

    void setUseTruth(bool on) { _use_truth = on; }

    void setUseCluster(bool on) { _use_cluster = on; }

    void setVertexProducer (std::string s) { _vertexProducer  = s; }

    /// set ROI size
    void setROI(double r) { _roi = r; }
    /// set amplitude cut for hits [hits < threshold do not pass hit-removal]
    void setHitAmpCut(double a) { _hit_amp = a; }

  protected:

    bool loadVertex(event_vertex *ev_vtx);

    bool _use_cluster;

    /// conversion factors for hits
    double _wire2cm, _time2cm;

    TTree* _tree;
    int _event, _run, _entry;
    double _qtot, _qtot_0, _qtot_1, _qtot_2;
    double _qremoved, _qremoved_0, _qremoved_1, _qremoved_2;
    int    _ntot, _ntot_0, _ntot_1, _ntot_2;
    int    _nremoved, _nremoved_0, _nremoved_1, _nremoved_2;
    double _frac, _frac_0, _frac_1, _frac_2;
    int _nc, _pi0;
    double _pi0E;
    // number of MCshwoers in event
    int _nshr;
    // total Edep by MCshowers
    double _edepshr;

    std::vector<int> _nclus_0_v, _nclus_1_v, _nclus_2_v;

    bool _use_truth;

    double _roi;
    double _hit_amp;

    std::vector<double> _vtx_w_cm, _vtx_t_cm;
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
