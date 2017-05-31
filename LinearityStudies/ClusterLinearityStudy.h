/**
 * \file ClusterLinearityStudy.h
 *
 * \ingroup LinearityStudies
 * 
 * \brief Class def header for a class ClusterLinearityStudy
 *
 * @author david
 */

/** \addtogroup LinearityStudies

    @{*/

#ifndef LARLITE_CLUSTERLINEARITYSTUDY_H
#define LARLITE_CLUSTERLINEARITYSTUDY_H

#include "Analysis/ana_base.h"

#include "TwoDimTools/Linearity.h"

#include "DataFormat/vertex.h"

namespace larlite {
  /**
     \class ClusterLinearityStudy
     User custom analysis class made by SHELL_USER_NAME
   */
  class ClusterLinearityStudy : public ana_base{
  
  public:

    /// Default constructor
    ClusterLinearityStudy();

    /// Default destructor
    virtual ~ClusterLinearityStudy(){}

    /** IMPLEMENT in ClusterLinearityStudy.cc!
        Initialization method to be called before the analysis event loop.
    */ 
    virtual bool initialize();

    /** IMPLEMENT in ClusterLinearityStudy.cc! 
        Analyze a data event-by-event  
    */
    virtual bool analyze(storage_manager* storage);

    /** IMPLEMENT in ClusterLinearityStudy.cc! 
        Finalize method to be called after all events processed.
    */
    virtual bool finalize();

    void SetClusterProducer(std::string s) { _clusterProducer = s; }
    void SetVertexProducer(std::string s) { _vertexProducer = s; }

  protected:

    bool loadVertex(event_vertex *ev_vtx);

    std::string _clusterProducer, _vertexProducer;

    TTree* _tree;
    int    _pl;
    int    _nhits;
    double _lin;
    double _ssv;
    double _local_lin_truncated;
    double _local_lin_avg;
    double _slope;
    double _dvtx_max, _dvtx_min;

    /// conversion factors for hits
    double _wire2cm, _time2cm;

    std::vector<double> _vtx_w_cm, _vtx_t_cm;
    
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
