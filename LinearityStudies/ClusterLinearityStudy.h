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

namespace larlite {
  /**
     \class ClusterLinearityStudy
     User custom analysis class made by SHELL_USER_NAME
   */
  class ClusterLinearityStudy : public ana_base{
  
  public:

    /// Default constructor
    ClusterLinearityStudy()
      : _tree(nullptr)
    { _name="ClusterLinearityStudy"; _fout=0;}

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

  protected:

    std::string _clusterProducer;

    TTree* _tree;
    int    _pl;
    int    _nhits;
    double _lin;
    double _local_lin_truncated;
    double _local_lin_avg;
    double _slope;

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
