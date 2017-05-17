/**
 * \file ROIRemoval.h
 *
 * \ingroup Clustering
 * 
 * \brief Class def header for a class ROIRemoval
 *
 * @author david caratelli
 */

/** \addtogroup Clustering

    @{*/

#ifndef LARLITE_ROIREMOVAL_H
#define LARLITE_ROIREMOVAL_H

#include "Analysis/ana_base.h"
#include "DataFormat/hit.h"

#include "HitRemovalBase.h"

#include "TwoDimTools/Linearity.h"

namespace larlite {
  /**
     \class ROIRemoval
     User custom analysis class made by SHELL_USER_NAME
   */
  class ROIRemoval : public ana_base, HitRemovalBase {
  
  public:

    /// Default constructor
    ROIRemoval();

    /// Default destructor
    virtual ~ROIRemoval(){}

    /** IMPLEMENT in ROIRemoval.cc!
        Initialization method to be called before the analysis event loop.
    */ 
    virtual bool initialize();

    /** IMPLEMENT in ROIRemoval.cc! 
        Analyze a data event-by-event  
    */
    virtual bool analyze(storage_manager* storage);

    /** IMPLEMENT in ROIRemoval.cc! 
        Finalize method to be called after all events processed.
    */
    virtual bool finalize();

    void setROI(double r) { _roi = r; }

    /// Set Producers
    void setClusterProducer(std::string s) { _clusProducer = s; }
    void setVertexProducer(std::string s) { _vertexProducer = s; }

  protected:
    
    /// Producers
    std::string _clusProducer;
    std::string _vertexProducer;

    // ROI box
    double _roi;
    
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
