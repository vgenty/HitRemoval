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

#include "HitRemovalBase.h"

#include "TwoDimTools/Linearity.h"

namespace larlite {
  /**
     \class ROIRemoval
     User custom analysis class made by SHELL_USER_NAME
   */
  class ROIRemoval : public HitRemovalBase {
  
  public:

    /// Default constructor
    ROIRemoval();

    /// Default destructor
    virtual ~ROIRemoval(){}

    virtual bool initialize();

    virtual bool analyze(storage_manager* storage);

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
