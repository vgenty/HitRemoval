/**
 * \file VertexDistanceCorrelation.h
 *
 * \ingroup Clustering
 * 
 * \brief Class def header for a class VertexDistanceCorrelation
 *
 * @author david caratelli
 */

/** \addtogroup Clustering

    @{*/

#ifndef LARLITE_VERTEXDISTANCECORRELATION_H
#define LARLITE_VERTEXDISTANCECORRELATION_H

#include "HitRemovalBase.h"

#include "TwoDimTools/Linearity.h"

namespace larlite {
  /**
     \class VertexDistanceCorrelation
     User custom analysis class made by SHELL_USER_NAME
   */
  class VertexDistanceCorrelation : public HitRemovalBase {
  
  public:

    /// Default constructor
    VertexDistanceCorrelation();

    /// Default destructor
    virtual ~VertexDistanceCorrelation(){}

    /** IMPLEMENT in VertexDistanceCorrelation.cc!
        Initialization method to be called before the analysis event loop.
    */ 
    virtual bool initialize();

    /** IMPLEMENT in VertexDistanceCorrelation.cc! 
        Analyze a data event-by-event  
    */
    virtual bool analyze(storage_manager* storage);

    /// set ROI size (cm away from vtx)
    void setROIRadius(double r) { _roi_radius = r; }

    /// set max linearity for cluster
    void setMaxLin(double l) { _lin_max = l; }

    /// minimum "far" and "close" distance
    void setDFarMin  (double d) { _dfar_min = d;   }
    void setDCloseMin(double d) { _dclose_min = d; }

    /// Set Producers
    void setClusterProducer(std::string s) { _clusProducer = s; }
    void setVertexProducer (std::string s) { _vtxProducer  = s; }

  protected:

    /// Producers
    std::string _clusProducer;
    std::string _vtxProducer;

    /// ROI radius
    double _roi_radius;
    /// max lin
    double _lin_max;
    /// min dfar / dclose
    double _dfar_min;
    double _dclose_min;
    
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
