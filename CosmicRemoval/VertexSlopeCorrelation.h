/**
 * \file VertexSlopeCorrelation.h
 *
 * \ingroup Clustering
 * 
 * \brief Class def header for a class VertexSlopeCorrelation
 *
 * @author david caratelli
 */

/** \addtogroup Clustering

    @{*/

#ifndef LARLITE_VERTEXSLOPECORRELATION_H
#define LARLITE_VERTEXSLOPECORRELATION_H

#include "Analysis/ana_base.h"
#include "DataFormat/hit.h"

#include "HitRemovalBase.h"

#include "TwoDimTools/Linearity.h"

namespace larlite {
  /**
     \class VertexSlopeCorrelation
     User custom analysis class made by SHELL_USER_NAME
   */
  class VertexSlopeCorrelation : public ana_base, HitRemovalBase {
  
  public:

    /// Default constructor
    VertexSlopeCorrelation();

    /// Default destructor
    virtual ~VertexSlopeCorrelation(){}

    /** IMPLEMENT in VertexSlopeCorrelation.cc!
        Initialization method to be called before the analysis event loop.
    */ 
    virtual bool initialize();

    /** IMPLEMENT in VertexSlopeCorrelation.cc! 
        Analyze a data event-by-event  
    */
    virtual bool analyze(storage_manager* storage);

    /** IMPLEMENT in VertexSlopeCorrelation.cc! 
        Finalize method to be called after all events processed.
    */
    virtual bool finalize();

    // set function to apply for cut
    void setCutFunction(double A, double xshift, double fact, double yshift)
    {_A = A; _xshift = xshift; _fact = fact; _yshift = yshift; }
    // set min number of hits for cluster to be considered
    void setMinNHits(int n) { _min_nhits = n; }
    // set ROI radius (cm)
    void setROIRadius(double r) { _roi_radius = r; }

    /// Set Producers
    void setClusterProducer(std::string s) { _clusProducer = s; }
    void setVertexProducer (std::string s) { _vtxProducer  = s; }

  protected:

    /// Producers
    std::string _clusProducer;
    std::string _vtxProducer;

    int _min_nhits;
    double _roi_radius;

    double _A;
    double _xshift;
    double _fact;
    double _yshift;

    double _ip;
    int    _nhits;

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
