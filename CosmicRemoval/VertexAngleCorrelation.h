/**
 * \file VertexAngleCorrelation.h
 *
 * \ingroup Clustering
 * 
 * \brief Class def header for a class VertexAngleCorrelation
 *
 * @author david caratelli
 */

/** \addtogroup Clustering

    @{*/

#ifndef LARLITE_VERTEXANGLECORRELATION_H
#define LARLITE_VERTEXANGLECORRELATION_H

#include "Analysis/ana_base.h"
#include "DataFormat/hit.h"

#include "HitRemovalBase.h"

#include "TwoDimTools/Linearity.h"

namespace larlite {
  /**
     \class VertexAngleCorrelation
     User custom analysis class made by SHELL_USER_NAME
   */
  class VertexAngleCorrelation : public ana_base, HitRemovalBase {
  
  public:

    /// Default constructor
    VertexAngleCorrelation();

    /// Default destructor
    virtual ~VertexAngleCorrelation(){}

    /** IMPLEMENT in VertexAngleCorrelation.cc!
        Initialization method to be called before the analysis event loop.
    */ 
    virtual bool initialize();

    /** IMPLEMENT in VertexAngleCorrelation.cc! 
        Analyze a data event-by-event  
    */
    virtual bool analyze(storage_manager* storage);

    /** IMPLEMENT in VertexAngleCorrelation.cc! 
        Finalize method to be called after all events processed.
    */
    virtual bool finalize();

    // set function to apply for cut
    void setCutFunction(double A, double xshift, double fact, double yshift)
    { _A = A; _xshift = xshift; _fact = fact; _yshift = yshift; }
    /// set maximum angle allowed
    void setMaxAngle(double a) { _angle_max = a; }
    
    /// Set Producers
    void setClusterProducer(std::string s) { _clusProducer = s; }
    void setVertexProducer (std::string s) { _vtxProducer  = s; }

  protected:

    /// Producers
    std::string _clusProducer;
    std::string _vtxProducer;

    double _A;
    double _xshift;
    double _fact;
    double _yshift;
    double _angle_max;

    int _nhits;
    double _ip, _angle;
    
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
