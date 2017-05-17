/**
 * \file HitRemovalBase.h
 *
 * \ingroup Base
 * 
 * \brief Class def header for a class HitRemovalBase
 *
 * @author david caratelli
 */

/** \addtogroup Base

    @{*/

#ifndef HITREMOVALBASE_H
#define HITREMOVALBASE_H

#include <iostream>

#include "LArUtil/GeometryHelper.h"
#include "LArUtil/Geometry.h"

#include "DataFormat/hit.h"
#include "DataFormat/vertex.h"
#include "DataFormat/cluster.h"

#include "TTree.h"

namespace larlite {

  /**
     \class HitRemovalBase
     User defined class HitRemovalBase ... these comments are used to generate
     doxygen documentation!
  */
  class HitRemovalBase{
    
  public:
    
    /// Default constructor
    HitRemovalBase();
    
    /// Default destructor
    ~HitRemovalBase(){}

    /// Verbosity setter
    void setVerbose(bool on) { _verbose = on; }


  protected:

    bool loadVertex(event_vertex *ev_vtx);

    /// vertex coordinates
    std::vector<double> _vtx_w_cm;
    std::vector<double> _vtx_t_cm;

    /// conversion factors for hits
    double _wire2cm, _time2cm;

    bool _verbose;

    TTree* _tree;
    
  };

}

#endif
/** @} */ // end of doxygen group 

