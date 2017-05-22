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

#ifndef LARLITE_HITREMOVALBASE_H
#define LARLITE_HITREMOVALBASE_H

#include "Analysis/ana_base.h"

#include "LArUtil/GeometryHelper.h"
#include "LArUtil/Geometry.h"

#include "DataFormat/hit.h"
#include "DataFormat/vertex.h"
#include "DataFormat/cluster.h"

#include "TwoDimTools/Linearity.h"

#include <TStopwatch.h>

#include "TTree.h"

namespace larlite {

  /**
     \class HitRemovalBase
     User defined class HitRemovalBase ... these comments are used to generate
     doxygen documentation!
  */
  class HitRemovalBase : public ana_base {
    
  public:
    
    /// Default constructor
    HitRemovalBase();
    
    /// Default destructor
    ~HitRemovalBase(){}

    bool finalize();

    /// Verbosity setter
    void setVerbose(bool on) { _verbose = on; }


  protected:

    TStopwatch _event_watch;
    double     _event_time;
    int        _event_num;

    std::vector<unsigned int> AvailableClusterIndices(const larlite::event_hit* ev_hit,
						      const std::vector< std::vector<unsigned int> >& clus_idx_v);

    double ImpactParameter(const twodimtools::Linearity& lin, const int& pl);

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

