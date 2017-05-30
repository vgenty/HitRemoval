#ifndef LARLITE_MAKEROI_CXX
#define LARLITE_MAKEROI_CXX

#include "MakeROI.h"

#include "DataFormat/vertex.h"
#include "DataFormat/roi.h"

namespace larlite {

  bool MakeROI::initialize() {

    
    return true;
  }
  
  bool MakeROI::analyze(storage_manager* storage) {
  
    auto ev_vtx = storage->get_data<event_vertex>(_vtx_producer);
    auto ev_roi = storage->get_data<event_roi>("ROI");

    storage->set_id(storage->run_id(),storage->subrun_id(),storage->event_id());
    
    auto const& vtx = ev_vtx->at(0);
    
    larlite::roi ROI(vtx,100.);

    ev_roi->emplace_back(ROI);
    
    return true;
  }

  bool MakeROI::finalize() {

    
    return true;
  }

}
#endif
