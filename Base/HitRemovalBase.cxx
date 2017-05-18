#ifndef HITREMOVALBASE_CXX
#define HITREMOVALBASE_CXX

#include "HitRemovalBase.h"

namespace larlite {

  HitRemovalBase::HitRemovalBase()
    : _tree(nullptr)
  {

    _vtx_w_cm = {0.,0.,0.};
    _vtx_t_cm = {0.,0.,0.};

    _wire2cm  = larutil::GeometryHelper::GetME()->WireToCm();
    _time2cm  = larutil::GeometryHelper::GetME()->TimeToCm();

  }
 
  bool HitRemovalBase::loadVertex(event_vertex* ev_vtx) {
    
    if (ev_vtx->size() != 1) return false;
    
    // get vertex position on each plane
    if ( (ev_vtx->size() == 1) ){
      auto const& vtx = ev_vtx->at(0);
      auto geoH = larutil::GeometryHelper::GetME();
      auto geom = larutil::Geometry::GetME();
      std::vector<double> xyz = {vtx.X(), vtx.Y(), vtx.Z()};
      for (size_t pl = 0; pl < 3; pl++){
	double *origin;
	origin = new double[3];
	geom->PlaneOriginVtx(pl,origin);
	auto const& pt = geoH->Point_3Dto2D(xyz,pl);
	_vtx_w_cm[pl] = pt.w;
	_vtx_t_cm[pl] = pt.t + 800 * _time2cm - origin[0];
      }
    }    

    return true;
  }

  std::vector<unsigned int> HitRemovalBase::AvailableClusterIndices(const larlite::event_hit* ev_hit,
								    const std::vector< std::vector<unsigned int> >& clus_idx_v)

  {

    std::vector<unsigned int> return_indices;
    
    // assumes that 1st GoodnessOfFit value is indicative of entire cluster

    for (size_t i=0; i < clus_idx_v.size(); i++) {

      auto idx_v = clus_idx_v.at(i);
    
      if (idx_v.size() == 0) continue;
      
      if (ev_hit->at(idx_v.at(0)).GoodnessOfFit() > 0)
	return_indices.push_back( i );
      
    }// for all clustrs
    
    return return_indices;
  }

}

#endif
