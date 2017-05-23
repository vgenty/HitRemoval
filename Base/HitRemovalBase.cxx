#ifndef HITREMOVALBASE_CXX
#define HITREMOVALBASE_CXX

#include "HitRemovalBase.h"
#include <iomanip>

namespace larlite {

  HitRemovalBase::HitRemovalBase()
    : _tree(nullptr)
  {

    _vtx_w_cm = {0.,0.,0.};
    _vtx_t_cm = {0.,0.,0.};

    _wire2cm  = larutil::GeometryHelper::GetME()->WireToCm();
    _time2cm  = larutil::GeometryHelper::GetME()->TimeToCm();

    _fout = 0;
    _name = "HitRemovalBase";

    _event_time = 0.;
    _event_num  = 0;

  }

  bool HitRemovalBase::finalize() {

    if (_fout) _fout->cd();
    if (_tree) _tree->Write();

    double time_per_event = _event_time / _event_num;
    std::cout << "Algo " << std::setw(25) << _name << " : "
	      << time_per_event * 1.e3 << " [ms/event]" << std::endl;

    return true;
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
    
    for (size_t i=0; i < clus_idx_v.size(); i++) {

      auto idx_v = clus_idx_v.at(i);

      if (idx_v.size() == 0) continue;

      bool removed = true;

      for (auto const& idx : idx_v) {
	if (ev_hit->at(idx).GoodnessOfFit() > 0) {
	  removed = false;
	  break;
	}
      }// for all hits in cluster

      if (removed == false)
	return_indices.push_back( i );
      
    }// for all clustrs
    
    return return_indices;
  }

  double HitRemovalBase::ImpactParameter(const twodimtools::Linearity& lin,
					 const int& pl) {
    
    // impact parameter to vertex:
    double x0 = _vtx_w_cm[pl];
    double y0 = _vtx_t_cm[pl];
    double IP    = ( - lin._slope * x0 + y0 - lin._intercept ) / sqrt( lin._slope * lin._slope + 1 );

    return IP;
  }


  BBox HitRemovalBase::GetBBox(const std::vector<unsigned int>& hit_idx_v,
			       larlite::event_hit* ev_hit) {

    // determine cluster's bbox
    double wmin, tmin, wmax, tmax;
    wmin = tmin = 10000.;
    wmax = tmax = 0.;
    
    for (auto const& hit_idx : hit_idx_v){
      double w = ev_hit->at(hit_idx).WireID().Wire  * _wire2cm;
      double t = ev_hit->at(hit_idx).PeakTime() * _time2cm;
      if (w > wmax) wmax = w;
      if (w < wmin) wmin = w;
      if (t > tmax) tmax = t;
      if (t < tmin) tmin = t;
    }// for all hits
    
    BBox box;
    box.wmin = wmin;
    box.wmax = wmax;
    box.tmin = tmin;
    box.tmax = tmax;

    return box;
  }


  bool HitRemovalBase::Intersect(const BBox& box, const double& radius, const int& pl) {

    BBox roi;
    
    roi.wmin = _vtx_w_cm[pl] - radius;
    roi.wmax = _vtx_w_cm[pl] + radius;
    roi.tmin = _vtx_t_cm[pl] - radius;
    roi.tmax = _vtx_t_cm[pl] + radius;

    if ( (box.wmax < roi.wmin) || (box.wmin > roi.wmax) ) return false;
    if ( (box.tmax < roi.tmin) || (box.tmin > roi.tmax) ) return false;

    return true;
  }

}

#endif
