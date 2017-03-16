#ifndef LARLITE_REMOVEDELTARAYS_CXX
#define LARLITE_REMOVEDELTARAYS_CXX

#include "RemoveDeltaRays.h"
#include "LArUtil/GeometryHelper.h"
#include "LArUtil/Geometry.h"
#include "DataFormat/cluster.h"
#include "DataFormat/vertex.h"

namespace larlite {

  RemoveDeltaRays::RemoveDeltaRays() {

    _name        = "RemoveDeltaRays";
    _fout        = 0;
    _verbose     = false;
    _clusProducer = "";
    _vtxProducer  = "";
    
    _vtx_w_cm = {0,0,0};
    _vtx_t_cm = {0,0,0};
    
  }

  bool RemoveDeltaRays::initialize() {

    _wire2cm  = larutil::GeometryHelper::GetME()->WireToCm();
    _time2cm  = larutil::GeometryHelper::GetME()->TimeToCm();
    
    return true;
  }
  
  bool RemoveDeltaRays::analyze(storage_manager* storage) {

    auto ev_vtx  = storage->get_data<event_vertex>(_vtxProducer);
    auto ev_clus = storage->get_data<event_cluster>(_clusProducer);

    larlite::event_hit* ev_hit = nullptr;
    auto const& ass_cluster_hit_v = storage->find_one_ass(ev_clus->id(), ev_hit, ev_clus->name());

    if (!ev_hit){
      std::cout << "No hits!" << std::endl;
      return false;
    }
    
    // get vertex position on each plane
    if ( (ev_vtx->size() == 1) ){
      auto const& vtx = ev_vtx->at(0);
      auto geoH = larutil::GeometryHelper::GetME();
      std::vector<double> xyz = {vtx.X(), vtx.Y(), vtx.Z()};
      for (size_t pl = 0; pl < 3; pl++){
	auto const& pt = geoH->Point_3Dto2D(xyz,pl);
	_vtx_w_cm[pl] = pt.w;
	_vtx_t_cm[pl] = pt.t + 800 * _time2cm;
      }
    }

    // loop trhough each cluster
    // find the "start" and "end" point of the cluster
    // if both outside of the ROI
    // remove
    // if above some thresdhold, remove cluster

    for (size_t i=0; i < ass_cluster_hit_v.size(); i++){

      // store output cluster hit indices
      std::vector<unsigned int> out_cluster_hit_idx_v;
      larlite::cluster out_clus;

      auto hit_idx_v = ass_cluster_hit_v[i];

      bool remove = false;

      int pl = ev_hit->at(hit_idx_v[0]).WireID().Plane;

      //if ( (pl !=1) || (hit_idx_v.size() != 48) ) continue;

      //std::cout << "Plane : " << pl << "\t N hits : " << hit_idx_v.size() << std::endl;
      //std::cout << "Vtx @ " << _vtx_w_cm[pl] << ", " << _vtx_t_cm[pl] << std::endl;
      
      // get coordinates of hits to calculate linearity
      std::vector<double> hit_w_v;
      std::vector<double> hit_t_v;
      
      // minimum distance to vtx
      // and coordinates for closest point
      double dmin = 1000000.;
      double wmin = 0;
      double tmin = 0;
      
      for (auto const& hit_idx : hit_idx_v){
	double w = ev_hit->at(hit_idx).WireID().Wire  * _wire2cm;
	double t = ev_hit->at(hit_idx).PeakTime() * _time2cm;
	hit_w_v.push_back( w );
	hit_t_v.push_back( t );
	double d = ( w - _vtx_w_cm[pl] ) * ( w - _vtx_w_cm[pl] ) + ( t - _vtx_t_cm[pl] ) * ( t - _vtx_t_cm[pl] );
	if (d < dmin) { dmin = d; wmin = w; tmin = t; }
      }

      dmin = sqrt(dmin);

      if (dmin > _roi_radius) remove = true;
      else {
	
	twodimtools::Linearity lin(hit_w_v,hit_t_v);
	
	// choose point on line from linearity fit
	double ptx = 1e4;
	double pty = ptx * lin._slope + lin._intercept;
	
	// find closest and furthest point in cluster to this point
	size_t pt_close_idx = -1;
	double d_close = 1000000.;
	size_t pt_far_idx = -1;
	double d_far = 0.;
	
	for (size_t j=0; j < hit_w_v.size(); j++) {
	  double xx = hit_w_v[j];
	  double yy = hit_t_v[j];
	  double dd = sqrt((xx-ptx)*(xx-ptx) + (yy-pty)*(yy-pty));
	  if (dd > d_far)   { d_far   = dd; pt_far_idx   = j; }
	  if (dd < d_close) { d_close = dd; pt_close_idx = j; }
	}

	if ( (pt_close_idx < 0) || (pt_far_idx < 0) ) continue;

	//std::cout << "compare point " << ptx << ", " << pty << std::endl;
	//std::cout << "\t\t near "    << hit_w_v[pt_close_idx] << ", " << hit_t_v[pt_close_idx] << std::endl;
	//std::cout << "\t\t far  "    << hit_w_v[pt_far_idx]   << ", " << hit_t_v[pt_far_idx]   << std::endl;

	// are near and far points outside of ROI?
	
	d_close = sqrt( (_vtx_w_cm[pl] - hit_w_v[pt_close_idx]) * (_vtx_w_cm[pl] - hit_w_v[pt_close_idx]) +
			(_vtx_t_cm[pl] - hit_t_v[pt_close_idx]) * (_vtx_t_cm[pl] - hit_t_v[pt_close_idx]) );
	
	d_far = sqrt( (_vtx_w_cm[pl] - hit_w_v[pt_far_idx]) * (_vtx_w_cm[pl] - hit_w_v[pt_far_idx]) +
		      (_vtx_t_cm[pl] - hit_t_v[pt_far_idx]) * (_vtx_t_cm[pl] - hit_t_v[pt_far_idx]) );

	double l_track = sqrt( pow(hit_w_v[pt_far_idx] - hit_w_v[pt_close_idx], 2) +
			       pow(hit_t_v[pt_far_idx] - hit_t_v[pt_close_idx], 2) );

	//std::cout << "from vtx, near : " << d_close << ", " << "\t far : " << d_far << std::endl;
	
	if ( (d_close > _roi_radius) && (d_far > _roi_radius) )
	  remove = true;

	
	// impact parameter to vertex:
	double x0 = _vtx_w_cm[pl];
	double y0 = _vtx_t_cm[pl];
	double IP = fabs( - lin._slope * x0 + y0 - lin._intercept ) / sqrt( lin._slope * lin._slope + 1 );

	// max IP allowed, but one end of the track must be outside of ROI
	if ( (IP > 30.) && ( (d_close > _roi_radius) || (d_far > _roi_radius) ) ) remove = true;

	// construct triangle OAB with O = vertex, A & B the two ends of the cluster
	// angle AOB is the angle on which to cut. If large cluster uncorrelated
	// with vertex.
	//std::cout << "O =  [" << x0 << ", " << y0 << std::endl;
	//std::cout << "A =  [" << hit_w_v[pt_close_idx] << ", " << hit_t_v[pt_close_idx] << std::endl;
	//std::cout << "B =  [" << hit_w_v[pt_far_idx]   << ", " << hit_t_v[pt_far_idx]   << std::endl;
	double OAx = hit_w_v[pt_close_idx] - x0;
	double OAy = hit_t_v[pt_close_idx] - y0;
	double OAm = sqrt(OAx*OAx + OAy*OAy); // magnitude
	double OBx = hit_w_v[pt_far_idx]   - x0;
	double OBy = hit_t_v[pt_far_idx]   - y0;
	double OBm = sqrt(OBx*OBx + OBy*OBy); // magnitude
	double cos = (OAy*OBy + OAx*OBx) / (OAm * OBm);
	//std::cout << "dot product is " << cos << std::endl;
	double angle = 180. * acos(cos) / 3.14;
	//std::cout << "angle is " << angle << std::endl;
	
	if ( (angle > 50) && (dmin > 4.) && ( (d_close > _roi_radius) || (d_far > _roi_radius) ) )
	  remove = true;

	if ( (l_track > _roi_radius) && ( (d_close > _roi_radius) || (d_far > _roi_radius) ) )
	  remove = true;
	
      }
      
      if (remove){
	for (auto const& hit_idx : hit_idx_v)
	  ev_hit->at(hit_idx).set_goodness(-1.0);
      }
      
    }// for all clusters

  
    return true;
  }

  bool RemoveDeltaRays::finalize() {

  
    return true;
  }

}
#endif
