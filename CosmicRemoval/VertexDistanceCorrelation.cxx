#ifndef LARLITE_VERTEXDISTANCECORRELATION_CXX
#define LARLITE_VERTEXDISTANCECORRELATION_CXX

#include "VertexDistanceCorrelation.h"
#include "LArUtil/GeometryHelper.h"
#include "LArUtil/Geometry.h"
#include "DataFormat/cluster.h"
#include "DataFormat/vertex.h"

namespace larlite {

  VertexDistanceCorrelation::VertexDistanceCorrelation() {

    _name        = "VertexDistanceCorrelation";
    _fout        = 0;
    _verbose     = false;
    _clusProducer = "";
    _vtxProducer  = "";
    
    _vtx_w_cm = {0,0,0};
    _vtx_t_cm = {0,0,0};
    
  }

  bool VertexDistanceCorrelation::initialize() {

    _wire2cm  = larutil::GeometryHelper::GetME()->WireToCm();
    _time2cm  = larutil::GeometryHelper::GetME()->TimeToCm();
    
    return true;
  }
  
  bool VertexDistanceCorrelation::analyze(storage_manager* storage) {

    auto ev_vtx  = storage->get_data<event_vertex>(_vtxProducer);
    auto ev_clus = storage->get_data<event_cluster>(_clusProducer);

    larlite::event_hit* ev_hit = nullptr;
    auto const& ass_cluster_hit_v = storage->find_one_ass(ev_clus->id(), ev_hit, ev_clus->name());

    if (!ev_hit){
      print(larlite::msg::kWARNING,__FUNCTION__,"no hits");
      return false;
    }

    if (loadVertex(ev_vtx) == false) {
      print(larlite::msg::kERROR,__FUNCTION__,"num. vertices != 1");
      return false;
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

      if (hit_idx_v.size() < 10) continue;

      bool remove = false;

      int pl = ev_hit->at(hit_idx_v[0]).WireID().Plane;

      //if ( (pl !=2) || (hit_idx_v.size() != 79) ) continue;

      if (_verbose){
	std::cout << "Plane : " << pl << "\t N hits : " << hit_idx_v.size() << std::endl;
	std::cout << "Vtx @ " << _vtx_w_cm[pl] << ", " << _vtx_t_cm[pl] << std::endl;
      }
      
      // get coordinates of hits to calculate linearity
      std::vector<double> hit_w_v;
      std::vector<double> hit_t_v;
      
      // minimum distance to vtx
      // and coordinates for closest point
      double dmin = 1000000.;
      double wmin = 0;
      double tmin = 0;

      // has this cluster already been removed? if so don't bother investigating
      if (ev_hit->at(0).GoodnessOfFit() < 0) continue;
      
      for (auto const& hit_idx : hit_idx_v){
	double w = ev_hit->at(hit_idx).WireID().Wire  * _wire2cm;
	double t = ev_hit->at(hit_idx).PeakTime() * _time2cm;
	hit_w_v.push_back( w );
	hit_t_v.push_back( t );
	double d = ( w - _vtx_w_cm[pl] ) * ( w - _vtx_w_cm[pl] ) + ( t - _vtx_t_cm[pl] ) * ( t - _vtx_t_cm[pl] );
	if (d < dmin) { dmin = d; wmin = w; tmin = t; }
      }

      dmin = sqrt(dmin);

      if (_verbose) std::cout << "dmin = " << dmin << std::endl;

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

	if (d_far < _roi_radius) continue;

	// are near and far points outside of ROI?
	
	d_close = sqrt( (_vtx_w_cm[pl] - hit_w_v[pt_close_idx]) * (_vtx_w_cm[pl] - hit_w_v[pt_close_idx]) +
			(_vtx_t_cm[pl] - hit_t_v[pt_close_idx]) * (_vtx_t_cm[pl] - hit_t_v[pt_close_idx]) );
	
	d_far = sqrt( (_vtx_w_cm[pl] - hit_w_v[pt_far_idx]) * (_vtx_w_cm[pl] - hit_w_v[pt_far_idx]) +
		      (_vtx_t_cm[pl] - hit_t_v[pt_far_idx]) * (_vtx_t_cm[pl] - hit_t_v[pt_far_idx]) );

	double dtmp;
	if (d_far < d_close) { dtmp = d_close; d_close = d_far; d_far = dtmp; }

	double l_track = sqrt( pow(hit_w_v[pt_far_idx] - hit_w_v[pt_close_idx], 2) +
			       pow(hit_t_v[pt_far_idx] - hit_t_v[pt_close_idx], 2) );

	if (_verbose) {
	  std::cout << "from vtx, near : " << d_close << ", " << "\t far : " << d_far << std::endl;
	  std::cout << "track length : " << l_track << std::endl;
	}
	
	if ( (d_close > _roi_radius) && (d_far > _roi_radius) )
	  remove = true;
	
	// impact parameter to vertex:
	double x0 = _vtx_w_cm[pl];
	double y0 = _vtx_t_cm[pl];

	// given the best-fit line and the vertex point
	// find the coordinate of the intersection between best-fit
	// and line perp. to best-fit passing through vertex
	double s1 = lin._slope;
	double m1 = lin._intercept;
	double s2 = -1./s1;
	double m2 = y0 + x0/s1;
	// incercept coordinates :
	double x1 = (m2-m1)/(s1-s2);
	double y1 = (s1*m2 - s2*m1) / (s1-s2);

	double IP = fabs( - lin._slope * x0 + y0 - lin._intercept ) / sqrt( lin._slope * lin._slope + 1 );

	if (_verbose) std::cout << "IP = " << IP << std::endl;

	// if angle between :
	// A) (x1,y1) -> point of closest approach
	// O) vtx
	// B) either end of the cluster line segment
	// is smaller than the angle between
	// B) O) and C) the other end of the cluster
	// then the IP point is ON the cluster line-segment
	// meaning that the line-segment passes close to the
	// vertex without actually being correlated

	double OAx = x1-x0;
	double OAy = y1-y0;
	double OAm = sqrt( OAx*OAx + OAy*OAy);
	double OBx = hit_w_v[pt_close_idx] - x0;
	double OBy = hit_t_v[pt_close_idx] - y0;
	double OBm = sqrt( OBx*OBx + OBy*OBy);
	double OCx = hit_w_v[pt_far_idx]   - x0;
	double OCy = hit_t_v[pt_far_idx]   - y0;
	double OCm = sqrt( OCx*OCx + OCy*OCy);

	double cosSegment   = (OBy*OCy + OBx*OCx) / (OBm * OCm);
	double angleSegment = 180. * fabs(acos(cosSegment)) / 3.14;
	double cosIP        = (OBy*OAy + OBx*OAx) / (OBm * OAm);
	double angleIP      = 180. * fabs(acos(cosIP)) / 3.14;

	if ( (angleIP < angleSegment) && (IP < 3 * d_close) && (d_close > 5.)) {
	  //     ( ( lin._local_lin_truncated_avg < _lin_max ) && ( (d_far > _dfar_min) && (d_close > _dclose_min) ) ) ) {
	  if (_verbose) { std::cout << "\t\t REMOVE" << std::endl; }
	  remove = true;
	}

      }
      
      if (remove){
	for (auto const& hit_idx : hit_idx_v)
	  ev_hit->at(hit_idx).set_goodness(-1.0);
      }
      
    }// for all clusters

  
    return true;
  }

  bool VertexDistanceCorrelation::finalize() {

  
    return true;
  }


  bool VertexDistanceCorrelation::loadVertex(event_vertex* ev_vtx) {
    
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

}
#endif
