#ifndef LARLITE_VERTEXANGLECORRELATION_CXX
#define LARLITE_VERTEXANGLECORRELATION_CXX

#include "VertexAngleCorrelation.h"
#include "LArUtil/GeometryHelper.h"
#include "LArUtil/Geometry.h"
#include "DataFormat/cluster.h"
#include "DataFormat/vertex.h"

namespace larlite {

  VertexAngleCorrelation::VertexAngleCorrelation() {

    _name        = "VertexAngleCorrelation";
    _fout        = 0;
    _verbose     = false;
    _clusProducer = "";
    _vtxProducer  = "";
    
  }

  bool VertexAngleCorrelation::initialize() {

    if (_tree) delete _tree;
    _tree = new TTree(_name.c_str(),_name.c_str());
    _tree->Branch("_nhits",&_nhits,"nhits/I");
    _tree->Branch("_ip",&_ip,"ip/D");
    _tree->Branch("_angle",&_angle,"angle/D");

    return true;
  }
  
  bool VertexAngleCorrelation::analyze(storage_manager* storage) {

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

    // select only clusters not previously removed
    auto const& clus_idx_v = AvailableClusterIndices(ev_hit, ass_cluster_hit_v);

    for (auto const& i : clus_idx_v) {

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
      
      for (auto const& hit_idx : hit_idx_v){
	double w = ev_hit->at(hit_idx).WireID().Wire  * _wire2cm;
	double t = ev_hit->at(hit_idx).PeakTime() * _time2cm;
	hit_w_v.push_back( w );
	hit_t_v.push_back( t );
	double d = ( w - _vtx_w_cm[pl] ) * ( w - _vtx_w_cm[pl] ) + ( t - _vtx_t_cm[pl] ) * ( t - _vtx_t_cm[pl] );
      }

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
      
      // impact parameter to vertex:
      double x0 = _vtx_w_cm[pl];
      double y0 = _vtx_t_cm[pl];
      
      double IP = fabs( - lin._slope * x0 + y0 - lin._intercept ) / sqrt( lin._slope * lin._slope + 1 );
      
      if (_verbose) std::cout << "IP = " << IP << std::endl;
      
      // measure angle between :
      // O) vertex
      // A) one end of cluster
      // B) other end of cluster
      
      double OAx = hit_w_v[pt_far_idx] - x0;
      double OAy = hit_t_v[pt_far_idx] - y0;
      double OAm = sqrt( OAx*OAx + OAy*OAy);
      double OBx = hit_w_v[pt_close_idx] - x0;
      double OBy = hit_t_v[pt_close_idx] - y0;
      double OBm = sqrt( OBx*OBx + OBy*OBy);
      
      double cosSegment   = (OBy*OAy + OBx*OAx) / (OBm * OAm);
      double angleSegment = 180. * fabs(acos(cosSegment)) / 3.14;
      
      _ip = IP;
      _angle = angleSegment;
      _nhits = hit_w_v.size();
      
      _tree->Fill();
      
      double curve = _A * exp( - (_angle + _xshift) / _fact ) + _yshift;
      
      if (IP > curve) remove = true;
      if (_angle > _angle_max) remove = true;
      
      if (remove){
	if (_verbose) { std::cout << "\t\t REMOVE clus with angle = " << _angle << " and IP = " << IP << std::endl; }
	for (auto const& hit_idx : hit_idx_v)
	  ev_hit->at(hit_idx).set_goodness(-1.0);
      }
      
    }// for all clusters
    
    
    return true;
  }
  
  bool VertexAngleCorrelation::finalize() {

    _fout->cd();
    _tree->Write();
  
    return true;
  }


}
#endif
