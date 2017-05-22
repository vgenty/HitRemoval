#ifndef LARLITE_VERTEXSLOPECORRELATION_CXX
#define LARLITE_VERTEXSLOPECORRELATION_CXX

#include "VertexSlopeCorrelation.h"

namespace larlite {

  VertexSlopeCorrelation::VertexSlopeCorrelation()
    : HitRemovalBase()
  {

    _name        = "VertexSlopeCorrelation";
    _clusProducer = "";
    _vtxProducer  = "";
    
  }

  bool VertexSlopeCorrelation::initialize() {

    if (_tree) delete _tree;
    _tree = new TTree(_name.c_str(),_name.c_str());
    _tree->Branch("_ip",&_ip,"ip/D");
    _tree->Branch("_nhits",&_nhits,"nhits/I");

    return true;
  }
  
  bool VertexSlopeCorrelation::analyze(storage_manager* storage) {

    _event_watch.Start();

    auto ev_vtx  = storage->get_data<event_vertex>(_vtxProducer);
    auto ev_clus = storage->get_data<event_cluster>(_clusProducer);

    larlite::event_hit* ev_hit = nullptr;
    auto const& ass_cluster_hit_v = storage->find_one_ass(ev_clus->id(), ev_hit, ev_clus->name());

    if (!ev_hit){
      print(larlite::msg::kERROR,__FUNCTION__,"no hits");
      return false;
    }

    if (loadVertex(ev_vtx) == false) {
      print(larlite::msg::kERROR,__FUNCTION__,"num. vertices != 1");
      return false;
    }

    // select only clusters not previously removed
    auto const& clus_idx_v = AvailableClusterIndices(ev_hit, ass_cluster_hit_v);

    if (_verbose) {
      std::cout << std::endl << " NEW EVENT : " << clus_idx_v.size() << " CLUSTERS "
		<< std::endl << std::endl;
    }

    for (auto const& i : clus_idx_v) {

      auto hit_idx_v = ass_cluster_hit_v[i];

      if (hit_idx_v.size() < _min_nhits) continue;

      bool remove = false;

      int pl = ev_hit->at(hit_idx_v[0]).WireID().Plane;

      if (_verbose) std::cout << "Nhits = " << hit_idx_v.size() << " @ Plane " << pl << " @ idx " << i << std::endl;

      // get coordinates of hits to calculate linearity
      std::vector<double> hit_w_v;
      std::vector<double> hit_t_v;

      double dmin = 100000.;
      double dmax = 0.;

      for (auto const& hit_idx : hit_idx_v){
	double w = ev_hit->at(hit_idx).WireID().Wire  * _wire2cm;
	double t = ev_hit->at(hit_idx).PeakTime() * _time2cm;
	hit_w_v.push_back( w );
	hit_t_v.push_back( t );
	double d = ( w - _vtx_w_cm[pl] ) * ( w - _vtx_w_cm[pl] ) + ( t - _vtx_t_cm[pl] ) * ( t - _vtx_t_cm[pl] );
	if (d < dmin) { dmin = d; }
	if (d > dmax) { dmax = d; }
      }
      
      dmin = sqrt(dmin);
      dmax = sqrt(dmax);

      // is the closest point further than the ROI radius?
      // if so ignore cluster
      if (dmin > _roi_radius) {
	if (_verbose) { std::cout << " REMOVE because Out of ROI " << std::endl; }
	remove = true;
      }

      else {
	twodimtools::Linearity lin(hit_w_v,hit_t_v);
	
	// impact parameter to vertex:
	double x0 = _vtx_w_cm[pl];
	double y0 = _vtx_t_cm[pl];
	double IP    = ( - lin._slope * x0 + y0 - lin._intercept ) / sqrt( lin._slope * lin._slope + 1 );
	double IPup = ( - (lin._slope + lin._slope_err) * x0 + y0 - lin._intercept ) / sqrt( (lin._slope + lin._slope_err) * (lin._slope + lin._slope_err) + 1 );
	double IPdn = ( - (lin._slope - lin._slope_err) * x0 + y0 - lin._intercept ) / sqrt( (lin._slope - lin._slope_err) * (lin._slope - lin._slope_err) + 1 );
	
	if (_verbose){
	  std::cout << "Plane : " << pl << "\t N hits : " << hit_idx_v.size() << std::endl;
	  std::cout << "Vtx @ " << _vtx_w_cm[pl] << ", " << _vtx_t_cm[pl] << std::endl;
	  std::cout << "IP = " << IP << "\t IPup = " << IPup << "\t IPdn = " << IPdn << std::endl;
	}
	
	// if all 3 the same sign: vertex external to cone swept by slope +/- uncertainty
	if ( ( (IP * IPup) < 0 ) || ( (IP * IPdn) < 0 ) || ( (IPup * IPdn) < 0 ) )
	  continue;
	
	IP = fabs(IP);
	IPup = fabs(IPup);
	IPdn = fabs(IPdn);
	
	double IPmin = IP;
	if (IPup < IPmin) IPmin = IPup;
	if (IPdn < IPmin) IPmin = IPdn;

	_ip = IPmin;
	_nhits = hit_w_v.size();
	_tree->Fill();

	double IPmaxExp = _A * exp(- (_nhits + _xshift) / _fact ) + _yshift;
	
	if (IPmin > IPmaxExp) remove = true;
      }// else statement
      
      if (remove){
	if (_verbose) { std::cout << "\t\t REMOVE!" << std::endl; }
	for (auto const& hit_idx : hit_idx_v)
	  ev_hit->at(hit_idx).set_goodness(-1.0);
      }
      
    }// for all clusters

    _event_time += _event_watch.RealTime();
    _event_num  += 1;
    
    return true;
  }
  
}
#endif
