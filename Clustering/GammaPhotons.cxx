#ifndef LARLITE_GAMMAPHOTONS_CXX
#define LARLITE_GAMMAPHOTONS_CXX

#include "GammaPhotons.h"
#include "LArUtil/GeometryHelper.h"
#include "LArUtil/Geometry.h"
#include "DataFormat/cluster.h"
#include "DataFormat/vertex.h"

#include "DataFormat/shower.h"
#include "DataFormat/pfpart.h"

#include "TwoDimTools/Linearity.h"

namespace larlite {

  GammaPhotons::GammaPhotons(){

    _name        = "GammaPhotons";
    _fout        = 0;
    _verbose     = false;
    _hitProducer = "";
    _out_clusterProducer = "rawclus";
    _radius      = 2.0;
    _cellSize    = 2;
    _max_rms     = 100;
    _max_clus_size = 100;
    _min_lin     = 0.1;
    _vtx_w_cm = {0,0,0};
    _vtx_t_cm = {0,0,0};
    _tick_min = 0;
    _tick_max = 9600;
    _roi_radius = 1000;

  }

  bool GammaPhotons::initialize() {

    _wire2cm  = larutil::GeometryHelper::GetME()->WireToCm();
    _time2cm  = larutil::GeometryHelper::GetME()->TimeToCm();

    std::cout << "********************************" << std::endl;
    std::cout << "Wire -> cm conversion : " << _wire2cm << std::endl;
    std::cout << "Time -> cm conversion : " << _time2cm << std::endl;
    std::cout << "********************************" << std::endl;

    return true;
  }
  
  bool GammaPhotons::analyze(storage_manager* storage) {

    auto evt_hits      = storage->get_data<event_hit>(_hitProducer);
    auto ev_clusters   = storage->get_data<event_cluster>(_out_clusterProducer);
    auto cluster_ass_v = storage->get_data<event_ass>(ev_clusters->name());
    auto ev_shower     = storage->get_data<event_shower>("showerreco");

    //set event ID through storage manager
    storage->set_id(storage->run_id(),storage->subrun_id(),storage->event_id());

    // grab associated pfparticles
    larlite::event_pfpart* ev_pfpart = nullptr;
    // grab associated clusters
    larlite::event_cluster* ev_cluster = nullptr;
    // grab associated hits
    larlite::event_hit* ev_hit         = nullptr;

    if (!ev_shower) { std::cout << "no showers...quit" << std::endl; return false; }
    auto const& ass_shr_pfpart_v  = storage->find_one_ass(ev_shower->id(), ev_pfpart,  ev_shower->name());
    if (!ev_pfpart) { std::cout << "no pfparts...quit" << std::endl; return false; }
    auto const& ass_pfpart_clus_v = storage->find_one_ass(ev_pfpart->id(), ev_cluster, ev_pfpart->name());
    if (!ev_cluster) { std::cout << "no clusters...quit" << std::endl; return false; }
    auto const& ass_clus_hit_v = storage->find_one_ass(ev_cluster->id(), ev_hit, ev_cluster->name());


    auto shower_hit_idx_v = getShowerHits(ass_shr_pfpart_v,ass_pfpart_clus_v,ass_clus_hit_v,ev_hit);

    if (shower_hit_idx_v.size() != ev_shower->size()) {
      std::cout << "different number of showers and hit vectors! exit." << std::endl;
      return false;
    }

    if (!evt_hits){
      std::cout << "No hits!" << std::endl;
      return false;
    }

    // a map to connect hit index wih a cluster index
    // each hit gets a cluster index
    // _clusterMap[hit_index] -> cluster_index
    std::map<size_t, size_t> _clusterMap;
    // a map to connect the cluster index with the vector of hit indices for that cluster
    // _clusters[index] -> vector of hit indices for that cluster
    std::map<size_t,std::vector<size_t> > _clusters;

    // keep track of largest cluster ID created
    size_t maxClusterID = 0;

    for (int pl=2; pl < 3; pl++){

      // hit map will only contain hits we want to use for clustering
      MakeHitMap(evt_hits,pl,shower_hit_idx_v);
      
      // iterator for hit cell map
      std::map<std::pair<int,int>, std::vector<size_t> >::iterator it;
      
      // loop through hits in each cell to find matches
      for (it = _hitMap.begin(); it != _hitMap.end(); it++){

	// pair = (i,j) indices of this cell in the _hitMap
	auto const& pair = it->first;
	
	// wire-space cell index
	// prepare a hit list of all neighboring cells
	// _________
	// |__|__|__|
	// |__|__|__|
	// |__|__|__|
	std::vector<size_t> cellhits = it->second;

	std::vector<size_t> neighborhits;
	getNeighboringHits(pair,neighborhits);

	for (size_t h1=0; h1 < cellhits.size(); h1++){
	  // has this hit been added to a cluster?
	  // if so not necessary to look at
	  auto const& hit1 = cellhits[h1];
	  // keep track if the hit will ever be matched to another
	  bool matched = false;
	  // if not find hits it should be clustered with and add it to the appropriate cluster
	  for (size_t h2=0; h2 < neighborhits.size(); h2++){
	    auto const& hit2 = neighborhits[h2];
	    if (hit1 == hit2) continue;
	    // are the hits compatible?
	    bool compat = HitsCompatible(evt_hits->at(hit1), evt_hits->at(hit2));
	    // should the hits go in the same cluster?
	    if (compat){
	      matched = true;
	      // if both hits have already been assigned to a cluster then we can merge the cluster indices!
	      if ( (_clusterMap.find(hit1) != _clusterMap.end()) and
		   (_clusterMap.find(hit2) != _clusterMap.end()) ){
		// if in the same cluster -> do nothing
		// if they are in different clusters:
		if (_clusterMap[hit1] != _clusterMap[hit2]){
		  auto idx1 = _clusterMap[hit1];
		  auto idx2 = _clusterMap[hit2];
		  // hit indices for 1st cluster:
		  auto hits1 = _clusters[idx1];
		  auto hits2 = _clusters[idx2];
		  // append hits2 to hits1
		  for (auto h : hits2){
		    hits1.push_back(h);
		    // also change the index that the hit goes to (idx1 instead of idx2)
		    _clusterMap[h] = idx1;
		  }
		  _clusters[idx1] = hits1;
		  // erase cluster @ index2
		  _clusters.erase(idx2);
		}// if they are in different clusters
	      }
	      // if compatible and the 2nd hit has been added to a cluster
	      // add hit1 to the same cluster
	      else if ( (_clusterMap.find(hit2) != _clusterMap.end()) and
			(_clusterMap.find(hit1) == _clusterMap.end()) ){
		auto clusIdx = _clusterMap[hit2];
		_clusterMap[hit1] = clusIdx;
		_clusters[clusIdx].push_back(hit1);
	      }
	      // otherwise, add both to a new cluster
	      else if ( (_clusterMap.find(hit1) != _clusterMap.end()) and
			(_clusterMap.find(hit2) == _clusterMap.end()) ){
		auto clusIdx = _clusterMap[hit1];
		_clusterMap[hit2] = clusIdx;
		_clusters[clusIdx].push_back(hit2);
	      }
	      // if neither has a cluster yet
	      else{
		// create a new cluster for this match
		_clusterMap[hit1] = maxClusterID;
		_clusterMap[hit2] = maxClusterID;
		std::vector<size_t> cl = {hit1,hit2};
		_clusters[maxClusterID] = cl;
		maxClusterID += 1;
	      }
	    }// if the two hits are compatible
	  }// 2nd loop through hits in the cell
	  // has this hit been matched? if not we still need to add it as its own cluster
	  if (matched == false){
	    _clusterMap[hit1] = maxClusterID;
	    _clusters[maxClusterID] = {hit1};
	    maxClusterID += 1;
	  }
	}// 1st loop through hits in the cell
      }// loop through all cells

    }// loop through all planes

    // make a vector for the clusters
    std::vector<std::vector<unsigned int> > cluster_vector;
    for (auto it = _clusters.begin(); it != _clusters.end(); it++){
      auto indices = it->second;
      // if there are enough indices, make a cluster
      if (indices.size() >= 1){
	std::vector<unsigned int> clus;
	for (auto idx : indices)
	  clus.push_back(idx);
	cluster_vector.push_back(clus);
      }// if there are 2 hits in cluster
    }
    
    // vector for assocaitions
    std::vector<std::vector<unsigned int> > _cluster_hit_ass;
    // for each cluster create a larlite::cluster
    for (size_t i=0; i < cluster_vector.size(); i++){
      
      if (cluster_vector[i].size() <= 0) continue;

      if (cluster_vector[i].size() > _max_clus_size) continue;
      
      larlite::cluster clus;
      clus.set_n_hits(cluster_vector[i].size());
      // grab wire/time bounds for this cluster
      float w_min = 9999;
      float t_min = 9999;
      float w_max = -1;
      float t_max = -1;
      float integral = 0;
      // grab list of hits with which to compute linearity
      std::vector<double> clus_wire_v;
      std::vector<double> clus_time_v;
      for (auto const& idx : cluster_vector[i]){
	auto const& hit = evt_hits->at(idx);
	float t = hit.PeakTime();
	if (t > t_max) t_max = t;
	if (t < t_min) t_min = t;
	float w = hit.WireID().Wire;
	if (w > w_max) w_max = w;
	if (w < w_min) w_min = w;
	integral += hit.Integral();
	clus_time_v.push_back( hit.PeakTime()    * _time2cm );
	clus_wire_v.push_back( hit.WireID().Wire * _wire2cm );
      }
      ::twodimtools::Linearity lin(clus_wire_v,clus_time_v);
      if (lin._local_lin_truncated_avg < _min_lin) continue;
      if ((lin._local_lin_truncated_avg == 1.0) and (clus_wire_v.size() >= 5)) continue;
      clus.set_start_wire(w_min,1.);
      clus.set_end_wire(w_max,1.);
      clus.set_start_tick(t_min,1.);
      clus.set_end_tick(t_max,1.);
      clus.set_integral(integral,0.,0.);
      clus.set_summedADC(integral, 0., 0.);
      clus.set_view( evt_hits->at( cluster_vector[i][0] ).View() );
      //clus.set_planeID( evt_hits->at( cluster_vector[i][0] ).WireID().Plane );
      //clus.set_planeID(2);
      // vector for associations
      ev_clusters->push_back(clus);
      _cluster_hit_ass.push_back(cluster_vector[i]);
      
    }// for all created clusters
    
    cluster_ass_v->set_association(ev_clusters->id(),product_id(data::kHit,evt_hits->name()),_cluster_hit_ass);    

    return true;
  }

  bool GammaPhotons::finalize() {

    return true;
  }


  std::vector< std::vector<size_t> > GammaPhotons::getShowerHits(const std::vector<std::vector<unsigned int> >& ass_shr_pfp_v,
								 const std::vector<std::vector<unsigned int> >& ass_pfp_cls_v,
								 const std::vector<std::vector<unsigned int> >& ass_cls_hit_v,
								 larlite::event_hit* ev_hit) {

    std::vector< std::vector<size_t> > shower_hits_v;

    for (auto const& shr_pfp : ass_shr_pfp_v) {

      if (shr_pfp.size() != 1) {
	std::cout << "shower associated with != PFP. Continue." << std::endl;
        continue;
      }

      auto const& pfp_idx = shr_pfp[0];

      // grab clusters associated with PFP                                                                                                                                                                                 
      auto const& clus_idx_v = ass_pfp_cls_v.at(pfp_idx);

      // loop through clusters, add hits which belong to the collection plane                                                                                                                                              
      for (auto const& clus_idx : clus_idx_v) {

        auto const& hit_idx_v = ass_cls_hit_v.at(clus_idx);

        // collection plane?                                                                                                                                                                                               
        if (ev_hit->at( hit_idx_v.at(0) ).WireID().Plane != 2) continue;

	std::vector<size_t> hits_v;

        for (auto const& hit_idx : hit_idx_v)
          hits_v.push_back( hit_idx );

        shower_hits_v.push_back( hits_v );

      }// for all cluster indices                                                                                                                                                                                          

    }// for all showers                                                                                                                                                                                                    

    return shower_hits_v;
  }


  // get all hits from neighboring cells
  void GammaPhotons::getNeighboringHits(const std::pair<int,int>& pair, std::vector<size_t>& hitIndices){
   
    auto const& i       = pair.first;
    // time-space cell index
    auto const& j       = pair.second;

    // _________
    // |__|__|__|
    // |__|XX|__|
    // |__|__|__|
    if (_hitMap.find(std::make_pair(i,j)) != _hitMap.end()){
      for (auto &h : _hitMap[std::make_pair(i,j)])
	hitIndices.push_back(h);
    }

    // now look at neighboring cells, if they exist
    // _________
    // |__|__|__|
    // |XX|__|__|
    // |__|__|__|
    if (_hitMap.find(std::make_pair(i-1,j)) != _hitMap.end()){
      for (auto &h : _hitMap[std::make_pair(i-1,j)])
	hitIndices.push_back(h);
    }
    // _________
    // |__|__|__|
    // |__|__|__|
    // |__|XX|__|
    if (_hitMap.find(std::make_pair(i,j-1)) != _hitMap.end()){
      for (auto &h : _hitMap[std::make_pair(i,j-1)])
	hitIndices.push_back(h);
    }
    // _________
    // |__|__|__|
    // |__|__|__|
    // |XX|__|__|
    if ( _hitMap.find(std::make_pair(i-1,j-1)) != _hitMap.end() ){
      for (auto &h : _hitMap[std::make_pair(i-1,j-1)])
	hitIndices.push_back(h);
    }
    // _________
    // |__|XX|__|
    // |__|__|__|
    // |__|__|__|
    if ( _hitMap.find(std::make_pair(i,j+1)) != _hitMap.end() ){
      for (auto &h : _hitMap[std::make_pair(i,j+1)])
	hitIndices.push_back(h);
    }
    // _________
    // |__|__|__|
    // |__|__|XX|
    // |__|__|__|
    if ( _hitMap.find(std::make_pair(i+1,j)) != _hitMap.end() ){
      for (auto &h : _hitMap[std::make_pair(i+1,j)])
	hitIndices.push_back(h);
    }
    // _________
    // |__|__|XX|
    // |__|__|__|
    // |__|__|__|
    if ( _hitMap.find(std::make_pair(i+1,j+1)) != _hitMap.end() ){
      for (auto &h : _hitMap[std::make_pair(i+1,j+1)])
	hitIndices.push_back(h);
    }
    // _________
    // |XX|__|__|
    // |__|__|__|
    // |__|__|__|
    if ( _hitMap.find(std::make_pair(i-1,j+1)) != _hitMap.end() ){
      for (auto &h : _hitMap[std::make_pair(i-1,j+1)])
	hitIndices.push_back(h);
    }
    // _________
    // |__|__|__|
    // |__|__|__|
    // |__|__|XX|
    if ( _hitMap.find(std::make_pair(i+1,j-1)) != _hitMap.end() ){
      for (auto &h : _hitMap[std::make_pair(i+1,j-1)])
	hitIndices.push_back(h);
    }
  }

  // if two hits are further apart then the set distance -> not compatible
  bool GammaPhotons::HitsCompatible(const hit& h1, const hit& h2){

    if (h1.WireID().Plane != h2.WireID().Plane)
      return false;

    double dt = ( h1.PeakTime() - h2.PeakTime() ) * _time2cm;
    //  if the hit time-ranges overlap, this distnce should be 0
    if (TimeOverlap(h1,h2,dt) == true)
      dt = 0;
    double dw = fabs( ((double)h1.Channel()-(double)h2.Channel())*_wire2cm );
    if (dw >  0.3) dw -= 0.3;
    //if (dw < -0.3) dw = 0;// += 0.3;
    double d = dt*dt + dw*dw;

    if (d > (_radius*_radius))
      return false;

    return true;
  }


  bool GammaPhotons::TimeOverlap(const larlite::hit& h1,
				    const larlite::hit& h2,
				    double& dmin) const
    {
      bool overlap = false;

      auto T1 = h1.PeakTime() * _time2cm; // time of first hit
      auto T2 = h2.PeakTime() * _time2cm;
      auto W1 = h1.RMS() * _time2cm;
      auto W2 = h2.RMS() * _time2cm;

      double d = dmin;
      
      if (T1 > T2) {

	if ( (T2+W2) > (T1-W1) ) return true;
	
	d = (T1-W1) - (T2+W2);
	if (d < dmin) dmin = d;
	
      }
      
      else {
	
	if ( (T1+W1) > (T2-W2) ) return true;
	
	d = (T2-W2) - (T1+W1);
	if (d < dmin) dmin = d;

      }
	
      return false;
    }

  void GammaPhotons::MakeHitMap(const event_hit* hitlist, int plane, const std::vector<std::vector<size_t> >& shr_hit_idx_v){

    _hitMap.clear();
    // temporary pair
    std::pair<int,int> tmpPair;

    for (auto const& shr_hit_v : shr_hit_idx_v) {

      for (auto const& hit_idx : shr_hit_v) {
	
	auto const& hit = hitlist->at(hit_idx);
	// skip if not of plane we want
	if (hit.View() != plane)
	  continue;

	// if RMS above threshold -> ignore
	if (hit.RMS() > _max_rms)
	  continue;
	
	// remove hits with time-tick < _tick_min or > _tick_max
	if ( (hit.PeakTime() < _tick_min) or (hit.PeakTime() > _tick_max) )
	  continue;
	
	double t = hit.PeakTime()*_time2cm;
	double w = hit.WireID().Wire*_wire2cm;
	
	// map is (i,j) -> hit list
	// i : ith bin in wire of some width
	// j : jth bin in time of some width
	int i = int(w/_cellSize);
	int j = int(t/_cellSize);
	tmpPair = std::make_pair(i,j);
	// does this entry exist in the map?
	// if yes -> append to vector
	// if no create new vector and add to map
	if (_hitMap.find(tmpPair) == _hitMap.end()){
	  std::vector<size_t> aaa = {hit_idx};
	  _hitMap[tmpPair] = aaa;
	}
	else
	  _hitMap[tmpPair].push_back(hit_idx);
      }// for all hits
    }// for all showers
    
    return;
  }
  
}
#endif
