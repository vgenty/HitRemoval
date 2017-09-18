#ifndef LARLITE_SPALLATIONCOVERAGE_CXX
#define LARLITE_SPALLATIONCOVERAGE_CXX

#include "SpallationCoverage.h"
#include "LArUtil/GeometryHelper.h"
#include "LArUtil/Geometry.h"
#include "DataFormat/cluster.h"
#include "DataFormat/vertex.h"

namespace larlite {

  SpallationCoverage::SpallationCoverage()
    : _tree(nullptr)
  {

    _name        = "SpallationCoverage";
    _fout        = 0;
    _verbose     = false;
    _hitProducer = "gaushit";
    _vtxProducer = "";
    _out_clusterProducer = "rawclus";
    _useVtx      = false;
    _radius      = 2.0;
    _cellSize    = 2;
    _vtx_radius  = 0;
    _max_rms     = 100;
    _vtx_w_cm = {0,0,0};
    _vtx_t_cm = {0,0,0};
    _tick_min = 0;
    _tick_max = 9600;
    _roi_radius = 1000;

    if (_tree) {delete _tree;}
    _tree = new TTree("tree","tree");
    _tree->Branch("_qcell",&_qcell,"qcell/D");
    _tree->Branch("_w",&_w,"w/D");
    _tree->Branch("_t",&_t,"t/D");

  }

  bool SpallationCoverage::initialize() {

    _wire2cm  = larutil::GeometryHelper::GetME()->WireToCm();
    _time2cm  = larutil::GeometryHelper::GetME()->TimeToCm();

    std::cout << "********************************" << std::endl;
    std::cout << "Wire -> cm conversion : " << _wire2cm << std::endl;
    std::cout << "Time -> cm conversion : " << _time2cm << std::endl;
    std::cout << "********************************" << std::endl;

    return true;
  }
  
  bool SpallationCoverage::analyze(storage_manager* storage) {

    auto evt_hits      = storage->get_data<event_hit>(_hitProducer);

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

    // hit map will only contain hits we want to use for clustering
    MakeHitMap(evt_hits,2);

    // loop throught TPC points and estimate the coverage
    for(int w=10; w < 93; w++) {
      for (int t=10; t < 25; t++) {
	
	_w = w * 10;
	_t = t * 10;

	int i = int(_w/_cellSize);
	int j = int(_t/_cellSize);
	
	_qcell = GetCellCharge(evt_hits, std::make_pair(i,j) );
	
	_tree->Fill();
	
      }// for all times
    }// for all wires
    
    return true;
  }

  bool SpallationCoverage::finalize() {

    if (_fout) 
      if (_tree) _tree->Write();

    return true;
  }

    /// get charge in cells
  double SpallationCoverage::GetCellCharge(const event_hit* ev_hit,
					   const std::pair<int,int>& cell) {

    double qtot = 0;

    auto const& hit_idx_v = _hitMap[cell];

    for (auto const& idx : hit_idx_v) {
      
      auto const& hit = ev_hit->at(idx);
      
      // is goodness of fit negative? if so ignore the hit
      if ( (hit.GoodnessOfFit()) < 0 )
	continue;
      
      qtot += ev_hit->at(idx).Integral();
      
    }// for hits in cell
    
    return qtot;
  }

  void SpallationCoverage::MakeHitMap(const event_hit* hitlist, int plane){
    
    _hitMap.clear();
    // temporary pair
    std::pair<int,int> tmpPair;

    
    for (size_t h=0; h < hitlist->size(); h++){
      
      auto const& hit = hitlist->at(h);
      // skip if not of plane we want
      if (hit.View() != plane)
	continue;

      // if RMS above threshold -> ignore
      if (hit.RMS() > _max_rms)
	continue;

      // is goodness of fit negative? if so ignore the hit
      if ( (hit.GoodnessOfFit()) < 0 )
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
	std::vector<size_t> aaa = {h};
	_hitMap[tmpPair] = aaa;
      }
      else
	_hitMap[tmpPair].push_back(h);
    }// for all hits

    return;
  }

}
#endif
