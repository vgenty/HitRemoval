/**
 * \file RemoveMuonHits.h
 *
 * \ingroup CalcEfficiency
 * 
 * \brief Class def header for a class RemoveMuonHits
 *
 * @author ariana Hackenburg
 */

/** \addtogroup CalcEfficiency

    @{*/

#ifndef LARLITE_REMOVEMUONHITS_H
#define LARLITE_REMOVEMUONHITS_H

#include "Analysis/ana_base.h"

namespace larlite {
  /**
     \class RemoveMuonHits
     User custom analysis class made by SHELL_USER_NAME
   */
  class RemoveMuonHits : public ana_base{
  
  public:

    /// Default constructor
    RemoveMuonHits(){ _name="RemoveMuonHits"; _fout=0; }

    /// Default destructor
    virtual ~RemoveMuonHits(){}

    virtual bool initialize();

    virtual bool analyze(storage_manager* storage);

    virtual bool finalize();

    void setVtxProducer( std::string vtx ) { _vtx_producer = vtx ; }

    void setTrkProducer( std::string trk ) { _trk_producer = trk ; }

    void setHitAssProducer( std::string hitass ) { _hit_ass_producer = hitass; }

  protected:

  int _event ;
  int e ;

  std::string _vtx_producer;
  std::string _trk_producer;
  std::string _hit_ass_producer; 
    
  };
}
#endif

//**************************************************************************
// 
// For Analysis framework documentation, read Manual.pdf here:
//
// http://microboone-docdb.fnal.gov:8080/cgi-bin/ShowDocument?docid=3183
//
//**************************************************************************

/** @} */ // end of doxygen group 
