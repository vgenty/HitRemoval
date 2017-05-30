/**
 * \file MakeROI.h
 *
 * \ingroup Base
 * 
 * \brief Class def header for a class MakeROI
 *
 * @author david
 */

/** \addtogroup Base

    @{*/

#ifndef LARLITE_MAKEROI_H
#define LARLITE_MAKEROI_H

#include "Analysis/ana_base.h"

namespace larlite {
  /**
     \class MakeROI
     User custom analysis class made by SHELL_USER_NAME
   */
  class MakeROI : public ana_base{
  
  public:

    /// Default constructor
    MakeROI(){ _name="MakeROI"; _fout=0;}

    /// Default destructor
    virtual ~MakeROI(){}

    virtual bool initialize();

    virtual bool analyze(storage_manager* storage);

    virtual bool finalize();

    void setVtxProducer(std::string s) { _vtx_producer = s; }

  protected:

    std::string _vtx_producer;
    
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
