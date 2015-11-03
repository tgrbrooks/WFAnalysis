/**
 * \file SimpleWFAna.h
 *
 * \ingroup SimpleWFAna
 * 
 * \brief Class def header for a class SimpleWFAna
 *
 * @author davidc1
 * adapted for triggger primitive studies by tom brooks
 */

/** \addtogroup SimpleWFAna

    @{*/

#ifndef SIMPLEWFANA_H
#define SIMPLEWFANA_H

#include "Analysis/ana_base.h"
#include <map>
#include <iostream>
#include <vector>
#include <algorithm>
#include <set>
#include <cmath>
#include "DataFormat/rawdigit.h"
#include "DataFormat/mctruth.h"
#include <TROOT.h>
#include <TH1.h>
#include <TH2.h>
#include <TCanvas.h>
#include <TGraph.h>
#include <TFrame.h>
#include <TAxis.h>
#include <TAttLine.h>
#include <TChain.h>
#include <TFile.h>
#include <string>
#include <fstream>
#include <map>
#include <time.h>

namespace larlite {
  /**
     \class SimpleWFAna
     User custom analysis class made by kterao
   */
  class SimpleWFAna : public ana_base{

  protected: //Pretty sure this doesn't need to be protected

    float _mean;

    int _evtN;

    int _hitNo, _hitNoU, _hitNoV, _hitNoY; 

    float _TDCstd, _TDCstdU, _TDCstdV, _TDCstdY;

    float _TDCiqr, _TDCiqrU, _TDCiqrV, _TDCiqrY;

    float _ADCamp, _ADCampU, _ADCampV, _ADCampY;

    float _WFint, _WFintU, _WFintV, _WFintY;

    // TDC standard deviations
    float stdTDC, UstdTDC, VstdTDC, YstdTDC;
    // Mean amplitudes
    double MampADC, UMampADC, VMampADC, YMampADC, SDampADC, USDampADC, VSDampADC, YSDampADC;
    // TDC interquartile ranges
    int iqrTDC, UiqrTDC, ViqrTDC, YiqrTDC;
    // Numbers of hits
    int _isHit, uHit ,vHit, yHit;
    // integrated ADC waveforms
    double intADC,UintADC,VintADC,YintADC;

    // wire number
    int wire;
    // event number
    int event;    

    double Tolerance;
    // file names for use in simplewfana.py
    std::string nm;
    std::string fnm;

    // Vectors for TDC times of hits
    std::vector<int> TDCvec, UTDCvec, VTDCvec, YTDCvec;
    // Vectors for ADC amplitudes of hits
    std::vector<double> ADCvec, UADCvec, VADCvec, YADCvec;

    TTree* _t_ch;

  public:

    /// Default constructor
    SimpleWFAna(double T,std::string name) { _name="SimpleWFAna"; _fout=0 ; Tolerance=T;nm=name;};

    /// Default destructor
    virtual ~SimpleWFAna(){};

    /** IMPLEMENT in SimpleWFAna.cxx
        Initialization method to be called before the analysis event loop.
    */ 
    virtual bool initialize();

    /** IMPLEMENT in SimpleWFAna.cxx 
        Analyze a data event-by-event  
    */
    virtual bool analyze(storage_manager* storage);

    /** IMPLEMENT in SimpleWFAna.cxx 
        Finalize method to be called after all events processed.
    */
    virtual bool finalize();

    // Accessor function for tolerance
    double GetT() {
      return Tolerance;
    }
    // Accessor functions for file names
    std::string GetName() {
      return nm;
    }

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
