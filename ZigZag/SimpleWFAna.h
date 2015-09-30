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

namespace larlite {
  /**
     \class SimpleWFAna
     User custom analysis class made by kterao
   */
  class SimpleWFAna : public ana_base{

  private:
    // Declare histograms
    TH1I  *h_HITS; 
    TH1I  *h_UHITS;
    TH1I  *h_VHITS;
    TH1I  *h_YHITS;
    TH2I  *h_UVHITS;
    TH1I  *h_HITvAMP;
    TH1D  *h_QCUT;
    TH1D  *h_QNCUT;

  protected:   

    float _mean;

    int _evtN;

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
    // TDC ranges
    int rangeTDC, urangeTDC, vrangeTDC, yrangeTDC;

    // wire number
    int wire;
    // event number
    int event;
    // number of removed events
    int removed;
    int removedu,removedCCQEu,removedNCQEu,removedCCREu,removedNCREu,removedCCDISu,removedNCDISu,removedCCCOu,removedNCCOu;
    int removedv,removedCCQEv,removedNCQEv,removedCCREv,removedNCREv,removedCCDISv,removedNCDISv,removedCCCOv,removedNCCOv;
    int removedy,removedCCQEy,removedNCQEy,removedCCREy,removedNCREy,removedCCDISy,removedNCDISy,removedCCCOy,removedNCCOy;
    int removeduv,removedCCQEuv,removedNCQEuv,removedCCREuv,removedNCREuv,removedCCDISuv,removedNCDISuv,removedCCCOuv,removedNCCOuv;
    int removeduy,removedCCQEuy,removedNCQEuy,removedCCREuy,removedNCREuy,removedCCDISuy,removedNCDISuy,removedCCCOuy,removedNCCOuy;
    int removedvy,removedCCQEvy,removedNCQEvy,removedCCREvy,removedNCREvy,removedCCDISvy,removedNCDISvy,removedCCCOvy,removedNCCOvy;
    int removedCCQE,removedNCQE,removedCCRE,removedNCRE,removedCCDIS,removedNCDIS,removedCCCO,removedNCCO;
    int CCQEno,CCREno,CCDISno,CCCOno,NCQEno,NCREno,NCDISno,NCCOno;
    // array of events
    std::vector<int> eventNo;
    // array of number of hits
    std::vector<double> hitNo, uhitNo, vhitNo, yhitNo;
    // array of standard deviations
    std::vector<double> sDev, usDev, vsDev, ysDev;
    // vector of interaction types
    std::vector<int> Type;
    // vector of neutrino energies
    std::vector<double> Qsq;

    double x;
    // file names for use in simplewfana.py
    std::string nm;
    std::string fnm;
    // cut limits to be set
    double Tmin, Tmax, Umin, Umax, Vmin, Vmax, Ymin, Ymax; 

    // Vectors for TDC times of hits
    std::vector<int> TDCvec, UTDCvec, VTDCvec, YTDCvec;
    // Vectors for ADC amplitudes of hits
    std::vector<double> ADCvec, UADCvec, VADCvec, YADCvec;

    int truthflag;

    int option;

  public:

    /// Default constructor
    SimpleWFAna(double T,std::string name,std::string fname,int op ) { _name="SimpleWFAna"; _fout=0 ; x=T;nm=name;fnm=fname; option=op; };

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
      return x;
    }
    // Accessor functions for file names
    std::string GetName() {
      return nm;
    }
    std::string GetfName() {
      return fnm;
    }
    // Accessor and mutator functions for cut values
    double GetTmin() {
      return Tmin;
    }
    void SetTmin(double t_min) {
      Tmin = t_min;
    }
    double GetTmax() {
      return Tmax;
    }
    void SetTmax(double t_max) {
      Tmax = t_max;
    }
    double GetUmin() {
      return Umin;
    }
    void SetUmin(double u_min) {
      Umin = u_min;
    }
    double GetUmax() {
      return Umax;
    }
    void SetUmax(double u_max) {
      Umax = u_max;
    }
    double GetVmin() {
      return Vmin;
    }
    void SetVmin(double v_min) {
      Vmin = v_min;
    }
    double GetVmax() {
      return Vmax;
    }
    void SetVmax(double v_max) {
      Vmax = v_max;
    }
    double GetYmin() {
      return Ymin;
    }
    void SetYmin(double y_min) {
      Ymin = y_min;
    }
    double GetYmax() {
      return Ymax;
    }
    void SetYmax(double y_max) {
      Ymax = y_max;
    }
    // Accessor and mutator function for interaction type
    std::vector<int> GetType() {
      return Type;
    }
    void SetType(std::vector<int> type) {
      Type = type;
    }
    // Accessor and mutator functions for energy
    std::vector<double> GetQsq() {
      return Qsq;
    }
    void SetQsq(std::vector<double> qsq) {
      Qsq = qsq;
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
