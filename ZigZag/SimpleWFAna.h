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
    TH2I  *h_UYHITS;
    TH2I  *h_VYHITS;
    TH1I  *h_HITvAMP;
    TH1D  *h_QCUT;
    TH1D  *h_QNCUT;
    TH1D  *h_TempHisto;
    TH1D  *h_nCutTempHisto;

  protected: //Pretty sure this doesn't need to be protected

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

    // wire number
    int wire;
    // event number
    int event;

    // number of removed events
    std::vector <int> Removedu;
    std::vector <int> Removedv;
    std::vector <int> Removedy;
    std::vector <int> Removeduv;
    std::vector <int> Removeduy;
    std::vector <int> Removedvy; 
    std::vector <int> Removeduvy;
    std::vector <int> RemovedType;
    std::vector <int> NeutrinoTypeNo;     
    
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

    double Tolerance;
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
    int plane;

    std::multimap<int,double> TypeEnergy;
    std::multimap<int,double> TypeEnergyBef;

  public:

    /// Default constructor
    SimpleWFAna(double T,std::string name,std::string fname,int op,int pln) { _name="SimpleWFAna"; _fout=0 ; Tolerance=T;nm=name;fnm=fname; option=op;plane=pln; };

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
