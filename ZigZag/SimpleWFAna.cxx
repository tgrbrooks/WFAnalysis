#ifndef SIMPLEWFANA_CXX
#define SIMPLEWFANA_CXX

#include "SimpleWFAna.h"

namespace larlite {

// ---- LOCAL FUNCTIONS ---- //

// Function that calculates number of hits on a wire assuming the beginning of one hit and end of another are separated by 150 TDCs
int hitPerWire(std::vector<short int> const & adcs, double T, double offset) {
    size_t l = 0;
    int hit = 0;
    // Determine if wire is hit or not
    while(l<adcs.size()){
        double x = adcs[l];
        if((x-offset)>T||(x-offset)<-T){
            // add one to hit counter
            hit ++;
            // skip ahead by 150 TDCs
            l = l + 150;
        }
        l++;
    }
    return hit;
}

// Function that returns the TDC of the start of each hit
std::vector<int> hitTDC(std::vector<short int> const & adcs, double T, double offset) {
    // Initialize TDC variable and vector to store
    size_t l = 0;
    std::vector<int> v_TDCs;
    // Determine if wire is hit or not
    while(l<adcs.size()){
        double x = adcs[l];
        if((x-offset)>T||(x-offset)<-T){
            // if wire is hit record TDC time of when it occurs
            int TDC = l;
            v_TDCs.push_back(TDC);
            l = l + 150;
        }
        l++;
    }
    return v_TDCs;
}

// Function that returns maximum height of each hit in ADC units
std::vector<double> ADCamp(std::vector<short int> const & adcs, double T, double offset) {
    size_t l = 0;
    std::vector<double> v_ADCamp;
    double max = 0;
    // Determine if wire is hit or not
    while(l<adcs.size()){
	double y = std::abs(adcs[l]-offset);
        if((y)>T){
            // If threshold is passed set max as this point
            max = y;
            // Loop over max size of peak and find the most extremal point
            for(size_t i = 1;i<150;i++){
                double yp = std::abs(adcs[l+i]-offset);
                if(yp>max){
	            max = yp;
                }
            }
            // Push on to vector of amplitudes
            v_ADCamp.push_back(max);
            l = l + 150;
        }
        l++;
    }
    return v_ADCamp;
}

// Function to return standard deviation of TDC
float TDCstd(std::vector<int> const & TDCvec) {
    // Calculate mean
    float stdTDC(0);
    float meanTDC(0);

    for(size_t f=0; f<TDCvec.size(); ++f) {
        meanTDC += ((float)TDCvec[f]);
    }
    meanTDC /= ((float)TDCvec.size());

    // Calculate standard deviation
    for(size_t k=0; k<TDCvec.size(); ++k)
      stdTDC += (TDCvec[k]-meanTDC)*(TDCvec[k]-meanTDC);
    stdTDC = sqrt( stdTDC / ((float)TDCvec.size()));
    return stdTDC;
}

// Function to calculate interquartile range of TDC
int TDCiqr(std::vector<int> TDCvec, int hitNo) {
    // order TDC's by value
    std::sort(TDCvec.begin(),TDCvec.end());
    int quart =(int) hitNo/4;
    int iqrTDC = (TDCvec[3*quart] - TDCvec[quart]);
    return iqrTDC;
}

// calculate mean ADC amplitude
double ampMean(std::vector<double> const & ADCvec) {
    // calculate mean
    double meanADC(0);
    for(size_t f=0; f<ADCvec.size(); ++f) {
        meanADC += ADCvec[f];
    }
    meanADC /= ((double)ADCvec.size());
    return meanADC;
}

// ---- MEMBER FUNCTIONS ---- //

// Called at beginning of event loop
bool SimpleWFAna::initialize() {

    // Set event number to 0
    _evtN = 0;

    _t_ch = new TTree("ch_tree","");
    _t_ch->Branch("evt",&_evtN,"evt/i");
    _t_ch->Branch("hitNo",&_hitNo,"hitNo/i");
    _t_ch->Branch("hitNoU",&_hitNoU,"hitNoU/i");
    _t_ch->Branch("hitNoV",&_hitNoV,"hitNoV/i");
    _t_ch->Branch("hitNoY",&_hitNoY,"hitNoY/i");
    _t_ch->Branch("TDCstd",&_TDCstd,"TDCstd/F");
    _t_ch->Branch("TDCstdU",&_TDCstdU,"TDCstdU/F");
    _t_ch->Branch("TDCstdV",&_TDCstdV,"TDCstdV/F");
    _t_ch->Branch("TDCstdY",&_TDCstdY,"TDCstdY/F");
    _t_ch->Branch("TDCiqr",&_TDCiqr,"TDCiqr/F");
    _t_ch->Branch("TDCiqrU",&_TDCiqrU,"TDCiqrU/F");
    _t_ch->Branch("TDCiqrV",&_TDCiqrV,"TDCiqrV/F");
    _t_ch->Branch("TDCiqrY",&_TDCiqrY,"TDCiqrY/F");
    _t_ch->Branch("ADCamp",&_ADCamp,"ADCamp/F");
    _t_ch->Branch("ADCampU",&_ADCampU,"ADCampU/F");
    _t_ch->Branch("ADCampV",&_ADCampV,"ADCampV/F");
    _t_ch->Branch("ADCampY",&_ADCampY,"ADCampY/F");
    _t_ch->Branch("WFint",&_WFint,"WFint/F");
    _t_ch->Branch("WFintU",&_WFintU,"WFintU/F");
    _t_ch->Branch("WFintV",&_WFintV,"WFintV/F");
    _t_ch->Branch("WFintY",&_WFintY,"WFintY/F");
    _t_ch->SetDirectory(0);

    return true;
  }
  
  // Called for every event
  bool SimpleWFAna::analyze(storage_manager* storage) {

    // Get rawdigit data
    auto wfs = storage->get_data<event_rawdigit>("daq");
    // Display error if rawdigit data not present
    if ( (!wfs) || (!wfs->size())){
    	print (msg::kERROR,__FUNCTION__,"RawDigit data product not found!");
    	return false;
    }

    // Initialize hit counter for event
    _isHit = 0;
    uHit = 0;
    vHit = 0;
    yHit = 0;

    // Initialize integration counter for event
    intADC = 0;
    UintADC = 0;
    VintADC = 0;
    YintADC = 0;

    // Clear TDC vectors
    TDCvec.clear();
    UTDCvec.clear();
    VTDCvec.clear();
    YTDCvec.clear();

    // Clear ADC vectors
    ADCvec.clear();
    UADCvec.clear();
    VADCvec.clear();
    YADCvec.clear();

    double offset = 0;
    double T = 0;

    // Loop over all wires in event
    for (size_t i=0; i < wfs->size(); i++){
        // get waveform from wire i
        auto const& wf = (*wfs).at(i);
        // Convert from analogue to digital
        auto const& adcs = wf.ADCs();
        // Calculate Mean
        for(size_t j=0; j<adcs.size(); ++j)
          offset += adcs[j];
        offset /= ((double)adcs.size());

        // Calculate offset from mean
        //offset = static_cast<double>(_mean);
        // Get tolerance from python script
        T = SimpleWFAna::GetT();

        // Calculate total number of hits
        // Add to hit counter for event
        // Change to hitCount() for number of hit wires
        _isHit = _isHit + hitPerWire(adcs,T,offset);
        // Separate into wire planes
        if(i<2400){uHit = uHit + hitPerWire(adcs,T,offset);}
        if(i>=2400&&i<4800){vHit = vHit + hitPerWire(adcs,T,offset);}
        if(i>=4800&&i<8256){yHit = yHit + hitPerWire(adcs,T,offset);}

        // Loop over TDC time
        for(size_t j=0;j<adcs.size()-1;++j){
            // Get offset adjusted modulus of two points
            double x1 = std::abs(adcs[j]-offset);
            double x2 = std::abs(adcs[j+1]-offset);
            // If both above tolerance add up area between them
            if(x1>T&&x2>T){
                intADC = intADC + 0.5*x1 + 0.5*x2;
                if(i<2400){UintADC = UintADC + 0.5*x1 + 0.5*x2;}
                if(i>=2400&&i<4800){VintADC = VintADC + 0.5*x1 + 0.5*x2;}
                if(i>=4800&&i<8256){YintADC = YintADC + 0.5*x1 + 0.5*x2;}
            }
        }


        // Get TDCs of hits and sort into planes
        std::vector<int> v_TDCs = hitTDC(adcs,T,offset);
        for (size_t f=0; f<v_TDCs.size(); ++f) {
            TDCvec.push_back(v_TDCs[f]);
            if(i<2400){UTDCvec.push_back(v_TDCs[f]);}
            if(i>=2400&&i<4800){VTDCvec.push_back(v_TDCs[f]);}
            if(i>=4800&&i<8256){YTDCvec.push_back(v_TDCs[f]);}
        }

        // Get ADC amplitudes of hits and sort into planes
        std::vector<double> v_ADCamp = ADCamp(adcs,T,offset);
        for (size_t f=0; f<v_ADCamp.size(); ++f) {
            ADCvec.push_back(v_ADCamp[f]);
            if(i<2400){UADCvec.push_back(v_ADCamp[f]);}
            if(i>=2400&&i<4800){VADCvec.push_back(v_ADCamp[f]);}
            if(i>=4800&&i<8256){YADCvec.push_back(v_ADCamp[f]);}
        }

  }

    // Calculate standard deviations of TDCs of hits in each plane
    stdTDC = TDCstd(TDCvec);
    UstdTDC = TDCstd(UTDCvec);
    VstdTDC = TDCstd(VTDCvec);
    YstdTDC = TDCstd(YTDCvec);
    // If not a number set to zero
    if(isnan(stdTDC)==1){stdTDC=0;}
    if(isnan(UstdTDC)==1){UstdTDC=0;}
    if(isnan(VstdTDC)==1){VstdTDC=0;}
    if(isnan(YstdTDC)==1){YstdTDC=0;}

    // If there are no hits set interquartile range to zero else calculte normally
    if(_isHit!=0){iqrTDC = TDCiqr(TDCvec,_isHit);}else{iqrTDC=0;}
    if(uHit!=0){UiqrTDC = TDCiqr(UTDCvec,uHit);}else{UiqrTDC=0;}
    if(vHit!=0){ViqrTDC = TDCiqr(VTDCvec,vHit);}else{ViqrTDC=0;}
    if(yHit!=0){YiqrTDC = TDCiqr(YTDCvec,yHit);}else{YiqrTDC=0;}

    // Calculate mean of amplitudes
    if(_isHit!=0){MampADC = ampMean(ADCvec);}else{MampADC=0;}
    if(uHit!=0){UMampADC = ampMean(UADCvec);}else{UMampADC=0;}
    if(vHit!=0){VMampADC = ampMean(VADCvec);}else{VMampADC=0;}
    if(yHit!=0){YMampADC = ampMean(YADCvec);}else{YMampADC=0;}

    // Fill arrays with whatever variable you want to cut - nothing beyond here should have to be changed
    // Work with _isHit/uHit, intADC, stdTDC, iqrTDC, MampADC
    _hitNo = _isHit;
    _hitNoU = uHit;
    _hitNoV = vHit;
    _hitNoY = yHit;

    _TDCstd = stdTDC;
    _TDCstdU = UstdTDC;
    _TDCstdV = VstdTDC;
    _TDCstdY = YstdTDC;

    _ADCamp = MampADC;
    _ADCampU = UMampADC;
    _ADCampV = VMampADC;
    _ADCampY = YMampADC;

    _WFint = intADC;
    _WFintU = UintADC;
    _WFintV = VintADC;
    _WFintY = YintADC;

    _TDCiqr = iqrTDC;
    _TDCiqrU = UiqrTDC;
    _TDCiqrV = ViqrTDC;
    _TDCiqrY = YiqrTDC;

    _t_ch->Fill();
std::cout<<MampADC<<std::endl;
    _evtN += 1;
    return true;
  }

  // Called at end of event loop
  bool SimpleWFAna::finalize() {

    // Save histograms to .root file
    if(_fout){
      _fout->cd();
      std::cout << "writing ch tree" << std::endl;
      _t_ch->Write();
    }
    //outfile->Close();

    return true;
  }

}
#endif
