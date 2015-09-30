#ifndef SIMPLEWFANA_CXX
#define SIMPLEWFANA_CXX

#include "SimpleWFAna.h"

namespace larlite {

// ---- LOCAL FUNCTIONS ---- //

// Function to draw digital waveforms
void drawWF( const std::vector<short int> ad, int size, double offset ) {
    // Variables to draw on graph
    size = 300;
    double x[size], y[size];
    int n = size;
    // Fill variables with digital wavform information
    // x in units of ticks = 0.5us, y in units of ADCs
    // To just draw part change size and l's on RHS of = to start point
    for(int l=0;l<n;l++) {
        x[l] = l+3500;
        y[l] = ad[l+3500]-offset;
    }
    // Create new canvas and graph and draw on waveform
    TCanvas *c1 = new TCanvas("c1","Waveform Graph",1100,400);
    TGraph *gr = new TGraph(size,x,y);
    gr->SetLineWidth(1);
    gr->SetLineColorAlpha(kBlue,1);
    gr->SetMarkerStyle(1);
    gr->SetMarkerSize(1);
    gr->GetXaxis()->SetTitle("Ticks /0.5us");
    gr->GetYaxis()->SetTitle("ADCs");
    gr->Draw();
    c1->Update();
    c1->GetFrame()->SetFillColor(21);
    c1->GetFrame()->SetBorderSize(12);
    c1->Modified();
}

// Function to determine of wire is hit or not
int hitCount(const std::vector<short int> adcs, double T, double offset) {
    int flag = 0;
    size_t l = 0;
    int hit = 0;
    // Determine if wire is hit or not
    while(flag == 0 && l<adcs.size()){
        double x = adcs[l];
        if((x-offset)>T||(x-offset)<-T){
            hit = hit + 1;
            flag = 1;
        }
        l++;
    }
    return hit;
}

// Function that calculates number of hits on a wire assuming the beginning of one hit and end of another are separated by 150 TDCs
int hitPerWire(const std::vector<short int> adcs, double T, double offset) {
    size_t l = 0;
    int hit = 0;
    // Determine if wire is hit or not
    while(l<adcs.size()){
        double x = adcs[l];
        if((x-offset)>T||(x-offset)<-T){
            // add one to hit counter
            hit = hit + 1;
            // skip ahead by 150 TDCs
            l = l + 150;
        }
        l++;
    }
    return hit;
}

// Function that returns the TDC of the start of each hit
std::vector<int> hitTDC(const std::vector<short int> adcs, double T, double offset) {
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
std::vector<double> ADCamp(const std::vector<short int> adcs, double T, double offset) {
    size_t l = 0;
    std::vector<double> v_ADCamp;
    double max = 0;
    // Determine if wire is hit or not
    while(l<adcs.size()){
        double y =pow(pow(adcs[l]-offset,2),0.5);
        if((y)>T){
            // If threshold is passed set max as this point
            max = y;
            // Loop over max size of peak and find the most extremal point
            for(size_t i = 1;i<150;i++){
                double yp = pow(pow(adcs[l+i]-offset,2),0.5);
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

// Function to return average
double wireAve(std::vector<double> hitNo) {
    // Calculate average - doesn't have to be hits
    double avHN = 0;
    for(int j=0; j<100; ++j)
    avHN += hitNo[j];
    avHN /= (100);
    return avHN;
}

// Function to return standard error
double wireStd(std::vector<double> hitNo, double avHN) {
    // Calculate standard error on the mean in number of hits
    double stDev = 0;
    for(int k=0; k<100; ++k)
      stDev += (hitNo[k]-avHN)*(hitNo[k]-avHN);
    stDev = sqrt( stDev /(100*99));
    return stDev;
}

// Function to return standard deviation of TDC
float TDCstd(std::vector<int> TDCvec) {
    // Calculate mean
    float meanTDC(0), stdTDC(0);
    for (size_t f=0; f<TDCvec.size(); ++f) {
        meanTDC += TDCvec[f];
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
    int TDCq1 = TDCvec[quart];
    int TDCq3 = TDCvec[3*quart];
    int iqrTDC = TDCq3 - TDCq1;
    return iqrTDC;
}
// calculate mean ADC amplitude
double ampMean(std::vector<double> ADCvec) {
    // calculate mean
    double meanADC(0);
    for(size_t f=0; f<ADCvec.size(); ++f) {
        meanADC += ADCvec[f];
    }
    meanADC /= ((double)ADCvec.size());
    return meanADC;
}
// calculate standard deviation of amplitude
double ampStd(std::vector<double> ADCvec,double meanADC) {
    // calculate standard deviation
    double stdADC(0);
    for(size_t k=0; k<ADCvec.size(); ++k)
      stdADC += (ADCvec[k]-meanADC)*(ADCvec[k]-meanADC);
    stdADC = sqrt( stdADC / ((double)ADCvec.size()));
    return stdADC;
}

// ---- MEMBER FUNCTIONS ---- //

// Called at beginning of event loop
bool SimpleWFAna::initialize() {

    // Set event number to 0
    _evtN = 0;

// ****** UNCOMMENT TO DRAW WAVEFORMS ****** //
    /*// Ask for event and wire number to display graph of digital waveform
    std::cout<<"Enter Event Number: ";
    std::cin>>event;std::cout<<std::endl;
    std::cout<<"Enter Wire Number: ";
    std::cin>>wire;std::cout<<std::endl;*/

    // Set counter of removed events to zero
    removed = 0;removedCCQE=0;removedNCQE=0;removedCCRE=0;removedNCRE=0;removedCCDIS=0;removedNCDIS=0;removedCCCO=0;removedNCCO=0;
    removedu = 0;removedCCQEu=0;removedNCQEu=0;removedCCREu=0;removedNCREu=0;removedCCDISu=0;removedNCDISu=0;removedCCCOu=0;removedNCCOu=0;
    removedv = 0;removedCCQEv=0;removedNCQEv=0;removedCCREv=0;removedNCREv=0;removedCCDISv=0;removedNCDISv=0;removedCCCOv=0;removedNCCOv=0;
    removedy = 0;removedCCQEy=0;removedNCQEy=0;removedCCREy=0;removedNCREy=0;removedCCDISy=0;removedNCDISy=0;removedCCCOy=0;removedNCCOy=0;
    removeduv = 0;removedCCQEuv=0;removedNCQEuv=0;removedCCREuv=0;removedNCREuv=0;removedCCDISuv=0;removedNCDISuv=0;removedCCCOuv=0;removedNCCOuv=0;
    removedvy = 0;removedCCQEvy=0;removedNCQEvy=0;removedCCREvy=0;removedNCREvy=0;removedCCDISvy=0;removedNCDISvy=0;removedCCCOvy=0;removedNCCOvy=0;
    removeduy = 0;removedCCQEuy=0;removedNCQEuy=0;removedCCREuy=0;removedNCREuy=0;removedCCDISuy=0;removedNCDISuy=0;removedCCCOuy=0;removedNCCOuy=0;
    // Set event type counters to zero - only used for last file
    CCQEno = 0;
    NCQEno = 0;
    CCREno = 0;
    NCREno = 0;
    CCDISno = 0;
    NCDISno = 0;
    CCCOno = 0;
    NCCOno = 0;
    // Book histograms
    h_HITS = new TH1I("h_HITS","",60,0,800000);
    h_UHITS = new TH1I("h_UHITS","",45,0,300000);
    h_VHITS = new TH1I("h_VHITS","",45,0,300000);
    h_YHITS = new TH1I("h_YHITS","",60,0,400000);
    h_UVHITS = new TH2I("h_UVHITS","",45,0,300000,60,0,400000);
    h_HITvAMP = new TH1I("h_HITvAMP","",80,0,200);
    h_QCUT = new TH1D("h_QCUT","",50,0,10);
    h_QNCUT = new TH1D("h_QNCUT","",50,0,10);

    // Count event types
    if(!Type.empty()){
        for(int i=0;i<100;i++){
            if(Type[i]==0){CCQEno=CCQEno+1;}
            if(Type[i]==1){NCQEno=NCQEno+1;}
            if(Type[i]==2){CCREno=CCREno+1;}
            if(Type[i]==3){NCREno=NCREno+1;}
            if(Type[i]==4){CCDISno=CCDISno+1;}
            if(Type[i]==5){NCDISno=NCDISno+1;}
            if(Type[i]==6){CCCOno=CCCOno+1;}
            if(Type[i]==7){NCCOno=NCCOno+1;}
	}
    }
    return true;
  }
  
  // Called for every event
  bool SimpleWFAna::analyze(storage_manager* storage) {
    
    // Set flag for mctruth info to 0
    int truthflag = 0;
    // Get mctruth data
    auto ev_mct = storage->get_data<event_mctruth>("generator");
    // Set truth data flag to 1 if data is present
    if ( (ev_mct) ){
    	truthflag = 1;
    }
    // If truth data is present record interaction type and energy of each event
    if(truthflag==1){
        auto const& mct = (*ev_mct).at(0);
        // Get neutrino info
        auto const& neut = (mct.GetNeutrino());
        // Get charged current or neutral current
        int ccnc = neut.CCNC();
        // Get interaction type
        int intMode = neut.Mode();
        // Store neutrino energy
        Qsq.push_back(neut.Nu().Trajectory()[0].E());
        if(intMode==0&&ccnc==0){Type.push_back(0);}//CCQE
        else if(intMode==0&&ccnc==1){Type.push_back(1);}//NCQE
        else if(intMode==1&&ccnc==0){Type.push_back(2);}//CCRE
        else if(intMode==1&&ccnc==1){Type.push_back(3);}//NCRE
        else if(intMode==2&&ccnc==0){Type.push_back(4);}//CCDIS
        else if(intMode==2&&ccnc==1){Type.push_back(5);}//NCDIS
        else if(intMode==3&&ccnc==0){Type.push_back(6);}//CCCO
        else if(intMode==3&&ccnc==1){Type.push_back(7);}//NCCO
        else{std::cout<<"ERROR IN MCTRUTH INFO"<<std::endl;return false;}
        _evtN += 1;
        return true;
    }

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

// ****** UNCOMMENT TO CUT ON INTEGRATED WAVEFORMS ****** //
    if(option==4){
    // Initialize integration counter for event
    intADC = 0;
    UintADC = 0;
    VintADC = 0;
    YintADC = 0;}

// ****** UNCOMMENT TO CUT ON TDC SPREAD ****** //
    if(option==2||option==5){
    // Clear TDC vectors
    TDCvec.clear();
    UTDCvec.clear();
    VTDCvec.clear();
    YTDCvec.clear();}

// ****** UNCOMMENT TO CUT ON ADC AMPLITUDE ******* //
    if(option==3){
    // Clear ADC vectors
    ADCvec.clear();
    UADCvec.clear();
    VADCvec.clear();
    YADCvec.clear();}

    // Loop over all wires in event
    for (size_t i=0; i < wfs->size(); i++){
        // get waveform from wire i
        auto const& wf = (*wfs).at(i);
        // Convert from analogue to digital
        auto const& adcs = wf.ADCs();

        // Calculate Mean
        for(size_t j=0; j<adcs.size(); ++j)
          _mean += adcs[j];
        _mean /= ((float)adcs.size());

        // Calculate offset from mean
        double offset = static_cast<double>(_mean);
        // Get tolerance from python script
        double T = SimpleWFAna::GetT();

        // Calculate total number of hits
        // Add to hit counter for event
        // Change to hitCount() for number of hit wires
        _isHit = _isHit + hitPerWire(adcs,T,offset);
        // Separate into wire planes
        if(i<2400){uHit = uHit + hitPerWire(adcs,T,offset);}
        if(i>=2400&&i<4800){vHit = vHit + hitPerWire(adcs,T,offset);}
        if(i>=4800&&i<8256){yHit = yHit + hitPerWire(adcs,T,offset);}

// ****** UNCOMMENT TO CUT ON INTEGRATED WAVEFORMS ****** //
    if(option==4){
    // Loop over TDC time
    for(size_t j=0;j<adcs.size()-1;++j){
        // Get offset adjusted modulus of two points
        double x1 = pow(pow(adcs[j]-offset,2),0.5);
        double x2 = pow(pow(adcs[j+1]-offset,2),0.5);
        // If both above tolerance add up area between them
        if(x1>T&&x2>T){
            intADC = intADC + 0.5*x1 + 0.5*x2;
            if(i<2400){UintADC = UintADC + 0.5*x1 + 0.5*x2;}
            if(i>=2400&&i<4800){VintADC = VintADC + 0.5*x1 + 0.5*x2;}
            if(i>=4800&&i<8256){YintADC = YintADC + 0.5*x1 + 0.5*x2;}
        }
    }}

// ****** UNCOMMENT TO CUT ON TDC SPREAD ****** //
    if(option==2||option==5){
    // Get TDCs of hits and sort into planes
    std::vector<int> v_TDCs = hitTDC(adcs,T,offset);
    for (size_t f=0; f<v_TDCs.size(); ++f) {
        TDCvec.push_back(v_TDCs[f]);
        if(i<2400){UTDCvec.push_back(v_TDCs[f]);}
        if(i>=2400&&i<4800){VTDCvec.push_back(v_TDCs[f]);}
        if(i>=4800&&i<8256){YTDCvec.push_back(v_TDCs[f]);}
    }}

// ****** UNCOMMENT TO CUT ON ADC AMPLITUDE ******* //
    if(option==3){
    // Get ADC amplitudes of hits and sort into planes
    std::vector<double> v_ADCamp = ADCamp(adcs,T,offset);
    for (size_t f=0; f<v_ADCamp.size(); ++f) {
        ADCvec.push_back(v_ADCamp[f]);
        if(i<2400){UADCvec.push_back(v_ADCamp[f]);}
        if(i>=2400&&i<4800){VADCvec.push_back(v_ADCamp[f]);}
        if(i>=4800&&i<8256){YADCvec.push_back(v_ADCamp[f]);}
    }}

// ****** UNCOMMENT TO DRAW WAVEFORMS ****** //
    //Draw digital waveform for channel i
    /*int w = i;
    if(w == wire && _evtN == event) {
        int size = 9600;
        drawWF(adcs,size,offset);
    }*/

  }

// ****** UNCOMMENT TO CUT ON TDC STANDARD DEVIATION ****** //
    if(option==2){
    // Calculate standard deviations of TDCs of hits in each plane
    stdTDC = TDCstd(TDCvec);
    UstdTDC = TDCstd(UTDCvec);
    VstdTDC = TDCstd(VTDCvec);
    YstdTDC = TDCstd(YTDCvec);
    // If not a number set to zero
    if(isnan(stdTDC)==1){stdTDC=0;}
    if(isnan(UstdTDC)==1){UstdTDC=0;}
    if(isnan(VstdTDC)==1){VstdTDC=0;}
    if(isnan(YstdTDC)==1){YstdTDC=0;}}

// ****** UNCOMMENT TO CUT ON TDC INTERQUARTILE RANGE ****** //
    if(option==5){
    // If there are no hits set interquartile range to zero else calculte normally
    if(_isHit!=0){iqrTDC = TDCiqr(TDCvec,_isHit);}else{iqrTDC=0;}
    if(uHit!=0){UiqrTDC = TDCiqr(UTDCvec,uHit);}else{UiqrTDC=0;}
    if(vHit!=0){ViqrTDC = TDCiqr(VTDCvec,vHit);}else{ViqrTDC=0;}
    if(yHit!=0){YiqrTDC = TDCiqr(YTDCvec,yHit);}else{YiqrTDC=0;}}

// ****** UNCOMMENT TO CUT ON ADC AMPLITUDE ******* //
    if(option==3){
    // Calculate mean and standard deviation of amplitudes
    if(_isHit!=0){MampADC = ampMean(ADCvec);SDampADC = ampStd(ADCvec,MampADC);}else{MampADC=0;SDampADC=0;}
    if(uHit!=0){UMampADC = ampMean(UADCvec);USDampADC = ampStd(UADCvec,UMampADC);}else{UMampADC=0;USDampADC=0;}
    if(vHit!=0){VMampADC = ampMean(VADCvec);VSDampADC = ampStd(VADCvec,VMampADC);}else{VMampADC=0;VSDampADC=0;}
    if(yHit!=0){YMampADC = ampMean(YADCvec);YSDampADC = ampStd(YADCvec,YMampADC);}else{YMampADC=0;YSDampADC=0;}}

    // Fill arrays of doubles with event number
    int en = _evtN;
    eventNo.push_back(_evtN);

    // Fill arrays with whatever variable you want to cut - nothing beyond here should have to be changed
    // Work with _isHit/uHit, intADC, stdTDC, iqrTDC, MampADC
    if(option==1){
    hitNo.push_back(_isHit);
    uhitNo.push_back(uHit);
    vhitNo.push_back(vHit);
    yhitNo.push_back(yHit);}
    if(option==2){
    hitNo.push_back(stdTDC);
    uhitNo.push_back(UstdTDC);
    vhitNo.push_back(VstdTDC);
    yhitNo.push_back(YstdTDC);}
    if(option==4){
    hitNo.push_back(MampADC);
    uhitNo.push_back(UMampADC);
    vhitNo.push_back(VMampADC);
    yhitNo.push_back(YMampADC);}
    if(option==4){
    hitNo.push_back(intADC);
    uhitNo.push_back(UintADC);
    vhitNo.push_back(VintADC);
    yhitNo.push_back(YintADC);}
    if(option==5){
    hitNo.push_back(iqrTDC);
    uhitNo.push_back(UiqrTDC);
    vhitNo.push_back(ViqrTDC);
    yhitNo.push_back(YiqrTDC);}

// ****** UNCOMMENT TO CALCULTE ADC AMPLITUDE STANDARD DEVIATIONS ****** //
    /**sDev.push_back(SDampADC);
    usDev.push_back(USDampADC);
    vsDev.push_back(VSDampADC);
    ysDev.push_back(YSDampADC);*/

    // Fill histograms
    // Different cuts will need different ranges
    // Change variables as before
    // Separate hits into the three wire planes
    if(option==1){
    h_HITS->Fill(_isHit);
    h_UHITS->Fill(uHit);
    h_VHITS->Fill(vHit);
    h_YHITS->Fill(yHit);
    h_UVHITS->Fill(uHit,vHit);}
    if(option==2){
    h_HITS->Fill(stdTDC);
    h_UHITS->Fill(UstdTDC);
    h_VHITS->Fill(VstdTDC);
    h_YHITS->Fill(YstdTDC);
    h_UVHITS->Fill(VstdTDC,YstdTDC);}
    if(option==3){
    h_HITS->Fill(MampADC);
    h_UHITS->Fill(UMampADC);
    h_VHITS->Fill(VMampADC);
    h_YHITS->Fill(YMampADC);
    h_UVHITS->Fill(VMampADC,YMampADC);}
    if(option==4){
    h_HITS->Fill(intADC);
    h_UHITS->Fill(UintADC);
    h_VHITS->Fill(VintADC);
    h_YHITS->Fill(YintADC);
    h_UVHITS->Fill(VintADC,YintADC);}
    if(option==5){
    h_HITS->Fill(iqrTDC);
    h_UHITS->Fill(UiqrTDC);
    h_VHITS->Fill(ViqrTDC);
    h_YHITS->Fill(YiqrTDC);
    h_UVHITS->Fill(ViqrTDC,YiqrTDC);}

    // Count number of removed events outside of some tolerance of hit wires
    // Total cut
    if(hitNo[en]>Tmax||hitNo[en]<Tmin){removed = removed + 1;
      if(!Type.empty()){if(Type[en]==0){removedCCQE=removedCCQE+1;}
      if(Type[en]==1){removedNCQE=removedNCQE+1;}
      if(Type[en]==2){removedCCRE=removedCCRE+1;}
      if(Type[en]==3){removedNCRE=removedNCRE+1;}
      if(Type[en]==4){removedCCDIS=removedCCDIS+1;}
      if(Type[en]==5){removedNCDIS=removedNCDIS+1;}
      if(Type[en]==6){removedCCCO=removedCCCO+1;}
      if(Type[en]==7){removedNCCO=removedNCCO+1;}}
    }
    // Cut on U plane
    if(uhitNo[en]>Umax||uhitNo[en]<Umin){removedu = removedu + 1;
      if(!Type.empty()){if(Type[en]==0){removedCCQEu=removedCCQEu+1;}
      if(Type[en]==1){removedNCQEu=removedNCQEu+1;}
      if(Type[en]==2){removedCCREu=removedCCREu+1;}
      if(Type[en]==3){removedNCREu=removedNCREu+1;}
      if(Type[en]==4){removedCCDISu=removedCCDISu+1;}
      if(Type[en]==5){removedNCDISu=removedNCDISu+1;}
      if(Type[en]==6){removedCCCOu=removedCCCOu+1;}
      if(Type[en]==7){removedNCCOu=removedNCCOu+1;}}
    }
    // Cut on V plane
    if(vhitNo[en]>Vmax||vhitNo[en]<Vmin){removedv = removedv + 1;
      if(!Type.empty()){if(Type[en]==0){removedCCQEv=removedCCQEv+1;}
      if(Type[en]==1){removedNCQEv=removedNCQEv+1;}
      if(Type[en]==2){removedCCREv=removedCCREv+1;}
      if(Type[en]==3){removedNCREv=removedNCREv+1;}
      if(Type[en]==4){removedCCDISv=removedCCDISv+1;}
      if(Type[en]==5){removedNCDISv=removedNCDISv+1;}
      if(Type[en]==6){removedCCCOv=removedCCCOv+1;}
      if(Type[en]==7){removedNCCOv=removedNCCOv+1;}}
    }
    // Cut on Y plane
    if(yhitNo[en]>Ymax||yhitNo[en]<Ymin){removedy = removedy + 1;
      if(!Type.empty()){if(Type[en]==0){removedCCQEy=removedCCQEy+1;}
      if(Type[en]==1){removedNCQEy=removedNCQEy+1;}
      if(Type[en]==2){removedCCREy=removedCCREy+1;}
      if(Type[en]==3){removedNCREy=removedNCREy+1;}
      if(Type[en]==4){removedCCDISy=removedCCDISy+1;}
      if(Type[en]==5){removedNCDISy=removedNCDISy+1;}
      if(Type[en]==6){removedCCCOy=removedCCCOy+1;}
      if(Type[en]==7){removedNCCOy=removedNCCOy+1;}}
    }
    // Cut on U and V planes also fill cutting histograms
    if((uhitNo[en]>Umax||uhitNo[en]<Umin)||(vhitNo[en]>Vmax||vhitNo[en]<Vmin)){removeduv = removeduv + 1;
      if(!Type.empty()){h_QCUT->Fill(Qsq[en]); 
          if(Type[en]==0){removedCCQEuv=removedCCQEuv+1;}
      if(Type[en]==1){removedNCQEuv=removedNCQEuv+1;}
      if(Type[en]==2){removedCCREuv=removedCCREuv+1;}
      if(Type[en]==3){removedNCREuv=removedNCREuv+1;}
      if(Type[en]==4){removedCCDISuv=removedCCDISuv+1;}
      if(Type[en]==5){removedNCDISuv=removedNCDISuv+1;}
      if(Type[en]==6){removedCCCOuv=removedCCCOuv+1;}
      if(Type[en]==7){removedNCCOuv=removedNCCOuv+1;}}
    }else {if(!Type.empty()){h_QNCUT->Fill(Qsq[en]);}}
    // Cut on U and Y planes
    if((uhitNo[en]>Umax||uhitNo[en]<Umin)||(yhitNo[en]>Ymax||yhitNo[en]<Ymin)){removeduy = removeduy + 1;
      if(!Type.empty()){if(Type[en]==0){removedCCQEuy=removedCCQEuy+1;}
      if(Type[en]==1){removedNCQEuy=removedNCQEuy+1;}
      if(Type[en]==2){removedCCREuy=removedCCREuy+1;}
      if(Type[en]==3){removedNCREuy=removedNCREuy+1;}
      if(Type[en]==4){removedCCDISuy=removedCCDISuy+1;}
      if(Type[en]==5){removedNCDISuy=removedNCDISuy+1;}
      if(Type[en]==6){removedCCCOuy=removedCCCOuy+1;}
      if(Type[en]==7){removedNCCOuy=removedNCCOuy+1;}}
    }
    // Cut on V and Y planes
    if((yhitNo[en]>Ymax||yhitNo[en]<Ymin)||(vhitNo[en]>Vmax||vhitNo[en]<Vmin)){removedvy = removedvy + 1;
      if(!Type.empty()){if(Type[en]==0){removedCCQEvy=removedCCQEvy+1;}
      if(Type[en]==1){removedNCQEvy=removedNCQEvy+1;}
      if(Type[en]==2){removedCCREvy=removedCCREvy+1;}
      if(Type[en]==3){removedNCREvy=removedNCREvy+1;}
      if(Type[en]==4){removedCCDISvy=removedCCDISvy+1;}
      if(Type[en]==5){removedNCDISvy=removedNCDISvy+1;}
      if(Type[en]==6){removedCCCOvy=removedCCCOvy+1;}
      if(Type[en]==7){removedNCCOvy=removedNCCOvy+1;}}
    }

    _evtN += 1;
    return true;
  }

  // Called at end of event loop
  bool SimpleWFAna::finalize() {

    // Check for truth run
    if ( hitNo.empty() ){
    delete h_HITS;
    delete h_UHITS;
    delete h_VHITS;
    delete h_YHITS;
    delete h_UVHITS;
    delete h_QCUT;
    delete h_QNCUT;
    delete h_HITvAMP;
    return true;}

    // Calculate min and max limits for each plane - only important for nnbar file
    Tmin = hitNo[0];
    Tmax = hitNo[0];
    Umin = uhitNo[0];
    Umax = uhitNo[0];
    Vmin = vhitNo[0];
    Vmax = vhitNo[0];
    Ymin = yhitNo[0];
    Ymax = yhitNo[0];
    for(int k=0; k<100; ++k){
      if(hitNo[k]<Tmin){Tmin = hitNo[k];}
      if(hitNo[k]>Tmax){Tmax = hitNo[k];}
      if(uhitNo[k]<Umin){Umin = uhitNo[k];}
      if(uhitNo[k]>Umax){Umax = uhitNo[k];}
      if(vhitNo[k]<Vmin){Vmin = vhitNo[k];}
      if(vhitNo[k]>Vmax){Vmax = vhitNo[k];}
      if(yhitNo[k]<Ymin){Ymin = yhitNo[k];}
      if(yhitNo[k]>Ymax){Ymax = yhitNo[k];}
    }

    // Get file names from python script
    std::string fileName = SimpleWFAna::GetName();
    const char * fName = fileName.c_str();
    std::string txtName = SimpleWFAna::GetfName();
    const char * tName = txtName.c_str();

    // Save histograms to .root file
    TFile* outfile = new TFile(fName, "RECREATE");
    h_HITS->Write();
    delete h_HITS;
    h_UHITS->Write();
    delete h_UHITS;
    h_VHITS->Write();
    delete h_VHITS;
    h_YHITS->Write();
    delete h_YHITS;
    h_UVHITS->Write();
    delete h_UVHITS;
    h_QCUT->Write();
    delete h_QCUT;
    h_QNCUT->Write();
    delete h_QNCUT;
    delete h_HITvAMP;
    outfile->Close();

    // Calculate averages and standard deviations
    double avHN = wireAve(hitNo);
    double stDev = wireStd(hitNo, avHN);
    double avHNu = wireAve(uhitNo);
    double stDevu = wireStd(uhitNo, avHNu);
    double avHNv = wireAve(vhitNo);
    double stDevv = wireStd(vhitNo, avHNv);
    double avHNy = wireAve(yhitNo);
    double stDevy = wireStd(yhitNo, avHNy);

// ****** UNCOMMENT FOR AVERAGE ADC AMPLITUDE STANDARD DEVIATIONS ****** //
    /*// CALCULATE AVERAGE STANDARD DEVIATIONS
    double avSDev = wireAve(sDev);
    double avSDevu = wireAve(usDev);
    double avSDevv = wireAve(vsDev);
    double avSDevy = wireAve(ysDev);*/

    // Write results to a .txt file
    std::ofstream myfile;
    myfile.open (tName);
    myfile<<"The average integrated charge was: "<<avHN<<" +/- "<<stDev<<std::endl;
  //  myfile<<"The average ADC standard deviation was: "<<avSDev<<std::endl;
    myfile<<"The average mean U integrated charge was: "<<avHNu<<" +/- "<<stDevu<<std::endl;
  //  myfile<<"The average U ADC standard deviation was: "<<avSDevu<<std::endl;
    myfile<<"The average mean V integrated charge was: "<<avHNv<<" +/- "<<stDevv<<std::endl;
  //  myfile<<"The average V ADC standard deviation was: "<<avSDevv<<std::endl;
    myfile<<"The average mean Y integrated charge was: "<<avHNy<<" +/- "<<stDevy<<std::endl;
  //  myfile<<"The average Y ADC standard deviation was: "<<avSDevy<<std::endl;
    myfile<<"Total number of events removed: "<<removed<<"% Total, "<<removedCCQE<<"/"<<CCQEno<<" CCQE, "<<removedNCQE<<"/"<<NCQEno<<" NCQE, "<<removedCCRE<<"/"<<CCREno<<" CCRE, "<<removedNCRE<<"/"<<NCREno<<" NCRE, "<<removedCCDIS<<"/"<<CCDISno<<" CCDIS, "<<removedNCDIS<<"/"<<NCDISno<<" NCDIS, "<<removedCCCO<<"/"<<CCCOno<<" CCCO, "<<removedNCCO<<"/"<<NCCOno<<" NCCO"<<std::endl<<std::endl;
    myfile<<"Number of events removed using u info: "<<removedu<<"% Total, "<<removedCCQEu<<"/"<<CCQEno<<" CCQE, "<<removedNCQEu<<"/"<<NCQEno<<" NCQE, "<<removedCCREu<<"/"<<CCREno<<" CCRE, "<<removedNCREu<<"/"<<NCREno<<" NCRE, "<<removedCCDISu<<"/"<<CCDISno<<" CCDIS, "<<removedNCDISu<<"/"<<NCDISno<<" NCDIS, "<<removedCCCOu<<"/"<<CCCOno<<" CCCO, "<<removedNCCOu<<"/"<<NCCOno<<" NCCO"<<std::endl<<std::endl;
    myfile<<"Number of events removed using v info: "<<removedv<<"% Total, "<<removedCCQEv<<"/"<<CCQEno<<" CCQE, "<<removedNCQEv<<"/"<<NCQEno<<" NCQE, "<<removedCCREv<<"/"<<CCREno<<" CCRE, "<<removedNCREv<<"/"<<NCREno<<" NCRE, "<<removedCCDISv<<"/"<<CCDISno<<" CCDIS, "<<removedNCDISv<<"/"<<NCDISno<<" NCDIS, "<<removedCCCOv<<"/"<<CCCOno<<" CCCO, "<<removedNCCOv<<"/"<<NCCOno<<" NCCO"<<std::endl<<std::endl;
    myfile<<"Number of events removed using y info: "<<removedy<<"% Total, "<<removedCCQEy<<"/"<<CCQEno<<" CCQE, "<<removedNCQEy<<"/"<<NCQEno<<" NCQE, "<<removedCCREy<<"/"<<CCREno<<" CCRE, "<<removedNCREy<<"/"<<NCREno<<" NCRE, "<<removedCCDISy<<"/"<<CCDISno<<" CCDIS, "<<removedNCDISy<<"/"<<NCDISno<<" NCDIS, "<<removedCCCOy<<"/"<<CCCOno<<" CCCO, "<<removedNCCOy<<"/"<<NCCOno<<" NCCO"<<std::endl<<std::endl;
    myfile<<"Number of events removed using u and v info: "<<removeduv<<"% Total, "<<removedCCQEuv<<"/"<<CCQEno<<" CCQE, "<<removedNCQEuv<<"/"<<NCQEno<<" NCQE, "<<removedCCREuv<<"/"<<CCREno<<" CCRE, "<<removedNCREuv<<"/"<<NCREno<<" NCRE, "<<removedCCDISuv<<"/"<<CCDISno<<" CCDIS, "<<removedNCDISuv<<"/"<<NCDISno<<" NCDIS, "<<removedCCCOuv<<"/"<<CCCOno<<" CCCO, "<<removedNCCOuv<<"/"<<NCCOno<<" NCCO"<<std::endl<<std::endl;
    myfile<<"Number of events removed using v and y info: "<<removedvy<<"% Total, "<<removedCCQEvy<<"/"<<CCQEno<<" CCQE, "<<removedNCQEvy<<"/"<<NCQEno<<" NCQE, "<<removedCCREvy<<"/"<<CCREno<<" CCRE, "<<removedNCREvy<<"/"<<NCREno<<" NCRE, "<<removedCCDISvy<<"/"<<CCDISno<<" CCDIS, "<<removedNCDISvy<<"/"<<NCDISno<<" NCDIS, "<<removedCCCOvy<<"/"<<CCCOno<<" CCCO, "<<removedNCCOvy<<"/"<<NCCOno<<" NCCO"<<std::endl<<std::endl;
    myfile<<"Number of events removed using y and u info: "<<removeduy<<"% Total, "<<removedCCQEuy<<"/"<<CCQEno<<" CCQE, "<<removedNCQEuy<<"/"<<NCQEno<<" NCQE, "<<removedCCREuy<<"/"<<CCREno<<" CCRE, "<<removedNCREuy<<"/"<<NCREno<<" NCRE, "<<removedCCDISuy<<"/"<<CCDISno<<" CCDIS, "<<removedNCDISuy<<"/"<<NCDISno<<" NCDIS, "<<removedCCCOuy<<"/"<<CCCOno<<" CCCO, "<<removedNCCOuy<<"/"<<NCCOno<<" NCCO"<<std::endl<<std::endl;
    myfile.close();

    eventNo.clear();
    hitNo.clear(), uhitNo.clear(), vhitNo.clear(), yhitNo.clear();
    //sDev.clear(), usDev.clear(), vsDev.clear(), ysDev.clear();

    return true;
  }

}
#endif
