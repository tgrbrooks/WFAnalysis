#ifndef SIMPLEWFANA_CXX
#define SIMPLEWFANA_CXX

#include "SimpleWFAna.h"

namespace larlite {

// ---- LOCAL FUNCTIONS ---- //

// Function to draw digital waveforms
void drawWF( std::vector<short int> const & ad, int size, double offset ) {
    // Variables to draw on graph
    size = 300;
    double x[size], y[size];
    // Fill variables with digital wavform information
    // x in units of ticks = 0.5us, y in units of ADCs
    // To just draw part change size and l's on RHS of = to start point
    for(int l=0;l<size;l++) {
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
int hitCount(std::vector<short int> const & adcs, double T, double offset) {
    int flag = 0;
    size_t l = 0;
    int hit = 0;
    // Determine if wire is hit or not
    while(flag == 0 && l<adcs.size()){
        double x = adcs[l];
        if((x-offset)>T||(x-offset)<-T){
            hit ++;
            flag = 1;
        }
        l++;
    }
    return hit;
}

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
	double y = abs(adcs[l]-offset);
        if((y)>T){
            // If threshold is passed set max as this point
            max = y;
            // Loop over max size of peak and find the most extremal point
            for(size_t i = 1;i<150;i++){
                double yp = abs(adcs[l+i]-offset);
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

// Function to return standard error
double StandErr(std::vector<double> const & cutVar, double avrg) {
    // Calculate standard error on the mean 
    double stErr = 0;
    for(size_t k=0; k<cutVar.size(); ++k)
    stErr += (cutVar[k]-avrg)*(cutVar[k]-avrg);
    stErr = sqrt( stErr / (cutVar.size()-1) );
    stErr = stErr/sqrt(cutVar.size());
    return stErr;
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

// Print to text file
void PrintText(std::ofstream & fileptr, std::vector<int> const & removedno, std::vector<int> const & nutypeno,int size) {
    fileptr<<removedno[0]/size<<"% Total, "<<removedno[1]<<"/"<<nutypeno[0]<<" CCQE, "
    <<removedno[2]<<"/"<<nutypeno[4]<<" NCQE, "<<removedno[3]<<"/"<<nutypeno[1]<<" CCRE, "
    <<removedno[4]<<"/"<<nutypeno[5]<<" NCRE, "<<removedno[5]<<"/"<<nutypeno[2]<<" CCDIS, "
    <<removedno[6]<<"/"<<nutypeno[6]<<" NCDIS, "<<removedno[7]<<"/"<<nutypeno[3]<<" CCCO, "
    <<removedno[8]<<"/"<<nutypeno[7]<<" NCCO"<<std::endl<<std::endl;
}

// ---- MEMBER FUNCTIONS ---- //

// Called at beginning of event loop
bool SimpleWFAna::initialize() {

    // Set event number to 0
    _evtN = 0;
    TH1::AddDirectory(kFALSE);

    // Initialize total time counters
    timeTav = 0; clockTav = 0;
    timeGetav = 0; clockGetav = 0;
    timeCutav = 0; clockCutav = 0;
    timeFillav = 0; clockFillav = 0;

    // ****** UNCOMMENT TO DRAW WAVEFORMS ****** //
    /*// Ask for event and wire number to display graph of digital waveform
    std::cout<<"Enter Event Number: ";
    std::cin>>event;std::cout<<std::endl;
    std::cout<<"Enter Wire Number: ";
    std::cin>>wire;std::cout<<std::endl;*/

    // Set event type counters to zero - only used for last file
    // NeutrinoTypeNo(8,0); // Neutrino Types CCQE, NCQE, CCRE, NCRE, CCDIS, NCDIS, CCCO and NCCO in that order.
    for(int i(0); i < 9; i++){
        Removedu.push_back(0); // Try using Removedu(9,0) instead?
        Removedv.push_back(0);
        Removedy.push_back(0);
        Removeduv.push_back(0);
        Removeduy.push_back(0);
        Removedvy.push_back(0);
        Removeduvy.push_back(0);
        RemovedType.push_back(0);
        if(i > 0) NeutrinoTypeNo.push_back(0);
    }

    // REMEMBER TO TEST AND CHANGE THESE LIMITS WHENEVER NEW DATA IS USED
    int tnbins=0; int unbins=0;int vnbins=0;int ynbins=0;int tmax=0;int umax=0;int vmax=0;int ymax=0;
    if(option==1){tnbins=70;unbins=40;vnbins=40;ynbins=60;tmax=7000;umax=2000;vmax=2000;ymax=3000;}
    if(option==2){tnbins=70;unbins=40;vnbins=40;ynbins=60;tmax=7000;umax=2000;vmax=2000;ymax=3000;}
    if(option==3){tnbins=60;unbins=60;vnbins=60;ynbins=60;tmax=120;umax=120;vmax=120;ymax=120;}
    if(option==4){tnbins=60;unbins=45;vnbins=45;ynbins=60;tmax=800000;umax=300000;vmax=300000;ymax=400000;}
    if(option==5){tnbins=70;unbins=40;vnbins=40;ynbins=60;tmax=7000;umax=2000;vmax=2000;ymax=3000;}

    // Book histograms
    h_HITS = new TH1I("h_HITS","",tnbins,0,tmax);
    h_UHITS = new TH1I("h_UHITS","",unbins,0,umax);
    h_VHITS = new TH1I("h_VHITS","",vnbins,0,vmax);
    h_YHITS = new TH1I("h_YHITS","",ynbins,0,ymax);
    h_UVHITS = new TH2I("h_UVHITS","",unbins,0,umax,vnbins,0,vmax);
    h_UYHITS = new TH2I("h_UYHITS","",unbins,0,umax,ynbins,0,ymax);
    h_VYHITS = new TH2I("h_VYHITS","",vnbins,0,vmax,ynbins,0,ymax);

    // Count event types
    if(!Type.empty()){
        for(size_t i=0;i<Type.size();i++){
            for(int j=0;j<8;j++){
                if(Type[i]==j) NeutrinoTypeNo[j] ++; // Loop over all neutrino types. See vector initialisation above.
            }
	}
    }
    return true;
  }
  
  // Called for every event
  bool SimpleWFAna::analyze(storage_manager* storage) {

    // Create initial total time and getting data time using both clock() and time()
    time_t timeT = time(NULL);
    time_t timeGet = time(NULL);
    clock_t clockT = clock();
    clock_t clockGet = clock();
  
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

    timeGet = time(NULL) - timeGet;
    clockGet = clock() - clockGet;
    // Create initial cut times
    time_t timeCut = time(NULL);
    clock_t clockCut = clock();

    // Initialize hit counter for event
    _isHit = 0;
    uHit = 0;
    vHit = 0;
    yHit = 0;

// ****** CUT ON INTEGRATED WAVEFORMS ****** //
    if(option==4){
        // Initialize integration counter for event
        intADC = 0;
        UintADC = 0;
        VintADC = 0;
        YintADC = 0;
    }

// ****** CUT ON TDC SPREAD ****** //
    if(option==2||option==5){
        // Clear TDC vectors
        TDCvec.clear();
        UTDCvec.clear();
        VTDCvec.clear();
        YTDCvec.clear();
    }

// ****** CUT ON ADC AMPLITUDE ******* //
    if(option==3){
        // Clear ADC vectors
        ADCvec.clear();
        UADCvec.clear();
        VADCvec.clear();
        YADCvec.clear();
    }

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
          _mean += adcs[j];
        _mean /= ((float)adcs.size());

        // Calculate offset from mean
        offset = static_cast<double>(_mean);
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

// ****** CUT ON INTEGRATED WAVEFORMS ****** //
    if(option==4){
        // Loop over TDC time
        for(size_t j=0;j<adcs.size()-1;++j){
            // Get offset adjusted modulus of two points
            double x1 = abs(adcs[j]-offset);
            double x2 = abs(adcs[j+1]-offset);
            // If both above tolerance add up area between them
            if(x1>T&&x2>T){
                intADC = intADC + 0.5*x1 + 0.5*x2;
                if(i<2400){UintADC = UintADC + 0.5*x1 + 0.5*x2;}
                if(i>=2400&&i<4800){VintADC = VintADC + 0.5*x1 + 0.5*x2;}
                if(i>=4800&&i<8256){YintADC = YintADC + 0.5*x1 + 0.5*x2;}
            }
        }
    }

// ****** CUT ON TDC SPREAD ****** //
    if(option==2||option==5){
        // Get TDCs of hits and sort into planes
        std::vector<int> v_TDCs = hitTDC(adcs,T,offset);
        for (size_t f=0; f<v_TDCs.size(); ++f) {
            TDCvec.push_back(v_TDCs[f]);
            if(i<2400){UTDCvec.push_back(v_TDCs[f]);}
            if(i>=2400&&i<4800){VTDCvec.push_back(v_TDCs[f]);}
            if(i>=4800&&i<8256){YTDCvec.push_back(v_TDCs[f]);}
        }
    }

// ****** CUT ON ADC AMPLITUDE ******* //
    if(option==3){
        // Get ADC amplitudes of hits and sort into planes
        std::vector<double> v_ADCamp = ADCamp(adcs,T,offset);
        for (size_t f=0; f<v_ADCamp.size(); ++f) {
            ADCvec.push_back(v_ADCamp[f]);
            if(i<2400){UADCvec.push_back(v_ADCamp[f]);}
            if(i>=2400&&i<4800){VADCvec.push_back(v_ADCamp[f]);}
            if(i>=4800&&i<8256){YADCvec.push_back(v_ADCamp[f]);}
        }
    }

// ****** UNCOMMENT TO DRAW WAVEFORMS ****** //
    //Draw digital waveform for channel i
    /*int w = i;
    if(w == wire && _evtN == event) {
        int size = 9600;
        drawWF(adcs,size,offset);
    }*/

  }

// ****** CUT ON TDC STANDARD DEVIATION ****** //
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
        if(isnan(YstdTDC)==1){YstdTDC=0;}
    }

// ****** CUT ON TDC INTERQUARTILE RANGE ****** //
    if(option==5){
        // If there are no hits set interquartile range to zero else calculte normally
        if(_isHit!=0){iqrTDC = TDCiqr(TDCvec,_isHit);}else{iqrTDC=0;}
        if(uHit!=0){UiqrTDC = TDCiqr(UTDCvec,uHit);}else{UiqrTDC=0;}
        if(vHit!=0){ViqrTDC = TDCiqr(VTDCvec,vHit);}else{ViqrTDC=0;}
        if(yHit!=0){YiqrTDC = TDCiqr(YTDCvec,yHit);}else{YiqrTDC=0;}
    }

// ****** CUT ON ADC AMPLITUDE ******* //
    if(option==3){
        // Calculate mean of amplitudes
        if(_isHit!=0){MampADC = ampMean(ADCvec);}else{MampADC=0;}
        if(uHit!=0){UMampADC = ampMean(UADCvec);}else{UMampADC=0;}
        if(vHit!=0){VMampADC = ampMean(VADCvec);}else{VMampADC=0;}
        if(yHit!=0){YMampADC = ampMean(YADCvec);}else{YMampADC=0;}
    }

    timeCut = time(NULL) - timeCut;
    clockCut = clock() - clockCut;
    // Create initial histogram filling end event removal time
    time_t timeFill = time(NULL);
    clock_t clockFill = clock();

    int en = _evtN;

    // Fill arrays with whatever variable you want to cut - nothing beyond here should have to be changed
    // Work with _isHit/uHit, intADC, stdTDC, iqrTDC, MampADC
    if(option==1){
        hitNo.push_back(_isHit);
        uhitNo.push_back(uHit);
        vhitNo.push_back(vHit);
        yhitNo.push_back(yHit);
    }
    if(option==2){
        hitNo.push_back(stdTDC);
        uhitNo.push_back(UstdTDC);
        vhitNo.push_back(VstdTDC);
        yhitNo.push_back(YstdTDC);
    }
    if(option==3){
        hitNo.push_back(MampADC);
        uhitNo.push_back(UMampADC);
        vhitNo.push_back(VMampADC);
        yhitNo.push_back(YMampADC);
    }
    if(option==4){
        hitNo.push_back(intADC);
        uhitNo.push_back(UintADC);
        vhitNo.push_back(VintADC);
        yhitNo.push_back(YintADC);
    }
    if(option==5){
        hitNo.push_back(iqrTDC);
        uhitNo.push_back(UiqrTDC);
        vhitNo.push_back(ViqrTDC);
        yhitNo.push_back(YiqrTDC);
    }

    // Fill histograms
    // Different cuts will need different ranges
    // Change variables as before
    // Separate hits into the three wire planes
    if(option==1){
        h_HITS->Fill(_isHit);
        h_UHITS->Fill(uHit);
        h_VHITS->Fill(vHit);
        h_YHITS->Fill(yHit);
        h_UVHITS->Fill(uHit,vHit);
        h_UYHITS->Fill(uHit,yHit);
        h_VYHITS->Fill(vHit,yHit);
    }
    if(option==2){
        h_HITS->Fill(stdTDC);
        h_UHITS->Fill(UstdTDC);
        h_VHITS->Fill(VstdTDC);
        h_YHITS->Fill(YstdTDC);
        h_UVHITS->Fill(UstdTDC,VstdTDC);
        h_UYHITS->Fill(UstdTDC,YstdTDC);
        h_VYHITS->Fill(VstdTDC,YstdTDC);
    }
    if(option==3){
        h_HITS->Fill(MampADC);
        h_UHITS->Fill(UMampADC);
        h_VHITS->Fill(VMampADC);
        h_YHITS->Fill(YMampADC);
        h_UVHITS->Fill(UMampADC,VMampADC);
        h_UYHITS->Fill(UMampADC,YMampADC);
        h_VYHITS->Fill(VMampADC,YMampADC);
    }
    if(option==4){
        h_HITS->Fill(intADC);
        h_UHITS->Fill(UintADC);
        h_VHITS->Fill(VintADC);
        h_YHITS->Fill(YintADC);
        h_UVHITS->Fill(UintADC,VintADC);
        h_UYHITS->Fill(UintADC,YintADC);
        h_VYHITS->Fill(VintADC,YintADC);
    }
    if(option==5){
        h_HITS->Fill(iqrTDC);
        h_UHITS->Fill(UiqrTDC);
        h_VHITS->Fill(ViqrTDC);
        h_YHITS->Fill(YiqrTDC);
        h_UVHITS->Fill(UiqrTDC,ViqrTDC);
        h_UYHITS->Fill(UiqrTDC,YiqrTDC);
        h_VYHITS->Fill(ViqrTDC,YiqrTDC);
    }

    if(!Type.empty()){
      for(int j=0; j < 8; j++){
        if(Type[en]==j){TypeEnergyBef.insert(std::pair<int,double>(Type[en],Qsq[en]));}
      }
    }

    // Count number of removed events outside of some tolerance of hit wires
    // Total cut
    if(hitNo[en]>Tmax||hitNo[en]<Tmin){
      RemovedType[0] ++;
      if(!Type.empty()){
        for(int j=0; j < 8; j ++){
          if(Type[en]==j){ RemovedType[j+1] ++;}
        }
      }
    }

    // Cut on U plane
    if(uhitNo[en]>Umax||uhitNo[en]<Umin){
      Removedu[0] ++;
      if(!Type.empty()){
        for(int j=0; j < 8; j++){
          if(Type[en]==j){Removedu[j+1] ++;}
        }
      }
    }else{
      if(!Type.empty()&&plane==1){
        for(int j=0; j < 8; j++){
          if(Type[en]==j){TypeEnergy.insert(std::pair<int,double>(Type[en],Qsq[en]));}
        }
      }}

    // Cut on V plane
    if(vhitNo[en]>Vmax||vhitNo[en]<Vmin){
      Removedv[0] ++;
      if(!Type.empty()){
        for(int j=0; j < 8; j++){
          if(Type[en]==j){Removedv[j+1] ++;}      
        }
      }
    }else{
      if(!Type.empty()&&plane==2){
        for(int j=0; j < 8; j++){
          if(Type[en]==j){TypeEnergy.insert(std::pair<int,double>(Type[en],Qsq[en]));}
        }
      }}

    // Cut on Y plane
    if(yhitNo[en]>Ymax||yhitNo[en]<Ymin){
      Removedy[0] ++;
      if(!Type.empty()){
        for(int j=0; j<8;j++){
          if(Type[en]==j){Removedy[j+1] ++;}
        }
      }
    }else{
      if(!Type.empty()&&plane==3){
        for(int j=0; j < 8; j++){
          if(Type[en]==j){TypeEnergy.insert(std::pair<int,double>(Type[en],Qsq[en]));}
        }
      }}

    // Cut on U and V planes also fill cutting histograms
    if((uhitNo[en]>Umax||uhitNo[en]<Umin)||(vhitNo[en]>Vmax||vhitNo[en]<Vmin)){
      Removeduv[0] ++;
      if(!Type.empty()){
        for(int j=0; j<8;j++){
          if(Type[en]==j){Removeduv[j+1] ++;}
        }
      }
    }else{
      if(!Type.empty()&&plane==4){
        for(int j=0; j < 8; j++){
          if(Type[en]==j){TypeEnergy.insert(std::pair<int,double>(Type[en],Qsq[en]));}
        }
      }}

    // Cut on U and Y planes
    if((uhitNo[en]>Umax||uhitNo[en]<Umin)||(yhitNo[en]>Ymax||yhitNo[en]<Ymin)){
      Removeduy[0] ++;
      if(!Type.empty()){
        for(int j=0; j<8;j++){
          if(Type[en]==j){Removeduy[j+1] ++;}
        }
      } 
    }else{
      if(!Type.empty()&&plane==5){
        for(int j=0; j < 8; j++){
          if(Type[en]==j){TypeEnergy.insert(std::pair<int,double>(Type[en],Qsq[en]));}
        }
      }}

    // Cut on V and Y planes
    if((yhitNo[en]>Ymax||yhitNo[en]<Ymin)||(vhitNo[en]>Vmax||vhitNo[en]<Vmin)){
      Removedvy[0] ++;
      if(!Type.empty()){
        for(int j=0; j<8;j++){
          if(Type[en]==j){Removedvy[j+1] ++;}
        }
      }   
    }else{
      if(!Type.empty()&&plane==6){
        for(int j=0; j < 8; j++){
          if(Type[en]==j){TypeEnergy.insert(std::pair<int,double>(Type[en],Qsq[en]));}
        }
      }}

    // Cut on U, V and Y planes
    if((uhitNo[en]>Umax||uhitNo[en]<Umin)||(yhitNo[en]>Ymax||yhitNo[en]<Ymin)||(vhitNo[en]>Vmax||vhitNo[en]<Vmin)){
      Removeduvy[0] ++;
      if(!Type.empty()){
        for(int j=0; j<8;j++){
          if(Type[en]==j){Removeduvy[j+1] ++;}
        }
      }   
    }else{
      if(!Type.empty()&&plane==7){
        for(int j=0; j < 8; j++){
          if(Type[en]==j){TypeEnergy.insert(std::pair<int,double>(Type[en],Qsq[en]));}
        }
      }}

    timeFill = time(NULL) - timeFill;
    timeT = time(NULL) - timeT;
    clockFill = clock() - clockFill;
    clockT = clock() - clockT;
    // Increment total time counters in seconds
    timeTav += (float)timeT;
    timeGetav += (float)timeGet;
    timeCutav += (float)timeCut;
    timeFillav += (float)timeFill;
    clockTav += ((float)clockT)/CLOCKS_PER_SEC;
    clockGetav += ((float)clockGet)/CLOCKS_PER_SEC;
    clockCutav += ((float)clockCut)/CLOCKS_PER_SEC;
    clockFillav += ((float)clockFill)/CLOCKS_PER_SEC;

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
    delete h_UYHITS;
    delete h_VYHITS;
    return true;}
    
    std::string StrArray[8] = {"h_QCCQE", "h_QNCQE", "h_QCCRE", "h_QNCRE", "h_QCCDIS", "h_QNCDIS", "h_QCCCO", "h_QNCCO"};
    std::string nCutStrArray[8] = {"h_NCUTQCCQE", "h_NCUTQNCQE", "h_NCUTQCCRE", "h_NCUTQNCRE", "h_NCUTQCCDIS", "h_NCUTQNCDIS", "h_NCUTQCCCO", "h_NCUTQNCCO"};

    // Calculate min and max limits for each plane - only important for nnbar file
    Tmin = hitNo[0];
    Tmax = hitNo[0];
    Umin = uhitNo[0];
    Umax = uhitNo[0];
    Vmin = vhitNo[0];
    Vmax = vhitNo[0];
    Ymin = yhitNo[0];
    Ymax = yhitNo[0];
    for(unsigned int k=0; k<hitNo.size(); ++k){
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
    h_UYHITS->Write();
    delete h_UYHITS;
    h_VYHITS->Write();
    delete h_VYHITS;
    // Write hitograms of neutrino type energy distributions before and after cut
    if(!TypeEnergyBef.empty()){
        for(int i=0;i<8;i++){
           h_TempHisto = new TH1D(StrArray[i].c_str(),"",50,0,10);
           // Loop over multimap of type and energy and find type, fill hist with energy and delete entry
           for(auto it=TypeEnergyBef.begin(); it!=TypeEnergyBef.end();){
             if((*it).first==i){
               h_TempHisto->Fill((*it).second);
               it = TypeEnergyBef.erase(it);
             } else {++it;}
           }
           h_nCutTempHisto = new TH1D(nCutStrArray[i].c_str(),"",50,0,10);
           // Do same for neutrinos that survived the cut
           for(auto it2=TypeEnergy.begin(); it2!=TypeEnergy.end();){
             if((*it2).first==i){
               h_nCutTempHisto->Fill((*it2).second);
               it2 = TypeEnergy.erase(it2);
             } else {++it2;}
           }
           // Clean up
           h_TempHisto->Write();
           h_nCutTempHisto->Write();
           delete h_TempHisto;
           delete h_nCutTempHisto;
        }
    }
    outfile->Close();

    // Calculate averages and standard deviations
    double avHN = ampMean(hitNo);
    double stErr = StandErr(hitNo,avHN);
    double avHNu = ampMean(uhitNo);
    double stErru = StandErr(uhitNo,avHNu);
    double avHNv = ampMean(vhitNo);
    double stErrv = StandErr(vhitNo,avHNv);
    double avHNy = ampMean(yhitNo);
    double stErry = StandErr(yhitNo,avHNy);

    // Write results to a .txt file
    std::ofstream myfile;
    myfile.open (tName);

    if(option==1){
        myfile<<"The average number of hits was: "<<avHN<<" +/- "<<stErr<<std::endl;
        myfile<<"The average number of U hits was: "<<avHNu<<" +/- "<<stErru<<std::endl;
        myfile<<"The average number of V hits was: "<<avHNv<<" +/- "<<stErrv<<std::endl;
        myfile<<"The average number of Y hits was: "<<avHNy<<" +/- "<<stErry<<std::endl;
    }

    if(option==2){
        myfile<<"The average TDC standard deviation was: "<<avHN<<" +/- "<<stErr<<std::endl;
        myfile<<"The average U TDC standard deviation was: "<<avHNu<<" +/- "<<stErru<<std::endl;
        myfile<<"The average V TDC standard deviation was: "<<avHNv<<" +/- "<<stErrv<<std::endl;
        myfile<<"The average Y TDC standard deviation was: "<<avHNy<<" +/- "<<stErry<<std::endl;
    }

    if(option==3){
        myfile<<"The average ADC amplitude was: "<<avHN<<" +/- "<<stErr<<std::endl;
        myfile<<"The average U ADC amplitude was: "<<avHNu<<" +/- "<<stErru<<std::endl;
        myfile<<"The average V ADC amplitude was: "<<avHNv<<" +/- "<<stErrv<<std::endl;
        myfile<<"The average Y ADC amplitude was: "<<avHNy<<" +/- "<<stErry<<std::endl;
    }

    if(option==4){
        myfile<<"The average integrated charge was: "<<avHN<<" +/- "<<stErr<<std::endl;
        myfile<<"The average U integrated charge was: "<<avHNu<<" +/- "<<stErru<<std::endl;
        myfile<<"The average V integrated charge was: "<<avHNv<<" +/- "<<stErrv<<std::endl;
        myfile<<"The average Y integrated charge was: "<<avHNy<<" +/- "<<stErry<<std::endl;
    }

    if(option==5){
        myfile<<"The average TDC interquartile range was: "<<avHN<<" +/- "<<stErr<<std::endl;
        myfile<<"The average U interquartile range was: "<<avHNu<<" +/- "<<stErru<<std::endl;
        myfile<<"The average V interquartile range was: "<<avHNv<<" +/- "<<stErrv<<std::endl;
        myfile<<"The average Y interquartile range was: "<<avHNy<<" +/- "<<stErry<<std::endl;
    }

    myfile<<std::endl<<"Cut Ranges:"<<std::endl<<"Total: Min = "<<Tmin<<"  Max = "<<Tmax<<std::endl<<"U Plane: Min = "<<Umin<<"  Max = "<<Umax<<std::endl;
    myfile<<"V Plane: Min = "<<Vmin<<"  Max = "<<Vmax<<std::endl<<"Y Plane: Min = "<<Ymin<<"  Max = "<<Ymax<<std::endl<<std::endl;

    int size = hitNo.size()/100;

    myfile<<"Total number of events removed: ";
    PrintText(myfile,RemovedType,NeutrinoTypeNo,size);

    myfile<<"Number of events removed using u info: ";
    PrintText(myfile,Removedu,NeutrinoTypeNo,size);

    myfile<<"Number of events removed using v info: ";
    PrintText(myfile,Removedv,NeutrinoTypeNo,size);

    myfile<<"Number of events removed using y info: ";
    PrintText(myfile,Removedy,NeutrinoTypeNo,size);

    myfile<<"Number of events removed using u and v info: ";
    PrintText(myfile,Removeduv,NeutrinoTypeNo,size);

    myfile<<"Number of events removed using v and y info: ";
    PrintText(myfile,Removedvy,NeutrinoTypeNo,size);

    myfile<<"Number of events removed using u and y info: ";
    PrintText(myfile,Removeduy,NeutrinoTypeNo,size);

    myfile<<"Number of events removed using u, v and y info: ";
    PrintText(myfile,Removeduvy,NeutrinoTypeNo,size);
    
    size = size*100;
    myfile<<"Timing Analysis:"<<std::endl<<"Total time per event= time( "<<timeTav/size<<"s ), clock( "<<clockTav/size<<"s )"<<std::endl;
    myfile<<"Time getting data = time( "<<timeGetav/size<<"s ), clock( "<<clockGetav/size<<"s )"<<std::endl;
    myfile<<"Time cutting = time( "<<timeCutav/size<<"s ), clock( "<<clockCutav/size<<"s )"<<std::endl;
    myfile<<"Time filling histograms and removing events = time( "<<timeTav/size<<"s ), clock( "<<clockTav/size<<"s )"<<std::endl;

    myfile.close();

    hitNo.clear(), uhitNo.clear(), vhitNo.clear(), yhitNo.clear();
    NeutrinoTypeNo.clear(); Removedu.clear(); Removedv.clear(); Removedy.clear(); 
    Removeduv.clear(); Removeduy.clear(); Removedvy.clear(); Removeduvy.clear(); RemovedType.clear(); 

    return true;
  }

}
#endif
