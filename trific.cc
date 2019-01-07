/**************************************************************
 * Work in progress Garfield++ simulation for TRIFIC.         *
 *                                                            *
 * Large chamber with vertically offset grids.  Meant to      *
 * perform more accurate simulation of window-to-first        *
 * cathode potential and electric field. Chamber directions   *
 * are not realistic!                                         *
 *                                                            *
 * A. Chester 16 Aug 2018                                     *
 *                                                            *
 **************************************************************/


#include <iostream>
#include <string>
#include <time.h>

#include <TApplication.h>
#include <TCanvas.h>
#include <TROOT.h>
#include <TMath.h>
#include <TRandom2.h>
#include <TFile.h>
#include <TH1D.h>
#include <TStyle.h>

#include "MediumMagboltz.hh"
#include "GeometrySimple.hh"
#include "SolidBox.hh"
#include "ComponentAnalyticField.hh"
#include "Sensor.hh"

#include "TrackHeed.hh"
#include "DriftLineRKF.hh"
#include "AvalancheMC.hh"

#include "ViewGeometry.hh"
#include "ViewCell.hh"
#include "ViewField.hh"
#include "ViewDrift.hh"

#include "Random.hh"

#include "FundamentalConstants.hh"

using namespace std;
using namespace Garfield;

// transfer function for convoluting electron signals
double transfer(double t) {
  const double a = 1.0;
  const double tau = 25.0; // ns
  return a*(t/tau)*TMath::Exp(1.0-(t/tau));
}

// get andode and cathode signals with time window and binning defined in sensor
TH1D GetGridElectronSignal(int nGr, Sensor* s, int nW, bool isCathode) {
  double tStart, binWidth;
  unsigned int nBins;
  s->GetTimeWindow(tStart, binWidth, nBins);
  double tStop = tStart + binWidth * double(nBins);

  TString hName,hTitle;
  if(isCathode) {
    hName  = TString::Format("hGrid%i_ElectronCathode", nGr);
    hTitle = TString::Format("Electron Signal Cathode Grid %i [ns];", nGr);
  }
  else {
    hName  = TString::Format("hGrid%i_ElectronAnode", nGr);
    hTitle = TString::Format("Electron Signal Anode Grid %i [ns];", nGr);
  }  
  
  TH1D hSignal(hName.Data(), hTitle.Data(), nBins, tStart, tStop);

  int wLow = nGr*nW;
  int wHigh = wLow + nW;
  TString wName;

  for(int i=wLow; i<wHigh; i++) {
    if(isCathode)
      wName = TString::Format("c%d",i);
    else 
      wName = TString::Format("a%d",i);
    // std::string whatIs;
    // if(isCathode) whatIs = "cathode";
    // else whatIs = "anode";
    // cout << "For " << whatIs << " grid " << nGr << " wLow = " << wLow << " wHigh = " << wHigh << " wire " << i << " named " << wName << endl;
    // getc(stdin);
    for(unsigned int j=1; j<=nBins; j++)
      hSignal.AddBinContent(j,s->GetElectronSignal(wName.Data(),j));
  }  
  return hSignal;
}

int main(int argc, char * argv[]) {

  string gasFileName = "cf4_NT_91.2Torr_nEF20_EFmin100.00_EFmax200000.gas";
  double gridAngle = 60.0;
  double anodeVoltage = 117.0; // V
  double cathodeVoltage = -117.0;
  string particleType = "alpha";
  double energy = 5.5e6; // 5.5 MeV in eV
  double particleTheta = 30.; // through center of grids 
  double driftPercentage = 0.05; // what % of e- to drift -- SLOW 

  if(argc != 9) {
    cout << "Usage: ./trific arg1 arg2 arg3 arg4 arg5 arg6 arg7" << endl;
    cout << "arg1: Magboltz gas file name" << endl;
    cout << "arg2: Grid angle w.r.t beam" << endl;
    cout << "arg3: Cathode [-] voltage in V" << endl;
    cout << "arg4: Anode [+] voltage in V" << endl;
    cout << "arg5: Particle type (alpha and 40Ar supported)" << endl;
    cout << "arg6: Particle kinetic energy in MeV" << endl;
    cout << "arg7: Particle theta w.r.t. beam in degrees" << endl;
    cout << "arg8: % of electrons to drift (0.0, ... , 100.0)" << endl;
    return(1);
  }

  for(int i=1; i<argc ; i++) {
    switch(i) {
    case 1: gasFileName = argv[1]; // Magboltz gas file
    case 2: gridAngle = atof(argv[2]); // Magboltz gas file
    case 3: cathodeVoltage = atof(argv[3]); // wire voltage settings in V
    case 4: anodeVoltage = atof(argv[4]);
    case 5: particleType = argv[5];
    case 6: energy = atof(argv[6])*1.0e6;   // eV
    case 7: particleTheta = atof(argv[7])*TMath::DegToRad(); // in radians
    case 8: driftPercentage = atof(argv[8])/100.0; // in radians
    }
  }

  if(gasFileName.length() == 0) {
    cout << "ERROR: String of length 0 passed as file name!" << endl;
    return(1);
  }
  if(gridAngle > 90 || gridAngle < 0) {
    cout << "ERROR: Grid angle must be in [0,90] degrees!" << endl;
    return(1);
  }
  if(anodeVoltage < 0) {
    cout << "ERROR: Anode voltage " << anodeVoltage << " is < 0!" << endl;
    return(1);
  }
  if(cathodeVoltage > 0) {
    cout << "ERROR: Cathode voltage " << cathodeVoltage << " is > 0!" << endl;
    return(1);
  }
  if(particleType != "alpha" && particleType != "Alpha" && 
     particleType != "ar40" && particleType != "40ar" && particleType != "40Ar" && particleType != "Ar40") {
    cout << "ERROR: Particle type " << particleType << " is not supported!" << endl;
    cout << "--> Please choose a supported particle types: 'alpha' or '40Ar'" << endl;
    cout << "\tor edit Track.cc to include your species of interest." << endl;
    return(1);
  }
  if(energy < 0) {
    cout << "ERROR: Particle energy " << energy << " is < 0!" << endl;
    return(1);
  }
  if(particleTheta > TMath::PiOver2() || particleTheta < -TMath::PiOver2()) {
    cout << "ERROR: Particle with angle " << particleTheta << " is moving upstream!" << endl;
    return(1);
  }
  if(driftPercentage < 0. || driftPercentage > 1.) {
    cout << "ERROR: Electron drift % " << driftPercentage*100. << " is not physical!" << endl;
    return(1);
  }

  // unsigned int theSeed = 1985810831;
  unsigned int theSeed = static_cast<unsigned int>(time(NULL));
  randomEngine.Seed(theSeed);

  MediumMagboltz* gas = new MediumMagboltz();
  if(gas->LoadGasFile(gasFileName) == false) {
    cout << "ERROR: Gas file not opened successfully!" << endl;
    cout << "Are you sure the file name is correct?" << endl;
    return(1);
  }
  else {
    cout << "Successfully opened Magboltz gas file:" << endl;
    cout << "---> " << gasFileName << endl;
  }

  // dimensions of TRIFIC chamber
  // approximated as a rectangular chamber
  // 30 deg. off normal = 60 deg. off beam axis
  // y,z: "radius"
  const double gridRotation = gridAngle*TMath::DegToRad(); // grid rotation angle w.r.t. beam axis
  const double chamberRadius = 7.341, chamberLength = 53.34; // chamber dimension in cm
  const double chamberVoltage = 0.; // chamber is ground
  const double gap = 1.3; // space between grids along beam axis
  // const double xgap = gap*TMath::Cos(gridRotation); // gap along x axis between wire grids in cm
  const double bracketToGrid = 0.44; // cm mylar window to first grid, braket to grid is 0.44 cm, total 1.94
  const int    nWires = 58; // number of wires/grid (yz-plane)
  //const int    nWires = 50; // number of wires/grid (xz-plane)
  //const int    nWires = 10; // number of wires/grid (58)
  const double dWire  = 0.0020; // wire diameter in cm = 20 um
  const double gridRadius = 5.775; // in cm (yz-plane)
  //const double gridRadius = 5.0; // in cm (xz-plane)

  // window as a bunch of wires
  // 100 um pitch ~ solid surface
  const int bwWires = 1662;
  //const int bwWires = 20;
  const double windowRadiusProfile = 3.5; // looking along beam axis
  const double windowRadius = windowRadiusProfile/TMath::Cos(TMath::PiOver2()-gridRotation); // cm, mylar window from 
  const double bracketRadius = 7.2/TMath::Cos(TMath::PiOver2()-gridRotation); // cm, window bracket
  const double windowToBracket = 1.5; // cm, bracket thickness
  const double bwXOffset = 4.5; // cm, offset from 0 in cm
  const double bwYOffset = 0.0; // cm, offset from 0 in cm
  const double px = TMath::Sin(TMath::PiOver2()-gridRotation)/TMath::Sin(gridRotation); // x coordinates for parallel wires
  const int npWires = 150; // number of parallel wires
  //const int npWires = 5; // number of parallel wires


  // for viewing, dimensions of cm
  double vcm = 0.5;
  double xvMin = -vcm,
    // xvMax = chamberLength + vcm,
    xvMax = chamberLength/3. + vcm,
    yvMin = -chamberRadius - vcm,
    yvMax = chamberRadius + vcm,
    zvMin = -chamberRadius - vcm,
    zvMax = chamberRadius + vcm;

  // create the TRIFIC geometry and gas
  GeometrySimple geo;
  SolidBox box(chamberLength/2.,0.,0.,
	       chamberLength/2.,chamberRadius,chamberRadius); // center x,y,z, half-length x,y,z
  // add the solid to the geometry, together with the gas inside.
  geo.AddSolid(&box, gas);

  // OpenGL driver error... seems to look OK before crashing
  // TString cgname = TString::Format("isopotential_flatchamber_Cathode%4.0fV_Anode%4.0fV",
  // 				   cathodeVoltage,anodeVoltage);
  // TCanvas cg(cgname.Data(),cgname.Data(),800,800);
  // ViewGeometry viewG;
  // viewG.SetGeometry(&geo);
  // viewG.SetCanvas(&cg);
  // cout << "Plotting geometry...\t";
  // viewG.Plot();
  // cout << "DONE!" << endl;
  // getc(stdin);
  
  cout << "Constructing TRIFIC chamber and grids... " << endl;
  ComponentAnalyticField driftCell;
  double windowVoltage = 0.0;
  // if(cathodeVoltage < 0)
  //   windowVoltage = 1.0; // for drawing tracks nicer

  // window and bracket "wires"
  for (int i=0; i<bwWires; i++) {
    // double yw = (2.*bracketRadius/static_cast<double>(nWires))*static_cast<double>(j)
    //  	-bracketRadius+gri%dRadius/static_cast<double>(nWires);
    double rw = (2.*bracketRadius/static_cast<double>(bwWires))*static_cast<double>(i)
      -bracketRadius+bracketRadius/static_cast<double>(bwWires);
    double xw = rw*TMath::Cos(gridRotation);
    if(rw > -windowRadius && rw < windowRadius)
      xw -= windowToBracket;
    double yw = rw*TMath::Sin(gridRotation);
    // TString bwName = TString::Format("bw%i",i);
    TString bwName = TString::Format("bw");
    driftCell.AddWire(bwXOffset+xw, bwYOffset+yw, dWire, windowVoltage, bwName.Data());
  }
  // parallel "wires" of window bracket
  for (int i=0; i<npWires; i++) {
    double xw = px*windowRadiusProfile-(windowToBracket/static_cast<double>(npWires))*static_cast<double>(i)
      -windowToBracket/static_cast<double>(npWires);
    double yw = windowRadiusProfile;
    // TString np1Name = TString::Format("np1_%i",i);
    TString np1Name = TString::Format("bw");
    driftCell.AddWire(bwXOffset+xw, bwYOffset+yw, dWire, windowVoltage, np1Name.Data());
  }
  for (int i=0; i<npWires; i++) {
    double xw = -px*windowRadiusProfile-(windowToBracket/static_cast<double>(npWires))*static_cast<double>(i)
      -windowToBracket/static_cast<double>(npWires);
    double yw = -windowRadiusProfile;
    // TString np2Name = TString::Format("np2_%i",i);
    TString np2Name = TString::Format("bw");
    driftCell.AddWire(bwXOffset+xw, bwYOffset+yw, dWire, windowVoltage, np2Name.Data());
  }

  // add wires for nGrids of wires at an angle w.r.t. beam axis. 
  // first grid is always cathode --> even grids are cathodes, odd grids are anodes
  const int nGrids = 21; // total number of grids
  int nC=0,nA=0; // number of anodes and cathodes
  for(int i=0; i<nGrids; i++) {
    // double xoffset = bracketToGrid*TMath::Cos(gridRotation) + static_cast<double>(i)*gap; 
    double xoffset = bwXOffset + bracketToGrid + static_cast<double>(i)*gap; 
    double yoffset = bwYOffset + 0.;
    // double yoffset = (bracketToGrid+1.3*static_cast<double>(i))*TMath::Sin(gridRotation); // y offset for rotated chamber
    for (int j=0; j<nWires; j++) {
      // double yw = (2.*gridRadius/static_cast<double>(nWires))*static_cast<double>(j)
      //  	-gridRadius+gridRadius/static_cast<double>(nWires);
      double rw = (2.*gridRadius/static_cast<double>(nWires))*static_cast<double>(j)
	-gridRadius+gridRadius/static_cast<double>(nWires);
      double xw = rw*TMath::Cos(gridRotation);
      double yw = rw*TMath::Sin(gridRotation);
      if(i%2 == 0) {
	// TString cathodeName = TString::Format("c%i",nC);
	 TString cathodeName = TString::Format("cathode");
	// driftCell.AddWire(xoffset, yoffset+yw, dWire, cathodeVoltage, cathodeName.Data());
	driftCell.AddWire(xoffset+xw, yoffset+yw, dWire, cathodeVoltage, cathodeName.Data());
	driftCell.AddReadout(cathodeName.Data());
	nC++;
      }
      else {
	// TString anodeName = TString::Format("a%i",nA);
	TString anodeName = TString::Format("anode");
	// driftCell.AddWire(xoffset, yoffset+yw, dWire, anodeVoltage, anodeName.Data());
	driftCell.AddWire(xoffset+xw, yoffset+yw, dWire, anodeVoltage, anodeName.Data());
	driftCell.AddReadout(anodeName.Data());
	nA++;
      }
    }
  }

  // TRIFIC chamber is a box looking perpendicular to the beam axis
  driftCell.AddPlaneX(0.0,chamberVoltage,"us_endcap");
  driftCell.AddPlaneX(chamberLength,chamberVoltage,"ds_endcap");
  driftCell.AddPlaneY(chamberRadius,chamberVoltage,"pos_radius");
  driftCell.AddPlaneY(-chamberRadius,chamberVoltage,"neg_radius");
  driftCell.SetGeometry(&geo);
  cout << "...DONE!" << endl;
  cout << "Grid width " << 2*gridRadius << " cm with " << nWires << " wires --> pitch = " 
       << 2.*gridRadius/static_cast<double>(nWires) << " cm." << endl;
  
  // make a sensor for the electrons
  cout << "Adding drift cell to sensor..." << endl;
  Sensor sensor;
  // sensor.DisableDebugging();
  // sensor.EnableDebugging();
  sensor.AddComponent(&driftCell);
  cout << "...DONE!" << endl;
  cout << "Request signal calculation on anode and cathode wires..." << endl;
  for(int i=0; i<nA; i++) {
    // TString wireName = TString::Format("a%d",i);
    TString wireName = TString::Format("anode");
    sensor.AddElectrode(&driftCell,wireName.Data());
  }
  for (int i=0; i<nC; i++) {
    // TString wireName = TString::Format("c%d",i);
    TString wireName = TString::Format("cathode");
    sensor.AddElectrode(&driftCell,wireName.Data());
  }
  cout << "DONE!" << endl;

  cout << "Setting signal calculation time limits..." << endl;
  //  Time in ns 
  const double tMin = 0.;
  const double tMax = 4000.; // = 4.0 us
  const double tStep = 10.; // 10 ns steps
  const int nTimeBins = static_cast<int>((tMax - tMin) / tStep);
  sensor.SetTimeWindow(0., tStep, nTimeBins);
  cout << "...DONE!" << endl;

  
  // // print out information on wires
  // // int nW = driftCell.GetNumberOfWires();
  // int nW = 40;
  // for (int i = 0; i < nW; ++i) {
  //   double xw, yw, dw, vw;
  //   std::string lbl;
  //   double lw, qw;
  //   int nw;
  //   driftCell.GetWire(i, xw, yw, dw, vw, lbl, lw, qw, nw);
  //   cout << i << "\t" << lbl << "\t" << xw << "\t" << yw << "\t" << vw << endl;
  // }
  // cout << "Press [Enter] to continue..." << endl;
  // getc(stdin);
  // cout << "Continuing..." << endl;

  // TApplication* theApp = new TApplication("theApp", &argc, argv);

  cout << "Plotting isopotential contour..." << endl;;
  TString cfname = TString::Format("isopotential_flatchamber_Cathode%.0fV_Anode%.0fV",
  				   cathodeVoltage,anodeVoltage);
  TCanvas cf(cfname.Data(), cfname.Data(), 800, 800);
  ViewField viewF;
  viewF.SetSensor(&sensor);
  viewF.SetArea(xvMin, yvMin, xvMax, yvMax);
  viewF.SetVoltageRange(1.05*cathodeVoltage, 1.05*anodeVoltage);
  viewF.SetCanvas(&cf);
  viewF.PlotContour("v");
  // viewF.Plot("v","cont1");
  cout << "...DONE!" << endl;
  // getc(stdin);

  cout << "Plotting electric field..." << endl;
  TString cename = TString::Format("electricfield_flatchamber_Cathode%.0fV_Anode%.0fV",
  				   cathodeVoltage,anodeVoltage);
  TCanvas ce(cename.Data(), cename.Data(), 800, 800);
  ViewField viewE;
  viewE.SetSensor(&sensor);
  viewE.SetArea(xvMin, yvMin, xvMax, yvMax);
  viewE.SetElectricFieldRange(0, 1.05*(anodeVoltage-cathodeVoltage));
  viewE.SetCanvas(&ce);
  viewE.PlotContour("e");
  // viewE.Plot("e","cont1");
  // viewE.PlotProfile(xvMin, yvMin, zvMin, xvMax, yvMax, zvMax, "e");
  cout << "...DONE" << endl;
  // ce.SetLogz();
   // ce.Draw();
   // theApp->Run(kTRUE);
  // getc(stdin);

  cout << "Creating canvas for plotting drift lines..." << endl;
  // visualization for drift lines to be computed later
  TString ccname = TString::Format("cell_drift_%sCathode%.0fV_Anode%.0fV",
  				   particleType.c_str(),cathodeVoltage,anodeVoltage);
  TCanvas cc(ccname.Data(), ccname.Data(), 800, 800);
  ViewCell viewC;
  viewC.SetComponent(&driftCell);
  viewC.SetArea(xvMin, yvMin, zvMin, xvMax, yvMax, zvMax);
  viewC.SetCanvas(&cc);
  ViewDrift viewD;
  viewD.SetArea(xvMin, yvMin, zvMin, xvMax, yvMax, zvMax);
  viewD.SetCanvas(&cc);
  cout << "...DONE!" << endl;
  // getc(stdin);

  // storage of electron drift lines
  cout << "Setting up electron tracking..." << endl;

  // Runge-Kutta-Fehlberg integration
  DriftLineRKF eDrift;
  eDrift.SetSensor(&sensor); // this is the evil part
  // const double maxStepSize=1.0; // cm
  // eDrift.SetMaximumStepSize(maxStepSize);
  // eDrift.EnableStepSizeLimit();
  // const unsigned int maxSteps = 10;
  // eDrift.SetMaxSteps(maxSteps);
  // const double accuracy = 1.;
  // eDrift.SetIntegrationAccuracy(accuracy);
  // eDrift.EnableDebugging();
  // eDrift.EnableVerbose();

  // Monte Carlo integration
  // AvalancheMC eDrift;
  // eDrift.SetSensor(&sensor);
  // eDrift.EnableSignalCalculation();
  // const double dist = 0.0002; // step distance in cm = 2 um
  // eDrift.SetDistanceSteps(dist);

  eDrift.EnablePlotting(&viewD);
  cout << "...DONE!" << endl;
  // getc(stdin);


  cout << "Initializing incident particle..." << endl;;
  // ionizing particle: [input energy] MeV [particle type]
  TrackHeed track;
  track.SetSensor(&sensor);
  track.SetParticle(particleType);
  track.SetKineticEnergy(energy);
  track.EnableElectricField();
  track.EnablePlotting(&viewD);

  // initialize ionizing particle position (in cm) and direction vector
  double azimuth = particleTheta, polar = TMath::PiOver2();
  double x0=bwXOffset-windowToBracket+0.1, y0, z0, t0, // real ion position + 0.1 cm
    dx = TMath::Cos(azimuth)*TMath::Sin(polar),
    dy = TMath::Sin(azimuth)*TMath::Sin(polar),
    dz = TMath::Cos(polar);
  y0 = z0 = t0 = 0.0;
  cout << "...DONE!" << endl;
  // getc(stdin);

  // start simulation
  sensor.ClearSignal();
  sensor.NewSignal();
  
  // create the track
  track.NewTrack(x0, y0, z0, t0, dx, dy, dz);
  cout << "\n=========== TrackHeed ==========="<<endl;
  cout << "Particle name        : " << track.GetName() << endl;
  cout << "Initial time         : t0 = " << t0 << endl;
  cout << "Initial position     : (" << x0 << "," << y0 << "," << z0 << ")" << endl;
  cout << "Initial direction    : (" << dx << "," << dy << "," << dz << ")" << endl;
  cout << "Mass [eV/c^2]        : " << track.GetMass() << "\tCharge [e]:" << track.GetCharge() << endl;
  cout << "Mass + eKin [eV/c^2] : " << track.GetEnergy() << "\tCharge [e]:" << track.GetCharge() << endl;
  cout << "Cluster Density      : " << track.GetClusterDensity() << " cm^-1" << endl;
  cout << "Stopping Power       : " << track.GetStoppingPower() << " eV/cm" << endl;
  cout << "Asymptotic W Value   : " << track.GetW() << " eV" << endl;
  cout << "Fano Factor          : " << track.GetFanoFactor() << endl;
  cout << "=================================\n" << endl;
  // getc(stdin);

  // position and time of the cluster and energy deposit
  // keep track of number of e-/cluster, clusters and total electrons
  double xcl,ycl,zcl,tcl,ecl,extra;
  // number of electrons in the cluster
  int ncl,
    // clusters counter
    nClusters=0,
    // total number of e-
    totNe=0,
    // total number of electrons drifted
    totNeDrift=0;
  
  TRandom2 randGen;
  randGen.SetSeed(theSeed);
  if(driftPercentage > 0) {
    cout << "BEGIN cluster generation, be patient..." << endl;
    while(track.GetCluster(xcl, ycl, zcl, tcl, ncl, ecl, extra)) {
      ++nClusters;
      totNe += ncl;

      // cout << "========================================================" << endl;
      // cout << "Cluster # " << nClusters << " @ (" << xcl << ", " << ycl << ", " << zcl
      // 	 << ") cm\ttime = " << tcl << " ns\t Edep = " << ecl << " eV" << endl; 	
      // cout << "Number of e- in the cluster: " << ncl << endl;
      // getc(stdin);

      // temporary variables for created electron
      // inital position, time, energy, direction
      double xe,ye,ze,te,EE,dxe,dye,dze; 	
      for(int ie = 0; ie < ncl; ++ie) {
	if(driftPercentage > 0.) {
	  double rndm = randGen.Rndm();
	
	  // track.GetElectron(ie,xe,ye,ze,te,EE,dxe,dye,dze);
	  // cout << "-------------------------------------------------------" << endl;
	  // cout << "Electron number " << ie << " @ (" << xe << ", " << ye << ", " << ze
	  // 	   << ") cm; time = " << te << " ns; E = " << EE << " eV;"
	  // 	   << " direction (" << dxe << ", " << dye << ", " << dze << ")" << endl;
	  // getc(stdin);
	
	  if(rndm < driftPercentage) {
	    // cout << "\r " << rndm  << " < " << driftPercentage << ": Drifting electron " << ie+1 << " of " << ncl << "..." << flush;
	    track.GetElectron(ie,xe,ye,ze,te,EE,dxe,dye,dze);
	    // only drift some % of electrons
	    // saves time for drawing
	    // cout << "Drifting..." << endl;
	    if(!eDrift.DriftElectron(xe,ye,ze,te)) {
	      // cout << "...Done! Continuing to next electron." << endl;
	      // continue;
	    }
	    totNeDrift++;
	  }
	  // else 
	  //   continue;
	}
      }
      cout << fixed;
      cout << setprecision(2);
      if(nClusters%100 == 0)
	cout << "\rDrifted " << totNeDrift << " of " << totNe << " (" << 100.*static_cast<double>(totNeDrift)/static_cast<double>(totNe) << "%) electrons in " << nClusters << " clusters..." << flush;

      // if(nClusters > 100) break;

    }

    cout << endl;
    cout << "...END cluster generation" << endl;
    cout << "TOTAL number of clusters: " << nClusters << endl;
    cout << "TOTAL number of drifted electrons: " << totNeDrift << endl;
    cout << "TOTAL number of electrons: " << totNe << endl;
  }

  cout << "Plotting drift lines..." << endl;
  bool is2D = true, drawAxes = true;
  viewD.Plot(is2D,drawAxes);
  viewC.Plot2d();
  cout << "...DONE!" << endl;
  // getc(stdin);
  // theApp->Run(kTRUE); 

  cout << "Saving data and exiting..." << endl;
  TString ffname = TString::Format("tiltTRIFIC_Cathode%.0fV_Anode%.0fV_%s%.0fMeV_%.0fDegrees_Drift%.2fPct_surf.root",
  				   cathodeVoltage,anodeVoltage,particleType.c_str(),energy/1.0e6,particleTheta*TMath::RadToDeg(),
				   driftPercentage*100.);
  TFile output(ffname.Data(),"RECREATE");
  cout << "Saving isopotential plot..." << endl;
  TString sname = TString::Format("%s.png",cfname.Data());
  cf.SaveAs(sname.Data());
  sname = TString::Format("%s.root",cfname.Data());
  cf.SaveAs(sname.Data());
  cf.Write();
  cout << "...DONE!" << endl;
  cout << "Saving electric field plot..." << endl;
  sname = TString::Format("%s.png",cename.Data());
  ce.SaveAs(sname.Data());
  sname = TString::Format("%s.root",cename.Data());
  ce.SaveAs(sname.Data());
  ce.Write();
  cout << "...DONE!" << endl;
  cout << "Saving driftline plot..." << endl;
  sname = TString::Format("%s.png",ccname.Data());
  cc.SaveAs(sname.Data());
  sname = TString::Format("%s.root",ccname.Data());
  cc.SaveAs(sname.Data());
  cc.Write();
  cout << "...DONE!" << endl;
  cout << "Calculating signals..." << endl;
  // sensor.SetTransferFunction(transfer);
  int cG = nC/nWires;
  int aG = nA/nWires;
  TH1D hist_scathode_ele[cG];
  TH1D hist_sanode_ele[aG];
  for(int i=0; i<cG; i++) {
    hist_scathode_ele[i] = GetGridElectronSignal(i, &sensor, nWires, true);
    hist_scathode_ele[i].Write();
  }
  for(int i=0; i<aG; i++) {
    hist_sanode_ele[i] = GetGridElectronSignal(i, &sensor, nWires, false);
    hist_sanode_ele[i].Write();
  }
  // sensor.ConvoluteSignal();
  cout << "...DONE!" << endl;
  output.Close();
  cout << "SUCCESS!" << endl;

  return 0;
}
