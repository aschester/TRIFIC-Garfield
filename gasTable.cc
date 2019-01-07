/****************************************************************
 *                                                              *
 * Program generates gas file for TRIFIC for given pressure     *
 * and reduced field range. This file is needed for simulating  *
 * response of TRIFIC detector and is application-specific.     *
 *                                                              *
 * Takes a while to run, be patient.                            *
 *                                                              *
 * A. Chester 7 June 2018                                       *
 *                                                              *
 ***************************************************************/

#include <iostream>
#include <cstdlib>

#include <TMath.h>
#include <TString.h>

#include "ComponentAnalyticField.hh"
#include "MediumMagboltz.hh"

using namespace Garfield;
using namespace std;

int makeGas(double Ptorr, int nEF, double EFmin, double EFmax) {

  MediumMagboltz* gas = new MediumMagboltz();
  gas->SetComposition("cf4",100.);
  gas->SetTemperature(293.15); // K
  gas->SetPressure(Ptorr); // Torr

  // flag for logarithmic spacing of field points
  bool useLog = true;
  if(nEF == 1) {
    useLog = false;
    EFmin = EFmax;
  }

  gas->SetFieldGrid(EFmin, EFmax, nEF, useLog);

  gas->SetMaxElectronEnergy(1.e3); // in eV - From Lars' alpha-g stuff

  gas->EnableDebugging();
  // gas->PrintGas();

  const int ncoll = 10; // number of collisions*10^7 for e- tracking
  const bool verbose = true;
  gas->GenerateGasTable(ncoll, verbose);

  TString filename = TString::Format("cf4_NT_%.1fTorr_nEF%.0i_EFmin%.2f_EFmax%.0f.gas",Ptorr,nEF,EFmin,EFmax);
  gas->WriteGasFile(filename.Data());

  return 0;
}

int main(int argc, char* argv[]) {

  double Ptorr = 91.2;    // gas pressure in Torr
  int nEF = 20;           // number of field points
  double EFmin = 100.;    // min reduced field in V/cm
  double EFmax = 200000.; // max reduced field in V/cm

  if(argc != 5) {
    cout << "Usage: ./trific arg1 arg2 arg3" << endl;
    cout << "arg1: CF4 gas pressure in Torr" << endl;
    cout << "arg2: Number of electric field points" << endl;
    cout << "arg3: Minimum reduced field in V/cm" << endl;
    cout << "arg4: Maximum reduced field in V/cm" << endl;
    exit(-1);
  }

  for(int i=1; i<argc ; i++) {
    switch(i) {
    case 1: Ptorr = atof(argv[i]); break; 
    case 2: nEF = atoi(argv[i]); break; 
    case 3: EFmin = atof(argv[i]); break; 
    case 4: EFmax = atof(argv[i]); break; 
    }
  }

  makeGas(Ptorr, nEF, EFmin, EFmax);

  return 0;
}
