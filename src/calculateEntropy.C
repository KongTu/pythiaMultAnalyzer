#include <memory>
#include <string>
#include <vector>
#include <iostream>
#include <fstream>
#include <math.h>
#include <map>
#include <sstream>
#include <iomanip>

#include "TFile.h"
#include "TH1.h"
#include "TCanvas.h"


using namespace std;


void calculateEntropy(){


	TFile* file = new TFile("../rootfiles/100k_CT10_pythia_multiplicity.root");
	TH1D* Ntrk_gen = (TH1D*) file->Get("Ntrk_gen");
	Ntrk_gen->Draw();

}