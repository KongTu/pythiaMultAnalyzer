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
#include "TMath.h"

using namespace std;


void calculateEntropy(){


	TFile* file = new TFile("../rootfiles/100k_CT10_pythia_multiplicity.root");
	TH1D* Ntrk_gen = (TH1D*) file->Get("Nch_gen");
	
	Ntrk_gen->Scale(1.0/Ntrk_gen->Integral());
	Ntrk_gen->Draw();

	double entropy = 0.0;
	for(int i = 0; i < Ntrk_gen->GetNbinsX(); i++){

		double pn = Ntrk_gen->GetBinContent(i+1);

		entropy += (-pn)*TMath::Log(pn);
	}

	cout << "entropy is " << entropy << endl;

}