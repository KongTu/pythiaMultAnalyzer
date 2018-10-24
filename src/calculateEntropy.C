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
	TH1D* Nch_gen_target = (TH1D*) file->Get("Nch_gen_target");
	TH1D* Nch_gen_current = (TH1D*) file->Get("Nch_gen_current");
	
	Nch_gen_target->Scale(1.0/Nch_gen_target->Integral());
	Nch_gen_current->Scale(1.0/Nch_gen_current->Integral());
	
	Nch_gen_current->Draw();

	double entropy_target = 0.0;
	for(int i = 0; i < Nch_gen_target->GetNbinsX(); i++){

		double pn = Nch_gen_target->GetBinContent(i+1);
		if( pn == 0. ) continue;
		entropy_target += (-pn)*TMath::Log(pn);
	}

	cout << "entropy_target is " << entropy_target << endl;

	double entropy_current = 0.0;
	for(int i = 0; i < Nch_gen_current->GetNbinsX(); i++){

		double pn = Nch_gen_current->GetBinContent(i+1);
		if( pn == 0. ) continue;
		entropy_current += (-pn)*TMath::Log(pn);
	}

	cout << "entropy_current is " << entropy_current << endl;

}