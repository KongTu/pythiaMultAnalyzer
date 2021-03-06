#include <iostream>
#include <vector>
#include <sstream>
#include <string>

#include <eicsmear/erhic/EventBase.h>
#include <eicsmear/erhic/EventMC.h>
#include <eicsmear/erhic/EventPythia.h>
#include <eicsmear/erhic/Particle.h>
#include <eicsmear/erhic/ParticleMC.h>
#include <eicsmear/erhic/Pid.h>



#include "TParticlePDG.h"
#include "TString.h"
#include "TF1.h"
#include "TH1.h"
#include "TH2.h"
#include "TH3.h"
#include "TMath.h"
#include "TTree.h"
#include "TChain.h"
#include "TFile.h"
#include "TCanvas.h"
#include "TSystem.h"
#include "TROOT.h"
#include "TGraph.h"
#include "TGraphErrors.h"
#include "TGraphAsymmErrors.h"
#include "TMultiGraph.h"
#include "TCanvas.h"
#include "TPad.h"
#include "TLegend.h"
#include "TLatex.h"
#include "TLine.h"
#include "TAxis.h"
#include "TGraph.h"
#include "TGraphErrors.h"
#include "TLorentzVector.h"
#include "TBranchElement.h"
#include "TLorentzRotation.h"

#define PI            3.1415926

#define MASS_JPSI 	  3.09688
#define MASS_PROTON   0.93827
#define MASS_NEUTRON  0.93957
#define MASS_DEUTERON 1.8751019
#define MASS_TRITON   2.7937167208086358
#define MASS_HE3      2.7937167208086358
#define MASS_ALPHA    3.7249556277448477
#define MASS_LI6      5.5874334416172715
#define MASS_C12      11.174866883234543
#define MASS_CA40     37.249556277448477
#define MASS_XE131    121.99229680864376
#define MASS_AU197    183.45406466643374
#define MASS_PB208    193.69769264273208

using namespace erhic;
using namespace std;

TLorentzRotation BoostToHCM(TLorentzVector const &eBeam_lab,
                            TLorentzVector const &pBeam_lab,
                            TLorentzVector const &eScat_lab) {
   TLorentzVector q_lab=eBeam_lab - eScat_lab;
   TLorentzVector p_plus_q=pBeam_lab + q_lab;
   // boost to HCM
   TLorentzRotation boost=TLorentzRotation(-1.0*p_plus_q.BoostVector());
   TLorentzVector pBoost=boost*pBeam_lab;
   TVector3 axis=pBoost.BoostVector();
   // rotate away y-coordinate
   boost.RotateZ(-axis.Phi());
   // rotate away x-coordinate
   boost.RotateY(M_PI-axis.Theta());
   return boost;
}

TH1D* etaStar_gen = new TH1D("etaStar_gen",";#eta",200,-10,10);
TH1D* eta_gen = new TH1D("eta_gen",";#eta",200,-10,10);
TH1D* pt_gen = new TH1D("pt_gen",";p_{T} (GeV/c)",200,0,20);


void pythiaMultAnalyzer(int nEvents, TString inputFilename ){

	TChain *tree = new TChain("EICTree");

	//tree->Add("/eicdata/eic0009/PYTHIA/ep/TREES/pythia.ep.27x920.5Mevents.1.RadCor=0.Q2-0.1.root" ); // Wild cards are allowed e.g. tree.Add("*.root" );
	
	// tree->Add("/gpfs/mnt/gpfs02/eic/ztu/condor/job1/output/pythia_PDF=10800_1.root" );
	// tree->Add("/gpfs/mnt/gpfs02/eic/ztu/condor/job1/output/pythia_PDF=10800_2.root" );
	// tree->Add("/gpfs/mnt/gpfs02/eic/ztu/condor/job1/output/pythia_PDF=10800_3.root" );
	// tree->Add("/gpfs/mnt/gpfs02/eic/ztu/condor/job1/output/pythia_PDF=10800_4.root" );
	// tree->Add("/gpfs/mnt/gpfs02/eic/ztu/condor/job1/output/pythia_PDF=10800_5.root" );
	// tree->Add("/gpfs/mnt/gpfs02/eic/ztu/condor/job1/output/pythia_PDF=10800_6.root" );
	// tree->Add("/gpfs/mnt/gpfs02/eic/ztu/condor/job1/output/pythia_PDF=10800_7.root" );
	// tree->Add("/gpfs/mnt/gpfs02/eic/ztu/condor/job1/output/pythia_PDF=10800_8.root" );
	// tree->Add("/gpfs/mnt/gpfs02/eic/ztu/condor/job1/output/pythia_PDF=10800_9.root" );
	// tree->Add("/gpfs/mnt/gpfs02/eic/ztu/condor/job1/output/pythia_PDF=10800_10.root" );

	//tree->Add("/eicdata/eic0009/PYTHIA/ep/TREES/pythia.ep.27x920.5Mevents.1.RadCor=0.Q2-0.1.root" ); // Wild cards are allowed e.g. tree.Add("*.root" );
	tree->Add("/gpfs/mnt/gpfs02/eic/ztu/condor/job2/output/pythia_PDF=60800_*.root" );
	

	EventPythia* event(NULL);// = new EventPythia;

	// EventBase* event(NULL);
	// EventBeagle* event_beagle(NULL);

	tree->SetBranchAddress("event", &event ); // Note &event, not event.

	TLorentzVector ebeam(0,0,-27.,27.);
	TLorentzVector pbeam(0,0,460,460);

	TH1D* Nch_gen_target[21];
	TH1D* Nch_gen_current[21];

	double x_region[21];
	for(int j = 0; j < 21; j++){

		x_region[j] = 1e-5+j*5e-5;
		Nch_gen_target[j] = new TH1D(Form("Nch_gen_target_%d",j),"",50,0,50);
		Nch_gen_current[j] = new TH1D(Form("Nch_gen_current_%d",j),"",50,0,50);
	}

	for(int i(0); i < nEvents; ++i ) {
      
		// Read the next entry from the tree.
		tree->GetEntry(i);

		//event information:
		double trueQ2 = event->GetTrueQ2();
		double trueW2 = event->GetTrueW2();
		double trueX = event->GetTrueX();
		double trueY = event->GetTrueY();
		double trueNu = event->GetTrueNu();
		double s_hat = event->GetHardS();
		double t_hat = event->t_hat;
		double u_hat = event->GetHardU();
		double photon_flux = event->GetPhotonFlux();
		int event_process = event->GetProcess();
		int nParticles = event->GetNTracks();
		int struck_nucleon = event->nucleon;
		
		if( event_process != 99 ) continue;
		if( trueQ2 < 10 || trueQ2 > 11 ) continue;
		if( trueX < 0.00001 || trueX > 0.0011 ) continue;
		
		int x_index = -1;

		for(int j = 0; j < 20; j++){

			if(trueX > x_region[j] && trueX < x_region[j+1]) x_index = j;
		}

		int nParticles_process = 0;
		int nParticles_process_current = 0;

		TLorentzVector scat_e;
		TLorentzVector part4v;
		TLorentzVector part4vStar;
		TLorentzRotation boost_REC_HCM;

		for(int j(0); j < nParticles; ++j ) {

			const erhic::ParticleMC* particle = event->GetTrack(j);

			TParticlePDG* info = particle->Id().Info();
			int charge = info->Charge();
			int pdg = particle->GetPdgCode();
			int status = particle->GetStatus();
			double pt = particle->GetPt();
			double eta = particle->GetEta();
			int index = particle->GetIndex();

			if( status != 1 ) continue;
			if( index == 3 ) scat_e = particle->Get4Vector();
			if( pdg == 22 || fabs(pdg) == 11 ) continue;
			if( charge == 0 ) continue;
			
			part4v = particle->Get4Vector();
			boost_REC_HCM=BoostToHCM(ebeam,pbeam,scat_e);
			part4vStar = boost_REC_HCM*part4v;

			pt_gen->Fill( part4vStar.Pt() );
			etaStar_gen->Fill( part4vStar.Eta() );
			eta_gen->Fill( eta );

			if( part4vStar.Pt() < 0.0 ) continue;
			
			if( part4vStar.Eta() < 0.0 ) {
				nParticles_process++;
			}
			if( part4vStar.Eta() > 0.0){
				nParticles_process_current++;
			}

		} // end of particle loop
		if( x_index >= 0 ){
			Nch_gen_target[x_index]->Fill( nParticles_process );
			Nch_gen_current[x_index]->Fill( nParticles_process_current );
		} 
		
		
	}

	TString outfilename;
	outfilename = "_pythia_multiplicity.root";


   	TFile output("../rootfiles/"+inputFilename+outfilename,"RECREATE");
   	for(int j = 0; j < 20; j++){
		Nch_gen_target[j]->Write();
	   	Nch_gen_current[j]->Write();
   	}
   	
   	pt_gen->Write();
   	etaStar_gen->Write();
   	eta_gen->Write();


}