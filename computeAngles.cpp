
#include "LHEF.h"
#include "TROOT.h"
#include "TH1.h"
#include "TH2.h"
#include "TFile.h"
#include "THStack.h"
#include "TCanvas.h"
#include "TLorentzVector.h"
#include <cmath>
#include "TTree.h"

using namespace std;

void computeAngles(TLorentzVector thep4H, TLorentzVector thep4Z1, TLorentzVector thep4M11, TLorentzVector thep4M12, TLorentzVector thep4Z2, TLorentzVector thep4M21, TLorentzVector thep4M22, double& costheta1, double& costheta2, double& Phi, double& costhetastar, double& Phi1);

void computeAnglesTTbar(TLorentzVector p4b1, TLorentzVector p4W1u, TLorentzVector p4W1d, TLorentzVector p4b2, TLorentzVector p4W2u, TLorentzVector p4W2d,
			double& costhetastar, double& costheta1, double& costheta2, double& costhetaW1, double& costhetaW2, 
			double& Phi, double& Phi1, double& PhiW1, double& PhiW2);

double deltaPhi (double phi1, double phi2)
{
    double deltaphi=fabs(phi1-phi2);
    if (deltaphi > 6.283185308) deltaphi -= 6.283185308;
    if (deltaphi > 3.141592654) deltaphi = 6.283185308-deltaphi;
    return deltaphi;
}


//! ========================================================================================


int main (int argc, char **argv) {
    
    TFile file(argv[2],"RECREATE");
    TTree* tree = new TTree("tree","tree");
    
    //float mWW, mWLep, mWHad; 
    float costheta1, costheta2, costhetastar, costhetaW1, costhetaW2;
    float phi, phi1, phiW1, phiW2;
    //float dEtajj, dPhijj, mjj;
    //int isSignal,isMuMinus;
    
    //tree->Branch("mWW",&mWW,"mWW/F");
    //tree->Branch("mWLep",&mWLep,"mWLep/F");
    //tree->Branch("mWHad",&mWHad,"mWHad/F");
    tree->Branch("costheta1",&costheta1,"costheta1/F");
    tree->Branch("costheta2",&costheta2,"costheta2/F");
    tree->Branch("costhetaW1",&costhetaW1,"costhetaW1/F");
    tree->Branch("costhetaW2",&costhetaW2,"costhetaW2/F");
    tree->Branch("costhetastar",&costhetastar,"costhetastar/F");
    tree->Branch("phi",&phi,"phi/F");
    tree->Branch("phi1",&phi1,"phi1/F");
    tree->Branch("phiW1",&phiW1,"phiW1/F");
    tree->Branch("phiW2",&phiW2,"phiW2/F");
    //tree->Branch("dEtajj",&dEtajj,"dEtajj/F");
    //tree->Branch("dPhijj",&dPhijj,"dPhijj/F");
    //tree->Branch("mjj",&mjj,"mjj/F");
    
    //tree->Branch("isSignal",&isSignal,"isSignal/I");
    //tree->Branch("isMuMinus",&isMuMinus,"isMuMinus/I");
    
    //PG loop over bkg
    //PG -------------
    
    int BKGnumber = 0 ;
    int BKGnumberWithObj = 0 ;
    int VBFBKGnumber = 0 ;
    
    int NSignal = 0;
    int NTotal = 0;
    
    std::ifstream ifsbkg (argv[1]) ;
    // Create the Reader object
    LHEF::Reader bkgReader (ifsbkg) ;
    
    //PG loop over BKG
    while ( bkgReader.readEvent () ) {
        ++BKGnumber;
        if (BKGnumber % 1000 == 0) 
	    std::cout << "BKG event " << BKGnumber << "\n" ;
        //if (BKGnumber > 500000) break;

        // --------------- Indices for final state particles  -------------------
        int i_b_1 = -1; // defined to be the pdg id = 5
        int i_Wu_1 = -1;        
        int i_Wd_1 = -1;
        int i_b_2 = -1;        
        int i_Wu_2 = -1;
        int i_Wd_2 = -1;        
        
        std::vector<int> leptons ;      
        std::vector<int> finalQuarks ;      
        std::vector<int> initialQuarks ;   
        std::vector<int> intermediates ;
        std::vector<int> tops ;
	std::vector<int> Ws ; 

        //PG loop over particles in the event
        for (int iPart = 0 ; iPart < bkgReader.hepeup.IDUP.size (); ++iPart){
            
            int mother1 = bkgReader.hepeup.MOTHUP.at(iPart).first - 1;
	    int mother1ID = -1;
	    if (mother1 != -1)
		mother1ID = bkgReader.hepeup.IDUP.at(mother1);
	    int myID = bkgReader.hepeup.IDUP.at(iPart);
        
	    if (BKGnumber % 1000 == 0)
	    {
		std::cout << "\t part type [" << iPart << "] " << myID
			  << "\t status " << bkgReader.hepeup.ISTUP.at(iPart)
			  << "\t mother index " << mother1;
		if(mother1!=-1)
		{
		    std::cout << "\t mother id " << mother1ID;
		}
		std::cout << "\n" ;
	    }
            //PG incoming particle          
            if (bkgReader.hepeup.ISTUP.at (iPart) == -1){
                initialQuarks.push_back (iPart) ;
            }

            //PG outgoing particles          
            if (bkgReader.hepeup.ISTUP.at (iPart) == 1){
		// find b quarks from top quark
		if( myID == 5 && mother1ID == 6)
		    i_b_1 = iPart;
		else if( myID == 5 && mother1ID == -6) // in case there is a three-body decay
		    i_Wd_2 = iPart;
		else if( myID == -5 && mother1ID == -6)
		    i_b_2 = iPart;
		else if( myID == -5 && mother1ID == 6)
		    i_Wd_1 = iPart; 
		// other quarks and leptons from W decay (could be three-body decay)
		else if( mother1ID == 24 || mother1ID == 6)
		{
		    if( myID == 2 || myID == 4 || myID == 12 || myID == 14 || myID == 16)
			i_Wu_1 = iPart;
		    else if( myID == -1 || myID == -3 || myID == -5 || myID == -11 || myID == -13 || myID == -15)
			i_Wd_1 = iPart;
		}
		// And W- decay products
		if( mother1ID == -24 || mother1ID == -6)
		{
		    if( myID == -2 || myID == -4 || myID == -12 || myID == -14 || myID == -16)
			i_Wu_2 = iPart;
		    else if( myID == 1 || myID == 3 || myID == 5 || myID == 11 || myID == 13 || myID == 15)
			i_Wd_2 = iPart;
		}

	       
                //PG leptons
                if (abs (bkgReader.hepeup.IDUP.at (iPart)) == 11 ||   //PG electron
                    abs (bkgReader.hepeup.IDUP.at (iPart)) == 13 ||   //PG muon
                    abs (bkgReader.hepeup.IDUP.at (iPart)) == 15 ||   //PG tau
                    abs (bkgReader.hepeup.IDUP.at (iPart)) == 12 ||   //PG neutrino
                    abs (bkgReader.hepeup.IDUP.at (iPart)) == 14 ||   //PG neutrino
                    abs (bkgReader.hepeup.IDUP.at (iPart)) == 16)     //PG neutrino                    
                    {
                    leptons.push_back (iPart) ;
                    } //PG leptons
                else
                    {
                    finalQuarks.push_back (iPart) ;
                    }
                
            } 
            
            //PG intermediates
            if (bkgReader.hepeup.ISTUP.at(iPart) == 2){
                intermediates.push_back (iPart) ;
            }
            
            //PG tops
            if (abs(bkgReader.hepeup.IDUP.at(iPart)) == 6){
                tops.push_back (iPart) ;
            }

            //PG Ws
            if (abs(bkgReader.hepeup.IDUP.at(iPart)) == 24){
                Ws.push_back (iPart) ;
            }

        } //PG loop over particles in the event

        TLorentzVector fs_b_1(
	    bkgReader.hepeup.PUP.at(i_b_1).at(0), //PG px
	    bkgReader.hepeup.PUP.at(i_b_1).at(1), //PG py
	    bkgReader.hepeup.PUP.at(i_b_1).at(2), //PG pz
	    bkgReader.hepeup.PUP.at(i_b_1).at(3) //PG E
	    ) ;
        TLorentzVector fs_Wu_1(
	    bkgReader.hepeup.PUP.at(i_Wu_1).at(0), //PG px
	    bkgReader.hepeup.PUP.at(i_Wu_1).at(1), //PG py
	    bkgReader.hepeup.PUP.at(i_Wu_1).at(2), //PG pz
	    bkgReader.hepeup.PUP.at(i_Wu_1).at(3) //PG E
	    ) ;
	
        TLorentzVector fs_Wd_1(
	    bkgReader.hepeup.PUP.at(i_Wd_1).at(0), //PG px
	    bkgReader.hepeup.PUP.at(i_Wd_1).at(1), //PG py
	    bkgReader.hepeup.PUP.at(i_Wd_1).at(2), //PG pz
	    bkgReader.hepeup.PUP.at(i_Wd_1).at(3) //PG E
	    ) ;
        TLorentzVector fs_b_2(
	    bkgReader.hepeup.PUP.at(i_b_2).at(0), //PG px
	    bkgReader.hepeup.PUP.at(i_b_2).at(1), //PG py
	    bkgReader.hepeup.PUP.at(i_b_2).at(2), //PG pz
	    bkgReader.hepeup.PUP.at(i_b_2).at(3) //PG E
	    ) ;
        TLorentzVector fs_Wu_2(
	    bkgReader.hepeup.PUP.at(i_Wu_2).at(0), //PG px
	    bkgReader.hepeup.PUP.at(i_Wu_2).at(1), //PG py
	    bkgReader.hepeup.PUP.at(i_Wu_2).at(2), //PG pz
	    bkgReader.hepeup.PUP.at(i_Wu_2).at(3) //PG E
	    ) ;
        TLorentzVector fs_Wd_2(
	    bkgReader.hepeup.PUP.at(i_Wd_2).at(0), //PG px
	    bkgReader.hepeup.PUP.at(i_Wd_2).at(1), //PG py
	    bkgReader.hepeup.PUP.at(i_Wd_2).at(2), //PG pz
	    bkgReader.hepeup.PUP.at(i_Wd_2).at(3) //PG E
	    ) ;
        
        //TLorentzVector p4_W1 = fs_Wu_1 + fs_Wd_1;
        //TLorentzVector p4_W2 = fs_Wu_2 + fs_Wd_2;        
        
        double a_costheta1, a_costheta2, a_costhetastar, a_costhetaW1, a_costhetaW2, a_Phi, a_Phi1, a_PhiW1, a_PhiW2;
        computeAnglesTTbar( fs_b_1, fs_Wu_1, fs_Wd_1, fs_b_2, fs_Wu_2, fs_Wd_2,
			    a_costhetastar, a_costheta1, a_costheta2, a_costhetaW1, a_costhetaW2, a_Phi, a_Phi1, a_PhiW1, a_PhiW2);

	costheta1 = (float) a_costheta1;                
        costheta2 = (float) a_costheta2;
	costhetaW1 = (float) a_costhetaW1;                
        costhetaW2 = (float) a_costhetaW2;
        phi = (float) a_Phi;
        costhetastar = (float) a_costhetastar;
        phi1 = (float) a_Phi1;
        phiW1 = (float) a_PhiW1;
        phiW2 = (float) a_PhiW2;

/*
       
        //mWW = (float) p4_WW.M();
        //mWLep = (float) p4_WLep.M();
        //mWHad = (float) p4_WHad.M();        

        dEtajj = (float) fabs( fs_Iqrk0.Eta() - fs_Iqrk1.Eta() );
        dPhijj = (float) deltaPhi(fs_Iqrk0.Phi(),fs_Iqrk1.Phi());     
        mjj = (float) (fs_Iqrk0 + fs_Iqrk1).M();

        isSignal = signalFlag;
*/        
        tree->Fill();
        
    }
    
    std::cout << "number of events = " << BKGnumber << std::endl;
    
    file.cd();
    tree->Write();
    file.Close();
    
    
// Now we are done.
return 0 ;
}

//////////////////////////////////
//// P A P E R   4 - V E C T O R   D E F I N I T I O N   O F   P H I   A N D   P H I 1
//////////////////////////////////
void computeAnglesTTbar(
  TLorentzVector p4b1,
  TLorentzVector p4W1u,
  TLorentzVector p4W1d,
  TLorentzVector p4b2, 
  TLorentzVector p4W2u,
  TLorentzVector p4W2d,
  double& costhetastar,
  double& costheta1,
  double& costheta2,
  double& costhetaW1,
  double& costhetaW2,
  double& Phi,
  double& Phi1,
  double& PhiW1,
  double& PhiW2
  ){

  //build W 4-vectors
  TLorentzVector p4W1 = p4W1u + p4W1d;
  TLorentzVector p4W2 = p4W2u + p4W2d;

  //build top 4-vectors
  TLorentzVector p4t1 = p4b1 + p4W1u + p4W1d;
  TLorentzVector p4t2 = p4b2 + p4W2u + p4W2d;

  // BEGIN THE CALCULATION

  // build H 4-vectors
  TLorentzVector p4H = p4t1 + p4t2;

  // -----------------------------------

  //// costhetastar  --> Nadja
  TVector3 boostX = -(p4H.BoostVector());
  TLorentzVector thep4t1inXFrame(p4t1);
  TLorentzVector thep4t2inXFrame(p4t2);
  thep4t1inXFrame.Boost(boostX);
  thep4t2inXFrame.Boost(boostX);
  TVector3 thet1X_p3 = TVector3(thep4t1inXFrame.X(), thep4t1inXFrame.Y(), thep4t1inXFrame.Z());
  TVector3 thet2X_p3 = TVector3(thep4t2inXFrame.X(), thep4t2inXFrame.Y(), thep4t2inXFrame.Z());
  costhetastar = thet1X_p3.CosTheta();

  TVector3 boostV1(0, 0, 0);
  TVector3 boostV2(0, 0, 0);
  //// --------------------------- costheta1
  //if (!(fabs(t1_bId)==5 || fabs(t1_uId)==2 || fabs(t1_dId)==1 || fabs(t1_uId)==4 || fabs(t1_dId)==3 )){ // can probably remove this for now?
    boostV1 = -(p4t1.BoostVector());
    if (boostV1.Mag()>=1.) {
      cout << "Warning: Mela::computeAngles: top1 boost with beta=1, scaling down" << endl;
      boostV1*=0.9999/boostV1.Mag();
    }
    TLorentzVector p4b1_BV1(p4b1);
    TLorentzVector p4W1u_BV1(p4W1u);
    TLorentzVector p4W1d_BV1(p4W1d);
    TLorentzVector p4b2_BV1(p4b2);
    TLorentzVector p4W2u_BV1(p4W2u);
    TLorentzVector p4W2d_BV1(p4W2d);
    p4b1_BV1.Boost(boostV1);
    p4W1u_BV1.Boost(boostV1);
    p4W1d_BV1.Boost(boostV1);
    p4b2_BV1.Boost(boostV1);
    p4W2u_BV1.Boost(boostV1);
    p4W2d_BV1.Boost(boostV1);

    TLorentzVector p4V2_BV1 = p4b2_BV1 + p4W2u_BV1 + p4W2d_BV1;
    //// costheta1
    costheta1 = -p4V2_BV1.Vect().Unit().Dot(p4b1_BV1.Vect().Unit());
    //}
    //else costheta1 = 0;

  //// --------------------------- costheta2
  //if (!(fabs(t1_bId)==5 || fabs(t1_uId)==2 || fabs(t1_dId)==1 || fabs(t1_uId)==4 || fabs(t1_dId)==3 )){
    boostV2 = -(p4t2.BoostVector());
    if (boostV2.Mag()>=1.) {
      cout << "Warning: Mela::computeAngles: top2 boost with beta=1, scaling down" << endl;
      boostV2*=0.9999/boostV2.Mag();
    }
    TLorentzVector p4b1_BV2(p4b1);
    TLorentzVector p4W1u_BV2(p4W1u);
    TLorentzVector p4W1d_BV2(p4W1d);
    TLorentzVector p4b2_BV2(p4b2);
    TLorentzVector p4W2u_BV2(p4W2u);
    TLorentzVector p4W2d_BV2(p4W2d);
    p4b1_BV2.Boost(boostV2);
    p4W1u_BV2.Boost(boostV2);
    p4W1d_BV2.Boost(boostV2);
    p4b2_BV2.Boost(boostV2);
    p4W2u_BV2.Boost(boostV2);
    p4W2d_BV2.Boost(boostV2);

    TLorentzVector p4V1_BV2 = p4b1_BV2 + p4W1u_BV2 + p4W1d_BV2;
    //// costheta2
    costheta2 = -p4V1_BV2.Vect().Unit().Dot(p4b2_BV2.Vect().Unit());
    //}
    //else costheta2 = 0;

  TVector3 boostV1_W(0, 0, 0);
  TVector3 boostV2_W(0, 0, 0);
  //// --------------------------- costhetaW1
  //if (!(fabs(t1_bId)==5 || fabs(t1_uId)==2 || fabs(t1_dId)==1 || fabs(t1_uId)==4 || fabs(t1_dId)==3 )){ // can probably remove this for now?
    boostV1_W = -(p4W1.BoostVector());
    if (boostV1_W.Mag()>=1.) {
      cout << "Warning: Mela::computeAngles: W1 boost with beta=1, scaling down" << endl;
      boostV1_W*=0.9999/boostV1_W.Mag();
    }
    TLorentzVector p4b1_BV1_W(p4b1);
    TLorentzVector p4W1u_BV1_W(p4W1u);
    TLorentzVector p4W1d_BV1_W(p4W1d);
    TLorentzVector p4b2_BV1_W(p4b2);
    TLorentzVector p4W2u_BV1_W(p4W2u);
    TLorentzVector p4W2d_BV1_W(p4W2d);
    p4b1_BV1_W.Boost(boostV1_W);
    p4W1u_BV1_W.Boost(boostV1_W);
    p4W1d_BV1_W.Boost(boostV1_W);
    p4b2_BV1_W.Boost(boostV1_W);
    p4W2u_BV1_W.Boost(boostV1_W);
    p4W2d_BV1_W.Boost(boostV1_W);

    TLorentzVector p4V2_BV1_W = p4b1_BV1_W;
    //// costheta1
    costhetaW1 = -p4V2_BV1_W.Vect().Unit().Dot(p4W1d_BV1_W.Vect().Unit());
    //}
    //else costhetaW1 = 0;

  //// --------------------------- costhetaW2
  //if (!(fabs(t1_bId)==5 || fabs(t1_uId)==2 || fabs(t1_dId)==1 || fabs(t1_uId)==4 || fabs(t1_dId)==3 )){
    boostV2_W = -(p4W2.BoostVector());
    if (boostV2_W.Mag()>=1.) {
      cout << "Warning: Mela::computeAngles: W2 boost with beta=1, scaling down" << endl;
      boostV2_W*=0.9999/boostV2_W.Mag();
    }
    TLorentzVector p4b1_BV2_W(p4b1);
    TLorentzVector p4W1u_BV2_W(p4W1u);
    TLorentzVector p4W1d_BV2_W(p4W1d);
    TLorentzVector p4b2_BV2_W(p4b2);
    TLorentzVector p4W2u_BV2_W(p4W2u);
    TLorentzVector p4W2d_BV2_W(p4W2d);
    p4b1_BV2_W.Boost(boostV2_W);
    p4W1u_BV2_W.Boost(boostV2_W);
    p4W1d_BV2_W.Boost(boostV2_W);
    p4b2_BV2_W.Boost(boostV2_W);
    p4W2u_BV2_W.Boost(boostV2_W);
    p4W2d_BV2_W.Boost(boostV2_W);

    TLorentzVector p4V1_BV2_W = p4b2_BV2_W;
    //// costheta2
    costhetaW2 = -p4V1_BV2_W.Vect().Unit().Dot(p4W2d_BV2_W.Vect().Unit());
    //}
    //else costhetaW2 = 0;

  //// --------------------------- Phi and Phi1 (old phistar1 - azimuthal production angle)
  TLorentzVector p4b1_BX(p4b1);
  TLorentzVector p4W1u_BX(p4W1u);
  TLorentzVector p4W1d_BX(p4W1d);
  TLorentzVector p4b2_BX(p4b2);
  TLorentzVector p4W2u_BX(p4W2u);
  TLorentzVector p4W2d_BX(p4W2d);

  p4b1_BX.Boost(boostX);
  p4W1u_BX.Boost(boostX);
  p4W1d_BX.Boost(boostX);
  p4b2_BX.Boost(boostX);
  p4W2u_BX.Boost(boostX);
  p4W2d_BX.Boost(boostX);
  TLorentzVector p4W1_BX = p4W1u_BX + p4W1d_BX;
  TLorentzVector p4W2_BX = p4W2u_BX + p4W2d_BX;
  TLorentzVector p4V1_BX = p4b1_BX + p4W1u_BX + p4W1d_BX;

  TVector3 beamAxis(0, 0, 1);
  TVector3 p3V1_BX = p4V1_BX.Vect().Unit();
  TVector3 normal1_BX = (p4b1_BX.Vect().Cross(p4W1_BX.Vect())).Unit();
  TVector3 normal2_BX = (p4b2_BX.Vect().Cross(p4W2_BX.Vect())).Unit();
  TVector3 normalSC_BX = (beamAxis.Cross(p3V1_BX)).Unit();


  //// Phi
  float tmpSgnPhi = p3V1_BX.Dot(normal1_BX.Cross(normal2_BX));
  float sgnPhi = 0;
  if (fabs(tmpSgnPhi)>0.) sgnPhi = tmpSgnPhi/fabs(tmpSgnPhi);
  float dot_BX12 = normal1_BX.Dot(normal2_BX);
  if (fabs(dot_BX12)>=1.) dot_BX12 *= 1./fabs(dot_BX12);
  Phi = sgnPhi * acos(-1.*dot_BX12);


  //// Phi1
  double tmpSgnPhi1 = p3V1_BX.Dot(normal1_BX.Cross(normalSC_BX));
  double sgnPhi1 = 0;
  if (fabs(tmpSgnPhi1)>0.) sgnPhi1 = tmpSgnPhi1/fabs(tmpSgnPhi1);
  double dot_BX1SC = normal1_BX.Dot(normalSC_BX);
  if (fabs(dot_BX1SC)>=1.) dot_BX1SC *= 1./fabs(dot_BX1SC);
  Phi1 = sgnPhi1 * acos(dot_BX1SC);

  if (isnan(costhetastar) || isnan(costheta1) || isnan(costheta2) || isnan(Phi) || isnan(Phi1)){
    cout << "WARNING: NaN in computeAngles: "
      << costhetastar << " "
      << costheta1  << " "
      << costheta2  << " "
      << Phi  << " "
      << Phi1  << " " << endl;
    cout << "   boostV1: " <<boostV1.Pt() << " " << boostV1.Eta() << " " << boostV1.Phi() << " " << boostV1.Mag() << endl;
    cout << "   boostV2: " <<boostV2.Pt() << " " << boostV2.Eta() << " " << boostV2.Phi() << " " << boostV2.Mag() << endl;
  }
}





void computeAngles(TLorentzVector thep4H, TLorentzVector thep4Z1, TLorentzVector thep4M11, TLorentzVector thep4M12, TLorentzVector thep4Z2, TLorentzVector thep4M21, TLorentzVector thep4M22, double& costheta1, double& costheta2, double& Phi, double& costhetastar, double& Phi1){
    
    ///////////////////////////////////////////////
    // check for z1/z2 convention, redefine all 4 vectors with convention
    ///////////////////////////////////////////////	
    TLorentzVector p4H, p4Z1, p4M11, p4M12, p4Z2, p4M21, p4M22;
    p4H = thep4H;
    
    p4Z1 = thep4Z1; p4M11 = thep4M11; p4M12 = thep4M12;
    p4Z2 = thep4Z2; p4M21 = thep4M21; p4M22 = thep4M22;
    //// costhetastar
	TVector3 boostX = -(thep4H.BoostVector());
	TLorentzVector thep4Z1inXFrame( p4Z1 );
	TLorentzVector thep4Z2inXFrame( p4Z2 );	
	thep4Z1inXFrame.Boost( boostX );
	thep4Z2inXFrame.Boost( boostX );
	TVector3 theZ1X_p3 = TVector3( thep4Z1inXFrame.X(), thep4Z1inXFrame.Y(), thep4Z1inXFrame.Z() );
	TVector3 theZ2X_p3 = TVector3( thep4Z2inXFrame.X(), thep4Z2inXFrame.Y(), thep4Z2inXFrame.Z() );    
    costhetastar = theZ1X_p3.CosTheta();
    
    //// --------------------------- costheta1
    TVector3 boostV1 = -(thep4Z1.BoostVector());
    TLorentzVector p4M11_BV1( p4M11 );
	TLorentzVector p4M12_BV1( p4M12 );	
    TLorentzVector p4M21_BV1( p4M21 );
	TLorentzVector p4M22_BV1( p4M22 );
    p4M11_BV1.Boost( boostV1 );
	p4M12_BV1.Boost( boostV1 );
	p4M21_BV1.Boost( boostV1 );
	p4M22_BV1.Boost( boostV1 );
    
    TLorentzVector p4V2_BV1 = p4M21_BV1 + p4M22_BV1;
    //// costheta1
    costheta1 = -p4V2_BV1.Vect().Dot( p4M11_BV1.Vect() )/p4V2_BV1.Vect().Mag()/p4M11_BV1.Vect().Mag();
    
    //// --------------------------- costheta2
    TVector3 boostV2 = -(thep4Z2.BoostVector());
    TLorentzVector p4M11_BV2( p4M11 );
	TLorentzVector p4M12_BV2( p4M12 );	
    TLorentzVector p4M21_BV2( p4M21 );
	TLorentzVector p4M22_BV2( p4M22 );
    p4M11_BV2.Boost( boostV2 );
	p4M12_BV2.Boost( boostV2 );
	p4M21_BV2.Boost( boostV2 );
	p4M22_BV2.Boost( boostV2 );
    
    TLorentzVector p4V1_BV2 = p4M11_BV2 + p4M12_BV2;
    //// costheta2
    costheta2 = -p4V1_BV2.Vect().Dot( p4M21_BV2.Vect() )/p4V1_BV2.Vect().Mag()/p4M21_BV2.Vect().Mag();
    
    //// --------------------------- Phi and Phi1
    //    TVector3 boostX = -(thep4H.BoostVector());
    TLorentzVector p4M11_BX( p4M11 );
	TLorentzVector p4M12_BX( p4M12 );	
    TLorentzVector p4M21_BX( p4M21 );
	TLorentzVector p4M22_BX( p4M22 );	
    
	p4M11_BX.Boost( boostX );
	p4M12_BX.Boost( boostX );
	p4M21_BX.Boost( boostX );
	p4M22_BX.Boost( boostX );
    
    TVector3 tmp1 = p4M11_BX.Vect().Cross( p4M12_BX.Vect() );
    TVector3 tmp2 = p4M21_BX.Vect().Cross( p4M22_BX.Vect() );    
    
    TVector3 normal1_BX( tmp1.X()/tmp1.Mag(), tmp1.Y()/tmp1.Mag(), tmp1.Z()/tmp1.Mag() ); 
    TVector3 normal2_BX( tmp2.X()/tmp2.Mag(), tmp2.Y()/tmp2.Mag(), tmp2.Z()/tmp2.Mag() ); 
    
    //// Phi
    TLorentzVector p4Z1_BX = p4M11_BX + p4M12_BX;    
    double tmpSgnPhi = p4Z1_BX.Vect().Dot( normal1_BX.Cross( normal2_BX) );
    double sgnPhi = tmpSgnPhi/fabs(tmpSgnPhi);
    Phi = sgnPhi * acos( -1.*normal1_BX.Dot( normal2_BX) );
    
    
    //////////////
    
    TVector3 beamAxis(0,0,1);
    TVector3 tmp3 = (p4M11_BX + p4M12_BX).Vect();
    
    TVector3 p3V1_BX( tmp3.X()/tmp3.Mag(), tmp3.Y()/tmp3.Mag(), tmp3.Z()/tmp3.Mag() );
    TVector3 tmp4 = beamAxis.Cross( p3V1_BX );
    TVector3 normalSC_BX( tmp4.X()/tmp4.Mag(), tmp4.Y()/tmp4.Mag(), tmp4.Z()/tmp4.Mag() );
    
    //// Phi1
    double tmpSgnPhi1 = p4Z1_BX.Vect().Dot( normal1_BX.Cross( normalSC_BX) );
    double sgnPhi1 = tmpSgnPhi1/fabs(tmpSgnPhi1);    
    Phi1 = sgnPhi1 * acos( normal1_BX.Dot( normalSC_BX) );    
    
    //    std::cout << "extractAngles: " << std::endl;
    //    std::cout << "costhetastar = " << costhetastar << ", costheta1 = " << costheta1 << ", costheta2 = " << costheta2 << std::endl;
    //    std::cout << "Phi = " << Phi << ", Phi1 = " << Phi1 << std::endl;    
    
}

