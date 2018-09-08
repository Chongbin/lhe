#include <iostream>
#include <fstream>
#include <string>
#include <sstream>
#include <vector>
#include <algorithm>
#include <iterator>
#include "particle.h"
#include "bucket_all.h"
#include "TFile.h"
#include "TCanvas.h"
#include "TH1.h"
#include "math.h"

using namespace std;

int main()
{
//  testbucket();
  //ifstream inFile("../tt_had_test_one.lhe");
  //ifstream inFile("../tt_had_test.lhe");
  //ifstream inFile("../tt_hadronic.lhe");
  ifstream inFile("../bbjjj.lhe");
  string line;
  bool event_flag = false; //switches on when finds an event
  bool event_meta = false; //event block readability switched off to skip the first event block line
  vector <finalstate::particle> evt;
  //mass
  TH1F htwmass("htwmass", "Mass of tw Buckets",150,0.0001,300); 
  TH1F htminmass("htminmass", "Mass of t- Buckets",150,0.0001,300); 
  TH1F ht0mass("ht0mass", "Mass of t0 Buckets",150,0.0001,300); 
  TH1F hXmass("hXmass", "Mass of the extra jets",110,-1,10); 
  // pT
  TH1F htwPt("htwPt", "Pt of tw Buckets",250,0,1200); 
  TH1F htminPt("htminPt", "Pt of t- Buckets",250,0,1200); 
  TH1F ht0Pt("ht0Pt", "Pt of t0 Buckets",250,0,1200); 
  TH1F hXPt("hXPt", "Pt of the extra jets",100,0,500); 
  // eta
  TH1F htweta("htweta", "#eta of tw Buckets",100,-10,10); 
  TH1F htmineta("htmineta", "#eta of t- Buckets",100,-10,10); 
  TH1F ht0eta("ht0eta", "#eta of t0 Buckets",100,-10,10);
  TH1F hXeta("hXeta", "#eta of the extra jets",100,-10,10);
  //W candidate mass
  TH1F hmW("hmW", "Mass of the (possible) W candidate",150,0.0001,300); 
  TH1F hmBucketPrim("hmBucketPrimitive", "Mass of the Entire Buckets before Recalculation",150,0,300); 
  TH1F hmratio("hmratio", "Mass Ratio Difference",120,-0.1,1.1); 
  

 
  int eventcounter = 0;
  int twcounter = 0;
  int tmincounter = 0;
  int t0counter = 0;
  int tXcounter = 0;
  while (getline(inFile, line)) 
  {
    if (line.find("<event>") != string::npos) 
    {
      event_flag = true;
      cout << "event: " << eventcounter << endl;
      eventcounter++;
      //cout << line << "\t" << event_flag << endl;
    }
    else if (line.find("</event>") != string::npos)
    {
      event_flag = false; //switch off the event block
      event_meta = false; //switch off the event block readability
      //cout << line << "\t" << event_flag << endl;
      finalstate::event ev1(evt);
      //cout << ">> " << ev1.nonbjet.size() << endl;
/*      for (int i = 0; i < evt.size(); ++i) {
        cout << evt[i].getPID() << "\t" << evt[i].getpX() << endl;
      }*/
      //discard events with less than two b jets
      if (ev1.bjet.size() == 2)
      {
        //bucket algo PASS1 (to find tw buckets)//
        vector <bucketAlgo::bucket> B;
        double bucketP1massMax = 200; //GeV
        double bucketP1massMin = 155; //GeV
        double firstP1Bucketwt = 100;
        B = bucketAlgo::doublebucket(ev1, bucketP1massMax, bucketP1massMin, "tw", firstP1Bucketwt);
        //cout << "B1: " << B[0].getBucketLabel() << "\tB2: " << B[1].getBucketLabel() << endl;
        //pass2 find t- buckets//
        vector <bucketAlgo::bucket> tmincand; //fill tmin candidates
        int telseindex; //tw or t0 bucket
        for (int i = 0; i < B.size(); ++i)
        {
          if (B[i].getBucketMass() > -1) {hmW.Fill(B[i].WcandMnum());}
	  cout << "W_cand_mass: " << B[i].WcandMnum() << endl;
	  cout << "ratio_min: " << B[i].WcandRatio() << endl;
          if (B[i].getBucketMass() > -1) {hmBucketPrim.Fill(B[i].getBucketMass());}
          if (B[i].getBucketMass() > -1) {hmratio.Fill(B[i].WcandRatio());}
          if (B[i].getBucketLabel() == "t-") {tmincand.push_back(B[i]);}
          else 
          {
            telseindex = i;
          } 
        }
        //redo both buckets for t-
        double bucketP2massMax = 155; //GeV
        double bucketP2massMin = 75; //GeV
        double firstP2Bucketwt = 1;
        //if (tmincand.size() != 0) {cout << "event: " << eventcounter << endl;}
        
        if (tmincand.size() == 2)
        {
          B = bucketAlgo::doublebucket(ev1, bucketP2massMax, bucketP2massMin, "t-", firstP2Bucketwt);
          //**//cout << "event: " << eventcounter << endl;
        }
        else if (tmincand.size() == 1)
        {
          B[1-telseindex] = bucketAlgo::singlebucket(ev1, B[telseindex], bucketP2massMax, bucketP2massMin);
        }
        //pT cut of 200 GeV on the final buckets while filling histograms//
        
        // fill histograms from the two buckets here
        for (int i = 0; i < B.size(); ++i)
        {
          /*cout << "label: " << B[i].getBucketLabel() << "\tmass: " << B[i].getBucketMass() << "\t[";
          vector<int> pidlist =  B[i].getPIDlist();
       
          for (int j = 0; j < pidlist.size(); ++j) {cout << pidlist[j] << ", ";}
          cout << "]" << endl;*/
          if (B[i].getBucketLabel() == "tw") 
          {
              if (B[i].getBucketMass() > -1) {
        	  htwmass.Fill(B[i].getBucketMass());
        	  htwPt.Fill(B[i].getBucketPt());
        	  htweta.Fill(B[i].getBucketEta()); }
        	  ++twcounter;
          }
          else if (B[i].getBucketLabel() == "t-")
          {       
              if (B[i].getBucketMass() > -1) {
      	          htminmass.Fill(B[i].getBucketMass());
        	  htminPt.Fill(B[i].getBucketPt());
        	  htmineta.Fill(B[i].getBucketEta());}
                  ++tmincounter;
          }
          else 
          {
		  //**//cout << "LL: " << B[i].getBucketLabel() << "\tmass: " << B[i].getBucketMass() << endl;
        	  if (B[i].getBucketMass() > -1) {
		  ht0mass.Fill(B[i].getBucketMass());
        	  ht0Pt.Fill(B[i].getBucketPt());
        	  ht0eta.Fill(B[i].getBucketEta());}
		  if (B[i].getBucketEta() == 0) 
		  {
		          vector<int> plll = B[i].getPIDlist();
			  cout << "bucket " << i << "<--{ " << B[i].getBucketLabel() << " : " << plll.size();
			  for (int lll = 0; lll < plll.size(); ++lll)
		          {
				  cout << plll[lll] << " , ";
			  }
       	                  //for (vector<int>::const_iterator l = plll.begin(); l != plll.end(); ++l)
	                    //  cout <<  *l << ", ";
			  cout << endl;
			  cout << "\tt0ETA: " << B[i].getBucketEta() << "\t event: " << eventcounter << endl;
	          }
                  ++t0counter;
          }
	//*^*//cout << "\t|\t";
	//*^*//cout << "bucket " << i << "<--{ " << B[i].getBucketLabel() << " : ";
	//*^*//vector<int> plll = B[i].getPIDlist();
       	//*^*//for (vector<int>::const_iterator l = plll.begin(); l != plll.end(); ++l)
	       //*^*//cout <<  *l << ", ";
        }
	//*^*//cout << endl;

          vector <finalstate::particle> Xtra = bucketAlgo::extra(ev1.EVT, B); // extra bucket
          cout << Xtra.size() << endl;
	  

          for (int nn=0; nn < Xtra.size(); ++nn)
          {
            ++tXcounter;
	    cout << "EXTRA " << Xtra[nn].getpX() << "\t" << Xtra[nn].getpY() << "\t" << Xtra[nn].getpZ() << "\t" << Xtra[nn].getE() << "\t" << Xtra[nn].getPID() << endl;
	    cout << "EXTRA MASS " << Xtra[nn].getM() << endl;
	    hXmass.Fill(Xtra[nn].getM());
            hXPt.Fill(Xtra[nn].getPt());
            hXeta.Fill(Xtra[nn].getEta());
          }
          
        //
      }
      evt.clear(); //insert event operations before clearing the vector 
    }
    else{
      if (event_flag) 
      {
        if (event_meta)
        {
          finalstate::particle p(line);
          if (p.getStatus() == 1) 
          {
            //cout << p.getPID() << "\t status: " << p.getStatus() << "\t" << p.getpX() << endl;
            evt.push_back(p);
          }
        }
        else
        {
         event_meta = true; //to make the rest of the event block readable
        }
      }
    }
  }
  htwmass.Draw();
  htminmass.Draw();
  ht0mass.Draw();
  cout << "tw: " << twcounter << "\tt-: " << tmincounter << "\tt0: " << t0counter << "\ttX: "<< tXcounter << endl;
  TFile f("test.root", "RECREATE");
  TCanvas c ("c", "c", 800, 600);
  //mass
  htwmass.GetXaxis()->SetTitle("Mass (GeV)");
  htwmass.Write();
  htwmass.Draw();
  c.Update();
  c.Print("cpp_reconstructed_top_mass_tw.eps");
  c.Clear();
  htminmass.GetXaxis()->SetTitle("Mass (GeV)");
  htminmass.Write();
  htminmass.Draw();
  c.Update();
  c.Print("cpp_reconstructed_top_mass_t_.eps");
  c.Clear();
  ht0mass.GetXaxis()->SetTitle("Mass (GeV)");
  ht0mass.Write();
  ht0mass.Draw();
  c.Update();
  c.Print("cpp_reconstructed_top_mass_t0.eps");
  c.Clear();
  hXmass.GetXaxis()->SetTitle("Mass (GeV)");
  hXmass.Write();
  hXmass.Draw();
  c.Update();
  c.Print("cpp_reconstructed_top_mass_x.eps");
  c.Clear();

  //Pt
  htwPt.GetXaxis()->SetTitle("Pt (GeV)");
  htwPt.Write();
  htwPt.Draw();
  c.Update();
  c.Print("cpp_reconstructed_top_pt_tw.eps");
  c.Clear();
  htminPt.GetXaxis()->SetTitle("Pt (GeV)");
  htminPt.Write();
  htminPt.Draw();
  c.Update();
  c.Print("cpp_reconstructed_top_pt_t_.eps");
  c.Clear();
  ht0Pt.GetXaxis()->SetTitle("Pt (GeV)");
  ht0Pt.Write();
  ht0Pt.Draw();
  c.Update();
  c.Print("cpp_reconstructed_top_pt_t0.eps");
  c.Clear();
  hXPt.GetXaxis()->SetTitle("Pt (GeV)");
  hXPt.Write();
  hXPt.Draw();
  c.Update();
  c.Print("cpp_reconstructed_top_pt_x.eps");
  c.Clear();

  //Eta
  htweta.GetXaxis()->SetTitle("#eta");
  htweta.Write();
  htweta.Draw();
  c.Update();
  c.Print("cpp_reconstructed_top_eta_tw.eps");
  c.Clear();
  htmineta.GetXaxis()->SetTitle("#eta");
  htmineta.Write();
  htmineta.Draw();
  c.Update();
  c.Print("cpp_reconstructed_top_eta_t_.eps");
  c.Clear();
  ht0eta.GetXaxis()->SetTitle("#eta");
  ht0eta.Write();
  ht0eta.Draw();
  c.Update();
  c.Print("cpp_reconstructed_top_eta_t0.eps");
  c.Clear();
  hXeta.GetXaxis()->SetTitle("#eta");
  hXeta.Write();
  hXeta.Draw();
  c.Update();
  c.Print("cpp_reconstructed_top_eta_x.eps");
  c.Clear();

  //W candidate
  hmW.GetXaxis()->SetTitle("Mass (GeV)");
  hmW.Write();
  hmW.Draw();
  c.Update();
  c.Print("cpp_reconstructed_top_m_jk.eps");
  c.Clear();
  hmBucketPrim.GetXaxis()->SetTitle("Mass (GeV)");
  hmBucketPrim.Write();
  hmBucketPrim.Draw();
  c.Update();
  c.Print("cpp_reconstructed_top_m_b.eps");
  c.Clear();
  hmratio.GetXaxis()->SetTitle("Mass Ratio");
  hmratio.Write();
  hmratio.Draw();
  c.Update();
  c.Print("cpp_reconstructed_top_ratio.eps");
  c.Clear();

  f.Close();

  inFile.close();
  return 0;
}
