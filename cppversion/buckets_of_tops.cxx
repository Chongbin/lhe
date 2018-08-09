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
#include "TH1.h"
#include "math.h"

using namespace std;

int main()
{
//  testbucket();
  //ifstream inFile("../tt_had_test_one.lhe");
  //ifstream inFile("../tt_had_test.lhe");
  ifstream inFile("../tt_hadronic.lhe");
  string line;
  bool event_flag = false; //switches on when finds an event
  bool event_meta = false; //event block readability switched off to skip the first event block line
  vector <finalstate::particle> evt;

  TH1F htwmass("htwmass", "Mass of tw Buckets",150,0,300); 
  int eventcounter = 0;
  int twcounter = 0;
  int tmincounter = 0;
  int t0counter = 0;
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
          cout << "event: " << eventcounter << endl;
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
        	  htwmass.Fill(B[i].getBucketMass());
        	  ++twcounter;
          }
          else if (B[i].getBucketLabel() == "t-"){++tmincounter;}
          else {++t0counter;}
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
  cout << "tw: " << twcounter << "\tt-: " << tmincounter << "\tt0: " << t0counter << endl;
  TFile f("test.root", "RECREATE");
  htwmass.Write();
  f.Close();


  inFile.close();
  return 0;
}
