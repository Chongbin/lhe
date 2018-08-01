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
using namespace std;

void testbucket()
{
    vector <int> IndexList;
    IndexList.push_back(1);
    IndexList.push_back(2);
    IndexList.push_back(4);
    IndexList.push_back(3);
    vector <vector <int> > p_set;
    p_set = bucketAlgo::pSet(IndexList);
    for (int i = 0; i < p_set.size(); ++i)
    {
      for (int j = 0; j < p_set[i].size(); ++j)
      {
        cout << p_set[i][j];
      }
      cout << "\t: ";
      vector <int> cset;
      cset = bucketAlgo::cSet(IndexList, p_set[i]);
      for (int k = 0; k < cset.size(); ++k)
      {
        cout << cset[k];
      }
      cset.clear();
      cout << "\n";
    }
}

double mTop = 173.5; //GeV

int main()
{
  testbucket();
  ifstream inFile("../tt_had_test.lhe");
  string line;
  bool event_flag = false; //switches on when finds an event
  bool event_meta = false; //event block readability switched off to skip the first event block line
  vector <finalstate::particle> evt;

  while (getline(inFile, line)) 
  {
    if (line.find("<event>") != string::npos) 
    {
      event_flag = true;
      //cout << line << "\t" << event_flag << endl;
    }
    else if (line.find("</event>") != string::npos)
    {
      event_flag = false; //switch off the event block
      event_meta = false; //switch off the event block readability
      cout << line << "\t" << event_flag << endl;
      finalstate::event ev1(evt);
      cout << ">> " << ev1.nonbjet.size() << endl;
      for (int i = 0; i < evt.size(); ++i) {
        cout << evt[i].getPID() << "\t" << evt[i].getpX() << endl;
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
	    cout << p.getPID() << "\t status: " << p.getStatus() << "\t" << p.getpX() << endl;
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

  inFile.close();
  TFile::Open("test.root", "RECREATE");
  return 0;
}
