#include <iostream>
#include <fstream>
#include <string>
#include <sstream>
#include <vector>
#include <algorithm>
#include <iterator>
#include  <cstdlib>
#include "math.h"

using namespace std;

//put event stuff inside finalstate namespace
namespace finalstate
{
  class particle 
  {
  private:
    double px;
    double py;
    double pz;
    double E;
    int pid;
    int status;
  public:
  
  ~particle() {} //destructor
  
  particle() 
  {
    px = 0;
    py = 0;
    pz = 0;
    E = 0; 
    pid = -9999;
    status = 0;
  }
  
  particle(string line) 
  {
      istringstream iss(line);
      vector<string> pinfo{istream_iterator<string>{iss},
        istream_iterator<string>{}};
      string pxstr(pinfo[6]);
      px = stof(pxstr);
      string pystr(pinfo[7]);
      py = stof(pystr);
      string pzstr(pinfo[8]);
      pz = stof(pzstr);
      string Estr(pinfo[9]);
      E = stof(Estr);
      string pidstr(pinfo[0]);
      pid = stoi(pidstr);
      string statusstr(pinfo[1]);
      status = stoi(statusstr);
  }
  
  double getpX() {return px;}

  double getpY() {return py;}

  double getpZ() {return pz;}

  double getE() {return E;}

  int getPID() {return pid;}

  int getStatus() {return status;}
  
  void setpX(double p1) {px = p1;}

  void setpY(double p2) {py = p2;}

  void setpZ(double p3) {pz = p3;}

  void setE(double p4) {E = p4;}
  
  double getM() {return sqrt(E*E - px*px - py*py - pz*pz);}

  double getPt() {return sqrt(px*py + py*py);}

  double getEta() {return atanh(pz/(sqrt(px*px + py*py + pz*pz)));}

  // overload + operator to add to particle objects
  particle operator+(particle p)
  {
    particle combined;
    combined.setpX(px + p.getpX());
    combined.setpY(py + p.getpY());
    combined.setpZ(pz + p.getpZ());
    combined.setE(E + p.getE());
    return combined;
  }

  // overload == operator to compare particle objects
  bool operator==(const particle& p2)
  {
    return ((px == p2.px) && (py == p2.py) && (pz == p2.pz) && (E == p2.E) && (pid == p2.pid) && (status == p2.status));
  }  

  // overload = operator to compare particle objects
  bool operator=(const particle& p2)
  {
    return ((px = p2.px) && (py = p2.py) && (pz = p2.pz) && (E = p2.E) && (pid = p2.pid) && (status = p2.status));
  }  

  };
  
  class event
  {
  private:
    
  public:
    vector <particle> bjet;
    vector <particle> nonbjet;
  
  ~event() {} //destructor
  
  event() {}
  
  event(vector <particle> evt) 
  {
    for (int i = 0; i < evt.size(); ++i) {
      int pid = evt[i].getPID();
      //cout << pid << "\t" << evt[i].getpX() << endl;
      if (abs(pid) == 5)
      {
        bjet.push_back(evt[i]);
      }
      else
      {
        nonbjet.push_back(evt[i]);
      }
    }
  }
  
  };
//namespace ends
}
