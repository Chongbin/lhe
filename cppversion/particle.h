#include <iostream>
#include <fstream>
#include <string>
#include <sstream>
#include <vector>
#include <algorithm>
#include <iterator>
#include  <cstdlib>
using namespace std;


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
      cout << pid << "\t" << evt[i].getpX() << endl;
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
}
