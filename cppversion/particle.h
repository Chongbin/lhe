#include <iostream>
#include <fstream>
#include <string>
#include <sstream>
#include <vector>
#include <algorithm>
#include <iterator>
using namespace std;



class particle 
{
private:
	  double px;
	  double py;
	  double pz;
	  double E;
	  int pid;
public:

~particle() {} //destructor

particle() 
{
	px = 0;
	py = 0;
	pz = 0;
	E = 0; 
	pid = -9999;
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
}

double getpX() {return px;}

double getpY() {return py;}

double getpZ() {return pz;}

double getE() {return E;}

int getPID() {return pid;}

};


