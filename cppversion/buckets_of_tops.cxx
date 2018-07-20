#include <iostream>
#include <fstream>
#include <string>
#include <sstream>
#include <vector>
#include <algorithm>
#include <iterator>
#include "particle.h"
using namespace std;

int main()
{
  ifstream inFile("../tt_had_test.lhe");
  string line;
  bool event_flag = false; //switches on when finds an event
  bool event_meta = false; //event block readability switched off to skip the first event block line

  while (std::getline(inFile, line)) 
  {
    if (line.find("<event>") != std::string::npos) 
    {
      event_flag = true;
      //cout << line << "\t" << event_flag << endl;
    }
    else if (line.find("</event>") != std::string::npos)
    {
      event_flag = false; //switch off the event block
      event_meta = false; //switch off the event block readability
      cout << line << "\t" << event_flag << endl;
    }
    else{
      if (event_flag) 
      {
        if (event_meta)
        {
  	particle p(line);
	cout << p.getpX() << endl;
        }
        else
        {
          event_meta = true; //to make the rest of the event block readable
        }
      }
    }
  }

  inFile.close();
  return 0;
}
