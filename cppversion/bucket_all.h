#include <stdio.h>
#include <math.h>
#include <string>
using namespace std;

//Forward declaration
namespace finalstate
{
  class particle;
  class finalstate;
}


//put top bucket algo stuff inside bucketAlgo namespace
namespace bucketAlgo
{
  //power set generator
  vector <vector <int> > pSet(vector <int> set, int offset=0)
  {
      /*set_size of power set of a set with set_size
        n is (2**n -1)*/
      unsigned int pow_set_size = pow(2, set.size());
      int counter, j;
      vector <vector <int> > pVec;
      /*Run from counter 000..1 until 111..1 (111..1 excluded)*/
      for(counter = 1; counter < pow_set_size-offset; counter++)
      {
        vector <int> v;
        for(j = 0; j < set.size(); j++)
        {
            /* Check if jth bit in the counter is set
               If set then pront jth element from set */
            if(counter & (1<<j))
              v.push_back(set[j]);
        }
        pVec.push_back(v); 
        v.clear();
      }
      return pVec;
  }
   
  //complementary set generator
  vector <int> cSet(vector <int> set, vector <int> subset)
  {
    vector <int> cVec;
    for (int i = 0; i < set.size(); ++i)
    {
       if (find(subset.begin(), subset.end(), set[i]) != subset.end())
       {
         continue;
       }
       //cout << set[i] <<endl;
       cVec.push_back(set[i]);
    }
    return cVec;
  }
  
  
  //class for bucket object
  class bucket
  {
  private:
    string bucket_label;
  public:

  ~bucket() {} //destructor

  bucket() 
  {
    bucket_label = "tx"; //tx means label unassigned
  }

  bucket(vector <finalstate::particle> nontopjets, finalstate::particle topjet)
  {
    return;
  } 

  string getBucketLabel() {return bucket_label;}

  };
}
