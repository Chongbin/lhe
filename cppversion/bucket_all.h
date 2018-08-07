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

  double mTop = 173.5; //GeV
  double mW = 80.4; //GeV



  //power set generator
  vector <vector <int> > pSet(vector <int> set, int offset=0)
  {
      //power set size : 2^n
      unsigned int pow_set_size = pow(2, set.size());
      int counter, j;
      vector <vector <int> > pVec;
      //Run counter from 000..1 and stop before 111..1 (111..1 excluded)
      for(counter = 1; counter < pow_set_size-offset; counter++)
      {
        vector <int> v;
        for(j = 0; j < set.size(); j++)
        {
            //if jth bit in the counter is 1 print the jth element
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
  };
  


  //class for bucket object
  class bucket
  {
  private:
    string bucket_label;
    double Mbucket, pTbucket, etabucket;
  public:
    vector <finalstate::particle> members;
    vector <finalstate::particle> nonbjets;
    finalstate::particle bjet;

  ~bucket() {} //destructor

  bucket() 
  {
    bucket_label = "tx"; //tx means label unassigned
    Mbucket = 0; //GeV
    pTbucket = 0; //GeV
  }

  bucket(vector <finalstate::particle> nonbjets, finalstate::particle bjet)
  {
    members.push_back(bjet);
    finalstate::particle b;
    b = b + bjet;
    for(int i =0; i < nonbjets.size(); ++i)
    {
      b = b + nonbjets[i];
      members.push_back(nonbjets[i]);
    }
    Mbucket = b.getM();
    pTbucket = b.getPt();
    etabucket = b.getEta();
  } 

  void setBucketLabel(string label)  {bucket_label=label;}

  string getBucketLabel() {return bucket_label;}

  double getBucketMass() {return Mbucket;}

  double getBucketPt() {return pTbucket;}

  double getBucketEta() {return etabucket;}

  vector<int> getPIDlist()
  {
    vector<int> pidlist;
    for (int i=0; i < members.size(); ++i)
    {
      pidlist.push_back(members[i].getPID());
    }
    return pidlist;
  }

  bool twflag() //true : tw ; false : not tw
  {
    bool flag = false; //defalt bucket not tw
    int sizeB = members.size();
    for (int i = 0; i < sizeB; ++i)
    {
      for (int j = 0; j < sizeB; ++j)
      {
        if (j > i)
        {
	  finalstate::particle temp = members[i] + members[j];
          if( abs((temp.getM()/Mbucket) - (mW/mTop)) < 0.15 )
          {
            flag = true; //it's a tw bucket
	    break;
          } 
        }
      }
    }
    return flag;
  }

  double twOptMetric() {return (Mbucket - mTop)*(Mbucket - mTop);} 

  double tminusOptMetric()
  {
    if (Mbucket > 155.0) {return 9999999.0;}
    else {return abs(Mbucket - 145.0);}
  }

  };




  //function to get two top buckets
  vector <bucketAlgo::bucket> doublebucket(finalstate::event ev, double MbucketMax, double MbucketMin, string target_label, double B1weight)
  {
    vector <bucketAlgo::bucket> B; //B1 and B2
    bucketAlgo::bucket B1, B2;
    int nonbjetsize = ev.nonbjet.size();
    //int bjetsize = ev.bjet.size();
    vector <vector <int> > nonbindexset1;
    vector <int> nonbset;
    for (int j = 0; j < nonbjetsize; ++j) {nonbset.push_back(j);}
    nonbindexset1 = pSet(nonbset); //no offset, bjetsize-1); //offset (bjetsize-1) to leave at least one jet for the other b-jet
    double Deltatw = pow(10,10); //arbit large number
    for (int i = 0; i < nonbindexset1.size(); ++i) //looping over all possible buckets
    {
      vector <vector <int> > nonbindexset2;
      vector <int> restjetset = cSet(nonbset, nonbindexset1[i]);
      /*switching of second bucket's freedom for now//nonbindexset2 = pSet(restjetset); // no offset as all the jets can now be used*/
      nonbindexset2.push_back(restjetset);
      vector <finalstate::particle> nonbA;
      for (int k = 0; k < nonbindexset1[i].size(); ++k)
      {
        nonbA.push_back(ev.nonbjet[nonbindexset1[i][k]]);
      }
      for (int i1 = 0; i1 < nonbindexset2.size(); ++i1)
      {
        vector <finalstate::particle> nonbB;
        for (int k1 = 0; k1 < nonbindexset2[i1].size(); ++k1)
        {
          nonbB.push_back(ev.nonbjet[nonbindexset2[i1][k1]]);
        }
	bucketAlgo::bucket Afirst(nonbA, ev.bjet[0]);
	double AfirstDistance = (target_label == "tw") ? Afirst.twOptMetric() : Afirst.tminusOptMetric();
        vector <int> pidlAfirst = Afirst.getPIDlist(); //
        bucketAlgo::bucket Asecond(nonbA, ev.bjet[1]);
	double AsecondDistance = (target_label == "tw") ? Asecond.twOptMetric() : Asecond.tminusOptMetric();
        bucketAlgo::bucket Bfirst(nonbB, ev.bjet[0]);
        vector <int> pidlAsecond = Asecond.getPIDlist(); //
	double BfirstDistance = (target_label == "tw") ? Bfirst.twOptMetric() : Bfirst.tminusOptMetric();
        bucketAlgo::bucket Bsecond(nonbB, ev.bjet[1]);
        vector <int> pidlBfirst = Bfirst.getPIDlist(); //
	double BsecondDistance = (target_label == "tw") ? Bsecond.twOptMetric() : Bsecond.tminusOptMetric();
        vector <int> pidlBsecond = Bsecond.getPIDlist(); //
        
        double del1 = (B1weight*AfirstDistance) + BsecondDistance;
        double del2 = (B1weight*BfirstDistance) + AsecondDistance;
        double del3 = AfirstDistance + (B1weight*BsecondDistance);
        double del4 = BfirstDistance + (B1weight*AsecondDistance);
        /*cout << del1 << "\t[";
        for (int j = 0; j < pidlAfirst.size(); ++j) {cout << pidlAfirst[j] << ", ";}
        cout << "]\t[";
        for (int j = 0; j < pidlBsecond.size(); ++j) {cout << pidlBsecond[j] << ", ";}
        cout << "]" << endl;
        cout << del2 << "\t[";
        for (int j = 0; j < pidlBfirst.size(); ++j) {cout << pidlBfirst[j] << ", ";}
        cout << "]\t[";
        for (int j = 0; j < pidlAsecond.size(); ++j) {cout << pidlAsecond[j] << ", ";}
        cout << "]" << endl;
        cout << del3 << "\t[";
        for (int j = 0; j < pidlBsecond.size(); ++j) {cout << pidlBsecond[j] << ", ";}
        cout << "]\t[";
        for (int j = 0; j < pidlAfirst.size(); ++j) {cout << pidlAfirst[j] << ", ";}
        cout << "]" << endl;
        cout << del4 << "\t[";
        for (int j = 0; j < pidlAsecond.size(); ++j) {cout << pidlAsecond[j] << ", ";}
        cout << "]\t[";
        for (int j = 0; j < pidlBfirst.size(); ++j) {cout << pidlBfirst[j] << ", ";}
        cout << "]" << endl;*/
        if ((del1 < del2) && (del1 < del3) && (del1 < del4))
        {
          if (del1 < Deltatw)
          {
            Deltatw = del1;
            B1 = Afirst;
            B2 = Bsecond;
          }
        }
        else if ((del2 < del1) && (del2 < del3) && (del2 < del4))
        {
          if (del2 < Deltatw)
          {
            Deltatw = del2;
            B1 = Bfirst;
            B2 = Asecond;
          }
        }
        else if ((del3 < del1) && (del3 < del2) && (del3 < del4))
        {
          if (del3 < Deltatw)
          {
            Deltatw = del3;
            B1 = Bsecond;
            B2 = Afirst;
          }
        }
        else 
        {
          if (del4 < Deltatw)
          {
            Deltatw = del4;
            B1 = Asecond;
            B2 = Bfirst;
          }
        }



	//cout << Deltatw << " @@@@@@@@ diff " << endl;
      }
    } // loop over all possible buckets ends
    cout << "del: " << Deltatw << endl;
    B.push_back(B1);
    B.push_back(B2);
    //cout << B.size() << "\t Bucketsize should be 2" << endl;
    for (int i = 0; i < B.size(); ++i)
    {
      string label; //label assignement
      double Bm = B[i].getBucketMass();
      //cout << "bucket mass: " << Bm << " : " << (Bm < MbucketMax) << endl;
      //cout << "bucket mass range: " << MbucketMin << " : " << MbucketMax << endl;
      if ((Bm < MbucketMax) && (Bm > MbucketMin))
      {
        if (target_label == "tw")
        {
          label = (B[i].twflag())?"tw":"t-";
        }
        else {label = "t-";}
      }
      else
      {
        label = "t0";
      }
      B[i].setBucketLabel(label);
    }

    return B;
  };

  //function to get one top bucket
  bucketAlgo::bucket singlebucket(finalstate::event ev, bucketAlgo::bucket twbucket, double MbucketMax, double MbucketMin)
  {
    //vector <finalstate::particle> twnonBjets = twbucket.nonbjets;
    finalstate::particle bucketBjet = (ev.bjet[0] == twbucket.bjet) ? ev.bjet[1] : ev.bjet[0]; //bjet to be used in the bucket
    vector <int> bucketnonbset;
    int Evnonbjetsize = ev.nonbjet.size();
    int twnonbjetsize = twbucket.nonbjets.size();
    for (int i = 0; i < Evnonbjetsize; ++i)
    {
      bool findtw = false;
      for (int j = 0; j < twnonbjetsize; ++j)
      {
        if (ev.nonbjet[i] == twbucket.nonbjets[j]) {findtw = true;} //particle part of tw bucket
      }
      if (findtw = false) {bucketnonbset.push_back(i);}
    }
    bucketAlgo::bucket finalbucket;
    vector <vector <int> > nonbindexset = pSet(bucketnonbset);  // no offset as all the jets can now be used
    double Deltatmin = pow(10,10); //arbit large number
    for (int i1 = 0; i1 < nonbindexset.size(); ++i1) //loop over all possible buckets
    {
      vector <finalstate::particle> nonB;
      for (int j1 = 0; j1 < nonbindexset[i1].size(); ++j1)
      {
        nonB.push_back(ev.nonbjet[nonbindexset[i1][j1]]);
      }
      bucketAlgo::bucket btemp(nonB, bucketBjet);
      double del = btemp.tminusOptMetric();
      if (Deltatmin > del)
      {
        Deltatmin = del;
        finalbucket = btemp;
      }
    }
    //label assignement
    double mfinalbucket = finalbucket.getBucketMass();
    string label = ((mfinalbucket < MbucketMax) && (mfinalbucket  > MbucketMin)) ? "t-" : "t0";
    finalbucket.setBucketLabel(label);
    return finalbucket;
  }

//namespace ends
}
