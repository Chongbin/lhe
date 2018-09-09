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

  vector <finalstate::particle> compvec(vector <finalstate::particle> EVT, vector <finalstate::particle> set) //another bug
  {
    vector <finalstate::particle> compset;
    for (int i = 0; i < EVT.size(); ++i)
    {
      bool findflag = false;
      for (int j = 0; j < set.size(); ++j)
      {
        findflag = findflag || (EVT[i]==set[j]); 
	//cout << "findflag: " << findflag  << "\tevtorder: " << EVT[i].getOrder() << "\tbucketpartOrder: " << set[j].getOrder() << endl;
      }
      if (!findflag) {compset.push_back(EVT[i]);}
    }
    return compset;
  };

  //power set generator
  vector <vector <int> > pSet(vector <int> set, int offset=0)
  {
      //power set size : 2^n
      unsigned int pow_set_size = pow(2, set.size());
      int counter;
      vector <vector <int> > pVec;
      //Run counter from 000..1 and stop before 111..1 (111..1 excluded)
      for(counter = 1; counter < pow_set_size-offset; counter++)
      {
        vector <int> v;
        for(int j = 0; j < set.size(); j++)
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
    double mpairnum, ratio_min;
  public:
    vector <finalstate::particle> members;
    vector <finalstate::particle> nonBJETS;
    finalstate::particle BJET;

  ~bucket() {} //destructor

  bucket()
  {
    bucket_label = "tx"; //tx means label unassigned
    Mbucket = -9999; //GeV
    pTbucket = -9999; //GeV
    etabucket = -9999;
    mpairnum = -9999; //GeV
  }

  bucket(vector <finalstate::particle> nonbjets, finalstate::particle bjet)
  {
    nonBJETS = nonbjets;
    BJET = bjet;
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

  vector<int> getOrderlist()
  {
    vector<int> orderlist;
    for (int i=0; i < members.size(); ++i)
    {
      orderlist.push_back(members[i].getOrder());
    }
    return orderlist;
  }

  bool twflag() //true : tw ; false : not tw
  {
    ratio_min = pow(10, 10); 
    bool flag = false; //defalt bucket not tw
    int sizeB = members.size();


    for (int i = 0; i < sizeB; ++i)
    {
      for (int j = 0; j < sizeB; ++j)
      {
        if (j < i) {continue;}
	//cout << "(i,j): " << i << "\t" << j << endl;
	//^//cout << "particle1: " << members[i].getOrder() << "\tpx: " << members[i].getpX() << "\tpy: " << members[i].getpY() << "\tpz: " << members[i].getpZ() << "\tE: " << members[i].getE() << "\tpid: " << members[i].getPID() << endl;
	//^//cout << "particle2: " << members[j].getOrder() << "\tpx: " << members[j].getpX() << "\tpy: " << members[j].getpY() << "\tpz: " << members[j].getpZ() << "\tE: " << members[j].getE() << "\tpid: " << members[j].getPID() << endl;
	finalstate::particle temp = members[i]+members[j];
        double dd = abs((temp.getM()/Mbucket) - (mW/mTop));
	//^//cout << "m_ij: " << temp.getM() << "\tmbucket: " << Mbucket << "\ttransient_mass_ratio: " << ratio_min << endl;
	//cout << "dd: " << dd << "\tratmin: " << ratio_min << endl;
        if (ratio_min > dd) //
        {
          ratio_min = dd;
          mpairnum = temp.getM();
	  //^//cout << "temp_mass_W: " << mpairnum << "\tupdating_mass_ratio: " << ratio_min << endl;
        }
      }
    }
    //^//cout << "final_mass_ratio: " << ratio_min << endl;
    if ( (ratio_min) < 0.15 ) {flag = true;}
    return flag;
  }

  double WcandMnum() {return mpairnum;}
  
  double WcandRatio() {return ratio_min;}

  double twOptMetric() {return (Mbucket - mTop)*(Mbucket - mTop);} 

  double tminusOptMetric()
  {
    if (Mbucket > 155.0) {return pow(10,10);}
    else {return abs(Mbucket - 145.0);}
  }

  };


  //function for the extra jets
  vector <finalstate::particle> extra(vector <finalstate::particle> EVT, vector <bucketAlgo::bucket> B)
  {
    vector <finalstate::particle> Xmembers = EVT;
    for (int i = 0; i < B.size(); ++i)
    {
      Xmembers = bucketAlgo::compvec(Xmembers, B[i].members);
      //cout << "iter" << i << ": ";
      //for (int j = 0; j < B[i].members.size(); ++j) { cout <<  B[i].members[j].getOrder() << " ,"; }
      //cout << endl;
    }
    return Xmembers;
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
    if (target_label == "tw")
    {
      nonbindexset1 = pSet(nonbset); //no offset for tw
    }
    else //for t- pair search 
    {
      for (int l1 = 0; l1 < nonbjetsize; ++l1)
      {
        vector <int> tempPar1;
        tempPar1.push_back(nonbset[l1]);
        nonbindexset1.push_back(tempPar1);
        tempPar1.clear();
      }
    }
    double Deltatw = pow(10,10); //arbit large number
    for (int i = 0; i < nonbindexset1.size(); ++i) //looping over all possible buckets
    {
      vector <vector <int> > nonbindexset2;
      vector <int> restjetset = cSet(nonbset, nonbindexset1[i]);
      /*switching of second bucket's freedom for now//nonbindexset2 = pSet(restjetset); // no offset as all the jets can now be used*/
      if (target_label == "tw") 
      {
        nonbindexset2.push_back(restjetset);
      }
      else //for t- pair search
      {
        for (int l2 = 0; l2 < restjetset.size(); ++l2)
        {
          vector <int> tempPar2;
          tempPar2.push_back(restjetset[l2]);
          nonbindexset2.push_back(tempPar2);
          tempPar2.clear();
        //**//cout << "--->" << nonbindexset1[i][0] << "\t" << nonbindexset2[l2][0] << endl;
        }
      }
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
        vector <int> pidlAfirst = Afirst.getOrderlist(); //
        bucketAlgo::bucket Asecond(nonbA, ev.bjet[1]);
	double AsecondDistance = (target_label == "tw") ? Asecond.twOptMetric() : Asecond.tminusOptMetric();
        vector <int> pidlAsecond = Asecond.getOrderlist(); //
        bucketAlgo::bucket Bfirst(nonbB, ev.bjet[0]);
	double BfirstDistance = (target_label == "tw") ? Bfirst.twOptMetric() : Bfirst.tminusOptMetric();
        vector <int> pidlBfirst = Bfirst.getOrderlist(); //
        bucketAlgo::bucket Bsecond(nonbB, ev.bjet[1]);
	double BsecondDistance = (target_label == "tw") ? Bsecond.twOptMetric() : Bsecond.tminusOptMetric();
        vector <int> pidlBsecond = Bsecond.getOrderlist(); //
        
        double del1 = (B1weight*AfirstDistance) + BsecondDistance;
        double del2 = (B1weight*BfirstDistance) + AsecondDistance;
        double del3 = AfirstDistance + (B1weight*BsecondDistance);
        double del4 = BfirstDistance + (B1weight*AsecondDistance);
        /*cout << "del1: " << del1 << "\t[";
        for (int j = 0; j < pidlAfirst.size(); ++j) {cout << pidlAfirst[j] << ", ";}
        cout << "]\t[";
        for (int j = 0; j < pidlBsecond.size(); ++j) {cout << pidlBsecond[j] << ", ";}
        cout << "]" << endl;
        cout << "del2: " << del2 << "\t[";
        for (int j = 0; j < pidlBfirst.size(); ++j) {cout << pidlBfirst[j] << ", ";}
        cout << "]\t[";
        for (int j = 0; j < pidlAsecond.size(); ++j) {cout << pidlAsecond[j] << ", ";}
        cout << "]" << endl;
        cout << "del3: " << del3 << "\t[";
        for (int j = 0; j < pidlBsecond.size(); ++j) {cout << pidlBsecond[j] << ", ";}
        cout << "]\t[";
        for (int j = 0; j < pidlAfirst.size(); ++j) {cout << pidlAfirst[j] << ", ";}
        cout << "]" << endl;
        cout << "del4: " << del4 << "\t[";
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



	//if (target_label == "t-") {cout << Deltatw << " @@@@@@@@ difftmp " << B1.members.size() << "\t" << B2.members.size() << "\tdel1: " <<  del1 << "\tdel2" << del2 << "\tdel3" << del3 << endl;}
      }
    } // loop over all possible buckets ends
    //if (target_label == "t-") {cout << "del: " << Deltatw << endl;}
    //cout << B.size() << "\t Bucketsize should be 2" << endl;
    B.push_back(B1);
    B.push_back(B2);
    //^//if (target_label == "tw") 
    //^//{
      //^//cout << "del: " << Deltatw << endl;
//      cout << "B1: [";
//      vector<int> b1pid = B[0].getPIDlist();
//      for (int j = 0; j < b1pid.size(); ++j) {cout << b1pid[j] << ",\t";}
//      cout << " ] " << endl;
//      cout << "B2: [";
//      vector<int> b2pid = B[1].getPIDlist();
//      for (int j = 0; j < b2pid.size(); ++j) {cout << b2pid[j] << ",\t";}
//      cout << " ] " << endl;

      //^//cout << "B1: [";
      //^//vector<int> b1pid = B[0].getOrderlist();
      //^//for (int j = 0; j < b1pid.size(); ++j) {cout << b1pid[j] << ",\t";}
      //^//cout << " ] " << endl;
      //^//cout << "B2: [";
      //^//vector<int> b2pid = B[1].getOrderlist();
      //^//for (int j = 0; j < b2pid.size(); ++j) {cout << b2pid[j] << ",\t";}
      //^//cout << " ] " << endl;
    //^//}
    for (int i = 0; i < B.size(); ++i)
    {
      string label; //label assignement , 
      double Bm = B[i].getBucketMass();
      //cout << "bucket mass: " << Bm << " : " << (Bm < MbucketMax) << endl;
      //cout << "bucket mass range: " << MbucketMin << " : " << MbucketMax << endl;
      if (target_label == "tw")
      {
        label = (B[i].twflag())?"tw":"t-";
      }
      else {label = "t-";}
      if ((Bm > MbucketMax) || (Bm < MbucketMin))
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

    bucketAlgo::bucket newbucket; 
    double Deltamin = pow(10,10); //arbit large number

    finalstate::particle bucketBjet = (ev.bjet[0] == twbucket.BJET) ? ev.bjet[1] : ev.bjet[0]; //bjet to be used in the bucket

    int Evnonbjetsize = ev.nonbjet.size();
    int twnonbjetsize = twbucket.nonBJETS.size();
    vector <int> plist = twbucket.getPIDlist();
    //**//cout << "[ ";
    //**//for (vector<int>::const_iterator l = plist.begin(); l != plist.end(); ++l)
	        //**//cout << *l << ", ";
    //**//cout << " ]" << endl;
    //cout << "***" << twbucket.nonBJETS.size() << endl;
    for (int i = 0; i < Evnonbjetsize; ++i)
    {
      bool findtw = false;
      for (int j = 0; j < twnonbjetsize; ++j)
      {
        bool findtemp = (ev.nonbjet[i]==twbucket.nonBJETS[j]); //particle part of tw bucket; can be added to extra
        findtw = findtw || findtemp;
        //cout << "found a match: " << findtw << endl;
	  //**//cout << "+ ev| px: " << ev.nonbjet[i].getpX() << "\tpy: " << ev.nonbjet[i].getpY() << "\tpz: " << ev.nonbjet[i].getpZ() << "\tE " << ev.nonbjet[i].getE() << "\tpid: " <<  ev.nonbjet[i].getPID() << "\tstatus: " << ev.nonbjet[i].getStatus() << endl;
  	  //**//cout << "+ twB| px: " << twbucket.nonBJETS[j].getpX() << "\tpy: " << twbucket.nonBJETS[j].getpY() << "\tpz: " << twbucket.nonBJETS[j].getpZ() << "\tE " << twbucket.nonBJETS[j].getE() << "\tpid: " <<  twbucket.nonBJETS[j].getPID() << "\tstatus: " << twbucket.nonBJETS[j].getStatus() << "\tdecision:" << findtemp << endl;
      }
      //**//cout << "!!!!decision: " << findtw << endl;
      if (!findtw)  //no match
      {
        vector <finalstate::particle> tempnb;
        tempnb.push_back(ev.nonbjet[i]);
        bucketAlgo::bucket tempbucket(tempnb, bucketBjet);
        if (Deltamin > tempbucket.tminusOptMetric())
        {
          Deltamin = tempbucket.tminusOptMetric();
          newbucket = tempbucket;
        }
	//**//cout << "~~~~~~~~~~~~~" << Deltamin << "\t" << ev.nonbjet[i].getPID() << endl;
      }
    }
    
    //label assignement
    double massbucket = newbucket.getBucketMass();
    string label = ((massbucket < MbucketMax) && (massbucket  > MbucketMin)) ? "t-" : "t0";
    newbucket.setBucketLabel(label);
    //**//cout << "M_ORI: " << newbucket.getBucketMass() << "\t" << newbucket.BJET.getM() << endl;
    return newbucket;
  };



//namespace ends
}
