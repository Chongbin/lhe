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

  int indexfinder(vector <finalstate::particle> EVT, finalstate::particle p) //(just to match the erroneous python version) (delete all of this in the final version) it is just to mimic a bug
  {
    int index = -99999;
    for (int i = 0; i < EVT.size(); ++i)
    {
      if (EVT[i]==p) {index = i;}
    }
    return index;
  };

  vector <finalstate::particle> compvec(vector <finalstate::particle> EVT, vector <finalstate::particle> set) //another bug
  {
    vector <finalstate::particle> compset;
    for (int i = 0; i < EVT.size(); ++i)
    {
      bool findflag = false;
      for (int j = 0; j < set.size(); ++j)
      {
        findflag = findflag || (EVT[i]==set[j]); 
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
    Mbucket = 0; //GeV
    pTbucket = 0; //GeV
    etabucket = 0;
    mpairnum = -1; //GeV
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

  bool twflag() //true : tw ; false : not tw
  {
    this->ratio_min = pow(10, 10); 
    bool flag = false; //defalt bucket not tw
    int sizeB = members.size();

    bool mflag = false; //I do labelling later (must be deleted later)
    if ((Mbucket < 200) && (Mbucket > 150)) {mflag = true;} // must be deleted


    for (int i = 0; i < sizeB; ++i)
    {
      for (int j = 0; j < sizeB; ++j)
      {
        if (j > i)
        {
          finalstate::particle temp = members[i] + members[j];
          double dd = abs((temp.getM()/Mbucket) - (mW/mTop));
          if (this->ratio_min > dd) 
          {
            this->ratio_min = dd;
            this->mpairnum = temp.getM();
          }
          cout << "----------L: " << this->mpairnum << endl;
          if( abs((temp.getM()/Mbucket) - (mW/mTop)) < 0.15 )
          {
            flag = flag || true; //it's a tw bucket
            if (mflag) {break;} //delete this line (not really the algo)
            //break; //uncomment this later
          } 
        }
      }
      if (flag && mflag) {break;} // modify this too (make if just flag and delete mflag)
    }
    cout << "----------R: " << this->mpairnum << endl;
    return flag;
  }

  double WcandMnum() {return this->mpairnum;}
  
  double WcandRatio() {return this->ratio_min;}

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
        cout << "--->" << nonbindexset1[i][0] << "\t" << nonbindexset2[l2][0] << endl;
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
    if (target_label == "t-") {cout << "del: " << Deltatw << endl;}
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


    bucketAlgo::bucket newbucket; 
    double Deltamin = pow(10,10); //arbit large number

    finalstate::particle bucketBjet = (ev.bjet[0] == twbucket.BJET) ? ev.bjet[1] : ev.bjet[0]; //bjet to be used in the bucket
    int bindex = indexfinder(ev.EVT, bucketBjet); //mimicing a bug

    int Evnonbjetsize = ev.nonbjet.size();
    int twnonbjetsize = twbucket.nonBJETS.size();
    vector <int> plist = twbucket.getPIDlist();
    cout << "[ ";
    for (vector<int>::const_iterator l = plist.begin(); l != plist.end(); ++l)
	        cout << *l << ", ";
    cout << " ]" << endl;
    //cout << "***" << twbucket.nonBJETS.size() << endl;
    vector <finalstate::particle> nonbalt = compvec(ev.nonbjet, twbucket.nonBJETS); //another bug included! 
    for (int i = 0; i < Evnonbjetsize; ++i)
    {
      bool findtw = false;
      if (indexfinder(ev.EVT, ev.nonbjet[i]) < bindex) {continue;} //this is just to mimic a bug
      for (int j = 0; j < twnonbjetsize; ++j)
      {
        bool findtemp = (ev.nonbjet[i]==twbucket.nonBJETS[j]); //particle part of tw bucket; can be added to extra
        findtw = findtw || findtemp;
        //cout << "found a match: " << findtw << endl;
	  cout << "+ ev| px: " << ev.nonbjet[i].getpX() << "\tpy: " << ev.nonbjet[i].getpY() << "\tpz: " << ev.nonbjet[i].getpZ() << "\tE " << ev.nonbjet[i].getE() << "\tpid: " <<  ev.nonbjet[i].getPID() << "\tstatus: " << ev.nonbjet[i].getStatus() << endl;
  	  cout << "+ twB| px: " << twbucket.nonBJETS[j].getpX() << "\tpy: " << twbucket.nonBJETS[j].getpY() << "\tpz: " << twbucket.nonBJETS[j].getpZ() << "\tE " << twbucket.nonBJETS[j].getE() << "\tpid: " <<  twbucket.nonBJETS[j].getPID() << "\tstatus: " << twbucket.nonBJETS[j].getStatus() << "\tdecision:" << findtemp << endl;
      }
      cout << "!!!!decision: " << findtw << endl;
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
	cout << "~~~~~~~~~~~~~" << Deltamin << "\t" << ev.nonbjet[i].getPID() << endl;
      }
    }
    
    if (newbucket.members.size() == 0) {bucketAlgo::bucket tempbugb(nonbalt, bucketBjet); newbucket = tempbugb;} //bug please delete later "(
    int pV = (newbucket.nonBJETS.size()!=0) ? newbucket.nonBJETS[0].getPID() : -9999;
    cout << "```````````````````" << Deltamin << "\t" << pV << endl;
    //label assignement
    double massbucket = newbucket.getBucketMass();
    string label = ((massbucket < MbucketMax) && (massbucket  > MbucketMin)) ? "t-" : "t0";
    newbucket.setBucketLabel(label);
    cout << "M_ORI: " << newbucket.getBucketMass() << "\t" << newbucket.BJET.getM() << endl;
    return newbucket;
  };



//namespace ends
}
