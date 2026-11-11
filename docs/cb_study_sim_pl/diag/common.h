// independent reimplementation (verifier); no code taken from plots/
#include <TFile.h>
#include <TNtuple.h>
#include <TBranch.h>
#include <cmath>
#include <cstdio>
#include <cstring>
#include <vector>
#include <string>
#include <map>
#include <unordered_map>
#include <set>
#include <algorithm>
#include <cstdint>
static const char* FN="/hdd3/przygoda/HADES/files/sim/smash_std_sim_epem.root";
static const double ME=0.51099895; // MeV
struct Leg{ float p,th,ph,id,simid,vx,vy,vz,par,proc; };
struct Row{ float ev,vrx,vry,vrz,isBest,trig; Leg a,b; };
struct Reader{
  TNtuple*t=nullptr; Row r; Long64_t n=0;
  void chk(const std::string&s){ if(!t->GetBranch(s.c_str())){printf("MISSING branch %s in %s\n",s.c_str(),t->GetName());exit(1);} }
  void bind(const std::string&s,float*p){ chk(s); t->SetBranchStatus(s.c_str(),1); t->SetBranchAddress(s.c_str(),p); }
  void bindLeg(const std::string&pre,Leg&l){
    bind(pre+"p",&l.p);bind(pre+"theta",&l.th);bind(pre+"phi",&l.ph);bind(pre+"id",&l.id);
    bind(pre+"sim_id",&l.simid);bind(pre+"sim_vertexx",&l.vx);bind(pre+"sim_vertexy",&l.vy);bind(pre+"sim_vertexz",&l.vz);
    bind(pre+"sim_parentid",&l.par);bind(pre+"sim_processid",&l.proc);
  }
  void open(TFile*f,const char*name,const char*p1,const char*p2){
    t=(TNtuple*)f->Get(name); if(!t){printf("no %s\n",name);exit(1);} n=t->GetEntries();
    t->SetBranchStatus("*",0);
    bind("event",&r.ev);bind("eVertReco_x",&r.vrx);bind("eVertReco_y",&r.vry);bind("eVertReco_z",&r.vrz);
    bind("isBest",&r.isBest);bind("trigbit",&r.trig);
    bindLeg(p1,r.a);bindLeg(p2,r.b);
  }
  inline void get(Long64_t i){ t->GetEntry(i); }
};
inline bool bitEq(float a,float b){ return std::memcmp(&a,&b,sizeof(float))==0; }
inline void kin(const Leg&l,double&nx,double&ny,double&nz){
  double t=l.th*M_PI/180., f=l.ph*M_PI/180.;
  nx=std::sin(t)*std::cos(f); ny=std::sin(t)*std::sin(f); nz=std::cos(t);
}
// mass in GeV, cos of opening angle
inline void massCos(const Leg&a,const Leg&b,double&M,double&cs){
  double n1[3],n2[3]; kin(a,n1[0],n1[1],n1[2]); kin(b,n2[0],n2[1],n2[2]);
  cs=n1[0]*n2[0]+n1[1]*n2[1]+n1[2]*n2[2];
  double p1=a.p,p2=b.p; double E1=std::sqrt(p1*p1+ME*ME),E2=std::sqrt(p2*p2+ME*ME);
  double M2=2*ME*ME+2*(E1*E2-p1*p2*cs);
  M=(M2>0?std::sqrt(M2):0.)/1000.;
}
inline std::vector<double> edges(){
  std::vector<double> e; for(int i=0;i<=20;i++) e.push_back(i*0.0025);
  e.push_back(0.08);e.push_back(0.11);e.push_back(0.14);
  for(int i=1;i<=33;i++) e.push_back(0.14+i*0.02);
  e.push_back(0.95);e.push_back(1.40);
  return e;
}
inline int findBin(const std::vector<double>&e,double M){
  if(M<e.front()||M>=e.back()) return -1;
  int i=std::upper_bound(e.begin(),e.end(),M)-e.begin()-1; return i;
}
