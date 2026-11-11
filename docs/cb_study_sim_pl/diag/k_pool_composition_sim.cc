// Diagnostic: k_mix when the pool holds only legs of the background population (truth-selected),
// versus the population used by the macro's pool. Compared with k_true per class.
#include "common.h"
#include <array>
struct K{ uint32_t ev,x,y,z; bool operator==(const K&o)const{return ev==o.ev&&x==o.x&&y==o.y&&z==o.z;} };
struct KH{ size_t operator()(const K&k)const{ size_t h=k.ev; h=h*1000003^k.x; h=h*1000003^k.y; h=h*1000003^k.z; return h;} };
inline uint32_t bits(float f){uint32_t u;memcpy(&u,&f,4);return u;}
struct LegD{ double p,E,nx,ny,nz; };
static LegD mk(float p,float th,float ph){ LegD l; l.p=p; l.E=std::sqrt(l.p*l.p+ME*ME); double t=th*M_PI/180.,f=ph*M_PI/180.; l.nx=sin(t)*cos(f);l.ny=sin(t)*sin(f);l.nz=cos(t); return l;}
inline double mass(const LegD&a,const LegD&b){ double cs=a.nx*b.nx+a.ny*b.ny+a.nz*b.nz; double m2=2*ME*ME+2*(a.E*b.E-a.p*b.p*cs); return (m2>0?sqrt(m2):0)/1000.; }
struct Ev{ std::vector<LegD> P,M; int n=0; bool bkg=false; bool ls=false; };
static std::vector<double> ED; static int NB;
struct H3{ std::vector<double> pm,pp,mm; void init(){pm.assign(NB,0);pp.assign(NB,0);mm.assign(NB,0);} };
static inline void addM(std::vector<double>&v,double M){ int b=findBin(ED,M); if(b>=0) v[b]+=1; }
static double kwin(const H3&h,double lo,double hi){ double a=0,g=0; for(int i=0;i<NB;i++){ double m=0.5*(ED[i]+ED[i+1]); if(m<lo||m>=hi) continue; a+=h.pm[i]; if(h.pp[i]>0&&h.mm[i]>0) g+=2*sqrt(h.pp[i]*h.mm[i]); } return g>0?a/g:0; }
static void run(const std::vector<const Ev*>&seq,int depth,H3&h){ h.init(); std::vector<const Ev*> pool(depth,nullptr); int head=0,filled=0;
  for(const Ev*e:seq){ for(int ib=0;ib<filled;ib++){ const Ev&B=*pool[ib];
      for(auto&a:e->P) for(auto&c:B.M) addM(h.pm,mass(a,c)); for(auto&a:e->M) for(auto&c:B.P) addM(h.pm,mass(a,c));
      for(auto&a:e->P) for(auto&c:B.P) addM(h.pp,mass(a,c)); for(auto&a:e->M) for(auto&c:B.M) addM(h.mm,mass(a,c)); }
    pool[head]=e; head=(head+1)%depth; if(filled<depth) filled++; } }
int main(){
  ED=edges(); NB=ED.size()-1; TFile*f=TFile::Open(FN);
  std::vector<Ev> os; os.reserve(8000000); std::unordered_map<K,uint32_t,KH> osmap; osmap.reserve(16000000);
  { Reader R; R.open(f,"EpEm_ID","ep_","em_"); K prev{0,0,0,0}; bool have=false; std::vector<Row> g;
    auto flush=[&](){ if(g.empty()) return; std::vector<std::array<float,3>> P,M; for(auto&r:g){ P.push_back({r.a.p,r.a.th,r.a.ph}); M.push_back({r.b.p,r.b.th,r.b.ph}); }
      std::sort(P.begin(),P.end());P.erase(std::unique(P.begin(),P.end()),P.end()); std::sort(M.begin(),M.end());M.erase(std::unique(M.begin(),M.end()),M.end());
      Ev e; for(auto&a:P)e.P.push_back(mk(a[0],a[1],a[2])); for(auto&a:M)e.M.push_back(mk(a[0],a[1],a[2])); e.n=P.size()+M.size();
      // event flagged as background if any written row is not genuine (one-particle-twice not counted as bkg)
      for(auto&r:g){ int sa=(int)std::lround(r.a.simid),sb=(int)std::lround(r.b.simid); bool lep=(sa==2||sa==3)&&(sb==2||sb==3); bool samev=bitEq(r.a.vx,r.b.vx)&&bitEq(r.a.vy,r.b.vy)&&bitEq(r.a.vz,r.b.vz)&&bitEq(r.a.par,r.b.par)&&bitEq(r.a.proc,r.b.proc); bool gen=lep&&samev&&sa!=sb; bool opt=lep&&samev&&sa==sb; if(!gen&&!opt) e.bkg=true; }
      K k{bits(g[0].ev),bits(g[0].vrx),bits(g[0].vry),bits(g[0].vrz)}; osmap[k]=os.size(); os.push_back(std::move(e)); g.clear(); };
    for(Long64_t i=0;i<R.n;i++){ R.get(i); const Row&r=R.r; if(!(r.vrz>-500)) continue; K k{bits(r.ev),bits(r.vrx),bits(r.vry),bits(r.vrz)}; if(have&&!(k==prev)) flush(); g.push_back(r); prev=k;have=true; } flush(); }
  std::vector<Ev> ls[2];
  for(int t=0;t<2;t++){ Reader R; R.open(f,t?"EmEm_ID":"EpEp_ID",t?"em1_":"ep1_",t?"em2_":"ep2_"); K prev{0,0,0,0}; bool have=false; std::vector<Row> g;
    auto flush=[&](){ if(g.empty()) return; K k{bits(g[0].ev),bits(g[0].vrx),bits(g[0].vry),bits(g[0].vrz)}; if(osmap.find(k)==osmap.end()){ std::vector<std::array<float,3>> S; for(auto&r:g){ S.push_back({r.a.p,r.a.th,r.a.ph}); S.push_back({r.b.p,r.b.th,r.b.ph}); } std::sort(S.begin(),S.end()); S.erase(std::unique(S.begin(),S.end()),S.end()); Ev e; e.n=S.size(); e.ls=true; for(auto&a:S){ (t==0?e.P:e.M).push_back(mk(a[0],a[1],a[2])); } ls[t].push_back(std::move(e)); } g.clear(); };
    for(Long64_t i=0;i<R.n;i++){ R.get(i); const Row&r=R.r; if(!(r.vrz>-500)) continue; K k{bits(r.ev),bits(r.vrx),bits(r.vry),bits(r.vrz)}; if(have&&!(k==prev)) flush(); g.push_back(r); prev=k;have=true; } flush(); }
  size_t nLs=ls[0].size()+ls[1].size(), nos=os.size();
  // uniform interleave of LS-only events into the OS sequence
  struct It{double pos; const Ev*e;}; std::vector<It> items; items.reserve(nos+nLs);
  for(size_t i=0;i<nos;i++) items.push_back({double(i)+0.5,&os[i]});
  for(int t=0;t<2;t++){ size_t nt=ls[t].size(); for(size_t j=0;j<nt;j++) items.push_back({(j+0.5)/nt*nos,&ls[t][j]}); }
  std::stable_sort(items.begin(),items.end(),[](const It&a,const It&b){return a.pos<b.pos;});
  auto mk_seq=[&](std::function<bool(const Ev&)> sel){ std::vector<const Ev*> s; for(auto&it:items) if(sel(*it.e)) s.push_back(it.e); return s; };
  auto report=[&](const char*name,std::function<bool(const Ev&)> sel){ auto s=mk_seq(sel); H3 h; run(s,10,h); double pm=0,pp=0,mm=0; for(int i=0;i<NB;i++){double m=0.5*(ED[i]+ED[i+1]); if(m<0.14||m>=0.30)continue; pm+=h.pm[i];pp+=h.pp[i];mm+=h.mm[i];} printf("%-62s events %9zu | k_mix [0.14,0.30)=%.4f  [0.30,0.80)=%.4f [0.05,0.14)=%.4f | window pairs pm %.0f pp %.0f mm %.0f\n",name,s.size(),kwin(h,0.14,0.30),kwin(h,0.30,0.80),kwin(h,0.05,0.14),pm,pp,mm); };
  report("all events (uniform interleave)",[](const Ev&){return true;});
  report("n=2 with genuine row only (no bkg, no LS-only)",[](const Ev&e){return !e.ls&&e.n==2&&!e.bkg;});
  report("n=2 with background row only",[](const Ev&e){return !e.ls&&e.n==2&&e.bkg;});
  report("LS-only events only",[](const Ev&e){return e.ls;});
  report("n=2 background row + LS-only events",[](const Ev&e){return (e.ls)||(e.n==2&&e.bkg);});
  report("n>=3 events only",[](const Ev&e){return e.n>=3;});
  report("n>=3 + LS-only (truth-free)",[](const Ev&e){return e.ls||e.n>=3;});
  report("n>=3 + n=2 background + LS-only (background population)",[](const Ev&e){return e.ls||e.n>=3||(e.n==2&&e.bkg);});
  return 0;
}
