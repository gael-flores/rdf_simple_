#include "ROOT/RDataFrame.hxx"
#include "ROOT/RVec.hxx"
#include "TCanvas.h"
#include "TH1D.h"
#include "TLatex.h"
#include "Math/Vector4D.h"
#include "TStyle.h"

using namespace ROOT;
using namespace ROOT::VecOps;
float Z_mass=91;

// Reconstruct the best Z-> ll candidate
RVec<size_t> best_z(RVecF pt, RVecF eta, RVecF phi, RVecF mass, RVecI charge)
{

  RVec<size_t> result;
  result.reserve(2);
  auto idx_cmb = ROOT::VecOps::Combinations(pt, 2);
  auto best_mass = -1;
  size_t best_i1=0;
  size_t best_i2=0;
  for (size_t i = 0; i < idx_cmb[0].size(); i++) {
    const auto i1 = idx_cmb[0][i];
    const auto i2 = idx_cmb[1][i];
    if (charge[i1] != charge[i2]) {
      ROOT::Math::PtEtaPhiMVector p1(pt[i1], eta[i1], phi[i1], mass[i1]);
      ROOT::Math::PtEtaPhiMVector p2(pt[i2], eta[i2], phi[i2], mass[i2]);
      const auto this_mass = (p1 + p2).M();
      if (std::abs(this_mass-Z_mass)<std::abs(best_mass-Z_mass)) {
	best_mass=Z_mass;
	best_i1=i1;
	best_i2=i2;
      }
    }
  }
  result.emplace_back(best_i1);
  result.emplace_back(best_i2);
  return result;
}

RVecF correct_gammaIso_for_muons_and_photons(const RVec<size_t>& idx,RVecF mpt, RVecF meta, RVecF mphi,RVecF gpt, RVecF geta, RVecF gphi,RVecF giso,RVec<bool> gID) {
  RVecF result;  
  result.reserve(gpt.size());
  //add the FSR
  
  for(size_t i=0;i<gpt.size();++i) {
    float iso=giso[i]*gpt[i];
    if (DeltaR(geta[i],meta[idx[0]],gphi[i],mphi[idx[0]])<0.3)
      iso=iso-mpt[idx[0]];
    if (DeltaR(geta[i],meta[idx[1]],gphi[i],mphi[idx[1]])<0.3)
      iso=iso-mpt[idx[1]];
    for(size_t j=0;j<gpt.size();++j) {
      if ((i==j)||gID[j]==0)
	continue;
      if (DeltaR(geta[i],geta[j],gphi[i],gphi[j])<0.3)
	iso=iso-gpt[j];
    }
    
    if(iso<0)
      iso=0;
    result.emplace_back(iso/gpt[i]);
  }
  return result;
}

RVecF correct_muoniso_for_photons(RVecF mpt, RVecF meta, RVecF mphi,RVecF miso,RVecF gpt, RVecF geta, RVecF gphi,RVec<bool> gFSR) {
  RVecF result;  
  result.reserve(mpt.size());
  //add the FSR
  
  for(size_t i=0;i<mpt.size();++i) {
    float iso=miso[i]*mpt[i];
    for(size_t j=0;j<gpt.size();++j) {
      if (gFSR[j]==0)
	continue;
      
      if (DeltaR(meta[i],geta[j],mphi[i],gphi[j])<0.3)
	iso=iso-gpt[j];
    }
    if(iso<0)
      iso=0;
    result.emplace_back(iso/mpt[i]);
  }
  return result;
}




RVec<bool> fsr_recovery(const RVec<size_t>& idx,RVecF mpt, RVecF meta, RVecF mphi, RVecF mmass,RVecF gpt, RVecF geta, RVecF gphi,RVec<bool> gID) {
  RVec<bool> result;  
  ROOT::Math::PtEtaPhiMVector p1(mpt[idx[0]], meta[idx[0]], mphi[idx[0]], mmass[idx[0]]);
  ROOT::Math::PtEtaPhiMVector p2(mpt[idx[1]], meta[idx[1]], mphi[idx[1]], mmass[idx[1]]);
  ROOT::Math::PtEtaPhiMVector Z = p1+p2;
  float m = Z.M();

   //add the FSR
   for(size_t i=0;i<gpt.size();++i) {
     bool keep=false;
     if ((DeltaR(geta[i],meta[idx[0]],gphi[i],mphi[idx[0]])<0.5 ||DeltaR(geta[i],meta[idx[1]],gphi[i],mphi[idx[1]])<0.5)&&gID[i]==1) {
       ROOT::Math::PtEtaPhiMVector gamma(gpt[i], geta[i], gphi[i],0);
       float m_fsr = (Z+gamma).M();
       if(std::abs(m_fsr-Z_mass)<std::abs(m-Z_mass)) {
	 keep=true;
	 m=m_fsr;
	 Z=Z+gamma;
       }
     }
     result.emplace_back(keep);
   }
  return result;
}





float ll_mass(const RVec<size_t>& idx,RVecF mpt, RVecF meta, RVecF mphi, RVecF mmass) {
  RVec<bool> result;  
  ROOT::Math::PtEtaPhiMVector p1(mpt[idx[0]], meta[idx[0]], mphi[idx[0]], mmass[idx[0]]);
  ROOT::Math::PtEtaPhiMVector p2(mpt[idx[1]], meta[idx[1]], mphi[idx[1]], mmass[idx[1]]);
  ROOT::Math::PtEtaPhiMVector Z = p1+p2;
  float m = Z.M();
  return m;
}

float calculate_llgamma_mass(const RVec<size_t>& idx,RVecF mpt, RVecF meta, RVecF mphi, RVecF mmass,RVecF gpt, RVecF geta, RVecF gphi, RVecF gFSR) {
  ROOT::Math::PtEtaPhiMVector p1(mpt[idx[0]], meta[idx[0]], mphi[idx[0]], mmass[idx[0]]);
  ROOT::Math::PtEtaPhiMVector p2(mpt[idx[1]], meta[idx[1]], mphi[idx[1]], mmass[idx[1]]);
  ROOT::Math::PtEtaPhiMVector Z = p1+p2;
  for (size_t i=0;i<gpt.size();++i) {
    ROOT::Math::PtEtaPhiMVector g(gpt[i], geta[i], gphi[i], 0);
    if(gFSR[i]) {
      Z=Z+g;
    }
  }
  return Z.M();
}

 
