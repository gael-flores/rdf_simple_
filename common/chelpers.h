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


using RVecB = ROOT::VecOps::RVec<bool>;
using RVecC = ROOT::VecOps::RVec<char>;
using RVecD = ROOT::VecOps::RVec<double>;
using RVecF = ROOT::VecOps::RVec<float>;
using RVecI = ROOT::VecOps::RVec<int>;
using RVecL = ROOT::VecOps::RVec<long int>;
using RVecLL = ROOT::VecOps::RVec<long long int>;
using RVecU = ROOT::VecOps::RVec<unsigned int>;
using RVecUL = ROOT::VecOps::RVec<unsigned long int>;
using RVecULL = ROOT::VecOps::RVec<unsigned long long int>;

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


RVec<bool> overlapClean(RVecF gphi, RVecF geta, RVecF lphi, RVecF leta){
  RVec<bool> out;
  out.reserve(gphi.size());
  for (size_t i = 0; i < gphi.size(); i++){
    bool overlap = false;
    for (size_t j = 0; j < lphi.size(); j++){
      if (pow(DeltaPhi(gphi[i], lphi[j]),2)/pow(0.5,2) + pow(abs(geta[i] - leta[j]), 2)/pow(0.4,2) < 1){
	overlap = true;
	break;
      }
    }
    out.emplace_back(overlap);
  }
  return out;
}

RVecF best_Z_info(RVecF pt, RVecF eta, RVecF phi, RVecF mass, const RVec<size_t>& idx){
 
  RVecF out;
  out.reserve(6);
  size_t id1 = idx[0];
  size_t id2 = idx[1];
  ROOT::Math::PtEtaPhiMVector p0(pt[id1], eta[id1], phi[id1], mass[id1]);
  ROOT::Math::PtEtaPhiMVector p1(pt[id2], eta[id2], phi[id2], mass[id2]);
  out.emplace_back((p0+p1).pt());
  out.emplace_back((p0+p1).eta());
  out.emplace_back((p0+p1).phi());
  out.emplace_back((p0+p1).mass());
  out.emplace_back(DeltaR(eta[id1], eta[id2], phi[id1], phi[id2]));
  out.emplace_back(DeltaPhi(phi[id1], phi[id2]));
  return out;
}

RVecF correct_gammaIso_for_muons(const RVec<size_t>& idx,RVecF mpt, RVecF meta, RVecF mphi,RVecF gpt, RVecF geta, RVecF gphi,RVecF giso,RVec<bool> fsr) {
  RVecF result;  
  result.reserve(gpt.size());
  //add the FSR
  
  for(size_t i=0;i<gpt.size();++i) {
    float iso=giso[i]*gpt[i];
    if (fsr[i]){
      if (DeltaR(geta[i],meta[idx[0]],gphi[i],mphi[idx[0]])<0.3)
	iso=iso-mpt[idx[0]];
      if (DeltaR(geta[i],meta[idx[1]],gphi[i],mphi[idx[1]])<0.3)
	iso=iso-mpt[idx[1]];
    }
   
    if(iso<0)
      iso=0;
    result.emplace_back(iso/gpt[i]);
  }
  return result;
}

RVecF correct_gammaIso_for_photons(const int idx1, const int idx2, RVecF gpt, RVecF geta, RVecF gphi, RVecF giso) {
  RVecF result = giso;  
  if (DeltaR(geta[idx1], geta[idx2], gphi[idx1], gphi[idx2]) < 0.3){
    float iso1 = giso[idx1]*gpt[idx1];
    float iso2 = giso[idx2]*gpt[idx2];
    result[idx1] = (iso1 - gpt[idx2])/gpt[idx1] > 0.0 ? (iso1 - gpt[idx2])/gpt[idx1] : 0.0;
    result[idx2] = (iso2 - gpt[idx1])/gpt[idx2] > 0.0 ? (iso2 - gpt[idx1])/gpt[idx2] : 0.0;
  }
  return result;
}


RVecF correct_gammaIso(const RVecF& gpt, const RVecF& geta, const RVecF& gphi, const RVecF& giso, const RVec<bool>& gID) {
  RVecF result;
  result.reserve(gpt.size());

  for (size_t i = 0; i < gpt.size(); ++i) {
    float iso = giso[i] * gpt[i];

    // Correction for other photons
    for (size_t j = 0; j < gpt.size(); ++j) {
      if ((i == j) || gID[j] == 0)
        continue;

      if (DeltaR(geta[i], geta[j], gphi[i], gphi[j]) < 0.3)
        iso = iso - gpt[j];
    }

    // Ensure the corrected isolation is non-negative
    if (iso < 0)
      iso = 0;

    result.emplace_back(iso / gpt[i]);
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
      
      if (DeltaR(meta[i],geta[j],mphi[i],gphi[j])<0.4)
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
     if (gID[i]==1) {
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

