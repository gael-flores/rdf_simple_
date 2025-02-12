#include "chelpers.h"
#include "../analysis/ddp_vertex.h"

float reweightCT(RVecF ct, const float ctOld, const float ctNew){
   float out = 1.0;
   //std::cout << "ct old: " << ctOld << "\n";
   //std::cout << "ct New: " << ctNew << "\n";
   if (ctOld == 0 || ctNew == 0)
      return out;
   for (size_t i = 0; i < ct.size(); i++){
      float w = (exp(-ct[i]/ctNew)/ctNew)/(exp(-ct[i]/ctOld)/ctOld);
      //std::cout << "ct: " << ct[i] << "\n";
      //std::cout << "weight: " << w << "\n";
      out = out * w;
   }
   //std::cout << "out: " << out;
   return out;
}

RVecF ptSmearUp(RVecF pt, RVecF eta, RVecF phi, RVecF scale){
   RVecF out;
   for (size_t i = 0; i < pt.size(); i++){
      ROOT::Math::PtEtaPhiMVector p4(pt[i], eta[i], phi[i], 0);
      float eCorr = p4.energy() + std::abs(scale[i]);
      float scale = eCorr / p4.energy();
      out.emplace_back(p4.pt() * scale);
   }
   return out;
}

RVecF ptSmearDown(RVecF pt, RVecF eta, RVecF phi, RVecF scale){
   RVecF out;
   for (size_t i = 0; i < pt.size(); i++){
      ROOT::Math::PtEtaPhiMVector p4(pt[i], eta[i], phi[i], 0);
      float eCorr = p4.energy() - std::abs(scale[i]);
      float scale = eCorr / p4.energy();
      out.emplace_back(p4.pt() * scale);
   }
   return out;
}

RVecF getDxy(RVecF pt, RVecF eta, RVecF phi, RVecF isEB, RVecF isEE, const float idx1, const float idx2, const float mass){
  RVecF out;
  out.reserve(2);
  VertexCalculator *calc = new VertexCalculator();
  std::vector<float> kin_fit = calc->getVertexInfo(pt[idx1], eta[idx1], phi[idx1], isEE[idx1], isEB[idx1], pt[idx2], eta[idx2], phi[idx2], isEE[idx2], isEB[idx2], mass);
  delete calc;
  ROOT::Math::PtEtaPhiMVector g1(pt[idx1], eta[idx1], phi[idx1], 0);
  ROOT::Math::PtEtaPhiMVector g2(pt[idx2], eta[idx2], phi[idx2], 0);
  out.emplace_back(kin_fit[6]);
  out.emplace_back((g1+g2).mass());
  return out;
}
