#include "math.h"
#include <Math/Vector4Dfwd.h>
#include <Math/GenVector/LorentzVector.h>
#include <Math/GenVector/PtEtaPhiM4D.h>
#include <Math/GenVector/PxPyPzM4D.h>
#include "TVector3.h"
#include "TMath.h"
class VertexCalculator {
 private:
  TVector3 vertex_;
  float pt_;
  float phi_;
  float d0_;
  float ip3d_;
  bool valid_;
  bool verbose_ = false;
  
 public:
  VertexCalculator() {}
  ~VertexCalculator() {} ;
  
  void setVertex(const TVector3& vertex){vertex_ = vertex;}
  void setPt(double pt){pt_ = pt;}  
  void setPhi(double phi){phi_ = phi;}
  void setD0(double d0){d0_ = d0;}
  void setIp3d(double ip3d){ip3d_ = ip3d;}
  void setValid(bool valid){valid_ = valid;}
  void setVerbose(bool verbose){verbose_ = verbose;}
  TVector3 vertex(){return vertex_;}    
  float pt(){return pt_;}
  float phi(){return phi_;}
  float d0(){return d0_;}
  float ip3d(){return ip3d_;}
  bool valid(){return valid_;}
  bool verbose(){return verbose_;}
    

    // Propagate prompt photon to ecal
TVector3 propToEcal(const double eta, const double phi, const int isEE, const int isEB){
  TVector3 out;
  float theta = 2*atan(exp(-eta));
  float ecal_z;
  if (theta == 0)
    ecal_z = 0;
  else
    ecal_z = 137./tan(theta);
  if (isEE){
    ecal_z = 324*ecal_z/abs(ecal_z);
    float r = ecal_z * tan(theta);
    float dx = cos(phi) * r;
    float dy = sin(phi) * r;
    out.SetXYZ(dx, dy, ecal_z);
  }
  else if (isEB){
    float r = 137.;
    float x = cos(phi) * r;
    float y = sin(phi) * r;
    out.SetXYZ(x, y, ecal_z);
  }
  return out;
}

float getRotAngle(const TVector3& v1, const TVector3& v2){
  TVector3 v3 = v1.Cross(v2);
  v3.SetMag(1.0);
  return acos(v3[2]);
}

TVector3 getRotAxis(const TVector3& v1, const TVector3& v2){
  if (v1[2]==0 && v2[2]==0)
    return TVector3(0,0,1.0);
  TVector3 v3 = v1.Cross(v2);
  v3.SetMag(1.0);
  TVector3 ez(0,0,1.0);
  return v3.Cross(ez);
}

double getPhi(const double mass, const double e1, const double e2){
  double cosphi = 1.0-mass*mass/(2.0*e1*e2);
  if (cosphi>=1 || cosphi<=-1)
    return M_PI;
  return acos(cosphi);
}

double getRadius(const TVector3& v1, const TVector3& v2, const double phi){
  double x1 = v1[0];
  double y1 = v1[1];
  double x2 = v2[0];
  double y2 = v2[1];
  return sqrt(pow((x1-x2)/2.,2)+pow((y1-y2)/2.,2))/sin(phi);
}

std::vector<std::pair<double, double>> getCenters(const TVector3& v1, const TVector3& v2, const double phi){
  double rad = getRadius(v1, v2, phi);
  double x1 = v1[0];
  double y1 = v1[1];
  double x2 = v2[0];
  double y2 = v2[1];
  double q = sqrt((x2-x1)*(x2-x1)+(y2-y1)*(y2-y1));
  double x3 = (x1+x2)/2.;
  double y3 = (y1+y2)/2.;
  double xp = x3 + sqrt(rad*rad-(q*q/4.))*(y1-y2)/q;
  double yp = y3 + sqrt(rad*rad-(q*q/4.))*(x2-x1)/q;
  double xm = x3 - sqrt(rad*rad-(q*q/4.))*(y1-y2)/q;
  double ym = y3 - sqrt(rad*rad-(q*q/4.))*(x2-x1)/q;
  std::pair<double, double> plus, minus;
  plus.first = xp;
  plus.second = yp;
  minus.first = xm;
  minus.second = ym;
  std::vector<std::pair<double, double>> out;
  out.push_back(plus);
  out.push_back(minus);
  return out;
}

TVector3 findIntersection(const double theta, const double radius, const double x0, const double y0, const TVector3& v1, const TVector3& v2, const double phi){
  TVector3 out(0,0,0);
  double a = 1.0+tan(theta)*tan(theta);
  double b = -2.*(x0+tan(theta)*y0);
  double c = pow(x0, 2) + pow(y0, 2) - pow(radius, 2);
  if ((pow(b, 2) - 4*a*c) < 0)
    return out;
  double xp = (-b + sqrt(pow(b, 2) - 4*a*c))/(2*a);
  double xm = (-b - sqrt(pow(b, 2) - 4*a*c))/(2*a);
  double yp = xp * tan(theta);
  double ym = xm * tan(theta);
  
  // Pick point that gives correct angle
  TVector3 p(xp, yp, 0.0);
  TVector3 m(xm, ym, 0.0);

  TVector3 p1 = p - v1;
  TVector3 p2 = p - v2;
  TVector3 m1 = m - v1;
  TVector3 m2 = m - v2;

  double phi1 = acos(p1*p2/(p1.Mag()*p2.Mag()));
  double phi2 = acos(m1*m2/(m1.Mag()*m2.Mag()));

  double dphi1 = fabs(phi1 - phi);
  double dphi2 = fabs(phi2 - phi);

  if (dphi1 <= dphi2){
    out.SetXYZ(xp, yp, 0);
  }else{
    out.SetXYZ(xm, ym, 0);
  }
  return out;
}


double getPt(const TVector3& vertex, const TVector3& v1, const TVector3& v2, const double e1, const double e2){
  TVector3 temp1 = v1-vertex;
  TVector3 temp2 = v2-vertex;
  double sinTheta1 = (vertex.Cross(temp1).Mag())/(vertex.Mag()*temp1.Mag())*copysign(1.0, vertex.Cross(temp1)[2]);
  double sinTheta2 = (vertex.Cross(temp2).Mag())/(vertex.Mag()*temp2.Mag())*copysign(1.0, vertex.Cross(temp2)[2]);
  double pt1 = e1 * sinTheta1;
  double pt2 = e2 * sinTheta2;
  return pt1 + pt2;
  
}


bool tooFar(const TVector3& vertex, const TVector3& v1, const TVector3& v2){
  double x,y;
  // Special cases for vertex at (0, y) or ecal clusters at (x, y1) (x, y2)
  if (vertex[0] == 0){
    if (v1[0]==v2[0])
      return fmax(v1[1], v2[1]) > vertex[1] && fmin(v1[1], v2[1]) < vertex[1];
    double m = (v2[1]-v1[1])/(v2[0]-v1[0]);
    y = v1[1]-m*v1[0];
    x = 0;
  } else {
    double mv = vertex[1]/vertex[0];
    if (v1[0]==v2[0]){
      x = v1[0];
      y = mv*x;
    } else {
      double mp = (v2[1]-v1[1])/(v2[0]-v1[0]);
      x = (v1[1]-mp*v1[0])/(mv-mp);
      y = mv*x;
    }
  }
  return (fmax(0, vertex[0]) > x && fmin(0, vertex[0]) < x && fmax(0, vertex[1]) > y && fmin(0, vertex[1]) < y);
}

int backwards(const TVector3& vertex, const TVector3& v1, const TVector3& v2){
  int out = 1;
  TVector3 midpoint;
  //  Find closest point to line connecting two ecal clusters
  if (v1[0]==v2[0])
    midpoint.SetXYZ(v1[0], 0, 0);
  else{
    double m = (v2[1]-v1[1])/(v2[0]-v1[0]);
    midpoint.SetX((v1[1]-m*v1[0])/(-1./m-m));
    midpoint.SetY(-1./m*(v1[1]-m*v1[0])/(-1./m-m));
  }

  if (midpoint[0]*vertex[0]+midpoint[1]*vertex[1] < 0)
    out = -1;
  if (verbose()){
    printf("Checking backwards: midpoint=(%f, %f) vertex=(%f,%f) \n", midpoint[0], midpoint[1],vertex[0],vertex[1]);
    printf("Midpoint * vertex = %f\n", midpoint*vertex);
    printf("scale=%d\n", out);
  }
  return out;
}

// Set the vertex, pt, and phi
void run(const TVector3& v1, const TVector3& v2, const double e1, const double e2, const double mass){


  double phi = M_PI;
  double pt = 999;
  TVector3 vertex(-999,-999,-999);

  setVertex(vertex);
  setPhi(phi);
  setPt(pt);
  setD0(-sqrt(vertex[0]*vertex[0]+vertex[1]*vertex[1]));
  setIp3d(-sqrt(vertex[0]*vertex[0]+vertex[1]*vertex[1]+vertex[2]*vertex[2]));
  setValid(false);

  TVector3 axis = getRotAxis(v1, v2);
  double theta = getRotAngle(v1, v2);
  phi = getPhi(mass, e1, e2);
  if (phi >= M_PI){
    if (verbose_)
      printf("Opening angle > pi, return false\n");
    return;
  }

  TVector3 temp1 = v1;
  TVector3 temp2 = v2;
  temp1.Rotate(theta, axis);
  temp2.Rotate(theta, axis);
  double radius = getRadius(temp1, temp2, phi);
  std::vector<std::pair<double, double>> centers = getCenters(temp1, temp2, phi);
  if (verbose()){
    printf("Ecal cluster 1 energy=%f (%f, %f)\nEcal cluster 2 energy=%f (%f, %f)\n",e1, temp1[0], temp1[1], e2, temp2[0], temp2[1]);
    printf("Photon phi=%f radius=%f\n", phi, radius);
  }

  for(uint i = 0; i < centers.size(); i++){
    pt = 999;
    std::pair<double, double> center = centers[i];
    if (verbose()){
      printf("   Checking circle with center (%f, %f)\n", center.first, center.second);
    }
    
    //Find ranges of phi values to check on locus of points
    double phi1 = atan2(temp1[1], temp1[0]);
    if (phi1 < 0)
      phi1 = M_PI*2 + phi1;
    double phi2 = atan2(temp2[1], temp2[0]);
    if (phi2 < 0)
      phi2 = M_PI*2 + phi2;
    double deltaPhi = phi2-phi1;
    //Bound deltaPhi between -pi,pi
    if (deltaPhi > M_PI){
      deltaPhi = M_PI*2 - deltaPhi;
      phi1 = atan2(temp2[1], temp2[0]);
      phi2 = atan2(temp1[1], temp1[0]);
    }
    if (deltaPhi < -M_PI){
      deltaPhi = M_PI*2 + deltaPhi;
      phi1 = atan2(temp1[1], temp1[0]);
      phi2 = atan2(temp2[1], temp2[0]);
    }
    
    int stepsPhi = 100;
    double stepSize = deltaPhi/stepsPhi;   
    for (int i=0; i<stepsPhi; i++){
      double angle = phi1 + i*stepSize;
      TVector3 temp = findIntersection(angle, radius, center.first, center.second, temp1, temp2, phi);
      double tempPt = getPt(temp, temp1, temp2, e1, e2);
      if (abs(tempPt) <= abs(pt)){
	pt = tempPt;
	vertex = temp;
      }
    }
    
    if (verbose())
      printf ("   Best fit vertex (%f, %f), pt = %f\n", vertex[0], vertex[1], pt);

    // Scale factor <0 means vertex is behind beamline
    int scale = backwards(vertex, temp1, temp2);
    bool unphysical = tooFar(vertex, temp1, temp2);
    if (verbose())
      printf("   checks: scale=%d toofar=%d\n", scale, unphysical);
    // Always pick physical vertex, otherwise only pick backwards if no other vertex found
    // The two circles cannot produce BOTH a physical and a backwards vertex
    if (!unphysical){
      if (verbose())
	printf("   Setting (%f, %f) as vertex\n", vertex[0], vertex[1]);
      vertex.Rotate(-theta, axis);
      if (verbose())
	printf("   Rotated vector around axis (%f, %f, %f) by angle %f\n", axis[0], axis[1], axis[2], -theta);
      setVertex(vertex);
      setPt(pt);
      setPhi(phi);
      setD0(scale*sqrt(vertex[0]*vertex[0]+vertex[1]*vertex[1]));
      setIp3d(scale*sqrt(vertex[0]*vertex[0]+vertex[1]*vertex[1]+vertex[2]*vertex[2]));
      setValid(true);
      if (verbose())
	printf("   New Vertex set at (%f, %f, %f), d0=%f, phi=%f\n", vertex[0], vertex[1], vertex[2], d0(), phi);
    }
  }
  return;
}

// With propagation, get the vertex within 3% of using precision calo clusters
void runPtEtaPhi(const double pt1, const double eta1, const double phi1, const int isEE1, const int isEB1, const double pt2, const double eta2, const double phi2, const int isEE2, const int isEB2, const double mass){
    TVector3 v1 = propToEcal(eta1, phi1, isEE1, isEB1);
    TVector3 v2 = propToEcal(eta2, phi2, isEE2, isEB2);
    ROOT::Math::PtEtaPhiMVector p4_1(pt1, eta1, phi1, 0);
    ROOT::Math::PtEtaPhiMVector p4_2(pt2, eta2, phi2, 0);
    run(v1, v2, p4_1.energy(), p4_2.energy(), mass);
    return;
}


// Call only after running to set p4 using vertex coords
ROOT::Math::PxPyPzMVector correctP4(const double pt, const double eta, const double phi, const int isEE, const int isEB){
  TVector3 v = vertex();
  TVector3 calo = propToEcal(eta, phi, isEE, isEB);
  TVector3 dR = calo - v;
  
  dR.SetMag(1.0);
  ROOT::Math::PtEtaPhiMVector p4(pt, eta, phi, 0);
  float e = p4.energy();
  ROOT::Math::PxPyPzMVector out(e*dR.x(), e*dR.y(), e*dR.z(), 0);
  return out;   
}

// Returns vector of [pt1, eta1, phi1, pt2, eta2, phi2, lxy, valid]
std::vector<float> getVertexInfo(const double pt1, const double eta1, const double phi1, const int isEE1, const int isEB1, const double pt2, const double eta2, const double phi2, const int isEE2, const int isEB2, const double mass){
  std::vector<float> out;
  runPtEtaPhi(pt1, eta1, phi1, isEE1, isEB1, pt2, eta2, phi2, isEE2, isEB2, mass);
  ROOT::Math::PxPyPzMVector p4_1 = correctP4(pt1, eta1, phi1, isEE1, isEB1);
  ROOT::Math::PxPyPzMVector p4_2 = correctP4(pt2, eta2, phi2, isEE2, isEB2);
  out.push_back(p4_1.pt());
  out.push_back(p4_1.eta());
  out.push_back(p4_1.phi());
  out.push_back(p4_2.pt());
  out.push_back(p4_2.eta());
  out.push_back(p4_2.phi());
  out.push_back(d0());
  out.push_back(valid());
  return out;
}

};

