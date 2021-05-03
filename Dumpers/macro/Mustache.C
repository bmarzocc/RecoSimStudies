#include "DataFormats/CaloRecHit/interface/CaloCluster.h"
#include "DataFormats/CaloRecHit/interface/CaloClusterFwd.h"
#include "DataFormats/EgammaReco/interface/SuperCluster.h"
#include "TMath.h"
#include "TVector2.h"
#include <cmath>
#include <vector>
using namespace std;

double Phi_mpi_pi(double x) {
     while (x >= TMath::Pi())
       x -= 2 * TMath::Pi();
     while (x < -TMath::Pi())
       x += 2 * TMath::Pi();
     return x;
}
 
bool inMustache(const float maxEta, const float maxPhi, 
		const float ClustE, const float ClusEta, 
		const float ClusPhi){

      //bool inMust=false;
      //float eta0 = maxEta;
      //float phi0 = maxPhi;      
      
      const auto log10ClustE = std::log10(ClustE);

      constexpr float log10EMin = -3.;
      constexpr float etaMin = 0.;

      if( log10ClustE<log10EMin || std::abs(ClusEta)<etaMin) return false;

      constexpr float sqrtLogClustETuning = 1.1;
      constexpr float pUp[3] = {-0.107537, 0.590969, -0.076494};
      constexpr float pLow[3] = {-0.0268843, 0.147742, -0.0191235};
      constexpr float w0Up[2] = {-0.00681785, -0.00239516};
      constexpr float w1Up[2] = {0.000699995, -0.00554331};
      constexpr float w0Low[2] = {-0.00681785, -0.00239516};
      constexpr float w1Low[2] = {0.000699995, -0.00554331};
      
      const float sineta0 = std::sin(maxEta);
      const float eta0xsineta0 = maxEta * sineta0;

      //2 parabolas (upper and lower)
      //of the form: y = a*x*x + b

      //b comes from a fit to the width
      //and has a slight dependence on E on the upper edge
      // this only works because of fine tuning :-D
      const float sqrt_log10_clustE = std::sqrt(log10ClustE + sqrtLogClustETuning);
      const float b_upper = w1Up[0] * eta0xsineta0 + w1Up[1] / sqrt_log10_clustE -
          0.5 * (w1Up[0] * eta0xsineta0 + w1Up[1] / sqrt_log10_clustE +
                 w0Up[0] * eta0xsineta0 + w0Up[1] / sqrt_log10_clustE);
      const float b_lower = w0Low[0] * eta0xsineta0 + w0Low[1] / sqrt_log10_clustE -
          0.5 * (w1Low[0] * eta0xsineta0 + w1Low[1] / sqrt_log10_clustE +
                 w0Low[0] * eta0xsineta0 + w0Low[1] / sqrt_log10_clustE);

      //the curvature comes from a parabolic
      //fit for many slices in eta given a
      //slice -0.1 < log10(Et) < 0.1
      const float curv_up = eta0xsineta0 * (pUp[0] * eta0xsineta0 + pUp[1]) + pUp[2];
      const float curv_low = eta0xsineta0 * (pLow[0] * eta0xsineta0 + pLow[1]) + pLow[2];

      //solving for the curviness given the width of this particular point
      const float a_upper = (1. / (4. * curv_up)) - std::abs(b_upper);
      const float a_lower = (1. / (4. * curv_low)) - std::abs(b_lower);

      const double dphi = Phi_mpi_pi(ClusPhi - maxPhi);
      const double dphi2 = dphi * dphi;
      // minimum offset is half a crystal width in either direction
      // because science.
      constexpr float half_crystal_width = 0.0087;
      const float upper_cut = (std::max((1. / (4. * a_upper)), 0.0) * dphi2 + std::max(b_upper, half_crystal_width)) + half_crystal_width;
      const float lower_cut = (std::max((1. / (4. * a_lower)), 0.0) * dphi2 + std::min(b_lower, -half_crystal_width));

      const float deta = (1 - 2 * (maxEta < 0)) * (ClusEta - maxEta);  // sign flip deta
      return (deta < upper_cut && deta > lower_cut);
}

bool inDynamicDPhiWindow(const float seedEta, const float seedPhi,
			 const float ClustE, const float ClusEta,
			 const float ClusPhi) {

      // from Rishi's fits 06 June 2013 in log base 10
      const double absSeedEta = std::abs(seedEta);
      const double logClustEt = std::log10(ClustE / std::cosh(ClusEta));
      const double clusDphi = std::abs(Phi_mpi_pi(seedPhi - ClusPhi));

      if( ClustE<0. || absSeedEta<0.) return false;

      double yoffset, scale, xoffset, width, saturation, cutoff, maxdphi;
      const int etaBin = ( (int)(absSeedEta >= 1.479) + 
			   (int)(absSeedEta >= 1.75)  +
			   (int)(absSeedEta >= 2.0)     );

      switch( etaBin ) {
      case 0: // EB
	yoffset = 0.0280506;
        scale = 0.946048;
        xoffset = -0.101172;
        width = 0.432767;
        saturation = 0.14;
        cutoff = 0.6;
	break;
      case 1: // 1.479 -> 1.75
	yoffset = 0.0497038;
        scale = 0.975707;
        xoffset = -0.18149;
        width = 0.431729;
        saturation = 0.14;
        cutoff = 0.55;
	break;
      case 2: // 1.75 -> 2.0
	yoffset = 0.05643;
        scale = 1.60429;
        xoffset = -0.642352;
        width = 0.458106;
        saturation = 0.12;
        cutoff = 0.45;
	break;
      case 3: // 2.0 and up
	yoffset = 0.0928887;
        scale = 1.22321;
        xoffset = -0.260256;
        width = 0.345852;
        saturation = 0.12;
        cutoff = 0.3;
	break;
      default:
	throw cms::Exception("InValidEtaBin")
	  << "Calculated invalid eta bin = " << etaBin 
	  << " in \"inDynamicDPhiWindow\"" << std::endl;
      }
      
      maxdphi = yoffset+scale/(1+std::exp((logClustEt-xoffset)*width));
      maxdphi = std::min(maxdphi,cutoff);
      maxdphi = std::max(maxdphi,saturation);
      
      return clusDphi < maxdphi;
}

