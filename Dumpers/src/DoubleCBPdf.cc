#include "RecoSimStudies/Dumpers/interface/DoubleCBPdf.h"

#include "RooRealVar.h"
#include "RooRealConstant.h"
using namespace RooFit;

//ClassImp(DoubleCBPdf) 

DoubleCBPdf::DoubleCBPdf(const char *name, const char *title, 
              RooAbsReal& _x,
              RooAbsReal& _mean,
              RooAbsReal& _width,
              RooAbsReal& _alpha1,
              RooAbsReal& _n1,
              RooAbsReal& _alpha2,
              RooAbsReal& _n2) 
  :
  RooAbsPdf(name,title), 
  x("x","x",this,_x),
  mean("mean","mean",this,_mean),
  width("width","width",this,_width),
  alpha1("alpha1","alpha1",this,_alpha1),
  n1("n1","n1",this,_n1),
  alpha2("alpha2","alpha2",this,_alpha2),
  n2("n2","n2",this,_n2)
{ 
} 


DoubleCBPdf::DoubleCBPdf(const DoubleCBPdf& other, const char* name) :  
  RooAbsPdf(other,name), 
  x("x",this,other.x),
  mean("mean",this,other.mean),
  width("width",this,other.width),
  alpha1("alpha1",this,other.alpha1),
  n1("n1",this,other.n1),
  alpha2("alpha2",this,other.alpha2),
  n2("n2",this,other.n2)
{ 
} 
 
double DoubleCBPdf::evaluate() const 
{ 
  double t = (x-mean)/width;
  double val = -99.;
  if(t>-alpha1 && t<alpha2){
     val = exp(-0.5*t*t);
  }else if(t<=-alpha1){
     double alpha1invn1 = alpha1/n1;
     val = exp(-0.5*alpha1*alpha1)*pow(1. - alpha1invn1*(alpha1+t), -n1);
  }else if(t>=alpha2){
     double alpha2invn2 = alpha2/n2;
     val = exp(-0.5*alpha2*alpha2)*pow(1. - alpha2invn2*(alpha2-t), -n2);        
  }
  //if(!std::isnormal(val)) {
  //   printf("bad val: x = %5f, t = %5f, mean = %5f, sigma = %5f, alpha1 = %5f, n1 = %5f, alpha2 = %5f, n2 = %5f\n",double(x), t, double(mean),double(width),double(alpha1),double(n1),double(alpha2), double(n2));
  //   printf("val = %5f\n",val);
  //}  
  return val;
} 
