#ifndef RecoSimStudies_Dumpers_DoubleCBPdf_H
#define RecoSimStudies_Dumpers_DoubleCBPdf_H

//from TrackingTools/TrackAssociator/test/DoubleCBPdf.h
// Developed by Wouter Hulsbergen

#include "RooAbsPdf.h"
#include "RooRealProxy.h"

class RooRealVar;
class RooAbsReal;

class DoubleCBPdf : public RooAbsPdf 
{
  public:
    explicit DoubleCBPdf(const char *name, const char *title,
                RooAbsReal& _x,
                RooAbsReal& _mean,
                RooAbsReal& _width,
                RooAbsReal& _alpha1,
                RooAbsReal& _n1,
                RooAbsReal& _alpha2,
                RooAbsReal& _n2);
  DoubleCBPdf(const DoubleCBPdf& other, const char* name=0) ;
  virtual TObject* clone(const char* newname) const { 
    return new DoubleCBPdf(*this,newname); }
  virtual ~DoubleCBPdf() { }
  
  protected:
    RooRealProxy x ;
    RooRealProxy mean;
    RooRealProxy width;
    RooRealProxy alpha1;
    RooRealProxy n1;
    RooRealProxy alpha2;
    RooRealProxy n2;
  
  double evaluate() const ;

};
#endif
