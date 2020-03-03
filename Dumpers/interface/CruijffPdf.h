#ifndef RecoSimStudies_Dumpers_CruijffPdf_H
#define RecoSimStudies_Dumpers_CruijffPdf_H

//from TrackingTools/TrackAssociator/test/CruijffPdf.h
// Developed by Wouter Hulsbergen

#include "RooAbsPdf.h"
#include "RooRealProxy.h"

class RooRealVar;
class RooAbsReal;

class CruijffPdf : public RooAbsPdf 
{
  public:
    explicit CruijffPdf(const char *name, const char *title, RooAbsReal& _m,
		RooAbsReal& _m0, 
		RooAbsReal& _sigmaL, RooAbsReal& _sigmaR,
		RooAbsReal& _alphaL, RooAbsReal& _alphaR);

  CruijffPdf(const CruijffPdf& other, const char* name = 0);
  virtual TObject* clone(const char* newname) const { 
    return new CruijffPdf(*this,newname); } 
  virtual ~CruijffPdf() { }

  protected:
    RooRealProxy m;
    RooRealProxy m0;
    RooRealProxy sigmaL;
    RooRealProxy sigmaR;
    RooRealProxy alphaL;
    RooRealProxy alphaR;

  double evaluate() const;

};

#endif
