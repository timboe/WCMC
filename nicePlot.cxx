// Dear emacs, this is -*- c++ -*-
// -------------------------------------------------------------
//  author: Tim Martin <Tim.Martin@cern.ch>
// -------------------------------------------------------------

#include "nicePlot.h"

#include <vector>
#include <cmath>
#include <set>
#include <sstream>
#include <algorithm>
#include <TF1.h>
#include <TProfile.h>
#include <TH1F.h>
#include <TH2F.h>
#include <TH1D.h>
#include <TH2D.h>
#include <TGraphAsymmErrors.h>
#include <TCanvas.h>
#include <TPad.h>
#include <TFile.h>
#include <TLatex.h>
#include <TLegend.h>
#include <TArrow.h>
#include <iostream>
#include <TStyle.h>
#include <TColor.h>
#include <TSystem.h>
#include <TLine.h>
#include <TError.h>

std::vector<nicePlot*> bookOutput::rp;
unsigned bookOutput::_break = 9999;

void err(TString e) {
  std::cout <<  "\033[1;31m ERROR:" << e << "\033[0m\n";
}

nicePlot::nicePlot() {
  reset();
}

nicePlot::nicePlot(TString _xTitle, TString _yTitle, TString _rTitle) {
  reset();
  init(_xTitle, _yTitle, _rTitle);
}

nicePlot::nicePlot(nicePlot& _copy) {
  reset();
  copy(_copy);
}

nicePlot::nicePlot(nicePlot* _copy) {
  reset();
  copy(*_copy);
}

void nicePlot::reset() {
  normAlreadDone = false;
  isReg = false;
  doLegend = true;
  line_width = 1;
  divideBinWidth = ratioOnly = includeBothInRatioError = killLowStatBins = false;
  legendX = legendY = 0.5;
  logX = logY = logZ = norm_ndc = false;
  xMin = xMax = yMin = yMax = rMin = rMax = 0;
  useAltColourScheme(0);
  n_mc = 0;
  stamp = "NONE";
  colab = "ATLAS";
  drawOp = "";
  plot2D = nullptr;
  norm_command = 0;
  norm_upper = norm_lower = -1;
  autoY = autoX = 0;
  rebin = 0;
  ratioIsDataStyle.clear();
  c = nullptr;
  fitMin = fitMax = 0;
  fit = "";
  mc_marker_size = 0;
  print_fit = true;
  ySpaceModLin = 1.2;
  ySpaceModLog = 2;
  ratioLineValue = 1;
  stretch = false;
  mc_stack = nullptr;
}

void nicePlot::init(TString _xTitle, TString _yTitle, TString _rTitle, bool quiet) {
  book();
  xTitle = _xTitle; yTitle = _yTitle; rTitle = _rTitle;
  if (!quiet) std::cout <<  "\033[1;36m Init:" << _xTitle << " | " << _yTitle << " | " << _rTitle << "\033[0m\n";
}

void nicePlot::book() {
  if (isReg == false) {
    bookOutput::get().book(this);
    isReg = !isReg;
  }
}

void nicePlot::unbook() {
  if (isReg == true) {
    bookOutput::get().unbook(this);
    isReg = !isReg;
  }
}

void nicePlot::useAltColourScheme(int _i) {
  static bool doOnce = true;
  int ci;
  if (doOnce == true) {
    doOnce = false;
    ci = 2555;
    new TColor(++ci, 22/256., 147/256., 167/256.); // blueLagoon
    new TColor(++ci, 200/256., 207/256., 2/256.); // springWind
    new TColor(++ci, 204/256., 12/256., 57/256.); // kaaskopPink
    new TColor(++ci, 230/256., 120/256., 30/256.); // orangeSakura
    new TColor(++ci, 167/256., 2/256., 103/256.); // mai
    new TColor(++ci, 14/256., 78/256., 173/256.); // skyblue
    new TColor(++ci, 131/256., 163/256., 0/256.); // kotak
    gStyle->SetHatchesLineWidth(2);
    //gStyle->SetHatchesSpacing(0.05);
  }
  switch (_i) {
    case 0: // default
    data_col[0] = kGreen-6;
    data_col[1] = kRed-6;
    data_col[2] = kBlue-6;
    data_col[3] = kAzure-6;
    data_col[4] = kMagenta-6;
    data_col[5] = kOrange-6;
       
    data_fill[0] = 1001;
    data_fill[1] = 3345;
    data_fill[2] = 3354;
    data_fill[3] = 1001;
    data_fill[4] = 3345;
    data_fill[5] = 3354;

    data_fill[0] = 1001;
    data_fill[1] = 1001; // All solid
    data_fill[2] = 1001;
    data_fill[3] = 1001;
    data_fill[4] = 1001; 
    data_fill[5] = 1001;

    data_sty[0] = 20;
    data_sty[1] = 22;
    data_sty[2] = 21;
    data_sty[3] = 23;
    data_sty[4] = 24;
    data_sty[5] = 25;

    mc_col[0] = kAzure-3;
    mc_col[1] = kBlue+1;
    mc_col[2] = kMagenta+1;
    mc_col[3] = kGreen+1;
    mc_col[4] = kAzure+3;
    mc_col[5] = kOrange-1;
    mc_col[6] = kMagenta+4;
    mc_sty[0] = 1;
    mc_sty[1] = 5;
    mc_sty[2] = 7;
    mc_sty[3] = 2;
    mc_sty[4] = 3;
    mc_sty[5] = 9;
    mc_sty[6] = 6;
    break;
    case 1: // A, B, C, D, E, F Comp. #2 Different data
    data_col[0] =  kGreen-6;
    data_col[1] = kRed-6;
    data_col[2] = kBlue-6;
    ci = 2555;
    mc_col[6] = ++ci;
    mc_col[1] = ++ci; // should stay 1
    mc_col[4] = ++ci; // should be 4
    mc_col[3] = ++ci; // should stay 3
    mc_col[2] = ++ci; // should be 2
    mc_col[5] = ++ci;
    mc_col[0] = ++ci; // should be 0
    mc_sty[0] = 1;
    mc_sty[1] = 1;
    mc_sty[2] = 1;
    mc_sty[3] = 1;
    mc_sty[4] = 1;
    mc_sty[5] = 1;
    mc_sty[6] = 1;
    break;
    case 2: // A, B, C, D, E, F Comp. #2 Different data
    data_col[0] =  kGreen-6;
    data_col[1] = kRed-6;
    data_col[2] = kBlue-6;
    ci = 2555;
    mc_col[3] = ++ci; // moves to 3
    mc_col[2] = ++ci; // moves to 2
    mc_col[0] = ++ci; // moves to 0
    mc_col[4] = ++ci; 
    mc_col[1] = ++ci; // moves to 1
    mc_col[5] = ++ci;
    mc_col[6] = ++ci;
    break;
    case 3: // A-D, B-F, C-G Comp# 1
    data_col[0] =  kGreen-6;
    data_col[1] = kRed-6;
    data_col[2] = kBlue-6;
    mc_col[0] = kRed+1;
    mc_col[3] = kRed-4;
    mc_col[1] = kBlue+1;
    mc_col[4] = kBlue-4;
    mc_col[2] = kGreen+1;
    mc_col[5] = kGreen-4;
    mc_col[6] = kBlack;
    break; // A-D, B-F, C-G Comp# 2
    case 4:
    data_col[0] =  kGreen-6;
    data_col[1] = kRed-6;
    data_col[2] = kBlue-6;
    mc_col[0] = kOrange-1;
    mc_col[3] = kOrange+4;
    mc_col[1] = kViolet-1;
    mc_col[4] = kViolet+4;
    mc_col[2] = kTeal-1;
    mc_col[5] = kTeal+4;
    mc_col[6] = kBlack;
    break;
  }
}

void nicePlot::useScaleColourScheme(int _i) {
  int _baseCol = mc_col[_i];
  mc_col[0] = _baseCol+3;
  mc_col[3] = _baseCol;
  mc_col[1] = _baseCol-3;
  mc_col[4] = _baseCol+1;
  mc_col[2] = _baseCol-2;
  mc_col[5] = _baseCol+2;
  mc_col[6] = _baseCol-1;
}

void nicePlot::copy( nicePlot& _copy ) {
  //cout << "INFO: Using copy constructor" << endl;
  book();
  divideBinWidth = _copy.divideBinWidth;
  includeBothInRatioError = _copy.includeBothInRatioError;
  ratioOnly = _copy.ratioOnly;
  xTitle = _copy.xTitle;
  yTitle = _copy.yTitle;
  rTitle = _copy.rTitle;
  doLegend = _copy.doLegend;
  legendX = _copy.legendX;
  legendY = _copy.legendY;
  logX = _copy.logX;
  logY = _copy.logY;
  logR = _copy.logR;
  logZ = _copy.logZ;
  line_width = _copy.line_width;
  for (int i=0; i < n_sty; ++i) {
    mc_col[i] = _copy.mc_col[i];
    mc_sty[i] = _copy.mc_sty[i];
    if (i >= 3) continue;
    data_col[i] =  _copy.data_col[i];
  }
  n_mc = _copy.n_mc;
  xMin = _copy.xMin;
  xMax = _copy.xMax;
  yMin = _copy.yMin;
  yMax = _copy.yMax;
  rMin = _copy.rMin;
  rMax = _copy.rMax;
  legendX = _copy.legendX;
  legendY = _copy.legendY;
  legendScale = _copy.legendScale;
  xTitle = _copy.xTitle;
  yTitle = _copy.yTitle;
  rTitle = _copy.rTitle;
  dataSyst = _copy.dataSyst;
  dataStat = _copy.dataStat;
  ratioIsDataStyle = _copy.ratioIsDataStyle;
  mc = _copy.mc;
  mc_do = _copy.mc_do;
  ratio = _copy.ratio;
  bin_label_x = _copy.bin_label_x;
  bin_label_y = _copy.bin_label_y;
  plotText = _copy.plotText;
  plotTextScale = _copy.plotTextScale;
  plotTextX = _copy.plotTextX;
  plotTextY = _copy.plotTextY;
  stamp = _copy.stamp;
  stampX = _copy.stampX;
  stampY = _copy.stampY;
  colab  = _copy.colab;
  plot2D = _copy.plot2D;
  norm_command = _copy.norm_command;
  norm_lower = _copy.norm_lower;
  norm_upper = _copy.norm_upper;
  norm_ndc = _copy.norm_ndc;
  autoY = _copy.autoY;
  autoX = _copy.autoX;
  lineX1 = _copy.lineX1;
  lineX2 = _copy.lineX2;
  lineY1 = _copy.lineY1;
  lineY2 = _copy.lineY2;
  lineStyle = _copy.lineStyle;
  lineWidth = _copy.lineWidth;
  linePoint = _copy.linePoint;
  rebin = _copy.rebin;
  c = _copy.c;
  drawOp = _copy.drawOp;
  killLowStatBins = _copy.killLowStatBins;
  fit = _copy.fit;
  fitMin = _copy.fitMin;
  fitMax = _copy.fitMax;
  mc_marker_size = _copy.mc_marker_size;
  print_fit = _copy.print_fit;
  ySpaceModLin = _copy.ySpaceModLin;
  ySpaceModLog = _copy.ySpaceModLog;
  ratioLineValue = _copy.ratioLineValue;
  stretch = _copy.stretch;
  mc_stack = _copy.mc_stack;
}

void nicePlot::setBounds(double _xMin, double _xMax, double _yMin, double _yMax, double _rMin, double _rMax) {
  xMin = _xMin; xMax = _xMax; 
  if (!(_yMin == 0. && _yMax == 0.)) { yMin = _yMin; yMax = _yMax; } 
  if (!(_rMin == 0. && _rMax == 0.)) { rMin = _rMin; rMax = _rMax; }
}

void nicePlot::setYBounds(double _yMin, double _yMax) {
  yMin = _yMin; yMax = _yMax;
}

void nicePlot::setRBounds(double _rMin, double _rMax) {
  rMin = _rMin; rMax = _rMax;
}

ObjType_t nicePlot::getObjectType(TObject* _o) {
  if ( _o->ClassName() == _aGraph.ClassName() ) {
    return kObjTypeGraph;
  } else if ( _o->ClassName() == _aFHist.ClassName() || _o->ClassName() == _aDHist.ClassName() || _o->ClassName() == _aProfile.ClassName() ) {
    return kObjTypeTH1;
  } else if ( _o->ClassName() == _aF2Hist.ClassName() || _o->ClassName() == _aD2Hist.ClassName()) {
    return kObjTypeTH2;
  } else return kUnsuportedType;
}

void nicePlot::addLine(double _x1, double _y1, double  _x2, double _y2, int _lineStyle, int _lineWidth, nicePlotArrow_t _point) {
  lineX1.push_back(_x1);
  lineX2.push_back(_x2);
  lineY1.push_back(_y1);
  lineY2.push_back(_y2);
  lineStyle.push_back(_lineStyle);
  lineWidth.push_back(_lineWidth);
  linePoint.push_back(_point);
}

nicePlot* nicePlot::addData(TGraphAsymmErrors* _gStat, TString _name, bool _ratioWithData) {
  if (_gStat == nullptr) { err("no data graph pointer");  return this; }
  if (dataSyst.size() >= n_sty) { err("Too much data"); return this; }
  TGraphAsymmErrors* _gSyst = (TGraphAsymmErrors*) _gStat->Clone();
  _gSyst->SetFillColor(data_col[  dataSyst.size() ]);
  _gSyst->SetFillStyle(data_fill[ dataSyst.size() ]);
  _gSyst->SetTitle(_name);

  _gStat->SetMarkerColor(0); // WC MC
  _gStat->SetLineColor(0); // WC MC
  _gStat->SetMarkerStyle(data_sty[  dataSyst.size() ]); // TIM M TEMP

  _gSyst->SetMarkerColor(0); // WC MC
  _gSyst->SetLineColor(0); // WC MC
  _gSyst->SetMarkerStyle(data_sty[  dataSyst.size() ]); // TIM M TEMP

  if (_ratioWithData == true) {
    if (dataSyst.size() == 0) {
      err("Cannot automagically do a ratio of this data to other data without the other data being added first");
      return this;
    }
    //If this is the first call, then also add Data/Data ratio
    if (ratio.size() == 0) {
      addRatio( dataSyst.back(), dataSyst.back(), true);
      addRatio( dataStat.back(), dataStat.back() );
    }
    addRatio( _gSyst, dataStat.back(), true );
    addRatio( _gStat, dataStat.back());
  }

  dataSyst.push_back(_gSyst);
  dataStat.push_back(_gStat);
  return this;
}

nicePlot* nicePlot::addData(TH1* _hStat, TString _name, float _norm, bool _ratioWithData) {
  if (_hStat == nullptr) { err("Null ptr"); return this; }
  TH1* _cloneStat = static_cast<TH1*>(_hStat->Clone());
  applyOptionsToHist(_cloneStat, _norm);
  return addData(new TGraphAsymmErrors(_cloneStat), _name, _ratioWithData);
}

nicePlot* nicePlot::addData(TFile* _f, TString _hist_name_in_file, TString _name, float _norm, bool _ratioWithData) {
  TObject* _syst = _f->Get(_hist_name_in_file);
  if ( _syst == nullptr) { err("Cannot find " + _hist_name_in_file); return this; }
  if      (getObjectType(_syst) == kObjTypeGraph) addData( static_cast<TGraphAsymmErrors*>( _syst ), _name, _ratioWithData );
  else if (getObjectType(_syst) == kObjTypeTH1)   addData( static_cast<TH1*>( _syst ), _name, _norm, _ratioWithData );
  else err("unsupported type");
  return this;
}

nicePlot* nicePlot::addDataSystematic(TGraphAsymmErrors* _gSyst, TGraphAsymmErrors* _gNominal, bool _symmetrise) {
  if (_gSyst == nullptr) {err("no data graph pointer");  return this; }
  if (dataSyst.size() == 0) { err("Add data before adding systematics"); return this; }
  graphAddSyst(_gSyst, _gNominal, _symmetrise);
  return this;
}

nicePlot* nicePlot::addDataSystematic(TH1* _hSyst, TH1* _hNominal, bool _symmetrise) {
  if (_hSyst == nullptr) { err("Null ptr"); return this; }
  TGraphAsymmErrors* _gNominal = nullptr;
  if (_hNominal != nullptr) _gNominal = new TGraphAsymmErrors(_hNominal);
  return addDataSystematic(new TGraphAsymmErrors(_hSyst), _gNominal, _symmetrise);
}

nicePlot* nicePlot::addDataSystematic(TFile* _f, TString _syst_name_in_file, bool _symmetrise) {
  TObject* _syst = _f->Get(_syst_name_in_file);
  if ( _syst == nullptr) { err("Cannot find " + _syst_name_in_file); return this; }
  if      (getObjectType(_syst) == kObjTypeGraph ) addDataSystematic( static_cast<TGraphAsymmErrors*>( _syst), nullptr, _symmetrise );
  else if (getObjectType(_syst) == kObjTypeTH1 )   addDataSystematic( static_cast<TH1*>( _syst), nullptr, _symmetrise );
  else err("unsuported type");
  return this;
}

nicePlot* nicePlot::addDataSystematic(TFile* _f, TString _syst_name_in_file, TString _nominal_name_in_file, bool _symmetrise) {
  TObject* _syst = _f->Get(_syst_name_in_file);
  TObject* _nominal = _f->Get(_nominal_name_in_file);
  if ( _syst == nullptr || _nominal == nullptr) { err("Cannot find " + _syst_name_in_file + " or " + _nominal_name_in_file); return this; }
  if (_syst->ClassName() != _nominal->ClassName()) { err("With addDataSystematic, both must be of the same type"); return this; }
  if      (getObjectType(_syst) == kObjTypeGraph ) addDataSystematic( static_cast<TGraphAsymmErrors*>(_syst), static_cast<TGraphAsymmErrors*>(_nominal), _symmetrise );
  else if (getObjectType(_syst) == kObjTypeTH1 )   addDataSystematic( static_cast<TH1*>(_syst), static_cast<TH1*>(_nominal), _symmetrise );
  else err("unsuported type");
  return this;
}

nicePlot* nicePlot::add2D(TH2* _h) {
  if (_h == nullptr) return this;
  if (plot2D != nullptr) { err("Already have a plot2D"); return this; }
  TH2* _clone = static_cast<TH2*>(_h->Clone());

  if (_clone->GetMaximum() > 1) {
    double _minval = 999.;
    for (int bin = 0; bin < _clone->GetNcells(); ++bin) {
      if (_clone->GetBinContent(bin) < _minval && fabs(_clone->GetBinContent(bin)) > 1e-15) _minval = _clone->GetBinContent(bin); // Set "empty" bins to the same as the min bin
    }
    for (int bin = 0; bin < _clone->GetNcells(); ++bin) {
      if (_clone->GetBinContent(bin) == 0) _clone->SetBinContent(bin, _minval); // Set "empty" bins to the same as the min bin
    }
  }

  if (norm_command == 1) _clone->Scale( 1./_clone->Integral() );
  plot2D = _clone;
  return this;
}

nicePlot* nicePlot::add2D(TFile* _f, TString _hist_name_in_file) {
  TObject* _o = _f->Get(_hist_name_in_file);
  if ( _o == nullptr) { err("Cannot find " + _hist_name_in_file); return this; }
  if ( getObjectType(_o) == kObjTypeTH2 ) add2D( static_cast<TH2*>( _o ) );
  else err("unsuported type");
  return this;
}

nicePlot* nicePlot::add2DEfficiency(TFile* _f, TString _hist_name_in_file_1, TString _hist_name_in_file_2) {
  TObject* _o1 = _f->Get(_hist_name_in_file_1);
  TObject* _o2 = _f->Get(_hist_name_in_file_2);
  if ( _o1 == nullptr ) { err("Cannot find " + _hist_name_in_file_1); return this; }
  if ( _o2 == nullptr ) { err("Cannot find " + _hist_name_in_file_2); return this; }
  if (getObjectType(_o1) != kObjTypeTH2 || getObjectType(_o2) != kObjTypeTH2) {
    err("either " + _hist_name_in_file_1 + " or " + _hist_name_in_file_2 + " are wrong type");
    return this;
  }

  TH2* _a = static_cast<TH2*>(_o1);
  TH2* _b = static_cast<TH2*>(_o2);
  TH2* _c = static_cast<TH2*>( _a->Clone() );
  _c->Divide(_b);
  // Bins missing entires in either hist have zero efficiency
  for (int bin = 0; bin < _a->GetNcells(); ++bin) {
    if (_a->GetBinContent(bin) == 0 || _b->GetBinContent(bin) == 0) _c->SetBinContent(bin, 1e-10);
  }
  return add2D(_c);
}

nicePlot* nicePlot::addStackMC(TH1* _h, TString _name) {
  if (_h == 0) return this;
  if (mc_stack == nullptr) {
    mc_stack = new THStack("", "");
  }

  TH1* _clone = static_cast<TH1*>(_h->Clone());
  applyOptionsToHist(_clone, 1);

  int _n = n_mc;
  if (_n >= n_sty) { err("Too many MC lines in " + _name); return this; }
  _clone->SetLineColor(mc_col[_n]);
  _clone->SetFillColor(mc_col[_n]);
  _clone->SetFillStyle(1001);
  _clone->SetLineWidth(line_width);
  _clone->SetMarkerColor(mc_col[_n]);
  _h->SetLineStyle(mc_sty[_n]);
  _clone->SetMarkerSize(mc_marker_size);
  _clone->SetTitle(_name);
  mc_stack->Add(_clone, "HIST");
  // mc_do.push_back(TString(""));
  ++n_mc;

  return this;
}

nicePlot* nicePlot::addMC(TGraphAsymmErrors* _g, TString _name, bool _ratioWithData) {
  if (_g == 0) return this;
  int _n = n_mc;
  if (_n >= n_sty) { err("Too many MC lines in " + _name); return this; }
  _g->SetLineColor(mc_col[_n]);
  _g->SetLineWidth(line_width);
  _g->SetMarkerColor(mc_col[_n]);
  _g->SetLineStyle(mc_sty[_n]);
  _g->SetMarkerSize(mc_marker_size);
  _g->SetTitle(_name);
  mc.push_back(_g);
  mc_do.push_back(TString(""));
  ++n_mc;

  // Kill low stat MC bins //TODO make this an option
  if (killLowStatBins) {
    for (int j = 0; j < _g->GetN(); ++j) {
      bool rm = false;
      if ( fabs(_g->GetY()[j]) > 1e-15 && (_g->GetEYhigh()[j] / _g->GetY()[j] > 0.5 || _g->GetEYlow()[j] / _g->GetY()[j] > 0.5) ) rm = true;
      if ( fabs(_g->GetEYhigh()[j]) < 1e-15 || fabs(_g->GetEYlow()[j]) < 1e-15 ) rm = true;
      if (rm) {
        _g->GetEYhigh()[j] = 0;
        _g->GetEYlow()[j] = 0;
        _g->GetY()[j] = 0;
      }
    }
  }

  if (_ratioWithData == true) {
    if (dataSyst.size() == 0) {
      err("Cannot automagically do a ratio of MC with data without data being added first");
      return this;
    }
    if (dataSyst.size() > 1) err("More than 1 data. Will do MC/Data ratio with last added data graph");
    //If this is the first call, then also add Data/Data ratio
    if (ratio.size() == 0) {
      addRatio( dataSyst.back(), dataSyst.back(), true );
      addRatio( dataStat.back(), dataStat.back() );
    }
    addRatio( _g, dataSyst.back() );
  }
  return this;
}

nicePlot* nicePlot::addMC(TH1* _h, TString _name, bool _ratioWithData, float _norm) {
  if (_h == 0) return this;
  TH1* _clone = static_cast<TH1*>(_h->Clone());
  applyOptionsToHist(_clone, _norm);
  return addMC(new TGraphAsymmErrors(_clone), _name, _ratioWithData);
}

nicePlot* nicePlot::addMC(TFile* _f, TString _hist_name_in_file, TString _name, bool _ratioWithData, float _norm) {
  TObject* _o = _f->Get(_hist_name_in_file);
  if ( _o == nullptr ) { err("Cannot find " + _hist_name_in_file); return this; }

  if (getObjectType(_o) == kObjTypeGraph) {
    return addMC( static_cast<TGraphAsymmErrors*>( _o ), _name, _ratioWithData );
  } else if (getObjectType(_o) == kObjTypeTH1) {
    return addMC( static_cast<TH1*>( _o ), _name, _ratioWithData, _norm );
  } else { 
    err("unsuported type");
    std::cerr << _o->ClassName() << std::endl;
  }
  return this;
}

TGraphAsymmErrors* nicePlot::graphDivide(TGraphAsymmErrors* _a, TGraphAsymmErrors* _b, bool _includeBothInRatioError) {
  TGraphAsymmErrors* _ratio = (TGraphAsymmErrors*) _a->Clone();
  if (_a->GetN() != _b->GetN()) {
    err(TString("different numbers of bins in the ratio! ") + _a->GetName() + " " + _b->GetName());
    return _ratio;
  }
  // Do this manually - much less can go wrong
  for (int i=0; i < _ratio->GetN(); i++) {
    double _rVal = _ratio->GetY()[i] / _b->GetY()[i];
    if (std::isinf(_rVal) || std::isnan(_rVal)) { // Prevent NaN or INF getting to the pad
      _ratio->GetY()[i]      = 0;
      _ratio->GetEYhigh()[i] = 0;
      _ratio->GetEYlow()[i]  = 0; 
    } else {
      // ERRORS - the error on the ratio reflects the fractional errors on the numerator and denominator
      const double _small = 1e-15;
      double _fErrHigh0 = 0, _fErrHigh1 = 0, _fErrLow0 = 1, _fErrLow1 = 0, _fErrHigh = 0, _fErrLow = 0;
      if (fabs(_a->GetY()[i]) > _small) { // Check not zero
        _fErrHigh0 = _a->GetEYhigh()[i] / _a->GetY()[i];
        _fErrLow0  = _a->GetEYlow()[i]  / _a->GetY()[i];
      }
      if (_includeBothInRatioError && fabs(_b->GetY()[i]) > _small ) { // Check not zero
        _fErrHigh1 = _b->GetEYhigh()[i] / _b->GetY()[i];
        _fErrLow1  = _b->GetEYlow()[i]  / _b->GetY()[i];
      }
      _fErrHigh = std::sqrt( (_fErrHigh0*_fErrHigh0) + (_fErrHigh1*_fErrHigh1) );
      _fErrLow  = std::sqrt( (_fErrLow0 *_fErrLow0 ) + (_fErrLow1 *_fErrLow1 ) );

      _ratio->GetY()[i]      /= _b->GetY()[i];
      _ratio->GetEYhigh()[i] = _ratio->GetY()[i] * _fErrHigh;
      _ratio->GetEYlow()[i]  = _ratio->GetY()[i] * _fErrLow;
    }
  }
  return _ratio;
}

void nicePlot::graphAddSyst(TGraphAsymmErrors* _gSyst, TGraphAsymmErrors* _gNominal, bool _symmetrise) {
  TGraphAsymmErrors* _g = dataSyst.back(); // What we're updating
  if (_gNominal == nullptr) _gNominal = _g; // Syst is w.r.t. the nominal data
  if (_gSyst->GetN() != _gNominal->GetN()) { err(TString("different numbers of bins in the data! ") + _gSyst->GetName() + " " + _gNominal->GetName()); return; }
  for (int i=0; i < _gNominal->GetN(); i++) {
    double _diff = _gSyst->GetY()[i] - _gNominal->GetY()[i];
    double _diffFrac = _diff / _gNominal->GetY()[i];
    if (std::isinf(_diff) || std::isnan(_diff)) { 
      err("Encountered a inf or nan when adding systematics");
      continue;
    } 
    if (_symmetrise || _diff > 0) {
      double _fracErr = _g->GetEYhigh()[i] / _g->GetY()[i];
      _fracErr = std::sqrt( (_fracErr * _fracErr) + (_diffFrac * _diffFrac) );
      _g->GetEYhigh()[i] = _g->GetY()[i] * _fracErr;
    }
    if (_symmetrise || _diff < 0) {
      double _fracErr = _g->GetEYlow()[i] / _g->GetY()[i];
      _fracErr = std::sqrt( (_fracErr * _fracErr) + (_diffFrac * _diffFrac) );
      _g->GetEYlow()[i] = _g->GetY()[i] * _fracErr;
    }
  }
}

nicePlot* nicePlot::addRatio(TH1* _mc, TH1* _data, int mcColourOffset) {
  if (_mc == nullptr || _data == nullptr) return this;
  TH1* _clone1 = static_cast<TH1*>(_mc->Clone());
  TH1* _clone2 = static_cast<TH1*>(_data->Clone());
  applyOptionsToHist(_clone1);
  applyOptionsToHist(_clone2);
  TGraphAsymmErrors* _mc_g = new TGraphAsymmErrors(_clone1);
  TGraphAsymmErrors* _data_g = new TGraphAsymmErrors(_clone2);
  int _n = n_mc;
  _n -= mcColourOffset;
  _mc_g->SetLineColor(mc_col[_n]);
  _mc_g->SetMarkerColor(0);
  _mc_g->SetLineStyle(mc_sty[_n]);
  _mc_g->SetLineWidth(line_width);
  _mc_g->SetMarkerSize(0);
  return addRatio(_mc_g, _data_g);
}

nicePlot* nicePlot::addRatio(TGraphAsymmErrors* _mc, TGraphAsymmErrors* _data, bool isDataStyle) {
  if (_mc == 0 || _data == 0) return this;
  ratioIsDataStyle.push_back(isDataStyle);
  ratio.push_back( graphDivide(_mc, _data, includeBothInRatioError) );
  return this;
}

nicePlot* nicePlot::addRatio(TFile* _f, TString _hist_name_in_file_1, TString _hist_name_in_file_2, int mcColourOffset) {
  return addRatio(_f, _hist_name_in_file_1, _f, _hist_name_in_file_2, mcColourOffset);
}

nicePlot* nicePlot::addRatio(TFile* _f1, TString _hist_name_in_file_1, TFile* _f2, TString _hist_name_in_file_2, int mcColourOffset) {
  TObject* _o1 = _f1->Get(_hist_name_in_file_1);
  TObject* _o2 = _f2->Get(_hist_name_in_file_2);
  if ( _o1 == nullptr ) { err("Cannot find " + _hist_name_in_file_1); return this; }
  if ( _o2 == nullptr ) { err("Cannot find " + _hist_name_in_file_2); return this; }

  if (_o1->ClassName() != _o2->ClassName()) { err("With manual ratio, both hists must currently be of the same type"); return this; }

  if (getObjectType(_o1) == kObjTypeGraph) {
    if (mcColourOffset != 0) err("mcColourOffset not supported yet for addRatio with two tgraphs");
    addRatio( static_cast<TGraphAsymmErrors*>( _o1 ), static_cast<TGraphAsymmErrors*>( _o2 ) );
  } else if (getObjectType(_o1) == kObjTypeTH1) {
    addRatio( static_cast<TH1*>( _o1 ), static_cast<TH1*>( _o2 ), mcColourOffset );
  } else err("unsuported type");
  return this;
}


nicePlot* nicePlot::addMCMCRatio(TFile* _f, TString _hist_name_in_file_1, TString _hist_name_in_file_2, TString _name, bool _ratioWithData) {
  return addMCMCRatio(_f, _hist_name_in_file_1, _f, _hist_name_in_file_2, _name, _ratioWithData);
}

// Add the ratio of TWO MC lines
nicePlot* nicePlot::addMCMCRatio(TFile* _f1, TString _hist_name_in_file_1, TFile* _f2, TString _hist_name_in_file_2, TString _name, bool _ratioWithData) {
  TObject* _o1 = _f1->Get(_hist_name_in_file_1);
  TObject* _o2 = _f2->Get(_hist_name_in_file_2);
  if ( _o1 == nullptr) { err("Cannot find " + _hist_name_in_file_1); return this;}
  if ( _o2 == nullptr) { err("Cannot find " + _hist_name_in_file_2); return this;}

  if ( _o1->ClassName() != _o2->ClassName() ) { err("type mismatch in addMCCMCRatio"); return this; }
  if (_name == "") _name = _hist_name_in_file_1;

  if ( getObjectType(_o1) == kObjTypeTH1 ) {
    TH1* clone1 = static_cast<TH1*>( _o1->Clone() );
    TH1* clone2 = static_cast<TH1*>( _o2->Clone() );
    applyOptionsToHist(clone1);
    applyOptionsToHist(clone2);
    clone1->Divide( clone2 );
    normAlreadDone = true;
    addMC(clone1, _name, _ratioWithData);
  } else if ( getObjectType(_o1) == kObjTypeGraph ) {
    TGraphAsymmErrors* clone = static_cast<TGraphAsymmErrors*>( _o1->Clone() );
    addMC(graphDivide(clone, static_cast<TGraphAsymmErrors*>(_o2), includeBothInRatioError), _name, _ratioWithData);
  } else err("unsuported type");
  return this;
}

nicePlot* nicePlot::addMCMCRatio(TGraphAsymmErrors* _a, TGraphAsymmErrors* _b, TString _name, bool _ratioWithData) {
  if ( _a == nullptr || _b == nullptr) { err("nullptr in addMCMCRatio"); return this; }
  return addMC(graphDivide(_a, _b, includeBothInRatioError), _name, _ratioWithData);
}

nicePlot* nicePlot::addDataDataRatio(TGraphAsymmErrors* _a, TGraphAsymmErrors* _b, TString _name) {
  if ( _a == nullptr || _b == nullptr) { err("nullptr in addDataDataRatio"); return this; }
  return addData(graphDivide(_a, _b, includeBothInRatioError), _name);
}

void nicePlot::applyOptionsToHist(TH1* _h, float _norm) {
  if (normAlreadDone) {
    normAlreadDone = false;
    return;
  }
  if (rebin > 0)         _h->Rebin(rebin);
  if (divideBinWidth)    _h->Scale( 1., "width");
  if (norm_command == 1) {
    int binLow = 1;
    int binHigh = _h->GetNbinsX();
    if(norm_lower != -1) binLow  = norm_ndc ? _h->FindBin(norm_lower) : norm_lower;
    if(norm_upper != -1) binHigh = norm_ndc ? _h->FindBin(norm_upper) : norm_upper;
    _h->Scale( 1./_h->Integral( binLow, binHigh ) );
  }
  if (_norm != 0)        _h->Scale( 1./_norm );
}

void nicePlot::addLable(double x, double y, TString text, double scale) {
  plotText.push_back(text);
  plotTextScale.push_back(scale);
  plotTextX.push_back(x);
  plotTextY.push_back(y);
}

void nicePlot::addBinLabel(TString text, TString axis) {
  if (axis == "x") bin_label_x.push_back(text);
  else if (axis == "y") bin_label_y.push_back(text);
  else err("wrong axis");
}

void nicePlot::setLegend(double x, double y, double scale) {
  legendX = x;
  legendY = y;
  legendScale = scale;
}

void nicePlot::setStamp(double x, double y, TString _stamp) {
  stampX = x;
  stampY = y;
  stamp = _stamp;
}

void nicePlot::setColab(double x, double y, TString _colab) {
  stampX = x;
  stampY = y;
  colab = _colab;
  if (stamp == TString("NONE")) stamp = "";
}

void nicePlot::scaleLastMC(double factor) {
  for (int i=0;i<mc.back()->GetN();i++) mc.back()->GetY()[i]      *= factor;
  for (int i=0;i<mc.back()->GetN();i++) mc.back()->GetEYhigh()[i] *= factor;
  for (int i=0;i<mc.back()->GetN();i++) mc.back()->GetEYlow()[i]  *= factor;
}

void nicePlot::scaleLastData(double factor) {
  for (int i=0;i<dataSyst.back()->GetN();i++)  {
    dataSyst.back()->GetY()[i] *= factor;
    dataSyst.back()->GetEYhigh()[i] *= factor;
    dataSyst.back()->GetEYlow()[i] *= factor;

    dataStat.back()->GetY()[i] *= factor;
    dataStat.back()->GetEYhigh()[i] *= factor;
    dataStat.back()->GetEYlow()[i] *= factor;
  }
}

void nicePlot::drawArrow(int i) {
  double _x2, _x1, _y2, _y1;
  _x2 = _x1 = (lineX1.at(i) + lineX2.at(i)) / 2.;
  _y2 = _y1 = (lineY1.at(i) + lineY2.at(i)) / 2.;
  switch (linePoint.at(i)) {
    default:
    case kNPArrowNone: return;
    case kNPArrowRight: _x2 = _x1 + ((xMax-xMin)/10.); break; // right
    case kNPArrowLeft:  _x2 = _x1 - ((xMax-xMin)/10.); break; // left
    case kNPArrowUp:    _y2 = _y1 + ((yMax-yMin)/10.); break; // up
    case kNPArrowDown:  _y2 = _y1 - ((yMax-yMin)/10.);  break; // down
  }
  TArrow* ar = new TArrow(_x1, _y1, _x2, _y2, 0.05, "|>");
  // std::cout << "x1:" << _x1 << " x2:" << _x2 << " y1:" << _y1 << " y2:" << _y2 << std::endl;
  ar->SetAngle(30);
  ar->SetLineWidth(lineWidth.at(i));
  ar->SetLineStyle(lineStyle.at(i));
  ar->SetFillStyle(0);
  if (plot2D != nullptr) { ar->SetLineColor(0); }
  ar->Draw();
}

TCanvas* nicePlot::getCanvas(TPad* _toDrawInto) {

  // if (dataStat.size() == 1) dataStat[0]->SetMarkerColor(kBlack);

  if (c != nullptr && _toDrawInto == nullptr) return c;

  RATIO = true;
  if (ratio.size() == 0) {
    RATIO = false;
  }

  // CANVAS
  if (_toDrawInto == nullptr) {
    c = new TCanvas();
    c->SetCanvasSize(800 * (stretch ? 2.5 : 1),600); // WC MC
    gPad->SetMargin(0,0,0,0);
    c->SetFillColor(1); // WC MC
    gPad->SetFillColor(1); // WC MC
    _toDrawInto = (TPad*) gPad; 
  }

  double _smallY = 1e-5, _smallR = 1e-5;
  if (logY || !autoY) _smallY = 0.;
  if (logR) _smallR = 0.;

  // AUTO X AXIS
  if (autoX && plot2D != nullptr) {
    xMin = plot2D->GetXaxis()->GetXmin();
    xMax = plot2D->GetXaxis()->GetXmax();
  } else if (autoX) {
    double minX = 9e99, maxX = -9e99, xRet, yRet, xErr;
    for (const TGraphAsymmErrors* g : dataSyst) { // DATA 
      for (int p=0; p < g->GetN(); ++p) {
        g->GetPoint(p, xRet, yRet);
        xErr = std::max(g->GetErrorXlow(p), g->GetErrorXhigh(p));
        if (xRet > maxX + xErr) maxX = xRet + xErr;
        if (xRet < minX - xErr) minX = xRet - xErr;
      }
    }
    for (const TGraphAsymmErrors* g : mc) { // MC 
      for (int p=0; p < g->GetN(); ++p) {
        g->GetPoint(p, xRet, yRet);
        xErr = std::max(g->GetErrorXlow(p), g->GetErrorXhigh(p));
        if (xRet > maxX + xErr) maxX = xRet + xErr;
        if (xRet < minX - xErr) minX = xRet - xErr;
      }
    }
    xMin = minX;
    xMax = maxX;
  }

  // AUTO Y AXIS
  if (autoY && plot2D != nullptr) {

    yMin = plot2D->GetYaxis()->GetXmin();
    yMax = plot2D->GetYaxis()->GetXmax();

    rMin = plot2D->GetMinimum();
    rMax = plot2D->GetMaximum();

    if (logZ && fabs(rMin) < 1e-15) {
      rMin = 1e-3;
    }

  } else if (autoY) {
    double minV = 9e99, maxV = -9e99, minVR = 9e99, maxVR = -9e99, xRet, yRet;

    if (yMin == 0 && yMax == 0) {
      for (unsigned int i=0; i < dataSyst.size(); ++i) { // DATA - Y
        for (int p=0; p < dataSyst[i]->GetN(); ++p) {
          dataSyst[i]->GetPoint(p, xRet, yRet);
          if (xRet < xMin || xRet > xMax) continue; // Out of range?
          if (yRet > maxV) maxV = yRet;
          if (yRet < minV && fabs(yRet) > 1e-15) minV = yRet;
        }
      }
      for (unsigned int i=0; i < mc.size(); ++i) { // MC - Y
        for (int p=0; p < mc[i]->GetN(); ++p) {
          mc[i]->GetPoint(p, xRet, yRet);
          if (xRet < xMin || xRet > xMax) continue; // Out of range?
          if (yRet > maxV) maxV = yRet;
          if (yRet < minV && fabs(yRet) > 1e-15) minV = yRet;
        }
      }
            
      if ( fabs(minV) < 1e-15 && logY == false) yMin = 1e-15;
      else if (logY) yMin = minV * 0.5;
      else           yMin = minV * 0.95;

      if (logY) yMax = maxV * ySpaceModLog;
      else      yMax = maxV * ySpaceModLin;
    }

    if (RATIO && rMin == 0 && rMax == 0) {
      for (unsigned int i=0; i < ratio.size(); ++i) { // RATIO - Y
        for (int p=0; p < ratio[i]->GetN(); ++p) {
          if (ratio[i]->GetPoint(p, xRet, yRet) == -1) continue; // Out of range?
          if (xRet <= xMin || xRet >= xMax) continue;
          if (!std::isinf((double)yRet) && yRet > maxVR) maxVR = yRet;
          if (!std::isinf((double)yRet) && fabs(yRet) > 1e-15 && yRet < minVR) minVR = yRet;
          if (ratioOnly) {
            if (!std::isinf((double)(1./yRet)) && 1./yRet > maxVR) maxVR = 1./yRet;
            if (!std::isinf((double)(1./yRet)) && fabs(1./yRet) > 1e-15 && yRet < minVR) minVR = 1./yRet;
          }
        }
      }
      rMin = minVR * 0.95;
      if (rMin > 1.) rMin = 0.5; // keep 1 on the axis
      rMax = maxVR * 1.05;
      // if (logR == false) logR = (rMax >= 10.); // do not force this
    } 

  } // END AUTO Y

  TPad* _p1 = nullptr;
  TPad* _p2 = nullptr;

  double xLow = 0., yLow = 0, xUp = 1., yUp = 1.;
  double yOffset = 0.35;
  if (ratioOnly) yOffset = 1.;

  double marginLeft = 0.162 * (stretch ? 1./2.5 : 1), marginRight = 0.04 * (stretch ? 1./2.5 : 1), marginTop = 0.06, marginBottom = 0.162;
  double labelSize = 0.05, titleOffset = (1.2 + titleOffsetMod) * (stretch ? 1./2.5 : 1);

  if (RATIO) { // WITH RATIO



    // TOP FRAME
    if (!ratioOnly) {
      _p1 = new TPad("","", xLow, yOffset, xUp, yUp);
      _p1->Draw();
      _p1->SetMargin(marginLeft, marginRight, 0, marginTop * (1./(1.-yOffset)));
      _p1->SetLogx(logX);
      _p1->SetLogy(logY);
      _p1->SetFillStyle(0);
      _p1->SetFillColor(1);
    }

    // BOT FRAME
    _p2 = new TPad("","",xLow, yLow, xUp, yOffset);
    _p2->Draw();
    _p2->SetMargin(marginLeft, marginRight, marginBottom * (1./yOffset), 0);
    if (ratioOnly) _p2->SetTopMargin(marginTop);
    _p2->SetLogx(logX);
    _p2->SetLogy(logR);
    _p2->SetFillStyle(0);
    _p2->SetFillColor(1);

    //std::cout <<" logR "<< logR << " | ";

    _p2->cd();
    _frame_bot = gPad->DrawFrame(xMin, rMin - _smallR, xMax, rMax - _smallR, "");
    _frame_bot->GetYaxis()->SetNdivisions(5,2,0,kTRUE);
    _frame_bot->GetXaxis()->SetMoreLogLabels(kTRUE);
    _frame_bot->GetXaxis()->SetTitle(xTitle);
    _frame_bot->GetYaxis()->SetTitle(rTitle);
    _frame_bot->GetXaxis()->SetLabelSize(labelSize * (1./yOffset));
    _frame_bot->GetXaxis()->SetTitleSize(labelSize * (1./yOffset));
    _frame_bot->GetYaxis()->SetLabelSize(labelSize * (1./yOffset));
    _frame_bot->GetYaxis()->SetTitleSize(labelSize * (1./yOffset));
    _frame_bot->GetYaxis()->SetTitleOffset(titleOffset * yOffset);



    _frame_bot->GetXaxis()->SetTitleColor(0);  // WC MC
    _frame_bot->GetYaxis()->SetTitleColor(0);

    _frame_bot->GetYaxis()->SetAxisColor(0);
    _frame_bot->GetXaxis()->SetAxisColor(0);

    _frame_bot->GetYaxis()->SetLabelColor(0);
    _frame_bot->GetXaxis()->SetLabelColor(0);

    // double extraY = 0;
    // if (extraYoffset > 0) extraY = extraYoffset/100.; // safe w super scripts

    if (!ratioOnly) {
      _p1->cd();
      _frame_top = gPad->DrawFrame(xMin, yMin + _smallY, xMax, yMax, "");
      _frame_top->GetYaxis()->SetTitle(yTitle);
      _frame_top->GetYaxis()->SetLabelSize(labelSize * (1./(1.-yOffset)));
      _frame_top->GetYaxis()->SetTitleSize(labelSize * (1./(1.-yOffset)));
      _frame_top->GetYaxis()->SetTitleOffset(titleOffset * (1.-yOffset));
      _frame_top->GetYaxis()->SetTitleColor(0);  // WC MC

      _frame_top->GetYaxis()->SetAxisColor(0);
      _frame_top->GetXaxis()->SetAxisColor(0);

      _frame_top->GetYaxis()->SetLabelColor(0);
     }

  } else if (plot2D != nullptr) { // 2D

    TPad* _p1 = new TPad("","",0, 0, 1, 1);
    _p1->SetFillStyle(0);
    _p1->Draw();
    _p1->SetLogx(logX);
    _p1->SetLogy(logY);
    _p1->SetLogz(logZ);
    _p1->SetMargin(marginLeft, marginLeft, marginBottom, marginTop);
    _p1->cd();

    plot2D->GetZaxis()->SetTitle(rTitle);
    plot2D->GetYaxis()->SetTitle(yTitle);
    plot2D->GetXaxis()->SetTitle(xTitle);

    plot2D->GetZaxis()->SetTitleColor(0);  // WC MC
    plot2D->GetYaxis()->SetTitleColor(0);
    plot2D->GetXaxis()->SetTitleColor(0);

    plot2D->GetZaxis()->SetAxisColor(0);
    plot2D->GetYaxis()->SetAxisColor(0);
    plot2D->GetXaxis()->SetAxisColor(0);

    plot2D->GetZaxis()->SetLabelColor(0);
    plot2D->GetYaxis()->SetLabelColor(0);
    plot2D->GetXaxis()->SetLabelColor(0);

    if ( bin_label_x.size() > 0 ) setBinLabels( plot2D, "x");
    if ( bin_label_y.size() > 0 ) setBinLabels( plot2D, "y");

    //gStyle->SetPalette(53); //53
    // kCool
    static bool doOnce = true;
    if (doOnce == true) {
      doOnce = false;
      Double_t stops[9] = { 0.0000, 0.1250, 0.2500, 0.3750, 0.5000, 0.6250, 0.7500, 0.8750, 1.0000};
      Double_t red[9]   = { 0/255., 51./255.,  43./255.,  33./255.,  28./255.,  35./255.,  74./255., 144./255., 246./255.};
      Double_t green[9] = { 0/255., 24./255.,  55./255.,  87./255., 118./255., 150./255., 180./255., 200./255., 222./255.};
      Double_t blue[9]  = { 0/255., 96./255., 112./255., 114./255., 112./255., 101./255.,  72./255.,  35./255.,   0./255.};
      TColor::CreateGradientColorTable(9, stops, red, green, blue, 255);
    }
    //
    plot2D->Draw("colz");
    plot2D->GetZaxis()->SetRangeUser(rMin, rMax);
    plot2D->GetXaxis()->SetMoreLogLabels(kTRUE);
    plot2D->GetXaxis()->SetLabelSize(labelSize);
    plot2D->GetXaxis()->SetTitleSize(labelSize);
    plot2D->GetYaxis()->SetLabelSize(labelSize);
    plot2D->GetYaxis()->SetTitleSize(labelSize);
    plot2D->GetYaxis()->SetTitleOffset(titleOffset);

  } else { // STANDARD - NO RATIO

    // TOP FRAME
    // TODO extra y offset here too
    TPad* _p1 = new TPad("","",0, 0, 1, 1);
    _p1->SetFillStyle(0);
    _p1->SetFillColor(1);
    _p1->Draw();
    _p1->SetLogx(logX);
    _p1->SetLogy(logY);
    _p1->SetMargin(marginLeft, marginRight, marginBottom, marginTop);
    _p1->cd();

    _frame_top = gPad->DrawFrame(xMin, yMin + _smallY, xMax, yMax, "");

    _frame_top->GetYaxis()->SetTitle(yTitle);
    _frame_top->GetXaxis()->SetTitle(xTitle);
    _frame_top->GetXaxis()->SetMoreLogLabels(kTRUE);
    //_frame_top->GetYaxis()->SetMoreLogLabels(kTRUE);
    _frame_top->GetXaxis()->SetLabelSize(labelSize);
    _frame_top->GetXaxis()->SetTitleSize(labelSize);
    _frame_top->GetYaxis()->SetLabelSize(labelSize);
    _frame_top->GetYaxis()->SetTitleSize(labelSize);
    _frame_top->GetYaxis()->SetTitleOffset(titleOffset);

    _frame_top->GetXaxis()->SetTitleColor(0);  // WC MC
    _frame_top->GetYaxis()->SetTitleColor(0);

    _frame_top->GetYaxis()->SetAxisColor(0);
    _frame_top->GetXaxis()->SetAxisColor(0);

    _frame_top->GetYaxis()->SetLabelColor(0);
    _frame_top->GetXaxis()->SetLabelColor(0);

  }


  if (!ratioOnly) {
    // DATA
    for (unsigned int i=0; i < dataSyst.size(); ++i) dataSyst[i]->Draw("same E2P" + drawOp);
    for (unsigned int i=0; i < dataStat.size(); ++i) dataStat[i]->Draw("same P" + drawOp);
    for (unsigned int i=0; i < mc.size();   ++i) {
      //cout << "drawing " << mc[i] << " " << mc[i]->GetTitle() << endl;
      mc[i]  ->Draw(TString("same P") + mc_do[i] + drawOp);
    }
    if (mc_stack != nullptr)  {
      mc_stack->Draw("same");
      // _frame_top = mc_stack->GetHistogram();
    }

    // FIT 
    if (fit != "") {
      float fit_y = 0.9;
      for (unsigned int i=0; i < mc.size(); ++i) {
        mc[i]->Fit(fit, "S", "", fitMin, fitMax);
        TF1* myfit = mc[i]->GetFunction(fit);
        if (myfit) {
          myfit->SetLineColor(mc_col[i]);
          myfit->SetLineWidth(1);
          myfit->SetLineStyle(3);
          if (print_fit) {
            for (int p=0; p < myfit->GetNpar(); ++p) {
              std::stringstream ss;
              ss << mc[i]->GetTitle() << " P" << p << " = " << myfit->GetParameter(p) << " #pm " << myfit->GetParError(p);
              addLable(0.8, fit_y, ss.str().c_str(), 0.2);
              fit_y -= 0.015;
            }
          }
        }
      }
    }

  }

  // LEGEND
  if (doLegend) {
    double _YSpaceUpper = 0.92;
    double _YEntriesSpacing = 0.04;
    if (RATIO) {
      _YSpaceUpper *= 0.98;
      _YEntriesSpacing = 0.05;
      legendY *= 0.9;
    }
    int _entries = dataSyst.size() + mc.size();
    if(mc_stack != nullptr) _entries += mc_stack->GetHists()->GetEntries();
    if ( legendY + (_YEntriesSpacing * _entries * legendScale) > _YSpaceUpper) {
      legendScale = ((_YSpaceUpper - legendY) / ( _entries * _YEntriesSpacing ));
      err("no room for legend ");
    }
    TLegend* L = new TLegend(legendX, legendY, 1, legendY + (_YEntriesSpacing * _entries * legendScale) , "","NDC");
    L->SetTextFont(43);
    L->SetTextColor(0); // WC MC
    L->SetTextSize(24. * legendScale);
    L->SetFillStyle(0);
    L->SetBorderSize(0);
    for (unsigned int i=0; i < dataSyst.size(); ++i) L->AddEntry(dataSyst[i],"","FPL");
    for (unsigned int i=0; i < mc.size(); ++i)       L->AddEntry(mc[i]  ,"",TString("L") + (mc[i]->GetLineWidth() ? TString("P") : TString("")) );
    if (mc_stack != nullptr) {
      TList* l = mc_stack->GetHists();
      for(const auto&& obj: *l) {
        TH1* h = (TH1*)obj;
        L->AddEntry(h, "", "LF");
      }
    }
    L->Draw("same");
  }

  if (RATIO) {

    // BOTTOM FRAME
    _p2->cd();

    // X AXIS BIN LABELS
    if ( bin_label_x.size() > 0 ) setBinLabels( _frame_bot, "x");

    TLine* _l = new TLine();
    _l->SetLineStyle(3);
    _l->DrawLine(xMin, ratioLineValue, xMax, ratioLineValue);

    // RATIOS
    for (unsigned int i=0; i < ratio.size(); ++i) {
      if (ratioIsDataStyle[i]) ratio[i]->Draw("same E2" + drawOp);
      else ratio[i]->Draw("same P" + drawOp);
    }

    // TIDY UP
    gPad->RedrawAxis();
    if (!ratioOnly) {
      _p1->cd();
      gPad->RedrawAxis();
    }

  } else {

    // X AXIS BIN LABELS
    if ( bin_label_x.size() > 0 ) {
      setBinLabels( _frame_top, "x" );
    }

    // TIDY UP
    gPad->RedrawAxis();

  }

  // LINES
  for (unsigned int i = 0; i < lineX1.size(); ++i) {
    TLine* _l = new TLine();
    _l->SetLineStyle( lineStyle.at(i) );
    _l->SetLineWidth( lineWidth.at(i) );
    if (plot2D != nullptr) _l->SetLineColor(0);
    bool _vertical = false;
    if (lineY1.at(i) == FLT_MAX && lineY2.at(i) == FLT_MAX) { // Set to pad
      lineY1.at(i) = yMin;
      lineY2.at(i) = yMax;
      _vertical = true;
    }
    _l->DrawLine( lineX1.at(i), lineY1.at(i), lineX2.at(i), lineY2.at(i) );
    drawArrow(i); 
    if (_vertical && RATIO) { // Draw in ratio box if vertical
      _p2->cd();
      lineY1.at(i) = rMin;
      lineY2.at(i) = rMax;
      _l->DrawLine( lineX1.at(i), lineY1.at(i), lineX2.at(i), lineY2.at(i) ); 
      _p1->cd();
    }
  }


  // TEXT
  _toDrawInto->cd();
  for (unsigned int i=0; i < plotText.size(); ++i) plotLable(plotTextX[i], plotTextY[i], plotText[i], plotTextScale[i]);
  if (stamp != TString("NONE")) ATLASLabel(stampX, stampY, stamp, 1, 1);

  return c;

}

void nicePlot::setBinLabels(TH1* _frame, TString axis) {
  TGraphAsymmErrors* _binning = nullptr;
  TH1* _binningStack = nullptr;
  if ( mc.size() > 0 ) _binning = mc.back();
  else if ( dataSyst.size() > 0 ) {
    _binning = dataSyst.back();
  } else if ( mc_stack != nullptr) { 
    _binningStack = (TH1*) (mc_stack->GetHists()->First());
  }
  else err("Cannot determin number of bins to apply labels to");
  if (axis == "x") {
    if (_binning != nullptr) {
      _frame->SetBins( _binning->GetN(), xMin, xMax );
    } else if (_binningStack != nullptr) { 
      _frame->SetBins( _binningStack->GetNbinsX(), xMin, xMax );
    }
    else err( "Cannot determin binning");
    _frame->GetXaxis()->SetLabelSize(0.08);
    for (unsigned i=0; i < bin_label_x.size(); ++i) {
      _frame->GetXaxis()->SetBinLabel(i+1, bin_label_x[i]);
    }
  } else if (axis == "y") {
    err("setBinLabels y not implimented");
  }
}

void nicePlot::setBinLabels(TH2* _th2, TString axis) {
  if (axis == "x") {
    for (unsigned i=0; i < bin_label_x.size(); ++i) {
      _th2->GetXaxis()->SetBinLabel(i+1, bin_label_x[i]);
    }
   } else if (axis == "y") {
    for (unsigned i=0; i < bin_label_y.size(); ++i) {
      _th2->GetYaxis()->SetBinLabel(i+1, bin_label_y[i]);
    }
  }
}

void nicePlot::plotLable(Double_t x,Double_t y, TString text, double scale) {
  TLatex l;
  l.SetNDC();
  l.SetTextFont(42);
  l.SetTextSize( l.GetTextSize() * scale);
  Color_t color = 1;
  if (true || plot2D != 0) color = 0; // WC MC
  l.SetTextColor(color);
  l.DrawLatex(x,y,text);
}

void nicePlot::ATLASLabel(Double_t x, Double_t y, const char* text, Color_t color, double scale)  {
  if (plot2D != nullptr) color = 0;
  TLatex l; //l.SetTextAlign(12); 
  l.SetNDC();
  l.SetTextFont(72);
  l.SetTextColor(color);
  l.SetTextSize( l.GetTextSize() * scale);

  double delx = 0.13 * scale;

  l.DrawLatex(x,y,colab);
  if (text) {
    TLatex p;
    p.SetNDC();
    p.SetTextFont(42);
    p.SetTextSize( p.GetTextSize() * scale);
    p.SetTextColor(color);
    p.DrawLatex(x+delx,y,text);
    //    p.DrawLatex(x,y,"#sqrt{s}=900GeV");
  }
}


void nicePlot::debug() {
  for (const TGraphAsymmErrors* g : mc) { // MC 
    for (int p=0; p < g->GetN(); ++p) {
      double xRet, yRet;
      g->GetPoint(p, xRet, yRet);
      std::cout << "p:" << p
        << " x:" << xRet << " + " << g->GetErrorXhigh(p) << " - " << g->GetErrorXlow(p)
        << " y:" << yRet << " + " << g->GetErrorYhigh(p) << " - " << g->GetErrorYlow(p) << std::endl;
    }
  }
}

void bookOutput::doMultipadOutput(TString _name, int _x, int _y) {
  if (_x * _y != (int)(rp.size()/_break)) { err(TString("Invalid dimensions for doMultipadOutput. "+std::to_string(_x)+" * "+std::to_string(_y)+" != "+std::to_string(rp.size())+" / "+std::to_string(_break))); return; }
  std::set<TCanvas*> _multiPads;
  std::cout <<  "\033[1;32m Doing " << _name << " Book MultiPad Output " << std::endl;

  for (unsigned _i = 0; _i < _break; ++_i) {
    TCanvas* _splitCanvas = new TCanvas();
    gPad->SetMargin(0,0,0,0);
    _splitCanvas->SetCanvasSize(_x*800 * (rp.at(0)->stretch ? 2.5 : 1), _y*600);
    _splitCanvas->Divide(_x, _y, 0, 0);
    static int _count = 1;
    for (unsigned _p = 0; _p < (unsigned)(_x*_y); ++_p) {
      TPad* _pad = (TPad*) _splitCanvas->cd(_p+1);
      _pad->SetFillColor(1); // WC MC
      rp.at(_i + (_p * _break))->drawInto(_pad);
    }
    std::cout << "|" << std::flush;
    TString _fname = _name + TString(".pdf");
    if (_break > 1) {
      if (_i == 0)             _fname += TString("(");
      else if (_i == _break-1) _fname += TString(")");
    }
    unsigned old = gErrorIgnoreLevel;
    TString _individualName = "img/" + _name + "_";
    _individualName += _count++;
    gErrorIgnoreLevel = 10000;
    _splitCanvas->Print(_fname,"pdf");
    _splitCanvas->Print(TString(_individualName + ".pdf"),"pdf");
    _splitCanvas->Print(TString(_individualName + ".png"),"png");
    gErrorIgnoreLevel = old;
    _multiPads.insert(_splitCanvas);
  }
  std::cout << "\033[0m\n";
}

void bookOutput::doBookOutput(TString _name) {
  unsigned max = _break;
  std::cout <<  "\033[1;32m Doing " << _name << " Book Output " << std::endl;
  if (rp.size() < _break) max = rp.size();
  int _count = 1;
  for (unsigned i=0; i< max; ++i) {
    std::cout << "|" << std::flush;
    // std::cout << rp.at(i)->xTitle << " "  << rp.at(i)->yTitle << std::endl;
    TString fname = _name + TString(".pdf");
    if (i == 0 && rp.size() > 1) {
      fname += TString("(");
      _count = 1;
    }
    else if (i == max-1 && rp.size() > 1) fname += TString(")");
    TString _individualName = "img/" + _name + "_";
    _individualName += _count++;
    unsigned old = gErrorIgnoreLevel;
    gErrorIgnoreLevel = 10000;
    rp.at(i)->getCanvas()->Print(fname,"pdf");
    rp.at(i)->getCanvas()->Print(TString(_individualName + ".png"),"png");
    rp.at(i)->getCanvas()->Print(TString(_individualName + ".pdf"),"pdf");
    rp.at(i)->getCanvas()->Print(TString(_individualName + ".root"),"root");
    gErrorIgnoreLevel = old;
  }
  std::cout << "\033[0m\n";
  rp.erase(rp.begin(), rp.begin()+max);
}

void bookOutput::clear() {
  rp.erase(rp.begin(), rp.end());
}

bookOutput& bookOutput::get() {
  static bookOutput instance; // Guaranteed to be destroyed.
  return instance;        // Instantiated on first use.
}

void bookOutput::setBreak(int _b) {
  _break = _b;
}

void bookOutput::setBreakNow() {
  if (_break == 9999) _break = rp.size();
}

void bookOutput::book(nicePlot* _rp) {
  rp.push_back( _rp );
}

void bookOutput::unbook(nicePlot* _rp) {
  const auto pos = std::find(rp.begin(), rp.end(), _rp);
  if (pos != rp.end()) rp.erase(pos);
}

void addDivider(TString _text1, TString _text2) {
  nicePlot* h = new nicePlot("","");
  h->setBounds(0, 1, 0, 1);
  h->addLable(0.2,0.5,_text1,2.);
  h->addLable(0.2,0.4,_text2,2.);
}
