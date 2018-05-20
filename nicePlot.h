
#include <vector>
#include <cmath>
#include <set>
#include <sstream>
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
#include <THStack.h>
#include <iostream>
#include <TStyle.h>
#include <TSystem.h>
#include <TLine.h>
#include <TError.h>

enum nicePlotArrow_t{kNPArrowNone, kNPArrowRight, kNPArrowLeft, kNPArrowUp, kNPArrowDown};

enum ObjType_t{kObjTypeGraph, kObjTypeTH1, kObjTypeTH2, kUnsuportedType};

class nicePlot;

void addDivider(TString _text1, TString _text2);

class bookOutput {
public:

  static bookOutput& get();
  static void doBookOutput(TString _name);
  static void doMultipadOutput(TString _name, int _x, int _y);
  static void setBreak(int _b);
  static void setBreakNow();
  static void book(nicePlot* _rp);
  static void unbook(nicePlot* _rp);
  static void clear();

private:
  static std::vector<nicePlot*> rp;
  static unsigned _break;

  bookOutput() {};
  bookOutput(bookOutput const&);
  void operator=(bookOutput const&);
};


class nicePlot {
public:

  nicePlot();
  nicePlot(TString _xTitle, TString _yTitle, TString _rTitle = "");
  nicePlot(nicePlot& _copy);
  nicePlot(nicePlot* _copy);
  void reset();
  void init(TString _xTitle, TString _yTitle, TString _rTitle = "", bool quiet = false);
  void book();
  void unbook();
  void useAltColourScheme(int _i = 0);
  void useScaleColourScheme(int _i = 0);
  void cancelNormaliseToOne() { norm_command = -1; }
  void normaliseToOne(int _lower = -1, int _upper = -1, bool _ndc = false) { norm_command = 1; norm_lower = _lower; norm_upper = _upper; norm_ndc = _ndc; }
  void copy( nicePlot& _copy );
  void divideByBinWidth(bool _divideByBinWidth = true) { divideBinWidth = _divideByBinWidth; }
  void setRebin(int r) { rebin = r; }
  void ATLASLabel(Double_t x, Double_t y, const char* text, Color_t color=1, double scale = 1.);
  void setBounds(double _xMin, double _xMax, double _yMin = 0., double _yMax = 0., double _rMin = 0., double _rMax = 0.);
  void setRBounds(double _rMin, double _rMax);
  void setYBounds(double _yMin, double _yMax);
  ObjType_t getObjectType(TObject* _o);
  static TGraphAsymmErrors* graphDivide(TGraphAsymmErrors* _a, TGraphAsymmErrors* _b, bool _includeBothInRatioError = false);
  void applyOptionsToHist(TH1* _h, float _norm = 0);
  void graphAddSyst(TGraphAsymmErrors* _gSyst, TGraphAsymmErrors* _gNominal = nullptr, bool _symmetrise = false);
  void addLine(double _x1, double _y1, double  _x2, double _y2, int _lineStyle = 1, int _lineWidth = 1, nicePlotArrow_t _point = kNPArrowNone);
  nicePlot* addDataDataRatio(TGraphAsymmErrors* _a, TGraphAsymmErrors* _b, TString _name);
  nicePlot* addData(TGraphAsymmErrors* _gStat, TString _name, bool _ratioWithData = false);
  nicePlot* addData(TH1* _hStat, TString _name, float _norm = 0., bool _ratioWithData = false);
  nicePlot* addData(TFile* _f, TString _hist_name_in_file, TString _name, float _norm = 0., bool _ratioWithData = false);
  nicePlot* addDataSystematic(TGraphAsymmErrors* _gSyst, TGraphAsymmErrors* _gNominal = nullptr, bool _symmetrise = true);
  nicePlot* addDataSystematic(TH1* _hSyst, TH1* _hNominal = nullptr, bool _symmetrise = true);
  nicePlot* addDataSystematic(TFile* _f, TString _syst_name_in_file, bool _symmetrise = true);
  nicePlot* addDataSystematic(TFile* _f, TString _syst_name_in_file, TString _nominal_name_in_file, bool _symmetrise = true);
  nicePlot* add2D(TH2* _h);
  nicePlot* add2D(TFile* _f, TString _hist_name_in_file);
  nicePlot* add2DEfficiency(TFile* _f, TString _hist_name_in_file_1, TString _hist_name_in_file_2);
  nicePlot* addMCMCRatio(TGraphAsymmErrors* _a, TGraphAsymmErrors* _b, TString _name, bool _ratioWithData = false);
  nicePlot* addMCMCRatio(TFile* _f, TString _hist_name_in_file_1, TString _hist_name_in_file_2, TString _name = "", bool _ratioWithData = false);
  nicePlot* addMCMCRatio(TFile* _f1, TString _hist_name_in_file_1, TFile* _f2, TString _hist_name_in_file_2, TString _name = "", bool _ratioWithData = false);
  nicePlot* addMC(TGraphAsymmErrors* _g, TString _name, bool _ratioWithData = false);
  nicePlot* addMC(TH1* _h, TString _name, bool _ratioWithData = false, float _norm = 0.);
  nicePlot* addMC(TFile* _f, TString _hist_name_in_file, TString _name, bool _ratioWithData = false, float _norm = 0.);
  nicePlot* addStackMC(TH1* _h, TString _name);
  nicePlot* addRatio(TGraphAsymmErrors* _mc, TGraphAsymmErrors* _data, bool isDataStyle = false);
  nicePlot* addRatio(TH1* _mc, TH1* _data, int mcColourOffset = 0);
  nicePlot* addRatio(TFile* _f, TString _hist_name_in_file_1, TString _hist_name_in_file_2, int mcColourOffset = 0);
  nicePlot* addRatio(TFile* _f1, TString _hist_name_in_file_1, TFile* _f2, TString _hist_name_in_file_2, int mcColourOffset = 0);
  void addLable(double x, double y, TString text, double scale = 1.);
  void addBinLabel(TString text, TString axis="x");
  void setLegend(double x, double y, double scale = 1);
  void setStamp(double x, double y, TString _stamp);
  void setColab(double x, double y, TString _colab);
  void setLogx(bool log) { logX = log; }
  void setLogy(bool log) { logY = log; }
  void setLogz(bool log) { logZ = log; }
  void setLogr(bool log) { logR = log; }
  void setLogxy(bool x = 0, bool y = 0) { logX = x; logY = y; }
  void setLogxyr(bool x = 0, bool y = 0, bool r = 0) { logX = x; logY = y; logR = r; }
  void setLogxyz(bool x = 0, bool y = 0, bool z = 0) { logX = x; logY = y; logZ = z; }
  void setAutoy(bool _autoY = 1) { autoY = _autoY; }
  void setAutox(bool _autoX = 1) { autoX = _autoX; }
  void setAutoxy(bool _autoX = 1, bool _autoY = 1) { autoX = _autoX; autoY = _autoY; }
  void scaleLastMC(double factor);
  TGraphAsymmErrors* getLastMC() { return mc.back(); }
  void addLastMCDrawOption(TString _do) { mc_do.back() = _do; }
  void setDrawOption(TString _do) { drawOp = _do; }
  int getNMC() { return mc.size(); }
  TGraphAsymmErrors* getMC(int _loc) { return mc.at(_loc); }
  void scaleLastData(double factor);
  void drawArrow(int i);
  void drawInto(TPad* _toDrawInto) { getCanvas(_toDrawInto); }
  TCanvas* getCanvas(TPad* _toDrawInto = nullptr);
  void setBinLabels(TH1* _frame, TString axis = "x");
  void setBinLabels(TH2* _th2, TString axis = "x");
  void plotLable(Double_t x,Double_t y, TString text, double scale);
  void setRatioOnly(bool _r) { ratioOnly = _r; }
  void setIncludeBothInRatioError(bool _i) { includeBothInRatioError = _i; }
  void debug();
  void setKillLowStatBins(bool _k = true) { killLowStatBins = _k; }
  void setLineWidth(int lw) { line_width = lw; }
  void setMCMarkerSize(int s) { mc_marker_size = s; }
  void setDoLegend(bool dl) { doLegend = dl; }
  void setFit(TString fitFn, double fMin = 0., double fMax = 0.) { fit = fitFn; fitMin = fMin; fitMax = fMax; }
  void setPrintFit(bool p) { print_fit = p; }
  void setTitleOffsetMod(double mod) { titleOffsetMod = mod; } // Tweek y axis
  void setYSpaceMod(double lin = 1.2, double log = 2.) { ySpaceModLin = lin; ySpaceModLog = log; }
  void setRatioLineValue(double v = 1.0) { ratioLineValue = v; }

  TH1* _frame_top;
  TH1* _frame_bot;

  //private:


  int n_mc;
  static const int n_sty = 7;
  int mc_col[n_sty];
  int mc_sty[n_sty];
  int data_col[n_sty];
  int data_fill[n_sty];
  int data_sty[n_sty];
  int line_width;
  int mc_marker_size;
  int norm_command, norm_lower, norm_upper;
  double xMin, xMax, yMin, yMax, rMin, rMax, legendX, legendY, legendScale, titleOffsetMod;
  double ySpaceModLin, ySpaceModLog, ratioLineValue;
  bool logX, logY, logZ, logR, norm_ndc, killLowStatBins;
  bool firstRatioIsData;
  bool includeBothInRatioError;
  bool RATIO;
  bool isReg;
  bool divideBinWidth;
  bool autoX, autoY;
  bool ratioOnly;
  bool doLegend;
  bool normAlreadDone;
  bool stretch;

  TString xTitle, yTitle, rTitle;

  std::vector<TGraphAsymmErrors*> dataSyst;
  std::vector<TGraphAsymmErrors*> dataStat;
  std::vector<TGraphAsymmErrors*> mc;
  std::vector<TGraphAsymmErrors*> ratio;
  std::vector<bool> ratioIsDataStyle;
  std::vector<TString> bin_label_x;
  std::vector<TString> bin_label_y;
  std::vector<TString> mc_do; //Draw Option
  TH2* plot2D;

  TGraphAsymmErrors _aGraph;
  TH1F _aFHist;
  TH1D _aDHist;
  TH2F _aF2Hist;
  TH2D _aD2Hist;
  TProfile _aProfile;
  TCanvas* c;
  THStack* mc_stack = nullptr;

  std::vector<TString> plotText;
  std::vector<double> plotTextScale;
  std::vector<double> plotTextX;
  std::vector<double> plotTextY;

  std::vector<double> lineX1;
  std::vector<double> lineX2;
  std::vector<double> lineY1;
  std::vector<double> lineY2;
  std::vector<int> lineStyle;
  std::vector<int> lineWidth;
  std::vector<nicePlotArrow_t> linePoint;

  TString stamp;
  TString colab;
  TString drawOp;
  double stampX, stampY;

  int rebin;
  int extraYoffset;

  TString fit;
  bool print_fit;
  double fitMin, fitMax;

};
