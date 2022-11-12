#include <map>
#include <string>
#include <iostream>
#include <fstream>
#include "nicePlot.cxx"
#include <sstream>
#include <TRandom3.h>
#include <TROOT.h>
#include <TH2.h>

std::vector<std::string> readLine(const std::string& line) {
  std::istringstream buf(line);
  std::istream_iterator<std::string> beg(buf), end;
  std::vector<std::string> results(beg, end);
  for (size_t i = 0; i < results.size(); ++i) {
    size_t location = results[i].find("-");
    if (location != std::string::npos) results[i].replace(location, 1, " ");
  }
  return results;
}

int main() {

  std::string line;
  std::ifstream odds("wc_2018_odds.txt");
  double winnings = 0;
  std::vector<double> oddsVec;

 TH1* oddsHist = new TH1D("", "", 18, 3.75, 12.75); 

  int game = 0;
  while ( getline(odds, line) ) {
    std::vector<std::string> results = readLine(line);
    if (results[0] == "#") {
      std::cout << line << std::endl;
      continue; // Comment
    } else if (results[0] == "##") break;
    oddsVec.push_back( stof(results[0]) );
    oddsHist->Fill(stof(results[0]));
    if (results[1] == "W") {
      winnings += stof(results[0]);
    }
    std:cout << "Game " << ++game << " | Odd:1 in " << results[0] << " (" << results[1] 
      << "). Winnings Status:" << winnings << (results[1] == "W" ? "         !!!" : "") << std::endl;
  }
  odds.close();

  std::cout << std::endl << "Total winnings:" << winnings << std::endl;
  
  TH1* probDist = new TH1D("", "", 180*2, 0., 180.); 
  TH1* probDistLow = new TH1D("", "", 180*2, 0., 180.); 
  TH1* probDistHgh = new TH1D("", "", 180*2, 0., 180.); 

 
  const double bias = 0.2;
  
  TRandom3 R;
  
  for (int round = 0; round < 100000000; ++round) {
    double winnings = 0, winningsLow = 0, winningsHgh = 0;
    for (const double odd : oddsVec) {
      if (R.Rndm() < 1. / odd) winnings += odd;
      if (R.Rndm() < 1. / (odd * (1. + bias))) winningsLow += odd;
      if (R.Rndm() < 1. / (odd * (1. - bias))) winningsHgh += odd;
    } 
    probDist->Fill(winnings - .01); 
    probDistLow->Fill(winningsLow - .01);
    probDistHgh->Fill(winningsHgh - .01);
    if (round % 10000000 == 0) std::cout << "Round " << round << " winnings " << winnings << std::endl;
  }


  nicePlot* npOdds = new nicePlot();
  npOdds->setAutox(1);
  npOdds->setYBounds(0, 25);
  npOdds->setLineWidth(4);
  npOdds->init("Odds (Exact score)", "Games");
  npOdds->setDoLegend(false);
  npOdds->addMC(oddsHist, "");
  npOdds->addMC(oddsHist, "");
  npOdds->addMC(oddsHist, ""); // Stack for pretty


  nicePlot* np = new nicePlot();
  np->setDrawOption("C");
  np->useAltColourScheme(1);
  np->addLine(oddsVec.size(), FLT_MAX, oddsVec.size(), FLT_MAX, 2, 3);
  np->addLine(73.5, FLT_MAX, 73.5, FLT_MAX, 2, 3);
  np->setFit("gaus", 0., 180.);
  np->setPrintFit(false);
  np->setAutoxy(1,1);
  np->setLineWidth(4);
  np->init("Return [GBP]", "1/N_{Rounds}");
  np->setLegend(.7, .7);
  np->normaliseToOne();
  np->addMC(probDistHgh, "Generous");
  np->addMC(probDist, "Fair");
  np->addMC(probDistLow, "Unfair");

  bookOutput::get().doBookOutput("toyBets");

}


int toyBets() {
  return main();
}

/*
| FCN=52401.2 FROM MIGRAD    STATUS=CONVERGED      88 CALLS          89 TOTAL
                     EDM=2.8995e-10    STRATEGY= 1      ERROR MATRIX ACCURATE
  EXT PARAMETER                                   STEP         FIRST
  NO.   NAME      VALUE            ERROR          SIZE      DERIVATIVE
   1  Constant     9.25655e-03   2.84259e-06   3.05878e-07   3.01205e+00
   2  Mean         8.02281e+01   1.74064e-02   1.66820e-03  -4.66960e-04
   3  Sigma        2.20319e+01   8.36931e-03   1.25516e-05   1.35906e-01
 FCN=49544.2 FROM MIGRAD    STATUS=CONVERGED      87 CALLS          88 TOTAL
                     EDM=2.49001e-08    STRATEGY= 1      ERROR MATRIX ACCURATE
  EXT PARAMETER                                   STEP         FIRST
  NO.   NAME      VALUE            ERROR          SIZE      DERIVATIVE
   1  Constant     1.01402e-02   3.37147e-06   3.42197e-07  -2.08050e+01
   2  Mean         6.37412e+01   1.88615e-02   1.66042e-03  -1.03735e-02
   3  Sigma        2.06054e+01   8.83148e-03   1.36570e-05   4.04369e-01
 FCN=54496.2 FROM MIGRAD    STATUS=CONVERGED      96 CALLS          97 TOTAL
                     EDM=4.8664e-09    STRATEGY= 1      ERROR MATRIX ACCURATE
  EXT PARAMETER                                   STEP         FIRST
  NO.   NAME      VALUE            ERROR          SIZE      DERIVATIVE
   1  Constant     1.08763e-02   4.13000e-06   3.96943e-07   2.19596e+01
   2  Mean         5.29768e+01   2.09294e-02   1.75128e-03   4.02066e-04
   3  Sigma        1.93999e+01   9.00018e-03   1.56167e-05  -3.61131e-01
*/