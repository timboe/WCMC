#include <TH1.h>
#include <TRandom3.h>
#include <fstream>
#include <iostream>

int main() {

  std::string line;
  std::ifstream odds("wc_2018_odds.txt");
  double v;
  std::vector<double> oddsVec;
  while (odds >> v) {
    oddsVec.push_back(v);
    std:cout << " | Odd:1 in " << v;
  }
  odds.close();
  
  TH1* probDist = new TH1D("", "", 128/2, -50., 130.); 
  TH1* probDistLow = new TH1D("", "", 128/2, -50., 130.); 
  TH1* probDistHgh = new TH1D("", "", 128/2, -50., 130.); 
  
  const double bias = 0.2;
  
  TRandom3 R;
  
  for (int round = 0; round < 100000000; ++round) {
    double winnings = oddsVec.size();
    winnings *= -1;
    double winningsLow = winnings;
    double winningsHgh = winnings;
    for (const double odd : oddsVec) {
      if (R.Rndm() < 1. / odd) winnings += odd;
      if (R.Rndm() < 1. / (odd * (1. + bias))) winningsLow += odd;
      if (R.Rndm() < 1. / (odd * (1. - bias))) winningsHgh += odd;
    } 
    probDist->Fill(winnings);
    probDistLow->Fill(winningsLow);
    probDistHgh->Fill(winningsHgh);
    if (round % 10000000 == 0) std::cout << "Round " << round << " winnings " << winnings << std::endl;
  }
  
  probDist->Draw("H");
  probDistLow->SetLineColor(kRed);
  probDistLow->Draw("SAMEH");
  probDistHgh->SetLineColor(kBlue);
  probDistHgh->Draw("SAMEH");

}


int toyBets() {
  return main();
}
