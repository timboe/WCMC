#include <map>
#include <string>
#include <iostream>
#include <fstream>
#include <sstream>
#include <TRandom3.h>
#include <TROOT.h>
#include <TH2.h>
#include "nicePlot.cxx"

class WCMC {  
  public:
    WCMC();
    void doMatch(const std::string& a, const std::string& b, const float low, const float high);
    void doGroup(const std::string& group, const float low, const float high);
    void addHistoric(int goalsA, int goalsB);
    void loadHistoricData();
    void addTeams();
    void addTeam(const std::string& t, const std::string& abreviation);
    void addGroups();
    void addGroup(const std::string& group, const std::string& A, const std::string& B, const std::string& C, const std::string& D);
    void resetTeamStatistics(const bool all);
    void resetLaterGroups();
    void runTraining(float& resultLow, float& resultHigh, const float startLow, const float stopLow, const float startHigh, const float stopHigh, const float step);
    const std::string getWinningTeam(const std::string group);
    void runFinal(const float goalinessLow, const float goalinessHigh);
    void recordStats(const std::string a, const std::string b, const int goalsA, const int goalsB);
    void execute();

    struct Team {
      Team() { m_points = 0; m_goalDiff = 0; m_goals = 0; }
      int m_rank;
      int m_points;
      int m_goalDiff;
      int m_goals;
      std::string m_abreviation;
    };

    std::map<std::string, Team > m_teams;
    std::map<std::string, std::vector<std::string>> m_groups;   
    TRandom3 R;
    TH1F* m_h_GoalsMC;
    TH1F* m_h_GoalsData;
    TH1F* m_h_GoalDiffMC;
    TH1F* m_h_GoalDiffData;
    std::map<std::string, TH2F*> m_h_matchResult;
    std::map<std::string, TH1F*> m_h_roundWinner;
    int m_trialsMax;
    bool m_matchPrint, m_matchStats, m_goalsScored;
    const std::vector<std::string> groups {"A", "B", "C", "D", "E", "F", "G", "H"};
    std::vector<std::string> m_teamsByRank;
    float m_bestChiG, m_bestChiGD;
};

void WCMC::doMatch(const std::string& a, const std::string& b, const float low, const float high) {
  const float reduction = m_teams.size() / high;

  float scoreA = low + ((m_teams.size() - m_teams[a].m_rank) / reduction);
  float scoreB = low + ((m_teams.size() - m_teams[b].m_rank) / reduction);
  while ( scoreA > low && scoreB > low) {
    scoreA -= R.Rndm() * low;
    scoreB -= R.Rndm() * low;
  }

  const int goalsA = R.Poisson(scoreA);
  const int goalsB = R.Poisson(scoreB);

  m_h_GoalsMC->Fill(goalsA + goalsB);
  m_h_GoalDiffMC->Fill( abs(goalsA - goalsB) );

  if (goalsA > goalsB) {
    m_teams[a].m_points += 3;
  } else if (goalsB > goalsA) {
    m_teams[b].m_points += 3;
  } else if (goalsA == goalsB) {
    m_teams[a].m_points += 1;
    m_teams[b].m_points += 1;
  }

  m_teams[a].m_goals += goalsA;
  m_teams[b].m_goals += goalsB;
  m_teams[a].m_goalDiff += goalsA - goalsB;
  m_teams[b].m_goalDiff += goalsB - goalsA;
  
  if (m_matchPrint) std::cout << a << ":" << goalsA << " - " << b << ":" << goalsB << " | "; 
  if (m_matchStats) recordStats(a, b, goalsA, goalsB);
  if (m_goalsScored) {
    m_h_roundWinner["5"]->Fill(m_teams[a].m_rank + 0.5, goalsA);
    m_h_roundWinner["5"]->Fill(m_teams[b].m_rank + 0.5, goalsB);
  }
}

void WCMC::recordStats(const std::string a, const std::string b, const int goalsA, const int goalsB) {
  const std::string key = a + "_" + b;
  std::map<std::string, TH2F*>::iterator it = m_h_matchResult.find(key);
  if (it == m_h_matchResult.end()) {
    m_h_matchResult[key] = new TH2F(key.c_str(), key.c_str(), 8, -.5, 7.5, 8, -.5, 7.5);
  }
  m_h_matchResult[key]->Fill(goalsA, goalsB);
}

void WCMC::doGroup(const std::string& group, const float low, const float high) {
  const vector<std::string>& teams = m_groups.at(group);
  for (unsigned i = 0; i < teams.size() - 1; ++i) {
    for (unsigned j = i + 1; j < teams.size(); ++j) {
      doMatch(teams.at(i), teams.at(j), low, high);
    }
  }
}

void WCMC::addTeam(const std::string& t, const std::string& abreviation) {
  const int pos = m_teams.size();
  m_teams[t] = Team();
  m_teams[t].m_rank = pos;
  m_teams[t].m_abreviation = abreviation;
  m_teamsByRank.push_back(t);
}

void WCMC::addHistoric(int goalsA, int goalsB) {
  m_h_GoalDiffData->Fill(abs( goalsA - goalsB ));
  m_h_GoalsData->Fill(goalsA + goalsB);
}

void WCMC::loadHistoricData() {  // 2014 WC
  std::string line;
  std::ifstream historic("wc_2014_results.txt");
  if (historic.is_open())  {
    while ( getline(historic, line) ) {
      std::istringstream buf(line);
      std::istream_iterator<std::string> beg(buf), end;
      std::vector<std::string> results(beg, end);
      addHistoric(std::stoi(results[0]), std::stoi(results[1]));
    }
    historic.close();
  }
  m_h_GoalDiffData->Scale( 1. / m_h_GoalDiffData->Integral() );
  m_h_GoalsData->Scale(1. / m_h_GoalsData->Integral() );
}

void WCMC::addTeams() {
  std::string line;
  std::ifstream teams("wc_2018_team_ranks.txt");
  if (teams.is_open())  {
    while ( getline(teams, line) ) {
      std::istringstream buf(line);
      std::istream_iterator<std::string> beg(buf), end;
      std::vector<std::string> results(beg, end);
      size_t location = results[0].find("-");
      if (location != std::string::npos) results[0].replace(location, 1, " ");
      addTeam(results[0], results[1]);
    }
    teams.close();
  }
  for (int i = 0; i < 6; ++i) m_h_roundWinner[std::to_string(i)] = new TH1F("", "", m_teams.size(), 0, m_teams.size()); // 5 is a special entry
}

void WCMC::addGroup(const std::string& group, const std::string& A, const std::string& B, const std::string& C, const std::string& D) {
  m_groups[group] = {A, B, C, D};
  for (unsigned i = 0; i < m_groups[group].size(); ++i) m_h_roundWinner[group + std::to_string(i)] = new TH1F("","", 4, -.5, 3.5);
}

void WCMC::resetTeamStatistics(const bool all) {
  for (auto& entry : m_teams) {
    entry.second.m_points = 0;
    if (all) {
      entry.second.m_goalDiff = 0;
      entry.second.m_goals = 0;
    }
  }
}

void WCMC::resetLaterGroups() {
  for (int i = 49; i < 65; ++i) {
    m_groups[std::to_string(i)].clear();
  }
}

void WCMC::addGroups() {
  addGroup("A", "Russia", "Saudi Ar.", "Egypt", "Uruguay");
  addGroup("B","Portugal", "Spain", "Morocco", "IR Iran");
  addGroup("C", "France", "Australia", "Peru", "Denmark");
  addGroup("D", "Argentina", "Iceland", "Croatia", "Nigeria");
  addGroup("E", "Brazil", "Switzerland", "Costa Rica", "Serbia");
  addGroup("F", "Germany", "Mexico", "Sweden", "Korea Rp.");
  addGroup("G", "Belgium", "Panama", "Tunisia", "England");
  addGroup("H", "Japan", "Senegal", "Colombia", "Poland");
}

WCMC::WCMC() {
  m_trialsMax = 100000;
  m_bestChiG = m_bestChiGD = -1;
    
  m_h_GoalsMC = new TH1F("MC",";Goals;Fraction",9,-0.5,8.5);
  m_h_GoalsData = new TH1F("Data",";Goals;",9,-0.5,8.5);
  m_h_GoalDiffMC = new TH1F("MC ",";Goals Difference;Fraction",9,-0.5,8.5);
  m_h_GoalDiffData = new TH1F("MC ",";Goals Difference;Fraction",9,-0.5,8.5);

  loadHistoricData();
  addTeams();
  addGroups();
}

void WCMC::runTraining(float& resultLow, float& resultHigh, const float startLow, const float stopLow, const float startHigh, const float stopHigh, const float step) {
  m_matchPrint = false;
  m_matchStats = false;
  m_goalsScored = false;
  float bestChi = 999;
  for (float trial_goalines_low = startLow; trial_goalines_low < stopLow; trial_goalines_low += step) {
    for (float trial_goalines_high = startHigh; trial_goalines_high > stopHigh; trial_goalines_high -= step) {
      if (trial_goalines_high <= trial_goalines_low) continue;
      int seed = 0;

      m_h_GoalsMC->Reset();
      m_h_GoalDiffMC->Reset();
      resetTeamStatistics(true);

      for (int trial = 0; trial < m_trialsMax; ++trial) {
        R.SetSeed(seed++);
        for (const std::string& group : groups)  doGroup(group, trial_goalines_low, trial_goalines_high);
      }

      m_h_GoalsMC->Scale( 1./m_h_GoalsMC->Integral() );
      m_h_GoalDiffMC->Scale( 1./m_h_GoalDiffMC->Integral() );

      float goodnessA = m_h_GoalsData->Chi2Test(m_h_GoalsMC, "NORM  UU CHI2/NDF");
      float goodnessB = m_h_GoalDiffData->Chi2Test(m_h_GoalDiffMC, "NORM  UU CHI2/NDF");

      if ( goodnessA + goodnessB < bestChi) {
        bestChi = goodnessA + goodnessB;
        m_bestChiG = goodnessA;
        m_bestChiGD = goodnessB;
        resultLow = trial_goalines_low;
        resultHigh = trial_goalines_high;
        std::cout << "--- Chi2 of:" << goodnessA + goodnessB << " for Low:" << resultLow << " High:" << resultHigh << endl;
      }
    }
  }
}

const std::string WCMC::getWinningTeam(const std::string group) {
  int winningPoints = -1, winningGD = -1, winningGoals = -1, winningRank = 999;
  std::string winningTeam;
  for (auto& team : m_groups[group]) {
    bool better = false;
    if ( m_teams[team].m_points > winningPoints ) better = true;
    else if ( m_teams[team].m_points == winningPoints) {
      if (m_teams[team].m_goalDiff > winningGD) better = true;
      else if (m_teams[team].m_goals > winningGoals) better = true;
      else if (m_teams[team].m_rank < winningRank) better = true; // This is not according to FIFA rules
    }
    if (better) { // TODO use a dummy Team for this
      winningPoints = m_teams[team].m_points;
      winningGD = m_teams[team].m_goalDiff;
      winningRank = m_teams[team].m_rank;
      winningGoals = m_teams[team].m_goals;
      winningTeam = team;
    }
  }
  return winningTeam;
}

void WCMC::runFinal(const float goalinessLow, const float goalinessHigh) {
  int seed = 0;
  m_h_GoalsMC->Reset();
  m_h_GoalDiffMC->Reset();
  m_goalsScored = true;
  for (int trial = 0; trial < m_trialsMax; ++ trial) {
    m_matchPrint = (trial == m_trialsMax-1);
    R.SetSeed(seed++);

    resetTeamStatistics(true);
    resetLaterGroups();
    m_matchStats = true;
    for (const std::string& group : groups)  {
      doGroup(group, goalinessLow, goalinessHigh);
      std::string teamPlace[4];
      for (int position = 0; position < 4; ++position) {
        teamPlace[position] = getWinningTeam(group);  
        m_teams[teamPlace[position]].m_points = -1; // Take out of action to get the next one
        m_h_roundWinner[group+std::to_string(position)]->Fill( std::distance(m_groups[group].begin(), std::find(m_groups[group].begin(), m_groups[group].end(), teamPlace[position])) );
      }
      if (m_matchPrint) std::cout << "Winner of group " << group << ":" << teamPlace[0] << ", runner up " << teamPlace[1] << std::endl;
      m_h_roundWinner["0"]->Fill( m_teams[teamPlace[0]].m_rank + 0.5 );
      m_h_roundWinner["0"]->Fill( m_teams[teamPlace[1]].m_rank + 0.5 );
      if (group == "A") {
        m_groups["49"].push_back(teamPlace[0]);
        m_groups["51"].push_back(teamPlace[1]);
      } else if (group == "B") {
        m_groups["51"].push_back(teamPlace[0]);
        m_groups["49"].push_back(teamPlace[1]);
      } else if (group == "C") {
        m_groups["50"].push_back(teamPlace[0]);
        m_groups["52"].push_back(teamPlace[1]);
      } else if (group == "D") {
        m_groups["52"].push_back(teamPlace[0]);
        m_groups["50"].push_back(teamPlace[1]);
      } else if (group == "E") {
        m_groups["53"].push_back(teamPlace[0]);
        m_groups["55"].push_back(teamPlace[1]);
      } else if (group == "F") {
        m_groups["55"].push_back(teamPlace[0]);
        m_groups["53"].push_back(teamPlace[1]);
      } else if (group == "G") {
        m_groups["54"].push_back(teamPlace[0]);
        m_groups["56"].push_back(teamPlace[1]);
      } else if (group == "H") {
        m_groups["56"].push_back(teamPlace[0]);
        m_groups["54"].push_back(teamPlace[1]);
      }
    }
    m_matchStats = false;

    resetTeamStatistics(false);
    for (int m = 49; m < 57; ++m) {
      doGroup(std::to_string(m), goalinessLow, goalinessHigh);
      const std::string winning = getWinningTeam(std::to_string(m));
      if (m_matchPrint) std::cout << "Winner of round " << m << ":" << winning << std::endl;
      m_h_roundWinner["1"]->Fill( m_teams[winning].m_rank + 0.5 ); // Many entries here, so we offset the axis ticks
      switch (m) {
        case 49: case 50: m_groups["57"].push_back(winning); break;
        case 51: case 52: m_groups["59"].push_back(winning); break;
        case 53: case 54: m_groups["58"].push_back(winning); break;
        case 55: case 56: m_groups["60"].push_back(winning); break;
      }
    }

    resetTeamStatistics(false);
    for (int m = 57; m < 61; ++m) {
      doGroup(std::to_string(m), goalinessLow, goalinessHigh);
      const std::string winning = getWinningTeam(std::to_string(m));
      if (m_matchPrint) std::cout << "Winner of QF match " << m << ":" << winning << std::endl;
      m_h_roundWinner["2"]->Fill( m_teams[winning].m_rank + 0.5 );
      switch (m) {
        case 57: case 58: m_groups["61"].push_back(winning); break;
        case 59: case 60: m_groups["62"].push_back(winning); break;
      }
    }

    resetTeamStatistics(false);
    doGroup("61", goalinessLow, goalinessHigh);
    const std::string finalistA = getWinningTeam("61");
    doGroup("62", goalinessLow, goalinessHigh);
    const std::string finalistB = getWinningTeam("62");
    m_h_roundWinner["3"]->Fill( m_teams[finalistA].m_rank + 0.5 ); 
    m_h_roundWinner["3"]->Fill( m_teams[finalistB].m_rank + 0.5 );

    resetTeamStatistics(false);
    m_groups["64"].push_back(finalistA);
    m_groups["64"].push_back(finalistB);
    doGroup("64", goalinessLow, goalinessHigh);
    const std::string winnerWinner = getWinningTeam("64");
    m_h_roundWinner["4"]->Fill( m_teams[winnerWinner].m_rank + 0.5 );
    if (m_matchPrint || trial % 1000 == 0) std::cout << "Trial:" << trial <<  " Winners of SFs " <<  finalistA << " & " << finalistB << ", WINNER WINNER:" << winnerWinner << std::endl << " ----------------- " << std::endl;
  }
  m_h_GoalsMC->Scale( 1./m_h_GoalsMC->Integral() );
  m_h_GoalDiffMC->Scale( 1./m_h_GoalDiffMC->Integral() );
}

void WCMC::execute() {
  float resultLowFine, resultHighFine;
  const bool reTrain = false;
  if (reTrain) {
    float resultLowCorse, resultHighCorse;
    runTraining(resultLowCorse, resultHighCorse, 0.1, 5.0, 5.0, 0.1, 0.1);
    runTraining(resultLowFine, resultHighFine, resultLowCorse - 0.5, resultLowCorse + 0.5, resultHighCorse + 0.5, resultHighCorse - 0.5, 0.01);
    std::cout << " ---->>>>> Tuned Low: "<< resultLowFine << " High: " << resultHighFine << "(Best chi2 G:" << m_bestChiG << ", GD:" << m_bestChiGD << ")" << std::endl;
  } else {
    resultLowFine = 1.52;
    resultHighFine = 1.79;
    m_bestChiG = 1.03463;
    m_bestChiGD = 0.616308;
  }
  runFinal(resultLowFine, resultHighFine);

  nicePlot* np_base = new nicePlot();
  np_base->setRBounds(1./m_trialsMax * 0.9, 2e-1);
  np_base->setLogz(true);
  np_base->normaliseToOne();
  const unsigned groupSize = m_groups["A"].size();
  for (unsigned i = 0; i < groupSize - 1; ++i) {
    for (unsigned j = i + 1; j < groupSize; ++j) {
      for (const std::string& group : groups) {
        const vector<std::string>& teams = m_groups.at(group);
        nicePlot* np = new nicePlot(np_base);
        np->init(teams.at(i) + " Goals", teams.at(j) + " Goals", "");
        TH2* h = m_h_matchResult[teams.at(i) + "_" + teams.at(j)];
        np->add2D(h);
        np->addLable(.5, .75, teams.at(i) + ": " + std::to_string( (int) (h->ProjectionX()->GetMaximumBin() - 1.) ) );
        np->addLable(.5, .80, teams.at(j) + ": " + std::to_string( (int) (h->ProjectionY()->GetMaximumBin() - 1.) ) );
        np->addLable(.5, .85, "Group: " + group);
      }
    }
  }
  bookOutput::setBreak(groups.size());
  bookOutput::get().doMultipadOutput("WCMC_GroupStage", 3, 2);
  bookOutput::clear();
  nicePlot* np_base_1d = new nicePlot();
  np_base_1d->setLineWidth(3);
  np_base_1d->setLegend(.7, .7);
  np_base_1d->setBounds(-.5, 8.5, 0.001, .85, 0., 2.);
  nicePlot* np_tuneGoals = new nicePlot(np_base_1d);
  np_tuneGoals->init("Total Goals", "Probability", "MC/Data");
  np_tuneGoals->addData(m_h_GoalsData, "Data 2014");
  np_tuneGoals->addMC(m_h_GoalsMC, "MC", true);
  np_tuneGoals->addLable(0.2,0.75, "Goaliness Low: " + std::to_string(resultLowFine));
  np_tuneGoals->addLable(0.2,0.80, "Goaliness High: " + std::to_string(resultHighFine));
  np_tuneGoals->addLable(0.2,0.85, "#chi^{2}/DoF: " + std::to_string(m_bestChiG));
  nicePlot* np_tuneGoalDiff = new nicePlot(np_base_1d);
  np_tuneGoalDiff->init("Goal Difference", "Probability", "MC/Data");
  np_tuneGoalDiff->addData(m_h_GoalDiffData, "Data 2014");
  np_tuneGoalDiff->addMC(m_h_GoalDiffMC, "MC", true);
  np_tuneGoalDiff->addLable(0.2,0.75, "Goaliness Low: " + std::to_string(resultLowFine));
  np_tuneGoalDiff->addLable(0.2,0.80, "Goaliness High: " + std::to_string(resultHighFine));
  np_tuneGoalDiff->addLable(0.2,0.85, "#chi^{2}/DoF: " + std::to_string(m_bestChiGD));
  bookOutput::setBreak(1);
  bookOutput::get().doMultipadOutput("WCMC_Tuning", 2, 1);
  bookOutput::clear();

  np_base_1d->useAltColourScheme(2);
  np_base_1d->setLineWidth(2);
  np_base_1d->normaliseToOne();
  np_base_1d->setLegend(.7, .75);
  np_base_1d->setBounds(-.5, 3.5, 0, 1.4);
  for (const std::string& group : groups) {
    nicePlot* np = new nicePlot(np_base_1d);
    np->init("Team", "Probability");
    np->addStackMC(m_h_roundWinner[group+"3"], "4th");
    np->addStackMC(m_h_roundWinner[group+"2"], "3rd");
    np->addStackMC(m_h_roundWinner[group+"1"], "Runner Up");
    np->addStackMC(m_h_roundWinner[group+"0"], "Winner");
    const vector<std::string>& teams = m_groups.at(group);
    for (const std::string& team : teams) np->addBinLabel(team);
    np->addLable(.25, .85, "Group " + group);
    np->addLable(.25, .80, "Winner: " + teams.at( m_h_roundWinner[group+"0"]->GetMaximumBin()-1 ) );
    np->addLable(.25, .75, "Runner Up: " + teams.at( m_h_roundWinner[group+"1"]->GetMaximumBin()-1 ) );
  }
  bookOutput::get().doMultipadOutput("WCMC_GroupResults", 2, 4);
  bookOutput::clear();

  np_base_1d->setBounds(0, m_teams.size());
  np_base_1d->setDoLegend(false);
  np_base_1d->stretch = true;
  np_base_1d->useAltColourScheme(1);
  np_base_1d->setLineWidth(6);
  for (const std::string& team : m_teamsByRank) np_base_1d->addBinLabel(m_teams[team].m_abreviation);
  int numberOfPassingTeams = 16;
  for (int i = 0; i < 5; ++i) {
    nicePlot* np_round = new nicePlot(np_base_1d);
    np_round->n_mc += i;
    np_round->init("Team", "Probability");
    np_round->setYBounds(0, 1. - (i * 0.20));
    np_round->addMC(m_h_roundWinner[std::to_string(i)], "");
    np_round->scaleLastMC(numberOfPassingTeams);
    numberOfPassingTeams /= 2;
    switch(i) {
      case 0: np_round->addLable(.7, .8, "Probability of passing the Group Stage"); break; 
      case 1: np_round->addLable(.7, .8, "Probability of passing the Round Of 16"); break; 
      case 2: np_round->addLable(.7, .8, "Probability of passing the Quarter Finals"); break; 
      case 3: np_round->addLable(.7, .8, "Probability of passing the Semi Finals"); break; 
      case 4: np_round->addLable(.7, .8, "Probability of Winning The Tournament"); break; 
    }
  }
  bookOutput::setBreak(5);
  bookOutput::get().doBookOutput("WCMC_KnockoutResults");
  bookOutput::clear();

  nicePlot* np_goals = new nicePlot(np_base_1d);
  np_goals->cancelNormaliseToOne();
  np_goals->init("Team", "Number of Goals");
  np_goals->addMC(m_h_roundWinner["5"], "");
  np_goals->scaleLastMC(1./m_trialsMax);
  np_goals->setYBounds(0, 12.);
  np_goals->addLable(.6, .8, "Mean Number of Goals Scored Per Team");
  bookOutput::setBreak(1);
  bookOutput::get().doBookOutput("WCMC_Goals");
}

int main() {
  gROOT->ProcessLine(".L AtlasStyle.C");
  gROOT->ProcessLine("SetAtlasStyle();");
  gErrorIgnoreLevel = 10000;
  WCMC wc2018;
  wc2018.execute();
}

int wcMC_2018() {
  return main();
}