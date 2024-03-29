#include <map>
#include <string>
#include <iostream>
#include <fstream>
#include <sstream>
#include <vector>
#include <iomanip>
#include <TRandom3.h>
#include <TROOT.h>
#include <TH2.h>
#include "nicePlot.cxx"

enum Mode {kFULL_TOURNAMENT, kAFTER_GROUP, kAFTER_16, kAFTER_QUARTER, kAFTER_SEMI};

class WCMC {  
  public:
    WCMC(const Mode mode);
    void doMatch(const std::string& a, const std::string& b, const float low, const float high);
    void doGroup(const std::string& group, const float low, const float high);
    void addHistoric(int goalsA, int goalsB, int year);
    void loadHistoricData();
    void addTeams();
    void addTeam(const std::string& t, const std::string& abreviation, const int rank);
    std::vector<std::string> readLine(const std::string& line);
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
      Team() { m_points = 0; m_goalDiff = 0; m_goals = 0; m_rank = 0; m_index = 0; }
      int m_rank;
      int m_index;
      int m_points;
      int m_goalDiff;
      int m_goals;
      std::string m_abreviation;
    };

    std::map<std::string, Team > m_teams;
    std::map<std::string, std::vector<std::string>> m_groups;   
    TRandom3 R;
    TH1F* m_h_GoalsMC;
    TH1F* m_h_GoalsData_Test;
    TH1F* m_h_GoalsData_Training;
    TH1F* m_h_GoalDiffMC;
    TH1F* m_h_GoalDiffData_Test;
    TH1F* m_h_GoalDiffData_Training;
    TH2F* m_h_trainCorse;
    TH2F* m_h_trainFine;
    std::map<std::string, TH2F*> m_h_matchResult;
    std::map<std::string, TH1F*> m_h_roundWinner;
    std::map<std::string, std::string> m_teamToAbrieviation;
    int m_trialsMax;
    int m_totalTeams;
    bool m_matchPrint, m_matchStats, m_goalsScored;
    std::vector<std::string> group_letters;
    std::vector<std::string> m_teamsByRank;
    std::vector<std::string> m_laterRoundTeams;
    float m_bestChiG_Test, m_bestChiGD_Test, m_bestChiG_Training, m_bestChiGD_Training;
    Mode m_mode; // Tournament progression
};

void WCMC::doMatch(const std::string& a, const std::string& b, const float low, const float high) {
  const float reduction = m_totalTeams / high;

  // std::cout << "     Match " << a << " vs " << b << std::endl;

  float scoreA = low + ((m_totalTeams - m_teams[a].m_rank) / reduction);
  float scoreB = low + ((m_totalTeams - m_teams[b].m_rank) / reduction);
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
    m_h_roundWinner["5"]->Fill(m_teams[a].m_index + 0.5, goalsA);
    m_h_roundWinner["5"]->Fill(m_teams[b].m_index + 0.5, goalsB);
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
  const std::vector<std::string>& teams = m_groups.at(group);
  // std::cout << " Do group " << group << ", size " << teams.size() << std::endl;
  for (unsigned i = 0; i < teams.size() - 1; ++i) {
    for (unsigned j = i + 1; j < teams.size(); ++j) {
      doMatch(teams.at(i), teams.at(j), low, high);
    }
  }
}

void WCMC::addTeam(const std::string& t, const std::string& abreviation, const int rank) {
  const int pos = m_teams.size();
  m_teams[t] = Team();
  m_teams[t].m_rank = rank;
  m_teams[t].m_index = pos;
  m_teams[t].m_abreviation = abreviation;
  std::cout << "Team " << t << " (" << abreviation << ") Rank:" << rank << " Index:" << pos << std::endl; 
  m_teamsByRank.push_back(t);
  m_teamToAbrieviation[t] = abreviation;
}

void WCMC::addHistoric(int goalsA, int goalsB, int year) {
  switch (year) {
    case 2014: m_h_GoalDiffData_Training->Fill(abs( goalsA - goalsB )); m_h_GoalsData_Training->Fill(goalsA + goalsB); break;
    case 2018: m_h_GoalDiffData_Training->Fill(abs( goalsA - goalsB )); m_h_GoalsData_Training->Fill(goalsA + goalsB); break;
    case 2022: m_h_GoalDiffData_Test->Fill(abs( goalsA - goalsB )); m_h_GoalsData_Test->Fill(goalsA + goalsB); break;
  }
}

void WCMC::loadHistoricData() {  // 2014 WC
  std::string line;
  std::ifstream historic2014("wc_2014_results.txt");
  while ( getline(historic2014, line) ) {
    std::vector<std::string> results = readLine(line);
    addHistoric(std::stoi(results[0]), std::stoi(results[1]), 2014);
  }
  std::ifstream historic2018("wc_2018_results.txt");
  while ( getline(historic2018, line) ) {
    std::vector<std::string> results = readLine(line);
    addHistoric(std::stoi(results[0]), std::stoi(results[1]), 2018);
  }
  m_h_GoalDiffData_Test->Scale( 1. / m_h_GoalDiffData_Test->Integral() );
  m_h_GoalsData_Test->Scale(1. / m_h_GoalsData_Test->Integral() );
  m_h_GoalDiffData_Training->Scale( 1. / m_h_GoalDiffData_Training->Integral() );
  m_h_GoalsData_Training->Scale(1. / m_h_GoalsData_Training->Integral() );
}

std::vector<std::string> WCMC::readLine(const std::string& line) {
  std::istringstream buf(line);
  std::istream_iterator<std::string> beg(buf), end;
  std::vector<std::string> results(beg, end);
  for (size_t i = 0; i < results.size(); ++i) {
    size_t location = results[i].find("-");
    if (location != std::string::npos) results[i].replace(location, 1, " ");
  }
  return results;
}

void WCMC::addTeams() {

  std::string line;
  if (m_mode > kFULL_TOURNAMENT) {
    std::ifstream pass;
    if      (m_mode == kAFTER_GROUP)   pass.open("wc_2022_pass_groups.txt");
    else if (m_mode == kAFTER_16)      pass.open("wc_2022_pass_16.txt");
    else if (m_mode == kAFTER_QUARTER) pass.open("wc_2022_pass_quarter.txt");
    else if (m_mode == kAFTER_SEMI)    pass.open("wc_2022_pass_semi.txt");

    while ( getline(pass, line) ) {
      std::vector<std::string> results = readLine(line);
      m_laterRoundTeams.push_back( results[0] );
      std::cout << "Passed stage " << (int)m_mode << ": '" << results[0] << "'" << std::endl;
    }
  }

  std::ifstream teams("wc_2022_team_ranks.txt");
  m_totalTeams = 0;
  while ( getline(teams, line) ) {
    std::vector<std::string> results = readLine(line);
    if (m_mode == kFULL_TOURNAMENT || std::count(m_laterRoundTeams.begin(), m_laterRoundTeams.end(), results[0]) != 0)  {
      addTeam(results[0], results[1], /*rank ==*/ m_totalTeams);
    } else std::cout << "  Dropping team: '" << results[0] << "'" << std::endl;
    ++m_totalTeams;
  }

  /*
  std::ifstream teams_updated("wc_2018_team_ranks_updated.txt");
  int rank = 0;
  while ( getline(teams_updated, line) ) {
    std::vector<std::string> results = readLine(line);
    if (m_mode == kFULL_TOURNAMENT || std::count(m_laterRoundTeams.begin(), m_laterRoundTeams.end(), results[0]) != 0)  {
      m_teams[results[0]].m_rank = rank;
    }
    ++rank;
  }
  */

  for (int i = 0; i < 6; ++i) m_h_roundWinner[std::to_string(i)] = new TH1F("", "", m_teams.size(), 0, m_teams.size()); // 5 is a special entry
}

void WCMC::addGroup(const std::string& group, const std::string& A, const std::string& B, const std::string& C, const std::string& D) {
  m_groups[group] = {A, B, C, D};
  group_letters.push_back(group);
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
  std::ifstream groups("wc_2022_groups.txt");
  std::string line;
  while ( getline(groups, line) ) {
    std::vector<std::string> r = readLine(line);
    addGroup(r[0], r[1], r[2], r[3], r[4]);
  }
}

WCMC::WCMC(const Mode mode) {
  m_trialsMax = 1000000;
  m_bestChiG_Test = m_bestChiGD_Test = m_bestChiG_Training = m_bestChiGD_Training = -1;

  m_mode = mode; 
    
  m_h_GoalsMC = new TH1F("MC",";Goals;Fraction",9,-0.5,8.5);
  m_h_GoalsData_Test = new TH1F("Data",";Goals;",9,-0.5,8.5);
  m_h_GoalsData_Training = new TH1F("Data",";Goals;",9,-0.5,8.5);
  m_h_GoalDiffMC = new TH1F("MC ",";Goals Difference;Fraction",9,-0.5,8.5);
  m_h_GoalDiffData_Test = new TH1F("MC ",";Goals Difference;Fraction",9,-0.5,8.5);
  m_h_GoalDiffData_Training = new TH1F("MC ",";Goals Difference;Fraction",9,-0.5,8.5);

  std::cout << "Loading Historic" << std::endl;
  loadHistoricData();
  std::cout << "Loading Teams" << std::endl;
  addTeams();
  if (m_mode == kFULL_TOURNAMENT) {
    std::cout << "Loading Groups" << std::endl;
    addGroups();
  }

  execute();
}

void WCMC::runTraining(float& resultLow, float& resultHigh, const float startLow, const float stopLow, const float startHigh, const float stopHigh, const float step) {
  m_matchPrint = false;
  m_matchStats = false;
  m_goalsScored = false;

  TH2F* hTrain;
  int nBins = (startHigh - stopHigh) / step; 
  if (step > 0.05) {
    m_h_trainCorse = new TH2F("TrainC", ";Low;High", nBins+1, startLow, stopLow, nBins+1, stopHigh, startHigh);
    std::cout << "New training CORSE " << nBins << ", " << startLow << " " << stopLow << ", " << startHigh << " " << stopHigh << std::endl; 
    hTrain = m_h_trainCorse;
  } else {
    m_h_trainFine = new TH2F("TrainF", ";Low;High", nBins, startLow, stopLow, nBins, stopHigh, startHigh);
    std::cout << "New training FINE " << nBins << ", " << startLow << " " << stopLow << ", " << startHigh << " " << stopHigh << std::endl; 
    hTrain = m_h_trainFine;
  }

  float bestChi = 999;
  for (float trial_goalines_low = startLow; trial_goalines_low < stopLow; trial_goalines_low += step) {
    for (float trial_goalines_high = startHigh; trial_goalines_high > stopHigh; trial_goalines_high -= step) {
      if (trial_goalines_high - trial_goalines_low < 1e-2) {
        //int b =  hTrain->FindBin(trial_goalines_low, trial_goalines_high);
        //hTrain->SetBinContent(b, 0);
        continue;
      }
      std::cout << std::setprecision(4) << "[" << trial_goalines_low << "," << trial_goalines_high << "] " << std::flush;
      int seed = 0;

      m_h_GoalsMC->Reset();
      m_h_GoalDiffMC->Reset();
      resetTeamStatistics(true);

      int multiplier = 10;
      int trials = 10000 * multiplier;
      if (hTrain == m_h_trainFine) trials = 1000 * multiplier;

      for (int trial = 0; trial < trials; ++trial) {
        R.SetSeed(seed++);
        for (const std::string& group : group_letters)  doGroup(group, trial_goalines_low, trial_goalines_high);
      }

      m_h_GoalsMC->Scale( 1./m_h_GoalsMC->Integral() );
      m_h_GoalDiffMC->Scale( 1./m_h_GoalDiffMC->Integral() );

      const float goodnessA = m_h_GoalsData_Training->Chi2Test(m_h_GoalsMC, "NORM  UU CHI2/NDF");
      const float goodnessB = m_h_GoalDiffData_Training->Chi2Test(m_h_GoalDiffMC, "NORM  UU CHI2/NDF");

      int b =  hTrain->FindBin(trial_goalines_low, trial_goalines_high);
      hTrain->SetBinContent(b, goodnessA+goodnessB);
      //std::cout << std::endl << "Fill bin " << trial_goalines_low << "," << trial_goalines_high << " = " << b << " to " << goodnessA+goodnessB << std::endl;

      if ( goodnessA + goodnessB < bestChi) {
        bestChi = goodnessA + goodnessB;
        //m_bestChiG_Test  = m_h_GoalsData_Test->Chi2Test(m_h_GoalsMC, "NORM  UU CHI2/NDF");
        //m_bestChiGD_Test = m_h_GoalDiffData_Test->Chi2Test(m_h_GoalDiffMC, "NORM  UU CHI2/NDF");
        m_bestChiG_Training = goodnessA;
        m_bestChiGD_Training = goodnessB;
        resultLow = trial_goalines_low;
        resultHigh = trial_goalines_high;
        std::cout << "--- Chi2 of:" << goodnessA + goodnessB << " for Low:" << resultLow << " High:" << resultHigh << std::endl;
      }
    }
  }
  std::cout << "chi2 when using the Training tuning dataset: G=" << m_bestChiG_Training << " GD=" << m_bestChiGD_Training << std::endl;
  std::cout << "chi2 against the Test dataset: G=" << m_bestChiG_Test << " GD=" << m_bestChiGD_Test << std::endl;
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
  bool firstEnglandWin = true;
  std::map<std::string, int> outcomes;
  std::map<std::string, int> outcomesToQuarter;
  std::map<std::string, int> outcomesToSemi;

  for (int trial = 0; trial < m_trialsMax; ++ trial) {
    m_matchPrint = (trial == m_trialsMax-1);
    R.SetSeed(seed++);

    std::set<std::string> dropOutAt16;
    std::set<std::string> dropOutAtQuarter;
    std::set<std::string> dropOutAtSemi;

    resetTeamStatistics(true);
    resetLaterGroups();

    if (m_mode == kFULL_TOURNAMENT) {
      m_matchStats = true;
      for (const std::string& group : group_letters)  {
        doGroup(group, goalinessLow, goalinessHigh);
        std::string teamPlace[4];
        for (int position = 0; position < 4; ++position) {
          teamPlace[position] = getWinningTeam(group);  
          m_teams[teamPlace[position]].m_points = -1; // Take out of action to get the next one
          m_h_roundWinner[group+std::to_string(position)]->Fill( std::distance(m_groups[group].begin(), std::find(m_groups[group].begin(), m_groups[group].end(), teamPlace[position])) );
        }
        if (m_matchPrint) std::cout << "Winner of group " << group << ":" << teamPlace[0] << ", runner up " << teamPlace[1] << std::endl;
        m_h_roundWinner["0"]->Fill( m_teams[teamPlace[0]].m_index + 0.5 );
        m_h_roundWinner["0"]->Fill( m_teams[teamPlace[1]].m_index + 0.5 );
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
    }
    m_matchStats = false;

    if (m_mode == kAFTER_GROUP) {
      m_matchStats = true;
      for ( size_t t = 0; t < m_laterRoundTeams.size(); ++t ) {
        switch (t) {
          case 0: case 3:   m_groups["49"].push_back(m_laterRoundTeams.at(t)); break;
          case 1: case 2:   m_groups["51"].push_back(m_laterRoundTeams.at(t)); break; 
          case 4: case 7:   m_groups["50"].push_back(m_laterRoundTeams.at(t)); break; 
          case 5: case 6:   m_groups["52"].push_back(m_laterRoundTeams.at(t)); break; 
          case 8: case 11:  m_groups["53"].push_back(m_laterRoundTeams.at(t)); break;
          case 9: case 10:  m_groups["55"].push_back(m_laterRoundTeams.at(t)); break;
          case 12: case 15: m_groups["54"].push_back(m_laterRoundTeams.at(t)); break;
          case 13: case 14: m_groups["56"].push_back(m_laterRoundTeams.at(t)); break;
        }
      }
    }

    if (m_mode < kAFTER_16) {
      resetTeamStatistics(false);
      for (int m = 49; m < 57; ++m) {
        doGroup(std::to_string(m), goalinessLow, goalinessHigh);
        const std::string winning = getWinningTeam(std::to_string(m));
        if (m_matchPrint) std::cout << "Winner of round " << m << ":" << winning << std::endl;
        m_h_roundWinner["1"]->Fill( m_teams[winning].m_index + 0.5 ); // Many entries here, so we offset the axis ticks
        m_teams[winning].m_points = -1; // Disable to get runner up
        const std::string dropOut = getWinningTeam(std::to_string(m));
        dropOutAt16.insert( m_teamToAbrieviation[dropOut] );
        switch (m) {
          case 49: case 50: m_groups["57"].push_back(winning); break;
          case 51: case 52: m_groups["59"].push_back(winning); break;
          case 53: case 54: m_groups["58"].push_back(winning); break;
          case 55: case 56: m_groups["60"].push_back(winning); break;
        }
      }
    }
    m_matchStats = false;

    if (m_mode == kAFTER_16) {
      m_matchStats = true;
      for ( size_t t = 0; t < m_laterRoundTeams.size(); ++t ) {
        switch (t) {
          case 0: case 1: m_groups["57"].push_back(m_laterRoundTeams.at(t)); break;
          case 2: case 3: m_groups["59"].push_back(m_laterRoundTeams.at(t)); break; 
          case 4: case 5: m_groups["58"].push_back(m_laterRoundTeams.at(t)); break; 
          case 6: case 7: m_groups["60"].push_back(m_laterRoundTeams.at(t)); break; 
        }
      }
    }

    if (m_mode < kAFTER_QUARTER) {
      resetTeamStatistics(false);
      for (int m = 57; m < 61; ++m) {
        doGroup(std::to_string(m), goalinessLow, goalinessHigh);
        const std::string winning = getWinningTeam(std::to_string(m));
        if (m_matchPrint) std::cout << "Winner of QF match " << m << ":" << winning << std::endl;
        m_h_roundWinner["2"]->Fill( m_teams[winning].m_index + 0.5 );
        m_teams[winning].m_points = -1; // Disable to get runner up
        const std::string dropOut = getWinningTeam(std::to_string(m));
        dropOutAtQuarter.insert( m_teamToAbrieviation[dropOut] );
        switch (m) {
          case 57: case 58: m_groups["61"].push_back(winning); break;
          case 59: case 60: m_groups["62"].push_back(winning); break;
        }
      }
    }
    m_matchStats = false;

    if (m_mode == kAFTER_QUARTER) {
      m_matchStats = true;
      for ( size_t t = 0; t < m_laterRoundTeams.size(); ++t ) {
        switch (t) {
          case 0: case 1: m_groups["61"].push_back(m_laterRoundTeams.at(t)); break;
          case 2: case 3: m_groups["62"].push_back(m_laterRoundTeams.at(t)); break; 
        }
      }
    }

    std::string finalistA, finalistB, runnerUpA, runnerUpB;
    if (m_mode < kAFTER_SEMI) {
      resetTeamStatistics(false);
      doGroup("61", goalinessLow, goalinessHigh);
      finalistA = getWinningTeam("61");
      m_teams[finalistA].m_points = -1; // Disable to get runner up
      runnerUpA = getWinningTeam("61");
      dropOutAtSemi.insert( m_teamToAbrieviation[runnerUpA] );
      m_groups["63"].push_back(runnerUpA); // Runner up
      doGroup("62", goalinessLow, goalinessHigh);
      finalistB = getWinningTeam("62");
      m_teams[finalistB].m_points = -1; // Disable to get runner up
      runnerUpB = getWinningTeam("62");
      dropOutAtSemi.insert( m_teamToAbrieviation[runnerUpB] );
      m_groups["63"].push_back(runnerUpB); // Runner up
      m_h_roundWinner["3"]->Fill( m_teams[finalistA].m_index + 0.5 ); 
      m_h_roundWinner["3"]->Fill( m_teams[finalistB].m_index + 0.5 );
    }
    m_matchStats = false;

    if (m_mode == kAFTER_SEMI) {
      m_matchStats = true;
      finalistA = m_laterRoundTeams.at(0);
      finalistB = m_laterRoundTeams.at(1);
      m_groups["63"].push_back(m_laterRoundTeams.at(2));
      m_groups["63"].push_back(m_laterRoundTeams.at(3));
    }

    resetTeamStatistics(false);
    m_groups["64"].push_back(finalistA);
    m_groups["64"].push_back(finalistB);
    doGroup("63", goalinessLow, goalinessHigh);
    const std::string thirdPlace = getWinningTeam("63");
    m_teams[thirdPlace].m_points = -1; // Disable to get 4th place
    const std::string fourthPlace = getWinningTeam("63");
    //
    doGroup("64", goalinessLow, goalinessHigh);
    const std::string winnerWinner = getWinningTeam("64");
    m_teams[winnerWinner].m_points = -1; // Disable to get 2th place
    const std::string secondPlace = getWinningTeam("64");
    m_h_roundWinner["4"]->Fill( m_teams[winnerWinner].m_index + 0.5 );
    if (m_matchPrint || trial % 10000 == 0) std::cout << "Trial:" << trial 
      << " 4th place:" << fourthPlace << " 3rd place:" << thirdPlace << ". Winners of SFs " <<  finalistA << " & " << finalistB 
      << ", WINNER WINNER:" << winnerWinner 
      << std::endl << " ----------------- " << std::endl;

    std::stringstream ss;
    ss << m_teamToAbrieviation[winnerWinner] << "/"
      << m_teamToAbrieviation[secondPlace] << "/";
    int i = 0;
    for (const std::string& s : dropOutAtSemi) {
      ss << s;
      if (++i < 2) ss << "_";
    }

    outcomesToSemi[ ss.str() ]++;

    ss << "/";
    i = 0;
    for (const std::string& s : dropOutAtQuarter) {
      ss << s;
      if (++i < 4) ss << "_";
    }

    outcomesToQuarter[ ss.str() ]++;

    ss << "/";
    i = 0;
    for (const std::string& s : dropOutAt16) {
      ss << s;
      if (++i < 8) ss << "_";
    }

    outcomes[ ss.str() ]++;

    if (firstEnglandWin && winnerWinner == "England") {
      std::cout << std::endl << std::endl << "1st England win on trial " << trial << " " << ss.str() << std::endl << std::endl;
      firstEnglandWin = false;
    }

  }

  m_h_GoalsMC->Scale( 1./m_h_GoalsMC->Integral() );
  m_h_GoalDiffMC->Scale( 1./m_h_GoalDiffMC->Integral() );

  if (m_mode == kAFTER_QUARTER) return;

  std::cout << "Outcomes to semi size is " << outcomesToSemi.size() << std::endl;

  // Find most current outcomes
  int iterations = 0;
  int print = 0;
  while (outcomesToSemi.size()) {
    // Find
    int highestScore = 0;
    for (auto const& [key, val] : outcomesToSemi) {
      if (val > highestScore) {
        highestScore = val;
      }
    }
    // Extract
    std::vector<std::string> outcomesWithScore;
    for (auto const& [key, val] : outcomesToSemi) {
      if (val == highestScore) {
        outcomesWithScore.push_back(key);
      }
    }
    // Erase & report
    bool doPrint = outcomesWithScore.size() <= 20 && ++print < 20;
    if (!doPrint && outcomesWithScore.size() > 1) std::cout << "Most common outcome (to semi) #" << ++iterations << ": with " << highestScore << " instances has " << outcomesWithScore.size() << " members" << std::endl;
    for (const std::string& s : outcomesWithScore) {
      outcomesToSemi.erase(s);
      if (doPrint) std::cout << "Most common outcome (to semi) #" << ++iterations << ": with " << highestScore << " instances = " << s << std::endl;
    }
  }

  if (m_mode == kAFTER_16) return;

  std::cout << "Outcomes to quarter size is " << outcomesToQuarter.size() << std::endl;

  // Find most current outcomes
  iterations = 0;
  print = 0;
  while (outcomesToQuarter.size()) {
    // Find
    int highestScore = 0;
    for (auto const& [key, val] : outcomesToQuarter) {
      if (val > highestScore) {
        highestScore = val;
      }
    }
    // Extract
    std::vector<std::string> outcomesWithScore;
    for (auto const& [key, val] : outcomesToQuarter) {
      if (val == highestScore) {
        outcomesWithScore.push_back(key);
      }
    }
    // Erase & report
    bool doPrint = outcomesWithScore.size() <= 20 && ++print < 20;
    if (!doPrint && outcomesWithScore.size() > 1) std::cout << "Most common outcome (to quarter) #" << ++iterations << ": with " << highestScore << " instances has " << outcomesWithScore.size() << " members" << std::endl;
    for (const std::string& s : outcomesWithScore) {
      outcomesToQuarter.erase(s);
      if (doPrint) std::cout << "Most common outcome (to quarter) #" << ++iterations << ": with " << highestScore << " instances = " << s << std::endl;
    }
  }

  if (m_mode == kAFTER_GROUP) return;

  std::cout << "Outcomes size is " << outcomes.size() << std::endl;

  // Find most current outcomes
  iterations = 0;
  print = 0;
  while (outcomes.size()) {
    // Find
    int highestScore = 0;
    for (auto const& [key, val] : outcomes) {
      if (val > highestScore) {
        highestScore = val;
      }
    }
    // Extract
    std::vector<std::string> outcomesWithScore;
    for (auto const& [key, val] : outcomes) {
      if (val == highestScore) {
        outcomesWithScore.push_back(key);
      }
    }
    // Erase & report
    bool doPrint = outcomesWithScore.size() <= 20 && ++print < 20;
    if (!doPrint && outcomesWithScore.size() > 1) std::cout << "Most common outcome #" << ++iterations << ": with " << highestScore << " instances has " << outcomesWithScore.size() << " members" << std::endl;
    for (const std::string& s : outcomesWithScore) {
      outcomes.erase(s);
      if (doPrint) std::cout << "Most common outcome #" << ++iterations << ": with " << highestScore << " instances = " << s << std::endl;
    }
  }




}

void WCMC::execute() {
  std::cout << "Execute with mode " << (int)m_mode << std::endl;
  float resultLowFine, resultHighFine;
  
  const bool reTrain = false;

  if (reTrain == true && m_mode != kFULL_TOURNAMENT) {
    std::cout << "Error. Can only train when m_mode = kFULL_TOURNAMENT";
    return;
  }

  if (reTrain) {
    float resultLowCorse, resultHighCorse;
    runTraining(resultLowCorse, resultHighCorse, 0.1, 5.0, 5.0, 0.1, /*step*/0.1);
    runTraining(resultLowFine, resultHighFine, resultLowCorse - 0.5, resultLowCorse + 0.5, resultHighCorse + 0.5, resultHighCorse - 0.5, /*step*/0.01);
    std::cout << " ---->>>>> Tuned Low: "<< resultLowFine << " High: " << resultHighFine << "(Best chi2 G:" << m_bestChiG_Training << ", GD:" << m_bestChiGD_Training << ")" << std::endl;
 
    nicePlot* npC = new nicePlot();
    npC->setLogz(true);
    npC->init("Goaliness Lower", "Goaliness Upper", "#chi^{2}/DoF");
    npC->setRBounds(1.5, m_h_trainCorse->GetMaximum() * 1.01);
    npC->add2D(m_h_trainCorse);

    nicePlot* npF = new nicePlot();
    npF->setLogz(false);
    npF->init("Goaliness Lower", "Goaliness Upper", "#chi^{2}/DoF");
    npF->setRBounds(2.45, 3.45);
    npF->add2D(m_h_trainFine);

    bookOutput::get().doBookOutput("WCMC_TuningGrid");
    bookOutput::clear();

  } else {
    // 2022
    resultLowFine = 1.53;
    resultHighFine = 1.54;
    m_bestChiG_Training = 1.594;
    m_bestChiGD_Training = 0.8898;

    // 2018
    //resultLowFine = 1.52;
    //resultHighFine = 1.79;
    //m_bestChiG_Training = 1.03463;
    //m_bestChiGD_Training = 0.616308;
  }
  runFinal(resultLowFine, resultHighFine);

  int numberOfPassingTeams = 16;
  for (int i=0; i < (int)m_mode; ++i) numberOfPassingTeams /= 2;
  unsigned start = 0, end = 0; // Used mid-tournament
  switch (m_mode) {
    case kFULL_TOURNAMENT: break;
    case kAFTER_GROUP:   start = 49; end = 57; break;
    case kAFTER_16:      start = 57; end = 61; break;
    case kAFTER_QUARTER: start = 61; end = 63; break; // Note game 63 is for 3rd place, we don't do this one
    case kAFTER_SEMI:    start = 63; end = 65; break;
  }

  bookOutput::clear();

  nicePlot* np_base = new nicePlot();
  np_base->setRBounds(1./m_trialsMax * 0.9, 2e-1);
  np_base->setLogz(true);
  np_base->normaliseToOne();
  if (m_mode == kFULL_TOURNAMENT) { // Group stage games
    const unsigned groupSize = m_groups["A"].size();
    for (unsigned i = 0; i < groupSize - 1; ++i) {
      for (unsigned j = i + 1; j < groupSize; ++j) {
        for (const std::string& group : group_letters) {
          const std::vector<std::string>& teams = m_groups.at(group);
          nicePlot* np = new nicePlot(np_base);
          np->init(teams.at(i) + " Goals", teams.at(j) + " Goals", "");
          TH2* h = m_h_matchResult[teams.at(i) + "_" + teams.at(j)];
          np->add2D(h);
          int maxX = -1, maxY = -1, maxZ = -1;
          h->GetBinXYZ(h->GetMaximumBin(), maxX, maxY, maxZ);
          np->addLable(.5, .75, teams.at(i) + ": " + std::to_string( maxX - 1 ));
          np->addLable(.5, .80, teams.at(j) + ": " + std::to_string( maxY - 1 ));
          np->addLable(.5, .85, "Group: " + group);
        }
      }
    }
    bookOutput::setBreak(group_letters.size());
    bookOutput::get().doMultipadOutput("WCMC_GroupStage", 3, 2);
    bookOutput::clear();
  } else if (m_mode >= kAFTER_GROUP) { // Knockout games
    for (unsigned i = start; i < end; ++i) {
      const std::string teamA = m_groups[std::to_string(i)].at(0);
      const std::string teamB = m_groups[std::to_string(i)].at(1);
      nicePlot* np = new nicePlot(np_base);
      np->init(teamA + " Goals", teamB + " Goals", "");
      TH2* h = m_h_matchResult[teamA + "_" + teamB];
      np->add2D(h);
      int maxX = -1, maxY = -1, maxZ = -1;
      h->GetBinXYZ(h->GetMaximumBin(), maxX, maxY, maxZ);
      np->addLable(.5, .75, teamA + ": " + std::to_string( maxX - 1 ));
      np->addLable(.5, .80, teamB + ": " + std::to_string( maxY - 1 ));
      np->addLable(.5, .85, "Game: " + std::to_string(i));
    }
    bookOutput::get().doBookOutput("WCMC_KnockoutGameScores_Mode" + std::to_string(m_mode));
    bookOutput::clear();    
  }

  nicePlot* np_base_1d = new nicePlot();
  np_base_1d->setLineWidth(3);
  np_base_1d->setLegend(.7, .7);
  
  np_base_1d->setBounds(-.5, 8.5, 0.001, .85, 0., 2.);
  nicePlot* np_tuneGoals = new nicePlot(np_base_1d);
  np_tuneGoals->init("Total Goals", "Probability", "MC/Data");
  np_tuneGoals->addData(m_h_GoalsData_Training, "Data 14+18", 0, false);
  np_tuneGoals->addMC(m_h_GoalsMC, "MC", true, 0.);
  //np_tuneGoals->addData(m_h_GoalsData_Test, "Data 2022", 0, false);
  np_tuneGoals->addLable(0.2,0.70, "Goaliness Low: " + std::to_string(resultLowFine));
  np_tuneGoals->addLable(0.2,0.75, "Goaliness High: " + std::to_string(resultHighFine));
  np_tuneGoals->addLable(0.2,0.80, "#chi^{2}/DoF 14+18: " + std::to_string(m_bestChiG_Training));
  //np_tuneGoals->addLable(0.2,0.85, "#chi^{2}/DoF 2022: " + std::to_string(m_bestChiG_Test));
  nicePlot* np_tuneGoalDiff = new nicePlot(np_base_1d);
  np_tuneGoalDiff->init("Goal Difference", "Probability", "MC/Data");
  np_tuneGoalDiff->addData(m_h_GoalDiffData_Training, "Data 14+18", 0, false);
  np_tuneGoalDiff->addMC(m_h_GoalDiffMC, "MC", true, 0);
  //np_tuneGoalDiff->addData(m_h_GoalDiffData_Test, "Data 2022", 0, false);
  np_tuneGoalDiff->addLable(0.2,0.70, "Goaliness Low: " + std::to_string(resultLowFine));
  np_tuneGoalDiff->addLable(0.2,0.75, "Goaliness High: " + std::to_string(resultHighFine));
  np_tuneGoalDiff->addLable(0.2,0.80, "#chi^{2}/DoF 14+18: " + std::to_string(m_bestChiGD_Training));
  //np_tuneGoalDiff->addLable(0.2,0.85, "#chi^{2}/DoF 2022: " + std::to_string(m_bestChiGD_Test));
  bookOutput::setBreak(1);
  bookOutput::get().doMultipadOutput("WCMC_Tuning", 2, 1);
  bookOutput::clear();

  if (m_mode == kFULL_TOURNAMENT) {
    np_base_1d->useAltColourScheme(2);
    np_base_1d->setLineWidth(2);
    np_base_1d->setLegend(.7, .75);
    np_base_1d->setBounds(-.5, 3.5, 0, 1.4);
    for (const std::string& group : group_letters) {
      nicePlot* np = new nicePlot(np_base_1d);
      np->init("Team", "Probability");
      np->normaliseToOne();
      np->addStackMC(m_h_roundWinner[group+"3"], "4th");
      np->addStackMC(m_h_roundWinner[group+"2"], "3rd");
      np->addStackMC(m_h_roundWinner[group+"1"], "Runner Up");
      np->addStackMC(m_h_roundWinner[group+"0"], "Winner");
      const std::vector<std::string>& teams = m_groups.at(group);
      for (const std::string& team : teams) np->addBinLabel(team);
      np->addLable(.25, .85, "Group " + group);
      np->addLable(.25, .80, "Winner: " + teams.at( m_h_roundWinner[group+"0"]->GetMaximumBin()-1 ) );
      np->addLable(.25, .75, "Runner Up: " + teams.at( m_h_roundWinner[group+"1"]->GetMaximumBin()-1 ) );
    }
    bookOutput::get().doMultipadOutput("WCMC_GroupResults", 2, 4);
    bookOutput::clear();
  }

  np_base_1d->setBounds(0, m_teams.size());
  np_base_1d->setDoLegend(false);
  double labelOffset = .3;
  if (m_mode < kAFTER_16) {
    np_base_1d->stretch = true;
    labelOffset = .7;
  }
  np_base_1d->useAltColourScheme(1);
  np_base_1d->setLineWidth(6);

  for (const std::string& team : m_teamsByRank) np_base_1d->addBinLabel(m_teams[team].m_abreviation);

  for (int i = (int)m_mode; i < 5; ++i) {
    nicePlot* np_round = new nicePlot(np_base_1d);
    np_round->normaliseToOne();
    np_round->n_mc += i;
    np_round->init("Team", "Probability");
    // np_round->setYBounds(0, 1. - (i * 0.90));
    
    TH1* hTemp = (TH1*) m_h_roundWinner[std::to_string(i)]->Clone();
    hTemp->Scale( (1./hTemp->Integral()) * numberOfPassingTeams );
    np_round->setYBounds(0, hTemp->GetMaximum() * 1.2);
    delete hTemp;

    np_round->addMC(m_h_roundWinner[std::to_string(i)], "");
    np_round->scaleLastMC(numberOfPassingTeams);
    numberOfPassingTeams /= 2;
    switch(i) {
      case 0: np_round->addLable(labelOffset, .8, "Probability of passing the Group Stage"); break; 
      case 1: np_round->addLable(labelOffset, .8, "Probability of passing the Round Of 16"); break; 
      case 2: np_round->addLable(labelOffset, .8, "Probability of passing the Quarter Finals"); break; 
      case 3: np_round->addLable(labelOffset, .8, "Probability of passing the Semi Finals"); break; 
      case 4: np_round->addLable(labelOffset, .8, "Probability of Winning The Tournament"); break; 
    }
  }
  bookOutput::setBreak(5 - (int)m_mode);
  bookOutput::get().doBookOutput("WCMC_KnockoutResults_Mode" + std::to_string((int)m_mode));
  bookOutput::clear();

  if (m_mode == kFULL_TOURNAMENT) {
    nicePlot* np_goals = new nicePlot(np_base_1d);
    np_goals->init("Team", "Number of Goals");
    np_goals->addMC(m_h_roundWinner["5"], "");
    np_goals->scaleLastMC(1./m_trialsMax);
    np_goals->setYBounds(0, 12.);
    np_goals->addLable(.6, .8, "Mean Number of Goals Scored Per Team");
    bookOutput::setBreak(1);
    bookOutput::get().doBookOutput("WCMC_Goals");
  }
}

int main() {
  gROOT->ProcessLine(".L AtlasStyle.C");
  gROOT->ProcessLine("SetAtlasStyle();");
  gErrorIgnoreLevel = 10000;
  //WCMC wc2022_a(kFULL_TOURNAMENT);
  //WCMC wc2018_b(kAFTER_GROUP);
  //WCMC wc2018_c(kAFTER_16);
  //WCMC wc2018_d(kAFTER_QUARTER);
  WCMC wc2018_e(kAFTER_SEMI);
}

int wcMC() {
  return main();
}