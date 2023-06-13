#ifndef ASSEMBLYTESTING_SIMULATED_ANNEALING_H
#define ASSEMBLYTESTING_SIMULATED_ANNEALING_H

#include <ctime>
#include <vector>
#include <string>
#include <cmath>
#include <numeric>
#include <algorithm>
#include <iostream>

using namespace std;

vector<int> simulated_annealing(double initial_temperature, double end_temperature,
        int graph_size, const vector<vector<int>> &graph);

#endif //ASSEMBLYTESTING_SIMULATED_ANNEALING_H
