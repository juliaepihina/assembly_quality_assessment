#ifndef ASSEMBLYTESTING_UTILITY_H
#define ASSEMBLYTESTING_UTILITY_H

#include "simulated_annealing.h"

struct cv_data {
    string graph_filename;
    string output_filename;
    int graph_size;
};

void print_help();

void parse_input(int argc, char *argv[], cv_data &input);

vector<vector<int>> load_graph(cv_data &input);

void save_state(vector<int> &state, string &output);

#endif //ASSEMBLYTESTING_UTILITY_H
