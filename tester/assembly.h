#ifndef ASSEMBLYTESTING_ASSEMBLY_H
#define ASSEMBLYTESTING_ASSEMBLY_H

#include "utility.h"

void fill_graph(const data &input, int graph_size, vector<vector<int>> &graph);

void build_dna(string &output_filename, string &result_dna, const int &graph_size,
               const vector<string> &reads_vec, const vector<vector<int>> &graph);

#endif //ASSEMBLYTESTING_ASSEMBLY_H
