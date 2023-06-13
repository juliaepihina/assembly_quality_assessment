#ifndef ASSEMBLYTESTING_SEQUENCING_H
#define ASSEMBLYTESTING_SEQUENCING_H

#include "utility.h"

void generate_sequence(int sequence_length, string &source_dna);

void replicate(string &source_dna, int n);

void reads_create(const data &input, int &graph_size, vector<string> &reads_vec);

#endif //ASSEMBLYTESTING_SEQUENCING_H
