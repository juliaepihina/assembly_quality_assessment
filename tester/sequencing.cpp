#include "sequencing.h"
#include <ctime>

void generate_sequence(int sequence_length, string &source_dna) {
    srand(time(nullptr));
    ofstream source(source_dna);
    string pattern = "atgc";
    string result_string;
    for(int i = 0; i < sequence_length; ++i) {
        result_string += pattern[rand()%4];
    }
    source << result_string;
    source.close();
}

void replicate(string &source_dna, int n) {
    ifstream source(source_dna);
    ofstream repl;
    for(int i = 0; i < n; ++i) {
        source.seekg(0, ifstream::beg);
        repl.open(source_dna + "_replica_" + to_string(i) + ".txt");
        repl << source.rdbuf();
        repl.close();
    }
    source.close();
}

void reads_create(const data &input, int &graph_size, vector<string> &reads_vec) {
    srand(time(nullptr));
    ifstream source(input.source_dna);
    ofstream reads(input.reads_filename);
    string in;
    source >> in;
    int delta = input.max_reads_length - input.min_reads_length;
    int current_index = 0;
    int current_random = 0;
    int ind = 0;
    while(current_index < in.length()) {
        current_random = input.min_reads_length + rand()%delta;
        reads_vec.emplace_back(in.substr(current_index, min(current_random, (int)in.length() - current_index)));
        reads <<  reads_vec[ind] << endl;
        current_index += current_random - rand()%current_random;
        ++ind;
    }
    graph_size = ind;
    reads.close();
    source.close();
}
