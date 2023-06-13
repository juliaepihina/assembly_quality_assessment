#include "assembly.h"

int intersect(const string& a, const string& b) {
    int length_first = a.length();
    int min_length = min((int)b.length(), length_first);
    int res = 0;
    while(true) {
        if (res >= min_length) {
            break;
        }
        if (a[length_first - res - 1] == b[res]) {
            ++res;
        } else {
            break;
        }
    }
    return res;
}

void fill_graph(const data &input, int graph_size, vector<vector<int>> &graph) {
    ifstream reads(input.reads_filename);
    string first, second;
    string temp;
    graph.resize(graph_size);
    for (int j = 0; j < graph_size; ++j) {
        graph[j].assign(graph_size, input.max_reads_length);
    }
    for(int i = 0; i < graph_size - 1; ++i) {
        if (!i) {
            getline(reads, temp);
            first = temp;
        }
        getline(reads, temp);
        second = temp;
        graph[i+1][i] = input.max_reads_length - intersect(first, second);
        graph[i][i+1] = graph[i+1][i];
        graph[i][i+1] = graph[i+1][i];
        first = second;
    }
}

void build_dna(string &output_filename, string &result_dna, const int &graph_size,
               const vector<string> &reads_vec, const vector<vector<int>> &graph) {
    vector<int> state(graph_size);
    ifstream output(output_filename);
    ofstream result(result_dna);
    for (int i = 0; i < graph_size; ++i) {
        output >> state[i];
    }
    for(int i = 0; i < graph_size; ++i){
        if (i != graph_size - 1) {
            result << reads_vec[i].substr(0, reads_vec[i].length() - graph[state[i]][state[i+1]]);
        } else {
            result << reads_vec[i];
        }
    }
    output.close();
    result.close();
}
