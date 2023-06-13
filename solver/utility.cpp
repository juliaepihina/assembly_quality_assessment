#include "utility.h"
#include <fstream>

void print_help() {
    cout << "-h - help about tool" << endl;
    cout << "-g [filename] - graph file" << endl;
    cout << "-o [filename] - output file" << endl;
}

void parse_input(int argc, char *argv[], cv_data &input) {
    int current = 1;
    while (current < argc) {
        if (string(argv[current]) == "-g") {
            input.graph_filename = string(argv[++current]);
        } else if (string(argv[current]) == "-o") {
            input.output_filename = string(argv[++current]);
        }
        ++current;
    }
}

vector<vector<int>> load_graph(cv_data &input) {
    ifstream file(input.graph_filename);
    file >> input.graph_size;
    vector<vector<int>> graph(input.graph_size, vector<int>(input.graph_size));
    for (int i = 0; i < input.graph_size; ++i) {
        for (int j = 0; j < input.graph_size; ++j) {
            file >> graph[i][j];
        }
    }
    file.close();
    return graph;
}

void save_state(vector<int> &state, string &output) {
    ofstream file(output);
    for (int i : state) {
        file << i << " ";
    }
    file.close();
}
