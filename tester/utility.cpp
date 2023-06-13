#include "utility.h"

void print_help() {
    cout << "-h - help about tool" << endl;
    cout << "-ar - remove all *.txt file in work directory" << endl;
    cout << "-r [length] - generate random sequence" << endl;
    cout << "-e [filename] - use existing sequence" << endl;
    cout << "-rl [min][max] - min and max reads length" << endl;
    cout << "-t [exe name] - name of tested program" << endl;
    cout << "-m [exe name] - name of meter" << endl;
}

void parse_input(int argc, char *argv[], data &input) {
    int current = 1;
    while (current < argc) {
        if (string(argv[current]) == "-r") {
            input.is_rand = true;
            input.source_dna = "source.txt";
            input.sequence_length = atoi(argv[++current]);
        } else if (string(argv[current]) == "-e") {
            input.is_rand = false;
            input.source_dna = string(argv[++current]);
        } else if (string(argv[current]) == "-rl") {
            input.min_reads_length = atoi(argv[++current]);
            input.max_reads_length = atoi(argv[++current]);
            if (input.max_reads_length < input.min_reads_length ||
                input.min_reads_length < 1 ||
                input.max_reads_length < 1) {
                cout << "-rl bad values" << endl;
                exit(0);
            }
        } else if (string(argv[current]) == "-t") {
            input.tested_program = string(argv[++current]);
        } else if (string(argv[current]) == "-m") {
            input.meter_program = string(argv[++current]);
        }
        ++current;
    }
    input.reads_filename = "reads.txt";
    input.graph_filename = "graph.txt";
    input.output_filename = "algo_output.txt";
    input.result_dna = "result.txt";
}

void save_graph(int graph_size, vector<vector<int>> &graph, const string& filename) {
    ofstream file(filename);
    file << graph_size << endl;
    for (auto & i : graph) {
        for (int j : i) {
            file << j << ' ';
        }
        file << endl;
    }
    file.close();
}

void auto_remove() {
#ifdef __linux__
    system("rm *.txt");
#elif _WIN32
    system("del *.txt");
#endif
}
