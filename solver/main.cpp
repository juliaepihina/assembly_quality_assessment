#include "utility.h"

int main(int argc, char *argv[]) {
    cv_data data = cv_data();
    if (argc == 5) {
        parse_input(argc, argv, data);
        vector<vector<int>> graph = load_graph(data);
        vector<int> state = simulated_annealing(100, 0.001, data.graph_size, graph);
        save_state(state, data.output_filename);
    } else if (argc == 2) {
        if (string(argv[1]) == "-h") {
            print_help();
        }
    } else {
        cout << "bad input in solver" << endl;
        print_help();
    }
}
