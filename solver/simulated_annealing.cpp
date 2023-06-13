#include "simulated_annealing.h"
#include <random>

double decrease_temperature(double initial_temperature, double i) {
    return initial_temperature * 0.1 / i;
}

double get_transition_probability(double delta_energies, double temp) {
    return exp(-delta_energies/temp);
}

bool is_transition(double probability) {
    srand(time(nullptr));
    double val = ((double)rand() / (RAND_MAX)) + 1;
    return val <= probability;
}

vector<int> generate_state_candidate(const vector<int> &seq, const int graph_size) {
    vector<int> res = seq;
    int size_seq = graph_size;
    int i = rand() % size_seq;
    int j = rand() % size_seq;
    i < j ? reverse(res.begin() + i, res.begin() + j) : reverse(res.begin() + j, res.begin() + i);
    return res;
}

int calculate_energy(vector<int> &state, const int graph_size, const vector<vector<int>> &graph) {
    int res = 0;
    for(int i = 0; i < graph_size - 1; ++i) {
        res += graph[state[i]][state[i+1]];
    }
    return res;
}

vector<int> simulated_annealing(double initial_temperature, double end_temperature,
        int graph_size, const vector<vector<int>> &graph) {
    vector<int> state(graph_size);
    iota(state.begin(), state.end(), 0);
    shuffle(state.begin(), state.end(), std::mt19937(std::random_device()()));

    int current_energy = calculate_energy(state, graph_size, graph);
    double temperature = initial_temperature;

    for(int i = 0; i < 100000; ++i) {
        vector<int> candidate_state = generate_state_candidate(state, graph_size);
        int candidate_energy = calculate_energy(candidate_state, graph_size, graph);

        if (candidate_energy < current_energy) {
            current_energy = candidate_energy;
            state = candidate_state;
        } else {
            double probability = get_transition_probability(candidate_energy - current_energy, temperature);
            if (is_transition(probability)){
                current_energy = candidate_energy;
                state = candidate_state;
            }
        }

        temperature = decrease_temperature(initial_temperature, 1);

        if(temperature < end_temperature) {
            break;
        }
    }

    return state;
}
