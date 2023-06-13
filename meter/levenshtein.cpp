#include "levenshtein.h"
#include <fstream>
#include <string>
#include <algorithm>
#include <vector>

using std::vector;
using std::ifstream;
using std::min;
using std::max;

int wagner_fischer(const string &first, const string &second) {
    int source_length;
    int result_length;
    source_length = first.length();
    result_length = second.length();
    vector<vector<int>> d_matrix(source_length + 1, vector<int>(result_length + 1));
    d_matrix[0][0] = 0;
    for(int i = 1; i <= source_length; ++i) {
        d_matrix[i][0] = i;
    }
    for(int i = 1; i <= result_length; ++i) {
        d_matrix[0][i] = i;
    }
    for(int i = 1; i <= source_length; ++i) {
        for(int j = 1; j <= result_length; ++j) {
            d_matrix[i][j] = min({d_matrix[i - 1][j] + 1, d_matrix[i][j - 1],
                                  d_matrix[i - 1][j - 1] + (first[i - 1] == second[j - 1] ? 0 : 1)});
        }
    }
    return d_matrix[source_length][result_length];
}

double editing_distance(const string &source_dna, const string &result_dna) {
    ifstream source(source_dna);
    ifstream result(result_dna);
    string in, out;
    source >> in;
    result >> out;
    result.close();
    source.close();
    return (1 - ((double)wagner_fischer(in, out) / (double)in.length())) * 100.0;
}
