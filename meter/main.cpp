#include "ssw_cpp.h"
#include <iostream>
#include <string.h>
#include <fstream>

using namespace std;
using namespace StripedSmithWaterman;

inline void calculate(char *argv[]);

int main(int argc, char *argv[]) {
    if (argc == 3) {
        // levenshtein
        // cout << "rate: " << editing_distance(string(argv[1]), string(argv[2])) << '%' << endl;
        calculate(argv);
    } else {
        cout << "bad input in meter" << endl;
    }
}

inline void calculate(char *argv[]) {
    ifstream file1(argv[1]);
    string ref, query;
    ifstream file2(argv[2]);
    file1 >> ref;
    file2 >> query;
    file1.close();
    file2.close();
    int32_t maskLen = strlen(query.c_str())/2;
    maskLen = maskLen < 15 ? 15 : maskLen;
    StripedSmithWaterman::Aligner aligner;
    StripedSmithWaterman::Filter filter;
    StripedSmithWaterman::Alignment alignment;
    aligner.Align(query.c_str(), ref.c_str(), ref.size(), filter, &alignment, maskLen);
    cout << (alignment.sw_score / 2) / ref.length() * 100 << "%" << endl;
}
