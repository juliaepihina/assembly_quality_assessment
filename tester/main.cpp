#include "sequencing.h"
#include "assembly.h"


#include <iostream>
#include <vector>
#include <numeric>
#include <map>
#include <algorithm>
#include <tuple>
#include <string>
#include <cstring>
#include <unordered_map>
#include <cctype>
#include <sstream>
using namespace std;

std::string read_all(std::ifstream& in) {
    std::ostringstream sstr;
    sstr << in.rdbuf();
    return sstr.str();
}


const int MAXN = 10000; // ������������ ����� �������������������
const int INF = 1e9; // ����������� ��������

int dp[MAXN][MAXN]; // DP �������

// ������� ��� ������� ���������� ��������-�����
tuple<int, string, string> needlemanWunsch(string seq1, string seq2, int gap_penalty, int mismatch_score, int match_score) {

    int n = seq1.length();
    int m = seq2.length();

    // �������������� ������� DP
    for (int i = 0; i <= n; i++) {
        for (int j = 0; j <= m; j++) {
            dp[i][j] = -INF;
        }
    }

    // ������� ������
    dp[0][0] = 0;
    for (int i = 1; i <= n; i++) {
        dp[i][0] = dp[i - 1][0] + gap_penalty;
    }
    for (int j = 1; j <= m; j++) {
        dp[0][j] = dp[0][j - 1] + gap_penalty;
    }

    // ��������� ������� DP
    for (int i = 1; i <= n; i++) {
        for (int j = 1; j <= m; j++) {
            int score = (seq1[i - 1] == seq2[j - 1]) ? match_score : mismatch_score; // ������������ ����
            dp[i][j] = max({ dp[i - 1][j - 1] + score, dp[i - 1][j] + gap_penalty, dp[i][j - 1] + gap_penalty }); // ��������� ������� DP
        }
    }

    // ������������ �����, ����� �������� ������������
    string alignedS1 = "";
    string alignedS2 = "";
    int i = n;
    int j = m;
    while (i > 0 || j > 0) {
        if (i > 0 && dp[i][j] == dp[i - 1][j] + gap_penalty) {
            alignedS1 = seq1[i - 1] + alignedS1;
            alignedS2 = "-" + alignedS2;
            i--;
        }
        else if (j > 0 && dp[i][j] == dp[i][j - 1] + gap_penalty) {
            alignedS1 = "-" + alignedS1;
            alignedS2 = seq2[j - 1] + alignedS2;
            j--;
        }
        else {
            alignedS1 = seq1[i - 1] + alignedS1;
            alignedS2 = seq2[j - 1] + alignedS2;
            i--;
            j--;
        }
    }

    return make_tuple(dp[n][m], alignedS1, alignedS2); 
}

vector<int> countFragmentsLengths(string s) {
    vector<int> lengths;
    int count = 0;
    for (int i = 0; i < s.length(); i++) {
        if (s[i] == '-') {
            lengths.push_back(count);
            count = 0;
        }
        else {
            count++;
        }
    }
    lengths.push_back(count);   // ��������� ����� ���������� ���������
    return lengths;
}

double calculate_N50(std::vector<int> contig_lengths) {
    double total_length = 0;
    for (int i = 0; i < contig_lengths.size(); i++) {
        total_length += contig_lengths[i];
    }
    std::sort(contig_lengths.begin(), contig_lengths.end());
    double half_length = total_length / 2;
    double current_length = 0;
    for (int i = contig_lengths.size() - 1; i >= 0; i--) {
        current_length += contig_lengths[i];
        if (current_length >= half_length) {
            return contig_lengths[i];
        }
    }
    return 0;
}

double calculate_L50(vector<int> contig_lengths)
{
    sort(contig_lengths.begin(), contig_lengths.end(), greater<int>());
    int total_sum = accumulate(contig_lengths.begin(), contig_lengths.end(), 0);
    int half_sum = total_sum / 2;
    int current_sum = 0;
    int L50 = 0;
    for (int value : contig_lengths) {
        current_sum += value;
        L50 += 1;
        if (current_sum >= half_sum) {
            break;
        }
    }
    return L50;
}


void drawHistogram(string dna, int k) {
    unordered_map<string, int> kmerCounts;
    int maxCount = 0;

    for (int i = 0; i <= dna.length() - k; i++) {
        string kmer = dna.substr(i, k);
        kmerCounts[kmer]++;
        maxCount = max(maxCount, kmerCounts[kmer]);
    }


    for (int i = maxCount; i > 0; i--) {
        cout << i << "| ";
        for (auto& kmerCount : kmerCounts) {
            if (kmerCount.second >= i) {
                cout << "* ";
            }
            else {
                cout << "  ";
            }
        }
        cout << endl;
    }

    cout << "  +";
    for (auto& kmerCount : kmerCounts) {
        cout << "--";
    }
    cout << endl << "   ";
    for (auto& kmerCount : kmerCounts) {
        cout << kmerCount.first << " ";
    }
    cout << endl;
}

double calculate_Q_value(string genome, int k, int pl, int pr) {
    unordered_map<string, int> kmer_freq;
    unordered_map<string, int> kmer_freq_complement;
    unordered_map<string, int> kmer_freq_assembly;

    int n = genome.length();
    for (int i = 0; i <= n - k; i++) {
        string kmer = genome.substr(i, k);
        string complement = "";
        for (int j = k - 1; j >= 0; j--) {
            if (kmer[j] == 'A') {
                complement += 'T';
            }
            else if (kmer[j] == 'T') {
                complement += 'A';
            }
            else if (kmer[j] == 'C') {
                complement += 'G';
            }
            else if (kmer[j] == 'G') {
                complement += 'C';
            }
        }

        if (kmer_freq.find(kmer) == kmer_freq.end()) {
            kmer_freq[kmer] = 1;
        }
        else {
            kmer_freq[kmer]++;
        }

        if (kmer_freq_complement.find(complement) == kmer_freq_complement.end()) {
            kmer_freq_complement[complement] = 1;
        }
        else {
            kmer_freq_complement[complement]++;
        }
    }

    for (auto itr = kmer_freq.begin(); itr != kmer_freq.end(); itr++) {
        if (itr->second == 1 && kmer_freq_complement[itr->first] == 1) {
            kmer_freq_assembly[itr->first] = 1;
        }
    }

    vector<int> kmer_freq_count;
    for (auto itr = kmer_freq_assembly.begin(); itr != kmer_freq_assembly.end(); itr++) {
        kmer_freq_count.push_back(itr->second);
    }
    sort(kmer_freq_count.begin(), kmer_freq_count.end());

    int kmer_freq_count_size = kmer_freq_count.size();
    double Q = 0;
    for (int i = 0; i < kmer_freq_count_size; i++) {
        if (i + 1 >= pl && i + 1 <= pr) {
            Q += count(kmer_freq_count.begin(), kmer_freq_count.end(), kmer_freq_count[i]);
        }
    }
    Q /= kmer_freq_count_size;

    return Q;
}

void grade(string sf, string rf) {

    std::ifstream f(sf);
    std::ifstream f1(rf);

    string seq1 = read_all(f);
    for(auto& i : seq1) i = toupper(i);
    string seq2 = read_all(f1);
    for(auto& i : seq2) i = toupper(i);
    int gap_penalty = -1;
    int match_score = 2;
    int mismatch_score = -1;
    tuple<int, string, string> result = needlemanWunsch(seq1, seq2, gap_penalty, mismatch_score, match_score);
    string res2 = get<2>(result);
    int count = 1;
    for (int i = 0; i < res2.length(); i++) {
        if (res2[i] == '-') {
            count++;
        }
    }

    //cout << "��������� ��������� ��������-�����: " << res2 << "\n";
    
    double N50_1 = calculate_N50(countFragmentsLengths(res2)) / seq2.length();
    cout << "N50: " << N50_1 << "\n";

    double L50_1 = 1 - calculate_L50(countFragmentsLengths(res2)) / count;
    cout << "L50: " << L50_1 << "\n";

    cout << "����������� \n";
    for (int k = 2; k < seq2.length(); k++) {
        drawHistogram(seq2, k);
    }

    int pl = 0;
    int pr = 10;
    int Q1;

    double Qmax = -1000000000;

    for(int i = 0; i < seq2.length() - 1; i++) {
        Q1 = calculate_Q_value(seq2, i, pl, pr);
        if(Q1 > Qmax)
            Qmax = Q1;
    }

    Q1 = Qmax;

    cout << "Q1 = " << Q1 << "\n";

    double grade = 0.5 * Q1 + 0.25 * N50_1 + 0.25 * L50_1;
    cout << "������: " << grade;
}


int main(int argc, char *argv[]) {
    data input = data();
    vector<vector<int>> graph;
    vector<string> reads_vec;
    int graph_size = 0;
    if (argc == 2) {
        if (string(argv[1]) == "-h") {
            print_help();
        } else if (string(argv[1]) == "-ar") {
            auto_remove();
        }
    } else if (argc == 10) {
        parse_input(argc, argv, input);
        if (input.is_rand) {
            generate_sequence(input.sequence_length, input.source_dna);
        }

        //std::cout << "Source DNA: " << input.source_dna;

        reads_create(input, graph_size, reads_vec);
        fill_graph(input, graph_size, graph);
        save_graph(graph_size, graph, input.graph_filename);
        system((input.tested_program + " -g " + input.graph_filename + " -o " + input.output_filename).c_str());
        build_dna(input.output_filename, input.result_dna, graph_size, reads_vec, graph);
        system((input.meter_program + " " + input.source_dna + " " + input.result_dna).c_str());

        //std::cout << "Result DNA: " << input.result_dna;
        grade(input.source_dna, input.result_dna);


    } else {
        cout << "bad input in tester" << endl;
        print_help();
    }
}
