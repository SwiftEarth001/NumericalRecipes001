#include <string>
#include <vector>
#include <fstream>
#include <iostream>
#include <sstream>
#include <chrono>

#include "../../../include/ludcmp.h"


int main(int argc, char *argv[])
{
    using std::chrono::high_resolution_clock;
    using std::chrono::duration_cast;
    using std::chrono::duration;
    using std::chrono::milliseconds;

    std::string floats_file[6];
    floats_file[0] = "..\\data\\sample0.txt";  // 10 x 10
    floats_file[1] = "..\\data\\sample1.txt";  // 25 x 25
    floats_file[2] = "..\\data\\sample2.txt";  // 50 x 50
    floats_file[3] = "..\\data\\sample3.txt";  // 100 x 100
    floats_file[4] = "..\\data\\sample4.txt";  // 250 x 250
    floats_file[5] = "..\\data\\sample5.txt";  // 500 x 500

    std::string output_file[6];
    output_file[0] = "results\\sample0_results.txt";
    output_file[1] = "results\\sample1_results.txt";
    output_file[2] = "results\\sample2_results.txt";
    output_file[3] = "results\\sample3_results.txt";
    output_file[4] = "results\\sample4_results.txt";
    output_file[5] = "results\\sample5_results.txt";

    int sizes[6] = { 10, 25, 50, 100, 250, 500 };

    // in case we want to clock the factorisation of
    //   just one of the files;
    int filenum = 0x0FFFFFFF;
    if (argc == 2) { filenum = std::stoi(argv[1]); }

    for (int m=0; m<6; m++) {
        // check if only a single file is specified;
        if ((filenum > -1) && (filenum < 6)) { m=filenum; }
        // get data from file;
        std::ifstream inFile(floats_file[m]);
        std::string line;
        std::vector< std::vector<float> > all_floats;
        while ( getline(inFile, line) ) {
            std::istringstream is(line);
            all_floats.push_back(
                std::vector<float>( std::istream_iterator<float>(is),
                                    std::istream_iterator<float>() ) 
                );
        }
        inFile.close();

        // put data to array;
        int N = sizes[m];
        double* randfloats = new double[(N)*(N)];
        for (int i=0; i < N; i++) {
            for (int j=0; j < N; j++) {
                randfloats[i*N+j] = all_floats[i][j];
            }
        }

        // put array to linbcg class;
        LUdcmp* model = new LUdcmp(N, N, randfloats);
        // time the factorisation;
        auto t1 = high_resolution_clock::now();
        model->complete_pivot_factorise();
        auto t2 = high_resolution_clock::now();
        // get the number of milliseconds as a double;
        duration<double, std::milli> ms_double = t2 - t1;

        // write results to file;
        std::ofstream test05file;
        test05file.open( output_file[m] );
        for (int i=0; i<N; i++) {
            for (int j=0; j<N; j++) {
                test05file << std::to_string(model->lu[i][j]) << " ";
            }
            test05file << "\n";
        }
        test05file << "\n";
        for (int i=0; i<N; i++) {
            test05file << std::to_string(model->P[i]) << " ";
        }
        test05file << "\n";
        for (int i=0; i<N; i++) {
            test05file << std::to_string(model->Q[i]) << " ";
        }
        test05file << "\n";
        test05file << "duration: " << std::to_string(ms_double.count()) << " ms\n";
        test05file.close();

        delete[] randfloats;
        delete model;
        // break if only a single file is specified;
        if ((filenum > -1) && (filenum < 6)) { break; }
    }
}