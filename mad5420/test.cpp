#include <string>
#include <vector>
#include <fstream>
#include <iostream>
#include <sstream>

#include "inclstsqrs.h"


int main(int argc, char *argv[])
{
    int numrows, numcols, update_count=1;
    std::vector< std::vector<double> > all_doubles;
    std::string line;

    IncLeastSqrs* stream_data_model;
    double* xvar;

    if ( argc >= 4 ) {
        if (argv[1][0] == '-') {
            std::ofstream outfile;
            outfile.open("outfile.txt");
            if (argv[1][1] == 'i') {
                numrows = std::stoi(argv[2]);
                numcols = std::stoi(argv[3]);
                if (argc >= 5) { update_count = std::stoi(argv[4]); }
                std::istream* strm_ptr = &std::cin;
                unsigned long long int id = 0;
                while (std::getline(*strm_ptr, line)) {
                    std::istringstream iss(line);
                    all_doubles.push_back(
                        std::vector<double>( 
                            std::istream_iterator<double>(iss),
                            std::istream_iterator<double>() 
                        )
                    );
                    ++id;
                    printf("%d %d %d %llu\n", numrows, numcols, update_count, id);
                    if (id == numcols-1) {
                        double** initmodel = new double*[numcols-1];
                        for (int i=0; i<numcols-1; i++) { 
                            initmodel[i] = new double[numcols]; 
                        }
                        for (int i=0; i<numcols-1; i++) {
                            for (int j=0; j<numcols; j++) {
                                initmodel[i][j] = all_doubles[i][j];
                            }
                        }
                        stream_data_model = new IncLeastSqrs(numrows, numcols, initmodel, numcols-1);
                        while (all_doubles.size() > 0) {
                            all_doubles.pop_back();
                        }
                    } else if (id >= numcols) {
                        stream_data_model->new_data(all_doubles[0]);
                        all_doubles.pop_back();
                        if ((id-numcols+1) % (update_count) == 0) {
                            stream_data_model->solve();
                            xvar = stream_data_model->get_least_sqr();
                            for (int k=0; k < (id<numrows ? id : numrows); k++) {
                                outfile << xvar[k] << " ";
                            }
                            outfile << std::endl;
                        }
                    }
                }
            }           
            if (argv[1][1] == 'f') {
                std::string floats_file = argv[2];
                numrows = std::stoi(argv[3]);
                numcols = std::stoi(argv[4]);
                if (argc >= 5) { update_count = std::stoi(argv[5]); }
                std::ifstream inFile(floats_file);
                unsigned long long int id = 0;
                while (std::getline(inFile, line)) {
                    std::istringstream iss(line);
                    all_doubles.push_back(
                        std::vector<double>( 
                            std::istream_iterator<double>(iss),
                            std::istream_iterator<double>() 
                        )
                    );
                    ++id;
                    if (id == numcols-1) {
                        double** initmodel = new double*[numrows];
                        for (int k=0; k<numcols; k++) { 
                            initmodel[k] = new double[numcols]; 
                        }
                        for (int i=0; i<numrows; i++) {
                            for (int j=0; j<numcols; j++) {
                                initmodel[i][j] = all_doubles[i][j];
                            }
                        }
                        stream_data_model = new IncLeastSqrs(numrows, numcols, initmodel, numrows);
                        while (all_doubles.size() > 0) {
                            all_doubles.pop_back();
                        }
                    } else if (id >= numcols) {
                        stream_data_model->new_data(all_doubles[0]);
                        all_doubles.pop_back();
                        if ((id-numcols) % (update_count) == 0) {
                            stream_data_model->solve();
                            xvar = stream_data_model->get_least_sqr();
                            for (int k=0; k < (id<numrows ? id : numrows); k++) {
                                outfile << xvar[k] << " ";
                            }
                            outfile << std::endl;
                        }
                    }
                }
            }
            outfile.close();
        }
    } // else ;
}