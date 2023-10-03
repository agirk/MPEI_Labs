#include <iostream>
#include <vector>
#include <time.h>
#include <vector>
#include "mpi.h"
using namespace std;

void consecutive(const int N) {
    const double S1 = 2.5;
    vector <double> A, B, C, Y;
    A.resize(N); B.resize(N); C.resize(N); Y.resize(N);

    for (int i = 0; i < N; i++) {
        A[i] = i;
        B[i] = i * 2;
        C[i] = (i + 1) / 2;
        Y[i] = (A[i] + B[i] * S1 + C[i]) * A[i];
    }

    for (int i = 0; i < N; i++) 
        cout << "  " << Y[i] << endl; 
}

int main(int argc, char** argv) {
    int ProcRank, size;
    vector <double> A, B, C, Y;
    vector <int> N_val = { 8, 64, 256, 1024, 4096 };
    const double S1 = 2.5;
    double timeStart, timeEnd;
    MPI_Status status;
    MPI_Init(&argc, &argv);
    MPI_Comm_rank(MPI_COMM_WORLD, &ProcRank);
    MPI_Comm_size(MPI_COMM_WORLD, &size);

    for (int k = 0; k < N_val.size(); k++) {
        MPI_Barrier(MPI_COMM_WORLD);
        A.resize(N_val[k]); B.resize(N_val[k]); C.resize(N_val[k]); Y.resize(N_val[k]);

        if (ProcRank == 0) {
            timeStart = MPI_Wtime();
            for (int circle = 0; circle < 1000; circle++) {
                for (int i = 0; i < N_val[k] / size; i++) {
                    //cout << "rank: " << ProcRank << "  i: " << i << endl;
                    A[i] = i;
                    B[i] = i * 2;
                    C[i] = (i + 1) / 2;
                    Y[i] = (A[i] + B[i] * S1 + C[i]) * A[i];
                }

                for (int i = 1; i < size; i++) { // Сбор результатов 
                    MPI_Recv(&Y[(int)i * N_val[k] / size], (int)N_val[k] / size, MPI_DOUBLE, i, 0, MPI_COMM_WORLD, &status);
                }
            }

            timeEnd = MPI_Wtime();

            cout << "Time: " << timeEnd - timeStart << endl << endl;// << "  Result:" << endl;
            //  for (int i = 0; i < N; i++)
            //      cout << "  " << Y[i] << endl;

              //cout << endl << endl <<  "  Consecutive programm" << endl;
              //consecutive(N); // Последовательная программа
        }
        else {
            for (int circle = 0; circle < 1000; circle++) {
                for (int i = ProcRank * N_val[k] / size; i < (ProcRank + 1) * N_val[k] / size; i++) {
                    //cout << "rank: " << ProcRank << "  i: " << i;
                    A[i] = i;
                    B[i] = i * 2;
                    C[i] = (i + 1) / 2;
                    Y[i] = (A[i] + B[i] * S1 + C[i]) * A[i];
                }
                MPI_Send(&Y[(int)ProcRank * N_val[k] / size], (int)N_val[k] / size, MPI_DOUBLE, 0, 0, MPI_COMM_WORLD);
            }
        }
    }

    MPI_Finalize();

    return 0;
}
