#include <iostream>
#include <vector>
#include <time.h>
#include <vector>
#include <chrono>
#include "mpi.h"
using namespace std;

void consecutive(const int N) {
    const double S1 = 1.0;
    vector <double> A, B, C, Y;
    A.resize(N); B.resize(N); C.resize(N); Y.resize(N);

    chrono::high_resolution_clock::time_point timeStart = chrono::high_resolution_clock::now();
    for (int i = 0; i < N; i++) {
        A[i] = i;
        B[i] = i * 2;
        C[i] = (i + 1) / 2;
        Y[i] = (A[i] + B[i] * S1 + C[i]) * A[i];
    }
    chrono::high_resolution_clock::time_point timeEnd = chrono::high_resolution_clock::now();
    chrono::duration<double, milli> milli_diff = timeEnd - timeStart;

    cout << "  Time (ms): " << milli_diff.count() << endl << endl << "  Result:" << endl;
    for (int i = 0; i < N; i++) 
        cout << "  " << Y[i] << endl;
   
}

int main(int argc, char** argv) {
    const int N = 12;
    
    int ProcRank, size;
    vector <double> A, B, C, Y;
    A.resize(N); B.resize(N); C.resize(N); Y.resize(N);
    const double S1 = 1.0;
    MPI_Status status;
    MPI_Init(&argc, &argv);
    MPI_Comm_rank(MPI_COMM_WORLD, &ProcRank);
    MPI_Comm_size(MPI_COMM_WORLD, &size);
    if (ProcRank == 0) {

        for (int i = 0; i < N / size; i++) {
            cout << "rank: " << ProcRank << "  i: " << i << endl;
            A[i] = i;
            B[i] = i * 2;
            C[i] = (i + 1) / 2;
            Y[i] = (A[i] + B[i] * S1 + C[i]) * A[i];
        }

        for (int i = 1; i < size; i++) // Сбор результатов 
            MPI_Recv(&Y[(int)i * N / size], (int)N / size, MPI_DOUBLE, i, 0, MPI_COMM_WORLD, &status);
       
        cout << endl << endl;
        for (int i = 0; i < N; i++)
            cout << "  " << Y[i] << endl;

        cout << endl << "  Consecutive programm" << endl;
        consecutive(N); // Последовательная программа
    }
    else {
        for (int i = ProcRank * N / size; i < (ProcRank + 1) * N / size; i++) {
            cout << "rank: " << ProcRank << "  i: " << i << endl;
            A[i] = i;
            B[i] = i * 2;
            C[i] = (i + 1) / 2;
            Y[i] = (A[i] + B[i] * S1 + C[i]) * A[i];
        }
        MPI_Send(&Y[(int)ProcRank * N / size], (int)N / size, MPI_DOUBLE, 0, 0, MPI_COMM_WORLD);
    }

    MPI_Finalize();

    return 0;
}
