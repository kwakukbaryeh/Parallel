#include <string.h>
#include <stdio.h>
#include <stdlib.h>
#include <mpi.h>

int main(int argc, char** argv) {
    // Initialize the MPI environment
    MPI_Init(&argc, &argv);

    // Get the number of processes
    int nproc;
    MPI_Comm_size(MPI_COMM_WORLD, &nproc);

    // Get the rank of the process
    int pid;
    MPI_Comm_rank(MPI_COMM_WORLD, &pid);

    // Define buffers
    int *in = (int *) malloc(sizeof(int) * nproc);

    // Allgather thingy
    MPI_Allgather(&pid, 1, MPI_INT, in, 1, MPI_INT, MPI_COMM_WORLD);

    // All processes print the received message
    
    printf("From proc %d: ", pid);
    for (int i = 0; i < nproc; i++)
        printf("%d", in[i]);
    printf("\n");

    // Finalize the MPI environment
    MPI_Finalize();
    return 0;
}