#include <iostream>
#include <cmath>
#include <stdlib.h>
#include <deque>
#include "mpi.h"
#define MASTER 0
#define WORKTAG 2
#define DIETAG 3
#define NEEDWORK 4
using namespace std;


void master(int arg1, int arg2);
void slave(int rank);

int main(int argc, char *argv[]){
  
  int size, rank, perJobs;
  long int arg1, arg2, totalJobs;

 //Getting the inputs
  if(argc != 3){
    printf("This program requires 2 int to be passed\n");
  }
  arg1 = strtol(argv[1], NULL, 10);// Pretend this is randSeedStart
  arg2 = strtol(argv[2], NULL, 10);// Pretend this is randSeedEnd

  MPI_Init(&argc, &argv);  
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);
  MPI_Comm_size(MPI_COMM_WORLD, &size);
 
  if(rank == MASTER) master(arg1, arg2);
  else slave(rank);
 
  MPI_Finalize();
  return 0;
}


void master(int arg1, int arg2){
  int totalJobs, rank, ntasks;
  deque<int>globalQueue;
  totalJobs = abs(arg2-arg1);

  MPI_Status status;
  MPI_Comm_size(MPI_COMM_WORLD, &ntasks);
 

//Fills the Global Queue
  int temp = arg1;                                    
  for(int i = 0; i <= totalJobs; i++, temp++){
    globalQueue.push_back(temp); 
  //  printf("%d, ", globalQueue[i]);
   } 
//   printf("\n");
    
    //intialize mass seed distribution
    for(rank = 1; rank < ntasks; rank++){
      if(!globalQueue.empty()){
        int buffer = globalQueue.front();
        globalQueue.pop_front();
        MPI_Send(&buffer, 1, MPI_INT,rank, 1, MPI_COMM_WORLD);
      }
    }




  while(!globalQueue.empty()){

   int buffer = globalQueue.front();
   globalQueue.pop_front();
   MPI_Recv(0, 0, MPI_INT, MPI_ANY_SOURCE, NEEDWORK,MPI_COMM_WORLD,&status);

   MPI_Send(&buffer, 1, MPI_INT, status.MPI_SOURCE, WORKTAG, MPI_COMM_WORLD);

  }

  for(rank = 1; rank < ntasks; rank++ ){
    MPI_Send(0,0,MPI_INT, rank, DIETAG, MPI_COMM_WORLD);
  } 

  
    
   return;
  }//end of master function


void slave(int rank){

   int randSeed;
   MPI_Status status;

   for(;;){

   MPI_Recv(&randSeed,1,MPI_INT,MASTER,MPI_ANY_TAG,MPI_COMM_WORLD,&status);
      
   if(status.MPI_TAG == DIETAG){
     return;
   }
   printf("processor %d has randSeed %d\n", rank, randSeed); 

   MPI_Send(0,0,MPI_INT,MASTER,NEEDWORK,MPI_COMM_WORLD);
   }

}
