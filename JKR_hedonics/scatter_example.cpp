#include <iostream>
#include <cmath>
#include <stdlib.h>
#include "mpi.h"
#define MASTER 0
using namespace std;


int main(int argc, char *argv[]){
  
  int size, rank, perJobs;
  long int arg1, arg2, totalJobs;

 //Getting the inputs
  if(argc != 3){
    printf("This program requires 2 int to be passed\n");
  }
  arg1 = strtol(argv[1], NULL, 10);// Pretend this is randSeedStart
  arg2 = strtol(argv[2], NULL, 10);// Pretend this is randSeedEnd
  //printf("%d, %d\n", arg1, arg2);

  totalJobs = abs(arg2 - arg1);

  int *globalArry = new int[totalJobs];

  MPI_Init(&argc, &argv);  
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);
  MPI_Comm_size(MPI_COMM_WORLD, &size);

//Determine the amount of randSeeds each Processor receives
  if(totalJobs<size){perJobs = 1;}
  else{perJobs = ceil((float)(totalJobs+1)/size);}

  printf("Jobs per Processors %d from %d / %d\n", perJobs, totalJobs+1, size);

  int *localArry = new int[perJobs];

//Fills the Global Arry 
  int temp = arg1;
  if(rank == MASTER){
    for(int i = 0; i <= totalJobs; i++, temp++){
      globalArry[i] = temp; 
      //printf("%d, ", globalArry[i]);
    } 
    printf("\n");
  }
  
  //Destributes the array among the processors
  MPI_Scatter(globalArry, perJobs, MPI_INT,   /*Send "perJobs" number of elements from the globalArry*/ 
              localArry, perJobs, MPI_INT, /* Receive "perJobs" number of elements and put them into localArry*/
              MASTER, MPI_COMM_WORLD); /* Identify the root node and the Comm group*/


//Prints the rank's localArry
 /* if(rank == 2){                 
    printf("Processor %d", rank);
    for(int k = 0; k < perJobs; k++){
      printf(" %d,", localArry[k]);
    }
    printf("\n");
  }*/                                        

 
//Got tired here but this is how I would make sure processors do not run the same randSeeds

int randSeed = arg1;
 int z = 0;
 for(randSeed = arg1; randSeed <= arg2; randSeed++){
   if(localArry[z] == randSeed && z <= perJobs){
     z++;
     printf("Processor %d is performing randSeed %d\n", rank, randSeed);
   }
 } 


  delete [] localArry;                                       
  MPI_Finalize();
  delete [] globalArry;
  return 0;
}
