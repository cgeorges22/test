
#include <math.h>
#include <stdlib.h>
#include <stdio.h>
#include <iostream>
#include <fstream>
#include <algorithm>
#include <string>  //JKR 10/7/16
#include <sstream> //JKR 10/7/16
//#include "hedonics-firm-v2.10.h"
//#include "hedonics-utilities-v2.10.h"
#include "kat_test.h"
//#include "kat_utilities.h"
#include "time.h"
#include "mpi.h"
#include <cassert> //7/4/09
#include <cstring>
#include <deque> //JKR 10/31/16
#include <cmath>
#define MASTER 0
#define WORKTAG 2
#define DIETAG 3
#define NEEDWORK 4
using namespace std;


double master(int numSearchTrials);
void slave(int rank);
void getSeeds(int& start, int& end);

int main(){
  
  int size, rank, perJobs, utility;

  MPI_Init(0, 0);  
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);
  MPI_Comm_size(MPI_COMM_WORLD, &size);
 
  if(rank == MASTER) simulation();
  else slave(rank);
 
  MPI_Finalize();

  return 0;

}

void slave(int rank){

   double utility;
   MPI_Status status;

   for(;;){

   MPI_Recv(&utility,1,MPI_DOUBLE,MASTER,MPI_ANY_TAG,MPI_COMM_WORLD,&status);
      
   if(status.MPI_TAG == DIETAG){
     return;
   }
  
   utility = SearchTrialLoop(utility);

   //printf("processor %d has randSeed %d\n", rank, randSeed); 

   MPI_Send(&utility,0,MPI_DOUBLE,MASTER,NEEDWORK,MPI_COMM_WORLD);
   }
}


void getSeeds(int& start, int& end){
  char line[256];
  char* variable;
  char* value;
  // Uses the file "input.txt"
  std::fstream input;
  input.open("input.txt", std::fstream::in);
  
  // Go through the input file line by line, tokenize the line into three strings
  input.getline(line, 256);
  variable = strtok(line, " ");

  while (strcmp(variable, "randSeedStart")) {
    input.getline(line, 256);
    variable = strtok(line, " ");
  }
  start = atoi(strtok(NULL, " "));
  input.getline(line, 256);
  strtok(line, " ");
  end = atoi(strtok(NULL, " ")); 
}

