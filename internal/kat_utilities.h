
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
//#include "kat_test.h"
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


double SearchTrialLoop(double utility){
	utility++;
	printf("here: %d\n", utility);
	return utility;
}


double master(int numSearchTrials){
  int totalJobs, rank, ntasks;
  deque<int>globalQueue;
  totalJobs = numSearchTrials;

  MPI_Status status;
  MPI_Comm_size(MPI_COMM_WORLD, &ntasks);
 

//Fills the Global Queue                            
  for(int i = 0; i <= totalJobs; i++){
    globalQueue.push_back(i); 
  //  printf("%d, ", globalQueue[i]);
   } 
//   printf("\n");


  double utility = 0, altUtility = 0;
    
    //intialize mass seed distribution
    for(rank = 1; rank < ntasks; rank++){
      if(!globalQueue.empty()){
        //int buffer = globalQueue.front();
        globalQueue.pop_front();
        MPI_Send(&utility, 1, MPI_DOUBLE,rank, 1, MPI_COMM_WORLD);
      }
    }
 
  while(!globalQueue.empty()){

   globalQueue. pop_front();
   MPI_Recv(&utility, 0, MPI_DOUBLE, MPI_ANY_SOURCE, NEEDWORK,MPI_COMM_WORLD,&status);

   if (utility > altUtility){
	altUtility = utility;
   }

   MPI_Send(&utility, 1, MPI_DOUBLE, status.MPI_SOURCE, WORKTAG, MPI_COMM_WORLD);

  }

  int waiting_on = numSearchTrials < ntasks ? numSearchTrials : ntasks;

  for(rank = 1; rank < waiting_on; rank++ ){
    MPI_Recv(&utility, 0, MPI_DOUBLE, MPI_ANY_SOURCE, NEEDWORK,MPI_COMM_WORLD,&status);	

   if (utility > altUtility){
	altUtility = utility;
   }

    MPI_Send(0,0,MPI_DOUBLE, rank, DIETAG, MPI_COMM_WORLD);
  } 

  
    
   return altUtility;
  }//end of master function

int func1(){

	int numSearchTrials = 10, max_u;

	// call master 
	max_u = master(numSearchTrials);

	printf("%d\n", max_u);
}

