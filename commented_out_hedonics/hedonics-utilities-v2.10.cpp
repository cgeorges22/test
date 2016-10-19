//procedures for firms >= v2.4.1 
//broken out of driver.cpp 5-15-07

#include <math.h>
#include <stdlib.h>
#include <stdio.h>
#include <iostream>
#include <fstream>
#include <algorithm> //8/13/08
#include "hedonics-firm-v2.10.h"
#include "hedonics-utilities-v2.10.h"

//Function Definitions

//initialize firms
void initializeFirms(firm firms[], const int firmNum, bool twoClasses, std::ofstream &testoutput, bool anotation, double probComplements, double probInnovate, double probIntegerConstrained, double intPriceDiscount, double wage, double salary, double markup) { //moved from main 10-4-10 //int round removed 6-14-15
        
    int i,j,k;
    for (i = 0; i < firmNum; i++) {
      firm newFirm;
      
      //set consumer good type // 5/11/16
      if (i < (int) firmNum/2) {
        newFirm.setType(1);
      } // for 1/2 firms of each type at start of sim
      else {
        newFirm.setType(2);
      }
      //set firm demand shares 
      if (twoClasses) {	
	if (newFirm.getType()==1) {
	  newFirm.setFinalDemandShare((2.0)/((double) firmNum));
	} //here for 1/2 and 1/2
	else {
	  newFirm.setFinalDemandShare((2.0)/((double) firmNum));
        }
      } //start with equal nominal demand shares by tpye 
      else {
        newFirm.setFinalDemandShare(1/((double) firmNum));
      }
 
      // expected nominal demand for final goods
      newFirm.setExptdFinalDemand(2.32558); //equilibrium is 2.32558 if h = 1.22558 if mu=2.0, gamma=0.1, and delta=0.7 -- which makes Ystar/n=1.0
      newFirm.setFinalDemand(2.32558); //7/30/08
      //expected real demand for intermediate goods
      newFirm.setExptdIntermedDemand(1.0); //is irrelevant under new timing
      newFirm.setWage(wage); //5/13/16
      newFirm.setMarkup(markup); //5/26/16
      newFirm.setPrice((newFirm.getMarkup()/(1-newFirm.getMarkup()*newFirm.getInputFactor()*intPriceDiscount))*(newFirm.getWage()/newFirm.getTechEff())); //9/9/10
      newFirm.setPrice0(newFirm.getPrice()); //added 6/8/09 
      newFirm.setSalary(salary); // 5/12/16
      for (j=0; j<numHedonicElements; j++) { //numHedonicElements is declared in firm.h
        if (randNumUni() > 1 - 0.1) {
          newFirm.setHedonicQualityElementJ(j,rand() % 2);
        } // random integer in {0,1}
        else {
          newFirm.setHedonicQualityElementJ(j,0);
        }
      } 

      //set up integer constraints
      if(randNumUni() > 1 - probIntegerConstrained) {
        newFirm.setIntegerConstrained(true);
      }
      for (k=0; k<maxNumComplements; k++) {
        newFirm.setComplementaryProductK(k,0);
      } //7/3/09 //maxNumComplements defined in firm.h
      for (k=0; k<maxNumComplements; k++) {
	for (j=0; j<numHedonicElements; j++) {
          newFirm.setComplementaryHedonicElementKJ(k,j,0);
        }
      }
      //set up complementarities
      if(randNumUni() > 1 - probComplements) {
	newFirm.setHasComplements(true);
	newFirm.setNumComplements(rand() % maxNumComplements);
	for(k=0; k<maxNumComplements; k++) {
          newFirm.setComplementaryProductK(k,rand() % firmNum);
        }
	for(k=0; k<maxNumComplements; k++) {
	  for(j=0;j<numHedonicElements;j++) {
	    if(randNumUni() > 1 - 0.1) {
	      newFirm.setComplementaryHedonicElementKJ(k,j,rand() % 2);
	    } // random integer in {0,1}
	    else {
	      newFirm.setComplementaryHedonicElementKJ(k,j,0);
	    }
	  }
	}
      }
      if (randNumUni() > 1 - probInnovate) {
        newFirm.setRandDInvestment(true);
      } //3/22/10
      firms[i] = newFirm;
    }
}
    
//Update firm price by markup over full MC
void updatePrice(firm firms[], const int firmNum, /*double inputFactor,*/ std::ofstream &testoutput, bool anotation, double intPriceDiscount) { //added inputFactor
    int i;
    double lastprice, markup;
    firm current, prev;
    
    lastprice = (firms[firmNum-1]).getPrice(); 
    
    for (i = firmNum - 1 ; i >= 0;  i--) { 
      double price;
      current = firms[i];  
         
      if(i == 0) //5-8-07
        price = (1+0.0)*current.getMarkup()*(current.getWage()/current.getTechEff() + current.getInputFactor()*lastprice*intPriceDiscount);
      else {
        prev = firms[i - 1];
        price = (1+0.0)*current.getMarkup()*(current.getWage()/current.getTechEff() + current.getInputFactor()*prev.getPrice()*intPriceDiscount);
      }
         
      current.setPrice(price);
      if (anotation && i==5) { 
        testoutput << "price of firm " << i+1 << " = " << current.getPrice() << "\n"; 
      }
      firms[i] = current;
    }
}

//production and sale of intermediate goods
void doIntermediateGoods(firm firms[], double &wageBill, const int firmNum, double overHeadGReq, std::ofstream &testoutput, int round, bool anotation, double intPriceDiscount, bool intermediateGoods) {
  
  double maxIntGoodProduction;
  int i;
  firm current, next, prev; 
  
  if (intermediateGoods) { //10-6-10    //bypasses entire routine if intermediateGoods = false         
                       
    //Firms [i] place orders for  I goods from neighbor firms [i-1] 
    for (i = 0; i < firmNum; i++) {
      double exptdMaxFinalProduction;
      current = firms[i];
      if(i == 0)
         prev = firms[firmNum - 1];
      else
         prev = firms[i - 1];
                
      if (anotation && i==5) {
        testoutput << "firm " << i+ 1 << " has capital of " << current.getCapital() << "\n";        
        testoutput << "firm " << i+ 1 << " expects nominal final goods demand of $" << current.getExptdFinalDemand() << " and thus real final demand of " 
		   << current.getExptdFinalDemand()/current.getPrice() << " units "  << "\n"; 
      }      
      //determine expected constraint on production of C goods (by own capital and own purchase of I goods from neighbor) (plus net revs from sale to other neighbor should be included here
      exptdMaxFinalProduction = ((current.getCapital()-current.getIntGoodOrder()*prev.getPrice()*intPriceDiscount)/(current.getWage()/current.getTechEff())); //5-08-07
      if (exptdMaxFinalProduction < 0.0) 
	exptdMaxFinalProduction = 0;

      //place order for I goods from next neighbor based on expected production of C goods (which is last period's actual demand)
      if (exptdMaxFinalProduction >= (current.getExptdFinalDemand()/current.getPrice()))
        current.setIntGoodOrder(overHeadGReq + current.getInputFactor()*(current.getExptdFinalDemand()/current.getPrice())); 
      else {
        current.setIntGoodOrder(overHeadGReq + current.getInputFactor()*exptdMaxFinalProduction); 
        if (anotation && i==5) { 
          testoutput << "capacity less than exp final demand, so restrict int demand" << "\n"; 
	}
      }
      if (anotation && i==5) { 
	testoutput << "firm " << i + 1 << " orders " << current.getIntGoodOrder() << " units from firm " << i <<  "\n"; 
      }
      firms[i] = current;           
    }
    
    //Firms [i] produce I goods for sale to neighbor firms [i+1].
    for(i = 0; i < firmNum; i++) {        
      current = firms[i];
      next = firms[(i + 1)%firmNum];
                
      current.setExptdIntermedDemand(next.getIntGoodOrder()); //can set this here or in next loop for alt timing
         
      //capacity to produce I goods for neighbor is limited by own capital
      maxIntGoodProduction = current.getCapital()/(current.getWage()/current.getTechEff());
      if (maxIntGoodProduction < 0) 
	maxIntGoodProduction = 0;

      //produce I goods for next neighbor, decrement capital by cost, increment ag wage bill
      if (maxIntGoodProduction > current.getExptdIntermedDemand()) { 
        current.setIntGoodProduction(current.getExptdIntermedDemand()); 
        current.setCapital(current.getCapital()-((current.getWage()/current.getTechEff())* current.getExptdIntermedDemand()));
	current.setProfit(current.getProfit() - ((current.getWage()/current.getTechEff())* current.getExptdIntermedDemand())); //3/22/10
	wageBill +=(current.getWage()/current.getTechEff())* current.getExptdIntermedDemand();
      }
      else {
        current.setIntGoodProduction(maxIntGoodProduction);
        wageBill += (current.getWage()/current.getTechEff())*maxIntGoodProduction ;
        current.setCapital(current.getCapital()-(current.getWage()/current.getTechEff())*maxIntGoodProduction );
        current.setProfit(current.getProfit() - (current.getWage()/current.getTechEff())*maxIntGoodProduction ); //3/22/10
        if (anotation && i==5) { 
	  testoutput << "capacity less than exp intermediate demand, so limit intermed production" << "\n";
        }
      }
      if (anotation && i==5) {
        testoutput << "firm " << i + 1 << " produces " << current.getIntGoodProduction() << " units for firm " << i + 2 <<  "\n"; 
        testoutput << "cost of production is " << current.getIntGoodProduction()*(current.getWage()/current.getTechEff()) << "\n"; 
        testoutput << "firm " << i+ 1 << " now has capital of " << current.getCapital() << "\n";  
      }
      firms[i] = current;
    }
  
    //neighbors [i+1] buy I goods 
    for (i = 0; i < firmNum; i++) {
      /*double afdUnits;*/
      current = firms[i];
      next = firms[(i + 1)%firmNum];

      //sell I goods to neighbor, update both firms capitals by sale,
      //this sale (and thus neighbor's capacity to produce C goods) can be constrained by my supply of I goods
      if (current.getIntGoodProduction() < next.getIntGoodOrder()) {
        next.setMaxFinalProduction((current.getIntGoodProduction() - overHeadGReq)/next.getInputFactor()); //take out OH // input factor now tied to firm 9/9/10
        next.setCapital(next.getCapital()-(current.getIntGoodProduction()*current.getPrice()*intPriceDiscount)); //5-08-07
	next.setProfit(next.getProfit() - (current.getIntGoodProduction()*current.getPrice()*intPriceDiscount)); //3/22/10
	current.setCapital(current.getCapital() + (current.getIntGoodProduction()*current.getPrice()*intPriceDiscount)); //5-08-07
	current.setProfit(current.getProfit() + (current.getIntGoodProduction()*current.getPrice()*intPriceDiscount)); //3/22/10
	if(anotation) {
          if (i==4) { 
            testoutput << "firm " << i + 2 << " buys " << current.getIntGoodProduction() << " units from firm " << i+1 << "\n"; 
            testoutput << "cost is " << (current.getIntGoodProduction()*current.getPrice()*intPriceDiscount) << "\n"; 
          }     
          if (i==5) { 
            testoutput << "can't meet neighbor intermed demand" << "\n"; 
            testoutput << "firm " << i + 1 << " sells " << current.getIntGoodProduction() << " units to firm " << i + 2 << "\n"; 
            testoutput << "revenues are " << (current.getIntGoodProduction()*current.getPrice()*intPriceDiscount) << "\n"; 
          }
        }
      }
      else {
        next.setMaxFinalProduction((next.getIntGoodOrder() - overHeadGReq)/next.getInputFactor()); //take out OH 
        next.setCapital(next.getCapital() - (next.getIntGoodOrder()*current.getPrice()*intPriceDiscount)); //5-07-08
	next.setProfit(next.getProfit() - (next.getIntGoodOrder()*current.getPrice()*intPriceDiscount)); //3/22/10
	current.setCapital(current.getCapital() + (next.getIntGoodOrder()*current.getPrice()*intPriceDiscount)); //5-08-07
	current.setProfit(current.getProfit() + (next.getIntGoodOrder()*current.getPrice()*intPriceDiscount)); //3/22/10
        if (anotation) {
          if (i==4) { 
            testoutput << "firm " << i + 2 << " buys " << next.getIntGoodOrder() << " units from firm " << i+1 << "\n"; 
            testoutput << "cost is " << (next.getIntGoodOrder()*current.getPrice()*intPriceDiscount) << "\n"; 
          }
          if (i==5) { 
            testoutput << "firm " << i + 1 << " sells " << next.getIntGoodOrder() << " units to firm " << i + 2 << "\n"; 
            testoutput << "revenues are " << (next.getIntGoodOrder()*current.getPrice()*intPriceDiscount) << "\n"; 
          }
        }
      }
      firms[i] = current;
      firms[(i + 1)%firmNum] = next;
      if(anotation && i==5) 
	{testoutput << "firm " << i+ 1 << " now has capital of " << current.getCapital() << "\n"; 
      }     
    }                 //end for
  }                   //end if (bypass)
}                     //end doInermediateGoods

//production of final (consumption) goods
void doProduction(firm firms[], double &wageBill, const int firmNum, double overHeadGReq, std::ofstream &testoutput, int round, bool anotation, double intPriceDiscount, bool intermediateGoods) {
                
    double maxFinalGoodProduction, unitLCost ; //rename unitPrice unitLCost (unit production labor cost)
    int i;
    firm current, next, prev; 

    //firm i produce C goods (for public demand); no sale of C goods yet
    for(i = 0; i < firmNum; i++) {

      current = firms[i];
      //capital constraint for purchasing labor after intermediate goods production, purchases, and sales
      unitLCost = current.getWage()/current.getTechEff();
      maxFinalGoodProduction = current.getCapital()/unitLCost; //5/13/16
      if (maxFinalGoodProduction < 0) 
        maxFinalGoodProduction = 0.0; 

      //constraint by more restrictive of (1) financial capacity to hire labor and (2) intermediate goods on hand
      if(maxFinalGoodProduction > current.getMaxFinalProduction() && intermediateGoods) 
        maxFinalGoodProduction = current.getMaxFinalProduction();
      if (maxFinalGoodProduction < 0) 
	maxFinalGoodProduction = 0.0; 
                  
      //produce C goods, decrement capital by labor cost, increment ag wage bill  
      if(maxFinalGoodProduction < current.getExptdFinalDemand()/current.getPrice()) { 
        current.setFinalProduction(maxFinalGoodProduction);
        wageBill += maxFinalGoodProduction*unitLCost; //changed unitPrice to unitLCost 5/13/16
        current.setCapital(current.getCapital()-(maxFinalGoodProduction*unitLCost));
	current.setProfit(current.getProfit() - (maxFinalGoodProduction*unitLCost)); //3/22/10 
        if (anotation && i==3 && current.getExptdFinalDemand()/current.getPrice() - maxFinalGoodProduction > 0.00001) { 
	  testoutput << "produce less than exp final demand due to constraint" << "\n"; 
	} 
      }
      else {
        current.setFinalProduction(current.getExptdFinalDemand()/current.getPrice()); 
        current.setCapital(current.getCapital() - ((current.getExptdFinalDemand()/current.getPrice())*unitLCost)); //changed unitPrice to unitLCost 5/13/16
	current.setProfit(current.getProfit() - ((current.getExptdFinalDemand()/current.getPrice())*unitLCost)); //3/22/10
        wageBill += (current.getExptdFinalDemand()/current.getPrice())*unitLCost;  
      }
      if (anotation && i==5) {
        testoutput << "firm " << i+1 << " has final demand share " << current.getFinalDemandShare() << "\n"; //5/17/16
        testoutput << "firm " << i+1 << " has expected real demand of " << current.getExptdFinalDemand()/current.getPrice() << "\n"; //5/17/16"
        testoutput << "firm " << i + 1 << " produces " << current.getFinalProduction() << " units of final good" << "\n"; 
        testoutput << "additional labor cost is " << current.getFinalProduction()*unitLCost << "\n"; //changed unitPrice to unitLCost 5/13/16
      }  
      firms[i] = current;
    }
}
  
//calculate Utility of the rep consumer given firm shares and total income  //6/7/09
//changed wageBill to income 5/13/16 to allow different consumer classes
double calcUtility(firm firms[], double income, const int firmNum, double shares[], std::ofstream &testoutput, int round, bool anotation, bool cesHedonics, double cesHedonicsElastGoods, double cesHedonicsElastChars, double complementIntConstraintFactor) {

    //shares[] will be a vector of firm shares // this is new
	
    int ii, jj, kk; //changed all from i,j,k 8/2/09
    double *repConsumption; //8/2/09 //real consumption bundle for the representative consumer //need pointer to put in heap
    repConsumption = new double[firmNum]; //8/2/09
    double *repHedonicConsumption; //8/2/09 //corersponding hedonic bundle
    repHedonicConsumption = new double[numHedonicElements]; //8/2/09
    double repUtility = 0.0; //3/26/10 //changed to 0 for CES utility
    double rhoOneInverse = 0.0; //6/30/11
    double rhoTwoInverse = 0.0; //7/7/11
	
    rhoOneInverse = 1.0/cesHedonicsElastGoods; // 6-30-11 //moved here 7/7/11
    rhoTwoInverse = 1.0/cesHedonicsElastChars; // 7/7/11
			
    //we will establish the utility of the consumption bundle in three steps: 
    //note this is all ag consumption, not per worker
    //a. calc real consumption bundle based on inputted firm shares 
    for (ii=0;ii<firmNum;ii++) {
      repConsumption[ii] = (income * shares[ii])/firms[ii].getPrice(); //change wageBill to income 5/13/16
    } 

    //b. calc hedonic consumption bundle corresponding to that goods bundle //use ces aggregator if cesHedonics
    //c. calc utility for that hedonic bundle //use ces aggregator
    for (jj=0; jj<numHedonicElements; jj++) { //for each characteristic
      double repHedCons = 0.0;
      for (ii=0; ii<firmNum; ii++) { //for each firm
        if (cesHedonics) { //added 6/24/09
          //first do direct characteristics
          double zz = 0; //added 11/15/10 to combine base and compelements before applying ces 
          zz += ((double) firms[ii].getHedonicQualityElementJ(jj)) * repConsumption[ii]; //1/15/10
          //then add complementarities
          if (firms[ii].getHasComplements()) {
            for (kk=0;kk<firms[ii].getNumComplements();kk++) {
              double overlap = 0;
              if (repConsumption[ii] < repConsumption[firms[ii].getComplementaryProductK(kk)]) 
                overlap = floor((1/complementIntConstraintFactor)*repConsumption[ii])*complementIntConstraintFactor; //added fractional indteger constraint 11/15/10
              else
                overlap = floor((1/complementIntConstraintFactor)*repConsumption[firms[ii].getComplementaryProductK(kk)])*complementIntConstraintFactor; //added fractional indteger constraint 11/15/10
              zz += ((double) firms[ii].getComplementaryHedonicElementKJ(kk,jj)) * overlap; //11/15/10
            } //end for
          } //end if hasComplements 
	  repHedCons += (pow(((double) zz), cesHedonicsElastGoods));		
	  // raise this to 1/rho2 below -- done 7/7/11      
        } //end if cesHedonics
        else {
          //first do direct characteristics  
          repHedCons += (((double) firms[ii].getHedonicQualityElementJ(jj)) * repConsumption[ii]);
          //then add complementarities
               
	  if (firms[ii].getHasComplements()) { //7/6/09
            for(kk=0;kk<firms[ii].getNumComplements();kk++) {
              double overlap = 0;
              if (repConsumption[ii] < repConsumption[firms[ii].getComplementaryProductK(kk)])
                overlap = floor(repConsumption[ii]);
              else
                overlap = floor(repConsumption[firms[ii].getComplementaryProductK(kk)]);
              repHedCons += (((double) firms[ii].getComplementaryHedonicElementKJ(kk,jj)) * overlap);
            }  //end for
          }     //end if hasComplements
        }        //end else
      }            //end for each firm
      if (cesHedonics) {
	repHedonicConsumption[jj] = (1/((double) firmNum))*pow(repHedCons, rhoOneInverse);
      } //7/7/11 added if(cesHedonics)
      else {
	repHedonicConsumption[jj] = repHedCons;
      } //original code
      repUtility += pow(repHedonicConsumption[jj], cesHedonicsElastChars); //3/26/10 //make this CES to avoid prob of 0 values
    }  //end for each characteristic
    repUtility = (1/((double) numHedonicElements))*pow(repUtility, rhoTwoInverse); // 7-9-11 //can simply comment this line out to keep simpler version w/dim returns //7/7/11
	
    delete[] repConsumption; //8/2/09
    repConsumption = NULL; //8/2/09
    delete[] repHedonicConsumption; //8/2/09
    repHedonicConsumption = NULL; //8/2/09
	
    return repUtility;

}

//set final goods demand shares by pseudo consumer search  //7/30/08 
	//return changed to double 6/02/09
	//upgraded to use calcUtility procedure 6/7/09
double doConsumerSearch(firm firms[], double wageBill, double salaryBill, const int firmNum, std::ofstream &testoutput, int round, int numMutationPairs, int numSearchTrials, double mutScale, bool anotation, bool cesHedonics, double cesHedonicsElastGoods, double cesHedonicsElastChars, double complementIntConstraintFactor) {
	
    //shares are initialized to equal shares in main()
	
    double utility = 1.0; //moved here from top of search loop 6/02/09 so that can return utility from the procedure
    double *shares; //8/2/09
    shares = new double[firmNum]; //8/2/09
    double *altShares; //8/2/09
    altShares = new double[firmNum]; //8/2/09
    double altUtility = 1.0;
    double income; //5/13/16 to be able to combine wage and salary bills once split up
    int firm1 = 0;
    int firm2 = 0;
	
    income = wageBill + salaryBill; //5/13/16
	
    //do consumer search
    if (round >= 1) {		
      int i,j;
				
      //start with shares and altShares at current shares
      for (i=0;i<firmNum;i++) {
        shares[i] = firms[i].getFinalDemandShare();
	altShares[i] = firms[i].getFinalDemandShare();
      }

      //now use calcUtility routine to get baseline utility
      //changed wageBill to income 5/13/16 to allow separated wageBill and salaryBill -- see above
      utility = calcUtility(firms, income, firmNum, shares, testoutput, round, anotation, cesHedonics, cesHedonicsElastGoods, cesHedonicsElastChars, complementIntConstraintFactor);		
		
      //now update the shares via a mutation loop: test mutation for each of numMutationPairs pairs of goods
      int s=0;
      int ss=0;
      for (ss=0; ss<numSearchTrials; ss++) { //11/15/10 //add here to have multiple tests as before
        int negCount = 0; //11/15/10
        for (s=0;s<numMutationPairs;s++) { //11/15/10 this is now nested below searchTrials - allows multiple mutation pairs before each test
          double mutFactor=0.0;
          int sss = 0; //7-15-11
	  int integerCase = 0; 
	  int firm3 = 0;
	  //pick two firms to consider changing shares
	  firm1 = rand() % firmNum;
	  firm2 = rand() % firmNum;
	  //and try for firm2 with positive share 7/11/11 //7-15-11
	  while ((firm2 == firm1 || altShares[firm2] == 0) && sss < 20) { //sss < 20 in case one firm has the entire market -- added () re && argument 6-14-15
            firm2 = rand() % firmNum;
            sss++; //7-15-11
          }
 	  //check on integer constraints //added 6/25/09
          if (!firms[firm1].getIntegerConstrained() && !firms[firm2].getIntegerConstrained()) {
	    integerCase = 1;
	  }
          if (firms[firm1].getIntegerConstrained() && !firms[firm2].getIntegerConstrained()) {
	    integerCase = 2;
	  }
          if (!firms[firm1].getIntegerConstrained() && firms[firm2].getIntegerConstrained()) {
	    integerCase = 3;
	  } 
          if (firms[firm1].getIntegerConstrained() && firms[firm2].getIntegerConstrained()) {
	    integerCase = 4;
	  } 
                 
          //pick a mutation factor m betwen 0 and mutScale for shares (mutScale < 1.0, should be scaled with numFirms)
	  if (integerCase == 1) { //neither firm is integer constrained
	    mutFactor = randNumUni()*mutScale; //rand number between 0 and mutscale 
	    //bound mutFactor so that shares are bounded by 0 and 1. but allow shares to be pushed to bound.
	    if (mutFactor > 1- altShares[firm1]) {
	      mutFactor = 1 - altShares[firm1];
	    }
	    if (mutFactor > altShares[firm2]) {
	      mutFactor = altShares[firm2];
	    } 
	    altShares[firm1] = altShares[firm1] + mutFactor;
	    altShares[firm2] = altShares[firm2] - mutFactor;
          }    //end if integerCase 1  
            
          //6/26/09 note that we are making decision to change SHARES lumpy but not so for
          //        changing q_i when p's change. really should have this for all changes in q_i's.
          //        Further, we are allowing violation of 0 sheres bound in cases 2 and 3 -- and then
          //        throwing out any changes for which there is a violation
           
          //this version will add one unit to firm1 and iteratively searches for a set of other firms to deduct shares from
          //currently just picks 10 other firms at random -- currently violates 0 lower bound and integer constraints
          //though no shares actually go negative as neg altShares cause neg utility ??? not true 
          //should set up a proppper search
            
          if (integerCase == 2) { //firm1 is integer constrained // consider adding 1 unit to firm1 and taking shares from all others 
            mutFactor = firms[firm1].getPrice() / income; //changed wageBill to income 5/13/16
            if (mutFactor < 1- altShares[firm1]) {
              altShares[firm1] = altShares[firm1] + mutFactor;
              for (j=0;j<10;j++) { //**** need to keep altShares from going negative here ***  
                firm3 = rand() % firmNum;
                altShares[firm3] = altShares[firm3] - mutFactor/((double) 10);
                if (altShares[firm3] < 0) {
		  ++negCount;
		} 
              } //end for
            } //end if mutFactor
          } //end if integerCase 2
                  
          if (integerCase == 3) { //firm2 is integer constrained // consider taking 1 unit from firm2 and adding shares from 10 others
            mutFactor = firms[firm2].getPrice() / income; //changed wageBill to income 5/13/16  
            if (mutFactor < altShares[firm2]) {  
              altShares[firm2] = altShares[firm2] - mutFactor;
              for (j=0;j<10;j++) { //**** need to keep altShares from going negative here ***  
                firm3 = rand() % firmNum;
                altShares[firm3] = altShares[firm3] + mutFactor/((double) 10);
                if (altShares[firm3] < 0) {
		  ++negCount;
		}
              } //end for  
            } //end if mutFactor
          } //end if integerCase 3
            
          if(integerCase == 4) { //both firm1 and firm2 are integer constrained
            //bag this for now //add it later  ******
          }
        }  //end for loop over mutationPairs 
			 
        //if statement kills comparrison if any altShares went negative
        if (negCount == 0) {
	  //changed wageBill to income to allow sepaarate wage and salary bills 5/131/6          
          altUtility = calcUtility(firms, income, firmNum, altShares, testoutput, round, anotation, cesHedonics, cesHedonicsElastGoods, cesHedonicsElastChars,complementIntConstraintFactor);
        }
        else {
	  altUtility = utility;
	}
      						
        //switch to altShares (and update utility) if altUtility > utility
        if (altUtility >= utility) { //11/15/10 switched from > to <= for drift (see Miller)
          for (j=0; j<firmNum; j++) {
	    firms[j].setFinalDemandShare(altShares[j]);
	  } //6/26/09 do for all, not just firm1 and firm2
	  utility = altUtility;
        }  //end if altUtility > utility          
        for (j=0; j<firmNum; j++) {
	  altShares[j] = firms[j].getFinalDemandShare();
	} //6/26/09 do for all, not just firm1 and firm2	
      }	//end mutation loop over search trials	     	
    } //end search 
	
    delete[] shares; //8/2/09
    shares = NULL;  //8/2/09
    delete[] altShares; //8/2/09
    altShares = NULL; //8/2/09   
    return utility; //added 6/02/09
} //end doConsumerSearch

//consumer search for class 1 (production workers) when two consumer/firm classes
double doConsumerSearchClass1(firm firms[], double wageBill, const int firmNum, std::ofstream &testoutput, int round, int numMutationPairs, int numSearchTrials, double mutScale, bool anotation, bool cesHedonics, double cesHedonicsElastGoods, double cesHedonicsElastChars, double complementIntConstraintFactor) {
	
    //shares are initialized in inititalizeFirms
	
    double utility = 1.0; //moved here from top of search loop 6/02/09 so that can return utility from the procedure
    double altUtility = 1.0;
    double *shares; //8/2/09
    shares = new double[firmNum]; //8/2/09
    double *altShares; //8/2/09
    altShares = new double[firmNum]; //8/2/09
	
    double income; //5/13/16 to be able to combine wage and salary bills once split up
    int firm1 = 0;
    int firm2 = 0;
    income = wageBill; //5/13/16
	
    //do consumer search
    if (round >= 1) {				
      int i,j;
				
      //start with shares and altShares at current shares for type 1 firms and set to zero for type 2 firms
      for (i=0; i<firmNum; i++) {
	if (firms[i].getType() == 1) {	//5/13/16
	  shares[i] = firms[i].getFinalDemandShare(); //5/13/16 note -- doubled for twoClasses case in initializeFirms
	  altShares[i] = firms[i].getFinalDemandShare();
	}
	else {	//5/13/16
	  shares[i] = 0;  //5/13/16
	  altShares[i] = 0;
	}		
      }

      //now use calcUtility routine to get baseline utility -- note that shares are really half shares at this point. but ok as ces.
      //changed wageBill to income 5/13/16 to allow separated wageBill and salaryBill -- see above
      utility = calcUtility(firms, income, firmNum, shares, testoutput, round, anotation, cesHedonics, cesHedonicsElastGoods, cesHedonicsElastChars, complementIntConstraintFactor);		

      if(round==1){printf("initial utility for production workers in round 1 is %f \n", utility);}
		
      //now update the shares via a mutation loop: test mutation for each of numMutationPairs pairs of goods
      int s=0;
      int ss=0;

      for (ss=0; ss<numSearchTrials; ss++) { //11/15/10 //add here to have multiple tests as before
        int negCount = 0; //11/15/10
        for(s=0; s<numMutationPairs; s++) { //11/15/10 this is now nested below searchTrials - allows multiple mutation pairs before each test
          double mutFactor = 0.0;
          int sss = 0; //7-15-11
	  int integerCase = 0;
	  int firm3 = 0;
	  //pick two firms to consider changing shares
	  firm1 = rand() % firmNum;
	  firm2 = rand() % firmNum;
	  //make sure they are right type //5/13/16
	  while (firms[firm1].getType()!=1) {
	    firm1 = rand() % firmNum;
	  }
	  while (firms[firm1].getType()!=1) {
	    firm2 = rand() % firmNum;
	  }
	  //make sure they are different and of the right type and try for positive share 7/11  6/15 5/16
	  while ((firm2 == firm1 || altShares[firm2] == 0 || firms[firm2].getType()!=1) && sss < 50) { //sss < 50 in case one firm has the entire market -- added () re && argument 6-14-15
            firm2 = rand() % firmNum;
            sss++; //7-15-11
          }
          //check on integer constraints //added 6/25/09
          if (!firms[firm1].getIntegerConstrained() && !firms[firm2].getIntegerConstrained()) {
	    integerCase = 1;
	  }
          if (firms[firm1].getIntegerConstrained() && !firms[firm2].getIntegerConstrained()) {
	    integerCase = 2;
	  }
          if (!firms[firm1].getIntegerConstrained() && firms[firm2].getIntegerConstrained()) {
	    integerCase = 3;
	  } 
          if (firms[firm1].getIntegerConstrained() && firms[firm2].getIntegerConstrained()) {
	    integerCase = 4;
	  } 
               
	  //make sure that both firms are type 1 //5/21/16  
          if (firms[firm1].getType()==1 && firms[firm2].getType()==1) { //5/21/16
            //pick a mutation factor m betwen 0 and mutScale for shares (mutScale < 1.0, should be scaled with numFirms)
	    if (integerCase == 1) { //neither firm is integer constrained
	      mutFactor = randNumUni() * mutScale; //rand number between 0 and mutscale 
              //determine altShares by adding m to one firm's share and subtracting from another
	      //bound mutFactor so that shares are bounded by 0 and 1. but allow shares to be pushed to bound.
	      if (mutFactor > 1- altShares[firm1]) {
		mutFactor = 1 - altShares[firm1];
	      }
	      if (mutFactor > altShares[firm2]) {
		mutFactor = altShares[firm2];
	      } 
	      altShares[firm1] = altShares[firm1] + mutFactor;
	      altShares[firm2] = altShares[firm2] - mutFactor;
            } //end if integerCase 1  
            
            //6/26/09 note that we are making decision to change SHARES lumpy but not so for
            //        changing q_i when p's change. really should have this for all changes in q_i's.
            //        Further, we are allowing violation of 0 sheres bound in cases 2 and 3 -- and then
            //        throwing out any changes for which there is a violation
           
            //this version will add one unit to firm1 and iteratively searches for a set of other firms to deduct shares from
            //currently just picks 10 other firms at random -- currently violates 0 lower bound and integer constraints
            //though no shares actually go negative as neg altShares cause neg utility ??? not true 
            //should set up a proppper search
            
            if (integerCase == 2) { //firm1 is integer constrained // consider adding 1 unit to firm1 and taking shares from all others 
              mutFactor = firms[firm1].getPrice() / income; //changed wageBill to income 5/13/16
              if (mutFactor < 1- altShares[firm1]) { 
                altShares[firm1] = altShares[firm1] + mutFactor;
                for (j=0;j<10;j++) { //**** need to keep altShares from going negative here ***  
                  firm3 = rand() % firmNum;
                  altShares[firm3] = altShares[firm3] - mutFactor/((double) 10);
                  if (altShares[firm3] < 0) {
		    ++negCount;
		  } 
                } //end for
              } //end if mutFactor
            } //end if integerCase 2
                  
       	    if (integerCase == 3) { //firm2 is integer constrained // consider taking 1 unit from firm2 and adding shares from 10 others
              mutFactor = firms[firm2].getPrice() / income; //changed wageBill to income 5/13/16  
              if(mutFactor < altShares[firm2]) {  
                altShares[firm2] = altShares[firm2] - mutFactor;
                for (j=0; j<10; j++) { //**** need to keep altShares from going negative here ***  
                  firm3 = rand() % firmNum;
                  altShares[firm3] = altShares[firm3] + mutFactor/((double) 10);
                  if (altShares[firm3] < 0) {
		    ++negCount;
		  }
                } //end for 
              } //end if mutFactor
            } //end if integerCase 3
            
            if (integerCase == 4) { //both firm1 and firm2 are integer constrained
              //bag this for now //add it later  ******
            }
          } //end if both firms type1 //5/21/16
        }  //end for loop over mutationPairs 
			
        //determine altUtility
        //double altUtility = 1.0;
        //if statement kills comparrison if any altShares went negative
        if (negCount == 0) {
          //changed wageBill to income to allow sepaarate wage and salary bills 5/131/6          
          altUtility = calcUtility(firms, income, firmNum, altShares, testoutput, round, anotation, cesHedonics, cesHedonicsElastGoods, cesHedonicsElastChars,complementIntConstraintFactor);
	}
        else {
	  altUtility = utility;
	}
	//switch to altShares (and update utility) if altUtility > utility
        if (altUtility >= utility) { //11/15/10 switched from > to <= for drift (see Miller)
	  for (j=0; j<firmNum; j++) {
	    if (firms[j].getType() == 1) { //added if type 5/18/16
	      firms[j].setFinalDemandShare(altShares[j]); 
	    } 
	  } //6/26/09 do for all, not just firm1 and firm2
	  utility = altUtility;
        }  //end if altUtility > utility          
        //reset altShares -- matters for firms who are tested and not updated 
        //altShares[firm1] = firms[firm1].getFinalDemandShare();
        //altShares[firm2] = firms[firm2].getFinalDemandShare();
        for (j = 0; j<firmNum; j++) {
          if (firms[j].getType() == 1){ //keep type 2 altShares zerio 5/24/16
            altShares[j] = firms[j].getFinalDemandShare();
	  } //6/26/09 do for all, not just firm1 and firm2	
	}
      }	//end mutation loop over search trials
    } //end search 
	
    delete[] shares; //8/2/09
    shares = NULL;  //8/2/09
    delete[] altShares; //8/2/09
    altShares = NULL; //8/2/09   
    return utility;    //added 6/02/09
} //end doConsumerSearchClass1

//consumer search for class 2 (OH workers) when two consumer/firm classes
double doConsumerSearchClass2(firm firms[], double salaryBill, const int firmNum, std::ofstream &testoutput, int round, int numMutationPairs, int numSearchTrials, double mutScale, bool anotation, bool cesHedonics, double cesHedonicsElastGoods, double cesHedonicsElastChars, double complementIntConstraintFactor) {
	
    //shares are initialized to equal shares in initializeFirms
	
    double utility = 1.0; //moved here from top of search loop 6/02/09 so that can return utility from the procedure
    double *shares; //8/2/09
    shares = new double[firmNum]; //8/2/09
    double *altShares; //8/2/09
    altShares = new double[firmNum]; //8/2/09
    double altUtility = 1.0;
    double income; //5/13/16 to be able to combine wage and salary bills once split up
    int firm1 = 0;
    int firm2 = 0;

    income = salaryBill; //5/13/16
	
    //do consumer search
    if (round >= 1) {				
      int i,j;
				
      //start with shares and altShares at current shares
      for (i=0; i<firmNum; i++) {
	if (firms[i].getType() == 2) {	//5/13/16
	  shares[i] = firms[i].getFinalDemandShare(); //5/13/16 note -- doubled for twoClasses case in initializeFirms
	  altShares[i] = firms[i].getFinalDemandShare();
	}
	else {	//5/13/16
	  shares[i] = 0; //5/13/16
	  altShares[i] = 0;
	}		
      }

      //now use calcUtility routine to get baseline utility -- note that shares are really half shares at this point. but ok as ces.
      //changed wageBill to income 5/13/16 to allow separated wageBill and salaryBill -- see above
      utility = calcUtility(firms, income, firmNum, shares, testoutput, round, anotation, cesHedonics, cesHedonicsElastGoods, cesHedonicsElastChars, complementIntConstraintFactor);		
      if (round==1) {
	printf("initial utility for overhead workers in round 1 is %f \n", utility);
      }

      //now update the shares via a mutation loop: test mutation for each of numMutationPairs pairs of goods
      int s=0;
      int ss=0;
      for (ss = 0; ss < numSearchTrials; ss++) { //11/15/10 //add here to have multiple tests as before
        int negCount = 0; //11/15/10
        for (s = 0; s < numMutationPairs; s++) { //11/15/10 this is now nested below searchTrials - allows multiple mutation pairs before each test
          double mutFactor = 0.0;
          int sss = 0; //7-15-11
          int integerCase = 0;
          int firm3 = 0;
          //pick two firms to consider changing shares
          firm1 = rand() % firmNum;
          firm2 = rand() % firmNum;
          //make sure they are right type //5/13/16
          while (firms[firm1].getType() != 2) {
            firm1 = rand() % firmNum;
	  }
	  while (firms[firm1].getType()!=2) {
	    firm2 = rand() % firmNum;
	  }
	  //make sure they are different and the right type and try for firm2 with pos share 7/11 7/15 5/16
	  while ((firm2 == firm1 || altShares[firm2] == 0 || firms[firm2].getType()!=2) && sss < 50) { //sss < 50 in case one firm has the entire market -- added () re && argument 6-14-15
            firm2 = rand() % firmNum;
            sss++; //7-15-11
          }
 	  //check on integer constraints //added 6/25/09
          if (!firms[firm1].getIntegerConstrained() && !firms[firm2].getIntegerConstrained()) {
	    integerCase = 1;
	  }
          if (firms[firm1].getIntegerConstrained() && !firms[firm2].getIntegerConstrained()) {
	    integerCase = 2;
	  }
          if (!firms[firm1].getIntegerConstrained() && firms[firm2].getIntegerConstrained()) {
	    integerCase = 3;
	  } 
          if (firms[firm1].getIntegerConstrained() && firms[firm2].getIntegerConstrained()) {
	    integerCase = 4;
	  } 
               
	  //make sure firms are both type 2 //5/21/16  
          if (firms[firm1].getType() == 2 && firms[firm2].getType() == 2) { //5/21/16
            //pick a mutation factor m betwen 0 and mutScale for shares (mutScale < 1.0, should be scaled with numFirms)
	    if (integerCase == 1) { //neither firm is integer constrained
	      mutFactor = randNumUni()*mutScale; //rand number between 0 and mutscale 

              //determine altShares by adding m to one firm's share and subtracting from another
   	      //bound mutFactor so that shares are bounded by 0 and 1. but allow shares to be pushed to bound.
	      if (mutFactor > 1- altShares[firm1]) {
	        mutFactor = 1 - altShares[firm1];
	      }
	      if (mutFactor > altShares[firm2]) {
		mutFactor = altShares[firm2];
	      } 
	      altShares[firm1] = altShares[firm1] + mutFactor;
	      altShares[firm2] = altShares[firm2] - mutFactor;
            } //end if integerCase 1  
            
            //6/26/09 note that we are making decision to change SHARES lumpy but not so for
            //        changing q_i when p's change. really should have this for all changes in q_i's.
            //        Further, we are allowing violation of 0 sheres bound in cases 2 and 3 -- and then
            //        throwing out any changes for which there is a violation
           
            //this version will add one unit to firm1 and iteratively searches for a set of other firms to deduct shares from
            //currently just picks 10 other firms at random -- currently violates 0 lower bound and integer constraints
            //though no shares actually go negative as neg altShares cause neg utility ??? not true 
            //should set up a proppper search
            
            if (integerCase == 2) { //firm1 is integer constrained // consider adding 1 unit to firm1 and taking shares from all others
              mutFactor = firms[firm1].getPrice() / income; //changed wageBill to income 5/13/16
              if (mutFactor < 1- altShares[firm1]) {
                altShares[firm1] = altShares[firm1] + mutFactor;
                for (j=0; j<10; j++) { //**** need to keep altShares from going negative here ***  
                  firm3 = rand() % firmNum;
                  altShares[firm3] = altShares[firm3] - mutFactor/((double) 10);
                  if (altShares[firm3] < 0) {
		    ++negCount;
		  }
                }  //end for
              }       //end if mutFactor
            }            //end if integerCase 2
                  
            if (integerCase == 3) { //firm2 is integer constrained // consider taking 1 unit from firm2 and adding shares from 10 others
              mutFactor = firms[firm2].getPrice() / income; //changed wageBill to income 5/13/16  
              if (mutFactor < altShares[firm2]) {  
                altShares[firm2] = altShares[firm2] - mutFactor;
		for (j = 0; j < 10;j++) { //**** need to keep altShares from going negative here ***  
                  firm3 = rand() % firmNum;
                  altShares[firm3] = altShares[firm3] + mutFactor/((double) 10);
                  if (altShares[firm3] < 0) {
		    ++negCount;
		  }
                }    //end for 
              }         //end if mutFactor
            }              //end if integerCase 3
            
            if(integerCase == 4) { //both firm1 and firm2 are integer constrained
              //bag this for now //add it later  ******
            }
          } //end if both firms type 2 //5/21/16
        }  //end for loop over mutationPairs 
			
        //if statement kills comparrison if any altShares went negative
        if (negCount == 0) {
	  //changed wageBill to income to allow sepaarate wage and salary bills 5/131/6          
          altUtility = calcUtility(firms, income, firmNum, altShares, testoutput, round, anotation, cesHedonics, cesHedonicsElastGoods, cesHedonicsElastChars,complementIntConstraintFactor);
        }
        else {
	  altUtility = utility;
	}

          //switch to altShares (and update utility) if altUtility > utility
        if (altUtility >= utility) { //11/15/10 switched from > to <= for drift (see Miller)
          for (j = 0; j < firmNum; j++) {
	    if (firms[j].getType() == 2) { //added if getType 5/18/16
	      firms[j].setFinalDemandShare(altShares[j]);
	    }
	  } //6/26/09 do for all, not just firm1 and firm2
	  utility = altUtility;
		       
        } //end if altUtility > utility          
        for (j = 0; j < firmNum; j++) {
          if (firms[j].getType() == 2) { //keep type 1 altShares zerio 5/24/16
            altShares[j] = firms[j].getFinalDemandShare();
	  } //6/26/09 do for all, not just firm1 and firm2	
        }		
      } //end mutation loop over search trials	
    } //end search 
	
    delete[] shares; //8/2/09
    shares = NULL;  //8/2/09
    delete[] altShares; //8/2/09
    altShares = NULL; //8/2/09   
    return utility;    //added 6/02/09
} //end doConsumerSearchClass2

				
//sell C goods and do GDP accounting //8/2/08 take out reference in &wageBill //6/8/09 add nomGDP, realGDP, gdpDeflator
//separate wageBill and salaryBill 5/13/16
double doFinalSales(firm firms[], bool twoClasses, double wageBill, double salaryBill, const int firmNum, int firmNum1, int firmNum2, double rpGain, std::ofstream &testoutput, int round, bool anotation, bool altPriceIndex, double &nomGDP, double &realGDP, double &gdpDeflator, double &avgProfit, double &avgProfitType1, double &avgProfitType2) { //added avgProfit 4/6/16
    int i;
    double priceSum, share, totQuant, totNomSales, totRealSales, priceIndex, profitSum, profitSum1, profitSum2, income; //added income 5/13/16
    priceSum = 0.0;
    totQuant = 0.0;
    totNomSales = 0.0; //added 6/8/09
    totRealSales = 0.0; //added 6/8/09
    priceIndex = 0.0; //added 6/8/09
    profitSum = 0.0; //added 4/6/16
    income = 0.0; //added 5/13/16
    profitSum1 = 0.0; //5/25/16
    profitSum2 = 0.0;
        
    for (i = 0; i < firmNum; i++) {
      double revenue, unitCost;
      firm current;
      current = firms[i];
      revenue = 0.0;
         
      unitCost = current.getWage()/current.getTechEff();
      share = current.getFinalDemandShare(); //5/17/16
      //if there are two consumer/firm classes, then segment consumers and firms, else don't -- note that DemandShares start are init 1/n not 1/2n
      if (twoClasses) { //5/13/16
	if (current.getType()==1) {
	  income = wageBill;
	}
	else {
	  income = salaryBill;
	}
      }
      else { //5/13/16
        income = wageBill + salaryBill; //5/13/16
      }
		 
      current.setFinalDemand(share * income); //actual final demand //7/30/08 explicitly created actual final demand rather than using updated expectedFinalDemand
      //replaced wageBill with income 5/13/16
      if (anotation && i == 5) {
        testoutput << "total income for this consumer group is " << income << "\n"; //5/17/16
	testoutput << "firm " << i+1 << " has demand share " << current.getFinalDemandShare() << "\n"; //5/17/16
	testoutput << "firm " << i + 1 << " has expected nominal final demand of $" <<  current.getExptdFinalDemand() << " and thus real expected final demand of " 
		   <<  current.getExptdFinalDemand()/current.getPrice() << " units" << "\n";
	testoutput << "firm " << i + 1 << " has actual nominal final demand of $" <<  current.getFinalDemand() << " and thus real final demand of " <<  current.getFinalDemand()/current.getPrice() 
		   << " units" << "\n"; //5/17/16
      } 

      if (current.getFinalDemand()/current.getPrice() >= current.getFinalProduction()) { //7/30/08
        revenue += current.getFinalProduction()*current.getPrice();
        if (anotation && i == 5){
          testoutput << "firm " << i + 1 << " sells " << current.getFinalProduction() << " units of final good" << "\n"; 
        }
      }
      else {
        revenue += current.getFinalDemand(); //7/30/08
        if (anotation && i == 5) {
	  testoutput << "firm " << i + 1 << " sells " << current.getFinalDemand()/current.getPrice() << " units of final good" << "\n"; 
	} //7/30/08
      }
      if(anotation && i==5){
	testoutput << "revenues are " <<  revenue << "\n"; 
      }
         
      current.setCapital(current.getCapital() + revenue);
      current.setProfit(current.getProfit() + revenue); //3/22/10
      if (anotation && i==5) {
	testoutput << "firm " << i+ 1 << " now has capital of " << current.getCapital() << "\n"; 
      }  
    
      if (current.getFinalProduction() > current.getFinalDemand()/current.getPrice()) {  //7/30/08
        totQuant += current.getFinalDemand()/current.getPrice();    //unadjusted for quality //7/30/08
        totNomSales += current.getFinalDemand(); //added 6/8/09 for nom gdp
        totRealSales += current.getFinalDemand() * (current.getPrice0()/current.getPrice()); //added 6/8/09 for real gdp
      }
      else {
        totQuant += current.getFinalProduction();     //unadjusted for quality 
        totNomSales += current.getFinalProduction() * current.getPrice(); //added 6/8/09 for nom gdp
        totRealSales += current.getFinalProduction() * current.getPrice0(); //added 6/8/09 for real gdp
      }         
      priceIndex += (current.getPrice()/current.getPrice0()) * current.getFinalDemandShare(); 
      //added 6/8/09 could consider using actual sales or production // oops -- this is using nom share not real share

      current.setExptdFinalDemand(current.getFinalDemand());  //update expected to actual nominal demand //7/30/08         
      firms[i] = current;
      double rp = 0.0; //3/22/10
      rp = (1-rpGain) * firms[i].getRecentProfits() + rpGain * firms[i].getProfit(); //3/22/10
      firms[i].setRecentProfits(rp); //3/22/10
		 
      profitSum += firms[i].getProfit(); //added 4/6/16
      if (firms[i].getType() == 1) {
	profitSum1 += firms[i].getProfit();
      }
      if (firms[i].getType() == 2) {
	profitSum2 += firms[i].getProfit();
      }		 
    }

    if (altPriceIndex){
      //this is for alt measure using current nominal shares as weights in price index
      nomGDP = totNomSales; //added 6/8/09
      gdpDeflator = priceIndex; //added 6/8/09
      realGDP = nomGDP/gdpDeflator; //added 6/8/09
    }
    else {
      //this is for std real gdp and gdp deflator measures
      nomGDP = totNomSales; //added 6/8/09
      realGDP = totRealSales; //added 6/8/09
      gdpDeflator = nomGDP/realGDP; //added 6/8/09  
    }
    
    if (anotation){
      testoutput << "totNomSales are " <<  totNomSales << "\n" ; 
      testoutput << "totRealSales are " <<  totRealSales << "\n" ;
      testoutput << "gdnDeflator is " <<  gdpDeflator << "\n" ;
      testoutput << "totQuant is " <<  totQuant << "\n" ;
    }
    
    avgProfit = profitSum/((double) firmNum); //added 4/6/16
    avgProfitType1 = profitSum1/((double) firmNum1); //5/25/16
    avgProfitType2 = profitSum2/((double) firmNum2); //5/25/16
    
    return totQuant;
}

//Returns a random number between 0-1 according to a uniform dist
double randNumUni() {
    return (1 - (double) rand() / ((double) RAND_MAX + 1)); 
}

//Returns a random number according to N(0,1) dist
double randNumNorm() {
    double x, y, ans;
    x = randNumUni();
    y = randNumUni();
    ans = (pow(((-2.0) * log(x)) , 0.5) * cos((2.0 * 3.14159) * y));   //M_PI
    return ans;
}

//Alters labor efficiency of firms stochastically
double updateTechEff(firm firms[], const int firmNum, std::ofstream &testoutput, int &round, double &avTechEff, double &strtTech, bool anotation, double probEffShock) {
    int i;
    double totTech, minTechEff, strtTechEff; 
    totTech = minTechEff = 0.0;
    strtTechEff = 0;
              
    for(i = 0; i < firmNum; i++) {
      firm current;
      double change, temp;
      current = firms[i];
      change = 1.0;
      if (randNumUni() >= 1 - probEffShock && round >=3) { //default is 0.1
        change = 1.00 + (randNumNorm() * current.getEffStdDev()); //default
        if (change <= 0.1) { //keep techEff positive
          change = 0.1;
          testoutput << " change less than 0.1 " << "\n";
        }
        current.setTechEff((current.getTechEff()) * change);
	firms[i] = current;
      }
      if (round == 3 && i == 5 && anotation) { 
	printf( "A shock is  %f\n", change );
      }
      if (i == 5 && anotation) {
	testoutput << "A shock factor for firm " << i+1 << " is " << change << "\n"; 
      }
         
      temp = current.getTechEff();
      totTech += temp;
      if (minTechEff == 0)
        minTechEff = temp;
      else if (minTechEff > temp)
        minTechEff = temp; 
    }
     
    avTechEff = totTech/((double)firmNum);
    //strtTech set to avg for now for record -- should remove as all done in distFund
    strtTech = avTechEff;
    strtTechEff = strtTech;
    return strtTechEff;     
}

//Alters labor efficiency of firms stochastically for representative agent/common shock case
double updateTechEffCommonShock(firm firms[], const int firmNum, std::ofstream &testoutput, int &round, double &avTechEff, double &strtTech, bool anotation) {
    int i;
    double totTech,  minTechEff,  strtTechEff, change; //moved change here for common shock
    totTech = minTechEff = 0.0;
    strtTechEff = 0;
    change = 1.0;
    
    if(round == 5) {
      change = 1.00 + (randNumNorm()*(firms[0]).getEffStdDev()); 
      if (change <= 0.1)  
        change = 0.1;
    }        
    if (anotation) {
      testoutput << "common A shock factor is " << change << "\n"; 
    }
                   
    for (i = 0; i < firmNum; i++) {
      firm current;
      double temp;
      current = firms[i];
      current.setTechEff((current.getTechEff()) * change);
      firms[i] = current;
      temp = current.getTechEff();
      totTech += temp;
      if (minTechEff == 0)
        minTechEff = temp;
      else if (minTechEff > temp)
        minTechEff = temp; 
    }
     
    avTechEff = totTech/((double)firmNum);
    //strtTech set to avg for now for record -- should remove as all done in distFund
    strtTech = avTechEff;
    strtTechEff = strtTech;
    //strtTechEff = avTechEff - ((avTechEff-minTechEff)/2);
    return strtTechEff;     
}

//removed legacy updateQual and updateQualCommonShock 8/31/10

//exogenous product innovation -- update hedonic qualities of firms' products
void updateHedonicQual(firm firms[], const int firmNum, std::ofstream &testoutput, int round, double probInnovate, double probUpTick, double probDownTick, bool anotation) {
    int i,j,k;
            
    for(i = 0; i < firmNum; i++) {     
      int value; 
      if (randNumUni() >= 1 - probInnovate  && round >= 3) { 
        for (j = 0; j < numHedonicElements; j++) {
	  value = firms[i].getHedonicQualityElementJ(j);
	  if (value > 0 || randNumUni() > 1- 10.0 * probUpTick) {
	    if (randNumUni() >= 1 - probUpTick) {
	      firms[i].setHedonicQualityElementJ(j,value + 1);
	    }
          }
	  if (randNumUni() >= 1 - probDownTick && value > 0) {
	    firms[i].setHedonicQualityElementJ(j,value - 1);
	  } //end if
	} //end for
			
	//now do for complements //7/19/09 not working correctly //must have fixed aug 2, as now tests correctly 3/19/10
	for (k = 0; k < maxNumComplements; k++) {
	  for (j = 0; j < numHedonicElements; j++) {
	    value = firms[i].getComplementaryHedonicElementKJ(k,j);
	    if (value > 0 || randNumUni() > 1- 10.0 * probUpTick) {
	      if (randNumUni() >= 1 - probUpTick) {
	        firms[i].setComplementaryHedonicElementKJ(k,j,value + 1);
	      }
	    }
	    if (randNumUni() >= 1 - probDownTick && value > 0) {
	      firms[i].setComplementaryHedonicElementKJ(k,j,value - 1);
	    } //end if
	  } //end for
	} //end for		
      } //end if
    } //end for 
}

//endogenous product innovation -- update hedonic qualities of firms' products for firms that are elibible to innovate conditional on recent RandD investment status
void doEndogProductInnovation(firm firms[], const int firmNum, std::ofstream &testoutput, int round, double probInnovate, double probUpdateZeroHQElement, double probUpTick, double probDownTick, bool multiplicativeMutation, double incrementFactor, bool anotation) { //added incrementFactor 7-11-11
    int i,j,k;        
    for(i = 0; i < firmNum; i++) { 
      int value; 
      //new increment terms to allow additive or multiplicative mutation 6-9-11
      double inc1, inc2, inc3; //incrementFactor //increment rate is now declared in main
      inc1 = 0.0; inc2 = 0.0; inc3 = 0.0;
      int increment = 1;
	if (randNumUni() >= 1 - probInnovate && randNumUni() <= firms[i].getRecentRandDInvestment()  && round >= 3) { 
	  //8/31/10 tied (stochastically) to recent rather than (deterministically to) current R&D investment

	  //first do for base characteristics	
	  for (j = 0; j < numHedonicElements; j++) {
	    value = firms[i].getHedonicQualityElementJ(j);
	    // additive or multiple increment new 6-9-11
	    inc1 = incrementFactor * ((double) firms[i].getHedonicQualityElementJ(j)); //multiplicative increment
	    inc2 = floor(inc1); // still integer constrianed mutation
	    if (inc2 > 1.0) {
	      inc3 = inc2;
	    }
	    else {
	      inc3 = 1.0;
	    } // but now use increment of 1 if mult increment below 2, else use mult increment
	    if (multiplicativeMutation) {
	      increment = (int) inc3;
	    } // if multiplicativeMutation turned on, then use this
	    else {
	      increment = 1;
	    } // else use additive increment 
	    // end setting up increment 6-9-11
	    if (value > 0 || randNumUni() > 1- probUpdateZeroHQElement) { // weak preferential attachment - made probUpdateHQElemnt separate param 11-9-15
	      if (randNumUni() >= 1 - probUpTick) {
		firms[i].setHedonicQualityElementJ(j,value + increment); //6-9-11 
	      }
            }
	    if (randNumUni() >= 1 - probDownTick && value > 0) {
	      firms[i].setHedonicQualityElementJ(j,value - increment); //6-9-11
	    } //end if
	  } //end for
			
	//now do for complements //7/19/09 fixed aug 2, tests correctly 3/19/10
	for (k = 0; k < maxNumComplements; k++) {
	  for (j = 0;j < numHedonicElements; j++) {
	    value = firms[i].getComplementaryHedonicElementKJ(k,j);
	    // increment setup added 6-9-11
	    inc1 = incrementFactor * ((double) firms[i].getComplementaryHedonicElementKJ(k,j));
	    inc2 = floor(inc1);
	    if (inc2 > 1.0) {
	      inc3 = inc2;
	    }
	    else {
	      inc3 = 1.0;
	    }
	    if (multiplicativeMutation) {
	      increment = (int) inc3;
	    }
	    else {
	      increment = 1;
	    }
	    // end increment setup 6-9-11
	    if (value > 0 || randNumUni() > 1- 10.0 * probUpTick) {
	      if (randNumUni() >= 1 - probUpTick) {
	        firms[i].setComplementaryHedonicElementKJ(k,j,value + increment); //increment is additive or multiplicative 6-9-11
	      }
	    }
	    if (randNumUni() >= 1 - probDownTick && value > 0) {
	      firms[i].setComplementaryHedonicElementKJ(k,j,value - increment); //6-9-11
	    } //end if
	  } //end for
	} //end for	
      } //end if
    } //end for    
}

//Update R&d investment status -- firms choose whether or not to invest in R&D via social learning with individual sampling of other firms, plus random mutation plus conditional on working capital  
void updateRandDInvestmentStatus(firm firms[], const int firmNum, double intensityOfChoice, double probRandDInvestMutation, bool intermediateGoods, double overHeadLReq, double overHeadRandDReq, double salary, double intPriceDiscount, std::ofstream &testoutput, int round, bool anotation) { //3/21/10
    int i,j;
    int firm1 = 0;
    int firm2 = 0;
		
    //firm i searches for two other firms, one with more and one will less recent RandD investment than i. Switches own behavior depending on relative profits of these two.
    for(i = 0; i < firmNum; i++) {
      int maxCompInnov = 5; //max number of comparison firms to search over // should put this in main // 1
      firm1 = rand() % firmNum;
      j = 0;
      while (firms[firm1].getRecentRandDInvestment() <= firms[i].getRecentRandDInvestment() && j<maxCompInnov) {
        firm1 = rand() % firmNum;
	j++;
      }
      //if can't find such a firm, set firm1 to firm i
      if (j == maxCompInnov) {
	firm1 = i;
      }
		
      //search for a firm2 with less innovation than firm i
      firm2 = rand() % firmNum;
      j = 0;
      while (firms[firm2].getRecentRandDInvestment() >= firms[i].getRecentRandDInvestment() && j<maxCompInnov) {
	firm2 = rand() % firmNum;
	j++;
      }
      //if can't find such a firm, set firm2 to firm i
      if (j == maxCompInnov) {
	firm2 = i;
      }
		
      //imitate the firm with the higher recent profit -- if that recent profit is greater than own
      if (firms[i].getRecentProfits() < firms[firm1].getRecentProfits() || firms[i].getRecentProfits() < firms[firm2].getRecentProfits()) {
	if (firms[firm1].getRecentProfits() > firms[firm2].getRecentProfits()) {
	  firms[i].setRandDInvestment(firms[firm1].getRandDInvestment());
	}
	if (firms[firm2].getRecentProfits() > firms[firm1].getRecentProfits()) {
	  firms[i].setRandDInvestment(firms[firm2].getRandDInvestment());
	}	
      }	//end if
		
      //add some randomness //3/22/10
      bool bb = firms[i].getRandDInvestment();
      if (randNumUni() > 1 - probRandDInvestMutation) { //0.01
	firms[i].setRandDInvestment(!bb);
      } 
		
      //updating of recentRandDInvestment, count of firms investing, and investment cost accounting moved to separate routine 8/31/10 
		
      //turn off RandD invesetment if working capital too low
      double expectedCost, ep, mlc, igc;
      firm current, prev; //just account for input demand not supply
      current = firms[i];
      if (i == 0) 
	prev = firms[firmNum - 1];
      else 
	prev = firms[i - 1];
      ep = current.getExptdFinalDemand()/current.getPrice(); //expected production
      mlc = current.getWage()/current.getTechEff(); //marginal labor cost
      igc = current.getIntGoodOrder()*prev.getPrice()*intPriceDiscount; //expected intermediate good purchase 
      if (intermediateGoods) 
	expectedCost = ep * mlc + igc + overHeadLReq * salary + overHeadRandDReq * salary; //expected costs this period with RandD
      else 
	expectedCost = ep * mlc + overHeadLReq * salary + overHeadRandDReq * salary; //expected costs this period with RandD
      if (current.getCapital() < 2 * expectedCost)
	current.setRandDInvestment(false);
      firms[i] = current;	
    } //end for
} //end updateEndogenousInnovation

//Update R&D investment status -- firms choose whether or not to invest in R&D via social learning by quasi global discrete choice, plus random mutation plus conditional on working capital //7/10
//firms who are using strategy with lower average profits switch with probability from logistic discrete choics
void updateRandDInvestmentStatusDC(firm firms[], const int firmNum, double intensityOfChoice, double weightRandDDC, double probRandDInvestMutation, bool intermediateGoods, double overHeadLReq, double overHeadRandDReq, double salary, double intPriceDiscount, double &medRecentRandDInvestment, double &avgRecentProfitAboveMedRandD, double &avgRecentProfitBelowMedRandD, double &avgRecentProfitType1, double &avgRecentProfitType2, std::ofstream &testoutput, int round, bool anotation) { //3/21/10

    int i,j;
    double discChoiceProb = 0.0; //this is the logistic prob for above med 
    double discChoiceProb1 = 0.0; //disc choice prob based on general profiability of R&D //5/28/16
    double discChoiceProb2 = 0.0; //disc choice prob based on rel profitability of own market //5/28/16
    int numAboveMed = 0;
    int numBelowMed = 0;
    int numType1 = 0; //5/25/16
    int numType2 = 0;
	
    double medRecentRandDInvest = 0.0; 	 //changed 4/6/16
    double avgRecentProfitAboveMed = 0.0; 
    double avgRecentProfitBelowMed = 0.0; 
    double avgRecentProfit1 = 0.0; //5/25/16
    double avgRecentProfit2 = 0.0;
    medRecentRandDInvestment = 0.0; //changed 4/6/16
    avgRecentProfitAboveMedRandD = 0.0; 
    avgRecentProfitBelowMedRandD = 0.0; 
    avgRecentProfitType1 = 0.0; //5/25/16
    avgRecentProfitType2 = 0.0;
	
    //calc median (or average) recentInnoation (for now using average **)
    for (i = 0; i<firmNum; i++) {
      medRecentRandDInvest += firms[i].getRecentRandDInvestment();
    }
    medRecentRandDInvestment = medRecentRandDInvest/((double) firmNum); //changed 4/6/16
	
    //calc average recentProfits for firms above medRecentRandDInvestment
    for (i = 0; i<firmNum; i++) {
      if (firms[i].getRecentRandDInvestment() > medRecentRandDInvestment) {
	avgRecentProfitAboveMed += firms[i].getRecentProfits();
	numAboveMed++;
      }
    }
    if (numAboveMed > 0) {
      avgRecentProfitAboveMedRandD = avgRecentProfitAboveMed/((double) numAboveMed);
    } //changed 4/6/16
	
    //calc average recentProfits for firms below medRecentRandDInvestment
    for (i = 0; i<firmNum; i++) {
      if (firms[i].getRecentRandDInvestment() < medRecentRandDInvestment) {
        avgRecentProfitBelowMed += firms[i].getRecentProfits();
	numBelowMed++;
      }
    }
    if (numBelowMed > 0) {
      avgRecentProfitBelowMedRandD = avgRecentProfitBelowMed/((double) numBelowMed);
    } //changed 4/6/16
	
    //calc average recentProfits for firms of types 1 and 2 //5/25/16
    for (i = 0; i<firmNum; i++) {
      if (firms[i].getType() == 1) {
	avgRecentProfit1 += firms[i].getRecentProfits();
	numType1++;
      }
      if (firms[i].getType() == 2) {
	avgRecentProfit2 += firms[i].getRecentProfits();
	numType2++;
      }
    }
    if (numType1 > 0) {
      avgRecentProfitType1 = avgRecentProfit1/((double) numType1);
    }
    if (numType2 > 0) {
      avgRecentProfitType2 = avgRecentProfit2/((double) numType2);
    }
	
    //calc discChoiceProbInnovate here -- disc choice prob with above med recent profit on top //changed discChoiceProb to discChoiceProb1 5/28/16	
    if (numAboveMed > 0 && numBelowMed > 0) {
      discChoiceProb1 = exp(intensityOfChoice * avgRecentProfitAboveMedRandD) / ( exp(intensityOfChoice * avgRecentProfitAboveMedRandD) + exp(intensityOfChoice * avgRecentProfitBelowMedRandD)); 
    }
    
    //calc discChoiceProbInnovate here -- disc choice prob with type 1 firms on top //5/28/16 // might want if n1 n2 > 0 ***
    discChoiceProb2 = exp(intensityOfChoice * avgRecentProfitType1) / ( exp(intensityOfChoice * avgRecentProfitType1) + exp(intensityOfChoice * avgRecentProfitType2));
                      
    //add probability that we focus on general R&D or own market profit in DC problem //5/28/16
    //double weightRandDDC = 0.5; //5/28/16 will move to input.txt later ***
    if (randNumUni() > 1 - weightRandDDC) { //5/28/16 condition on relative overall R&D profitability
      //switch those on lower profit end with probability discChoiceProbInnovate or 1-discChoicProbInnovate
      if (numAboveMed > 0 && numBelowMed > 0) {
       	if (avgRecentProfitAboveMedRandD > avgRecentProfitBelowMedRandD) { //changed 4/6/16
	  for (i = 0; i<firmNum; i++) {
            bool bb = firms[i].getRandDInvestment();
            if (firms[i].getRecentRandDInvestment() < medRecentRandDInvestment) {
              if (randNumUni() > (2 - 2 * discChoiceProb1)) {
                firms[i].setRandDInvestment(true);
              }
            }
          } 
        }	
        if (avgRecentProfitAboveMedRandD < avgRecentProfitBelowMedRandD) { //changed 4/6/16
	  for (i = 0; i < firmNum; i++) {
            bool bb = firms[i].getRandDInvestment();
            if (firms[i].getRecentRandDInvestment() > medRecentRandDInvestment) {
              if (randNumUni() > (2 * discChoiceProb1)) {
                firms[i].setRandDInvestment(false);
              }
            }
          } 
        } 
      }
    }
    else {  //else condition on own market relative profit  //5/28/16
      if (avgRecentProfitType1 > avgRecentProfitType2) { //5/28/16
        for (i = 0; i<firmNum; i++) {
	  if (randNumUni() > (2 - 2 * discChoiceProb1)){
            if (firms[i].getType()==1 && firms[i].getRandDInvestment()==false) { //only switch if not investing
              firms[i].setRandDInvestment(true);
            }
            if (firms[i].getType()==2 && firms[i].getRandDInvestment()==true) { //only switch if not investing
	      firms[i].setRandDInvestment(false);
            }
          }
        } 
      }	
      if (avgRecentProfitType1 < avgRecentProfitType2) { //changed 4/6/16
        for (i = 0; i<firmNum; i++) {
	  if (randNumUni() > (2 * discChoiceProb1)) {
            if (firms[i].getType() == 1 && firms[i].getRandDInvestment() == true) {
              firms[i].setRandDInvestment(false);
            }
            if (firms[i].getType() == 2 && firms[i].getRandDInvestment() == false) {
              if (randNumUni() > (2 * discChoiceProb1)) {
                firms[i].setRandDInvestment(true);
              }
            } 
          } 
	}
      }
    } //end else -- condition on own market relative profit
	
    if((round % 1000) == 0){
      printf("average recent profit below and above median recnt R&D activity in round %d are %f and %f \n", round, avgRecentProfitBelowMedRandD, avgRecentProfitAboveMedRandD );
    } //4/6/16 testing 

	
	
    //add some random choice turn off RandD invesetment if working capital too low
    for (i = 0; i<firmNum; i++) {
		
      //add some randomness //3/22/10
      bool bb = firms[i].getRandDInvestment();
      if (randNumUni() > 1 - probRandDInvestMutation) { 
        firms[i].setRandDInvestment(!bb);
      } 
		
      //turn off RandD invesetment if working capital too low	//10-6-10
      double expectedCost, ep, mlc, igc;
      firm current, prev; //just account for input demand not supply
      current = firms[i];
      if (i == 0)
	 prev = firms[firmNum - 1];
      else
	 prev = firms[i - 1];
      ep = current.getExptdFinalDemand()/current.getPrice(); //expected production
      mlc = current.getWage()/current.getTechEff();          //marginal labor cost
      igc = current.getIntGoodOrder()*prev.getPrice()*intPriceDiscount;  //expected intermediate good purchase 
      if (intermediateGoods) 
	expectedCost = ep * mlc + igc + overHeadLReq * salary + overHeadRandDReq * salary; //expected costs this period with RandD
      else 
	expectedCost = ep * mlc + overHeadLReq * salary + overHeadRandDReq * salary;                  //expected costs this period with RandD
      if (current.getCapital() < 2*expectedCost)
	current.setRandDInvestment(false);
      firms[i] = current;    		
      } //end for	
    //  updating of recentRandDInvestment, count of firms investing, and investment cost accounting moved to separate routine 8/31/10 
} //end updateEndogenousInnovation

//increment wageBill and firm profits for fixed overhead labor and current RandD overhead labor and increment recentRandDInvestment //10-7-10
void doOHandRandDAccounting(firm firms[], const int firmNum, double riGain, double overHeadLReq, double overHeadRandDReq, double salary, double &wageBill, double &salaryBill, int &numRandD, int &numRandD1, int &numRandD2, std::ofstream &testoutput, int round, bool anotation) { //3/21/10  
    int i;
    for (i = 0; i<firmNum; i++) {
      //all firms incur fixed overHead labor cost // 10/4/10
      firms[i].setProfit(firms[i].getProfit() - overHeadLReq * salary); //10-6-10
      firms[i].setCapital(firms[i].getCapital() - overHeadLReq * salary); //6-1-11
      salaryBill += overHeadLReq * salary; //5/13/16 put OH salaries in salaryBill rather than wageBill
            
      //firms with innovation status true incur labor cost and update count of firms innovating
      if (firms[i].getRandDInvestment()) { //3/23/10
        firms[i].setProfit(firms[i].getProfit() - overHeadRandDReq * salary); //10-6-10
        firms[i].setCapital(firms[i].getCapital() - overHeadRandDReq * salary); //6-1-11
        salaryBill += overHeadRandDReq * salary; //5/13/16 put OH salaries in salaryBill rather than wageBill
        numRandD++; //7/9/10
        if (firms[i].getType() == 1){
	  numRandD1++;
        } //5/29/16 market 1
        if (firms[i].getType() == 2){
	  numRandD2++;
        } //5/29/16 market 2
        if (anotation && i==5) {
	  testoutput << "firm " << i+1 << " invests in R&D " << "\n";
        } //5/17/16
      }
           
      //update recentRandDInvestment with constant gain riGain
      double ri = (1 - riGain) * firms[i].getRecentRandDInvestment() + riGain * firms[i].getRandDInvestment();
      firms[i].setRecentRandDInvestment(ri); 		    
    } //end for
} //end doRandDAccounting

//restart when capital becomes negative and also randomly -- for representative agent/common shocks case // need to redo this for hedonics ****
int doRestartsCommonShock(firm firms[], const int firmNum, double strtTech, double wageBill, double overHeadGReq, double intPriceDiscount, bool imitation, double probRandomRestart, double wage, double salary, double markup) {
    
    firm current, prev, comparison; 
    int i, j, rndrestrt; 
    rndrestrt = 0;
    
    return rndrestrt;
}

//restart firms when capital becomes negative and also randomly -- here for hedonic case //8/12/08  -- //dramatically simplified 6-27-11
int doRestartsHedonics(firm firms[], const int firmNum, double strtTech, double wageBill, double salaryBill, double overHeadGReq, double intPriceDiscount, bool imitation, double probRandomRestart, double probIntegerConstrained, double probComplements, double wage, double salary, double markup) {
    
    firm current, prev, comparison; 
    int i, j, k, rndrestrt; 
    rndrestrt = 0;
    
    for (i = 0; i < firmNum; i++) {
      double price, share;
      current = firms[i];
      if (i == 0)
        prev = firms[firmNum - 1];
      else
        prev = firms[i - 1];
         
      j = rand() % firmNum;   
      comparison = firms[j]; 
                             
      if (current.getCapital() < 0.0  || randNumUni() < probRandomRestart) { //with random failure 0.0005
        firm newFirm;
        if (i==5) {
	  printf("firm 6 has restarted \n" );
	}
    
        if (imitation  && comparison.getTechEff() >= 0.5) { //combined 7/29/09
	  newFirm = comparison;
	} //this is new 6-27-11 //cut out all of the specific copying		 
        else {
          newFirm.setTechEff(1.0); //1.0
          for (j=0; j<numHedonicElements; j++) {
            if (randNumUni() > 1 - 0.1) {
	      newFirm.setHedonicQualityElementJ(j,rand() % 2);
            }
            else {
	      newFirm.setHedonicQualityElementJ(j,0);
	    }
	  }
	  if (randNumUni() > 1 - probIntegerConstrained) {
	    newFirm.setIntegerConstrained(true);
          }
	  for (k = 0; k < maxNumComplements; k++) {
	    newFirm.setComplementaryProductK(k,0);
	  } //7/3/09
	  for(k = 0; k < maxNumComplements; k++) {
	    for(j = 0; j < numHedonicElements; j++) {
	      newFirm.setComplementaryHedonicElementKJ(k,j,0);
            }
	  }
	  if (randNumUni() > 1 - probComplements) {
	    newFirm.setHasComplements(true);
	    newFirm.setNumComplements(rand() % maxNumComplements);
	    for (k = 0; k < maxNumComplements; k++) {
	      newFirm.setComplementaryProductK(k,rand() % firmNum);
            }
	    for (k = 0; k < maxNumComplements; k++) {
	      for(j = 0; j < numHedonicElements; j++) {
	        if (randNumUni() > 1 - 0.1) {
	          newFirm.setComplementaryHedonicElementKJ(k,j,rand() % 2);
                } // random integer in {0,1}
	        else {
	          newFirm.setComplementaryHedonicElementKJ(k,j,0);
                }
	      }
	    } 
	  }				   
        }
                
        newFirm.setType(current.getType()); //5/13/16
        newFirm.setMarkup(markup); //5/26/16
        price = newFirm.getMarkup() * (newFirm.getWage()/newFirm.getTechEff() + newFirm.getInputFactor() * prev.getPrice() * intPriceDiscount); //9/9/10
        newFirm.setPrice(price);
        newFirm.setPrice0(current.getPrice0()); //added   6/9/09
        newFirm.setWage(wage); //5/26/16
        newFirm.setSalary(salary); //5/12/16
        share = current.getFinalDemandShare(); //8/12/08
        newFirm.setFinalDemandShare(share); //8/12/08 
        newFirm.setExptdFinalDemand(share * (wageBill + salaryBill)); // added salaryBill 5/13/16
        prev.setExptdIntermedDemand(overHeadGReq + newFirm.getInputFactor()*(current.getExptdFinalDemand()/newFirm.getPrice())); //changed wageBill to current.getExptdFinalDemand() 5/13/16
        newFirm.setExptdIntermedDemand(current.getExptdIntermedDemand());
        newFirm.setIntegerConstrained(current.getIntegerConstrained());
        newFirm.setCapital(100.0); //7-12-11
        newFirm.setRandDInvestment(current.getRandDInvestment()); //7-15-11
        newFirm.setRecentRandDInvestment(current.getRecentRandDInvestment());  //7-15-11
        newFirm.setProfit(current.getProfit()); //7-15-11
        newFirm.setRecentProfits(current.getRecentProfits());  //7-15-11
				
        current = newFirm;
        rndrestrt++;
      }
         
      if (i == 0)
        firms[firmNum - 1] = prev;
      else
        firms[i - 1] = prev; 
      firms[i] = current;
    }
    return rndrestrt;
}

//restart firms when capital becomes negative and also randomly -- here for hedonic case //8/12/08  -- //dramatically simplified 6-27-11
// and now for endog n1 n2 at entry by comparison of avg recent profit in each market
int doRestartsHedonicsEndogN1N2(firm firms[], const int firmNum, int &firmNum1, int &firmNum2, double strtTech, double wageBill, double salaryBill, double overHeadGReq, double intPriceDiscount, bool imitation, double probRandomRestart, double probIntegerConstrained, double probComplements, double wage, double salary, double markup, double avgRecentProfitType1, double avgRecentProfitType2, double probSwitchMarkets, double probRandomSwitchMarkets, int ownMarketBias) {

    firm current, prev, comparison; 
    int i, j, k, ii, jj, kk,  rndrestrt; 
    rndrestrt = 0;
    
    for (i = 0; i < firmNum; i++) {
      double price, share;
      current = firms[i];
      if (i == 0)
        prev = firms[firmNum - 1];
      else
        prev = firms[i - 1];
         
      if (current.getCapital() < 0.0  || randNumUni() < probRandomRestart) { //with random failure 0.0005
        firm newFirm;

        if(i==5) {
	  printf("firm 6 has restarted \n" );
	  printf("avgRecentProfitType1 and avgRecentProfitType2 are %f and %f \n ", avgRecentProfitType1, avgRecentProfitType2);
	}
				
	// decide whether to switch markets // added 6/8/16
	bool switch1to2 = false; 
	bool switch2to1 = false; 

	if ((current.getType()==1 && avgRecentProfitType1 < avgRecentProfitType2 && randNumUni() < probSwitchMarkets) || randNumUni() < probRandomSwitchMarkets) {
	  current.setType(2);
	  switch1to2 = true; 
	}	
					
	if ((current.getType()==2 && avgRecentProfitType1 > avgRecentProfitType2 && randNumUni() < probSwitchMarkets) || randNumUni() < probRandomSwitchMarkets) {
	  current.setType(1);
	  switch2to1 = true;
	}
				
	//recount n1 and n2
	ii = kk = jj = 0;
	for (ii = 0; ii < firmNum; ii++) {
	  if (firms[ii].getType() == 1) 
	    jj++;
          else
	     kk++;
	}
	firmNum1 = jj;
	firmNum2 = kk;
								
	//pick a potential firm to immitate //updated 6/8/16
	//bias this toward own market by adjustable factor
	int sss = 0; 
	j = rand() % firmNum;
	while (firms[j].getType() != current.getType() && sss<ownMarketBias) {
	  j = rand() % firmNum;   //moved here 6/8/16
	  sss++;
	}

   	comparison = firms[j]; 
    
        if (imitation  && comparison.getTechEff() >= 0.5) { //combined 7/29/09
	  newFirm = comparison;
	} //this is new 6-27-11 //cut out all of the specific copying				 
        else {
          newFirm.setTechEff(1.0); //1.0
	  for (j = 0; j < numHedonicElements; j++) {
            if (randNumUni() > 1 - 0.1) {
	      newFirm.setHedonicQualityElementJ(j,rand() % 2);
	    }
            else {
	      newFirm.setHedonicQualityElementJ(j,0);
	    }
	  }
	  if (randNumUni() > 1 - probIntegerConstrained) {
	    newFirm.setIntegerConstrained(true);
	  }
	  for (k = 0; k < maxNumComplements; k++) {
	    newFirm.setComplementaryProductK(k,0);
	  } //7/3/09
	  for(k = 0; k < maxNumComplements; k++) {
	    for(j = 0; j < numHedonicElements; j++) {
	      newFirm.setComplementaryHedonicElementKJ(k,j,0);
	    }
	  }
	  if (randNumUni() > 1 - probComplements) {
	    newFirm.setHasComplements(true);
	    newFirm.setNumComplements(rand() % maxNumComplements);
	    for (k = 0; k < maxNumComplements; k++) {
	      newFirm.setComplementaryProductK(k,rand() % firmNum);
	    }
	    for(k = 0; k < maxNumComplements; k++) {
	      for(j = 0; j < numHedonicElements; j++) {
		if (randNumUni() > 1 - 0.1) {
		  newFirm.setComplementaryHedonicElementKJ(k,j,rand() % 2);
		} // random integer in {0,1}
		else {
		  newFirm.setComplementaryHedonicElementKJ(k,j,0);
		}
	      }
	    }
	  }				   
        }
                
        newFirm.setType(current.getType()); //5/13/16
        newFirm.setMarkup(markup); //5/26/16
        price = newFirm.getMarkup()*(newFirm.getWage()/newFirm.getTechEff() + newFirm.getInputFactor()*prev.getPrice()*intPriceDiscount); //9/9/10
        newFirm.setPrice(price);
        newFirm.setPrice0(current.getPrice0()); //added   6/9/09
        newFirm.setWage(wage); //5/26/16
        newFirm.setSalary(salary); //5/12/16
			   	
	if (!switch1to2 && !switch2to1) {
	  share = current.getFinalDemandShare(); //8/12/08	
	}
	if (switch1to2) {	
	  share = 1/(double)firmNum2;
	}	
	if (switch2to1) {
	  share = 1/(double)firmNum1; 
	}
				
	newFirm.setFinalDemandShare(share); //8/12/08				
				
        newFirm.setExptdFinalDemand(share * (wageBill + salaryBill)); // added salaryBill 5/13/16
        prev.setExptdIntermedDemand(overHeadGReq + newFirm.getInputFactor()*(current.getExptdFinalDemand()/newFirm.getPrice())); //changed wageBill to current.getExptdFinalDemand() 5/13/16
	newFirm.setExptdIntermedDemand(current.getExptdIntermedDemand());
	newFirm.setIntegerConstrained(current.getIntegerConstrained());
	newFirm.setCapital(100.0); //7-12-11
	newFirm.setRandDInvestment(current.getRandDInvestment()); //7-15-11
	newFirm.setRecentRandDInvestment(current.getRecentRandDInvestment());  //7-15-11
	newFirm.setProfit(current.getProfit()); //7-15-11
        newFirm.setRecentProfits(current.getRecentProfits());  //7-15-11
				
	current = newFirm;
        rndrestrt++;

     	if (i == 0)
          firms[firmNum - 1] = prev;
        else
          firms[i - 1] = prev;
    
        firms[i] = current; 
      }	 
    } 
    
    //renormalize shares to sum to one in each sector //6/8/16
    //should do intermed demand, etc. as well 
    double sumShares1=0.0;
    double sumShares2=0.0;
    for (i = 0; i < firmNum; i++){
      if (firms[i].getType() == 1)
	sumShares1 = sumShares1 + firms[i].getFinalDemandShare();
      else 
	sumShares2 = sumShares2 + firms[i].getFinalDemandShare();
    }
    for(i = 0; i < firmNum; i++){
      if (firms[i].getType() == 1) {
	firms[i].setFinalDemandShare(firms[i].getFinalDemandShare()/sumShares1);
	firms[i].setExptdFinalDemand(firms[i].getFinalDemandShare() * wageBill);	
      }
      else {
	firms[i].setFinalDemandShare(firms[i].getFinalDemandShare()/sumShares2);
	firms[i].setExptdFinalDemand(firms[i].getFinalDemandShare() * salaryBill);
      }
    }
    return rndrestrt;
}


