/*
 *  Marginal.cpp 
 *
 *  Copyright (C) 2007, 2008 Lydia Steiner and Thomas Kahle
 *  <{steiner,kahle}@mis.mpg.de>
 *
 * ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
 *
 *  This program is free software; you can redistribute it and/or modify
 *  it under the terms of the GNU General Public License as published by
 *  the Free Software Foundation; either version 2 of the License, or (at
 *  your option) any later version.
 *
 *  This program is distributed in the hope that it will be useful, but
 *  WITHOUT ANY WARRANTY; without even the implied warranty of
 *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
 *  General Public License for more details.
 *
 *  You should have received a copy of the GNU General Public License along
 *  with this program; if not, write to the Free Software Foundation, Inc.,
 *  59 Temple Place, Suite 330, Boston, MA 02111-1307 USA.
 *
 * ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
 */
#include <cstdlib>
#include <stdlib.h>
#include <math.h>
#include <iostream>
#include <fstream>
#include <string>
#include <vector>

#include "Marginal.h"
#include "functions.h"

using namespace std;

/*
 * PARAMS:
 * setA - set A which should be used to create this marginal
 * omegaA - "omega_A" which should be used to create this marginal
 * pEmp - samples ditribution
 * N - length of setA
 * alphabetSize - size of the alphabets 
 */
Marginal::Marginal(bool setA[], unsigned int omegaA, vector<double>* pEmp,unsigned int N,unsigned int alphabetSize){
          //initialize variables
          alpha = 0.0;
          belongToMe.push_back(0);
          //vector needed for calculate elements of belongToMe
          vector<unsigned int> help;
          //for all postions in current A
          for(unsigned int i = 0; i < N; i++){
              //if there is not a one          
              if(!setA[i]){
                 //than for all elements j in alphabet
                 for(unsigned int j = 0; j < alphabetSize; j++){
                     //for all elements k already in belongToMe
                     for(unsigned int k = 0; k < belongToMe.size(); k++){
                         //append j to k and save result in help
                         help.push_back(belongToMe[k] + j * static_cast<unsigned int>(pow(static_cast<double>(alphabetSize),static_cast<double>(i))));//intPower(alphabetSize,i));
                     }//rof
                 }//rof
                 //update belongToMe with help
                 belongToMe = help;
                 //delete all from help
                 help.clear();
              }//fi
              else{
                    //if there is a one
                    //then get the letter which should stand at this position
                   unsigned int factor = omegaA % alphabetSize;
                   omegaA /= alphabetSize;
                   //append it to all elements in belongToMe
                   for(unsigned int j = 0; j < belongToMe.size(); j++){
                       belongToMe[j] += (factor * static_cast<unsigned int>(pow(static_cast<double>(alphabetSize),static_cast<double>(i))));//intPower(alphabetSize, i));
                   }//rof
              }//esle
          }//rof
          //calculate alphab by summing up the value from the sample distribution for all elements in belongToMe
          for(unsigned int i = 0; i < belongToMe.size(); i++){
              unsigned int index = belongToMe[i];                                
              alpha += (*pEmp)[index];
          }//rof
}


/*
 * PARAMS:
 * pk - distribution
 * RETURN:
 * norm for given distribution 
 */
double Marginal::getNorm(vector<double>* pk){
                 double norm = 0.0;
                 //sum up all ditribution values of th elements in belongToMe
                 for(unsigned int i = 0; i < belongToMe.size(); i++){
                     norm += (*pk)[belongToMe[i]];
                 }//rof
                 
                 return norm;
}

/*
 * RETURN:
 * alpha of this marginal
 */
double Marginal::getAlpha(){
                 return alpha;
}

/*
 * RETURN:
 * iterator set to the begin of the vector belongToMe
 */
vector<unsigned int>::iterator Marginal::getWordsBegin(){
                                         return belongToMe.begin();
}

/*
 * RETURN:
 * iterator set to the end of the vector belongToMe
 */
vector<unsigned int>::iterator Marginal::getWordsEnd(){
                                         return belongToMe.end();
}
