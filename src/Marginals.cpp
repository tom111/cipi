/*
 *  Marginals.cpp 
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
#include <list>

#include "Marginals.h"
#include "functions.h"
#include "Marginal.h"

using namespace std;

/*
 * PARAMS:
 * setsA - all A's for that marginals should be created
 * N - length of one sample
 * pEmp - distribution of the samples
 * order - current order 
 * alphabetSize - numbers of elementsin the alphabet
 */
Marginals::Marginals(vector<unsigned int> setsA,unsigned int N, vector<double>* pEmp, unsigned int order,unsigned int alphabetSize){
           //initialize varibales
           //vector represent "Omega_A"
           vector<unsigned int> omegaAs;
           //vector needed to calculate "Omega_A"
           vector<unsigned int> omegaAsNew;
           //number of elements in "Omege_A"
           unsigned int omegaSize = static_cast<unsigned int>(pow(static_cast<double>(alphabetSize),static_cast<double>(order)));//intPower(alphabetSize,order);
           //initialize vectors
           omegaAs.reserve(omegaSize);
           omegaAsNew.reserve(omegaSize);
           omegaAs.push_back(0);
           unsigned int exp = 0;
           //for all length i of element in Omega_A
           for(unsigned int i = 0; i < order; i++){
               //for all elements j in alphabet
               for(unsigned int j = 0; j < alphabetSize; j++){
                   //for all elemts k already in "Omega_A"
                   for(unsigned int k = 0; k < omegaAs.size(); k++){
                       //append j to k and save new element
                       unsigned int pushback = omegaAs.at(k) + j * static_cast<unsigned int>(pow(static_cast<double>(alphabetSize),static_cast<double>(exp)));//intPower(alphabetSize, exp);
                       omegaAsNew.push_back(pushback);
                   }//rof
               }//rof
               exp++;
               omegaAs = omegaAsNew;
               omegaAsNew.clear();
           }//rof
           //for all A's 
           for(unsigned int h = 0; h < setsA.size(); h++){
               //make A binary
               unsigned int setA = setsA[h];
               bool binary[N];
               makeBinary(binary, N, setA);
               //for all elements of "Omagea_A"
               for(unsigned int i = 0; i < omegaAs.size(); i++){
                   //create and save marginal
                   Marginal marginal(binary, omegaAs.at(i), pEmp, N, alphabetSize);
                   marginals.push_back(marginal);
               }
           }//rof
}

/*
 * PARAMS:
 * pk - distribution for one  iterative scaling run
 * RETURN:
 * updated distribution
 */
vector<double> Marginals::doScaling(vector<double> pk){
                          //get end of marginals
                          mEnd = marginals.end();
                          //for all marginals
                          for (list<Marginal>::iterator i = marginals.begin(); i != mEnd; ++i){
                               //calculate c
                               double c = (*i).getAlpha() / (*i).getNorm(&pk);
                               //get end of elements belong to marginal
                               vector<unsigned int>::iterator end = (*i).getWordsEnd();
                               //for all element of the maginal
                               for(vector<unsigned int>::iterator j = (*i).getWordsBegin(); j != end; ++j){
                                    //if not already null update with c
                                    if(pk[*j] != 0.0){
                                       pk[*j] *= c;
                                    }//fi
                                }//rof
                          }//rof                    
                          
                          return pk;
}
