#ifndef MARGINAL_H
#define MARGINAL_H

#include <cstdlib>
#include <stdlib.h>
#include <math.h>
#include <iostream>
#include <fstream>
#include <string>
#include <vector>

#include "functions.h"



class Marginal{

//alpha of Marginal
double alpha;
//vector for collecting all samples belong to this marginal
std::vector<unsigned int> belongToMe;

public:
       //create a marginal 
       Marginal(bool setA[], unsigned int omegaA, std::vector<double>* pEmp,unsigned int N,unsigned int alphabetSize);
       //calculate norm for the given ditribution
       double getNorm(std::vector<double>* pk);
       //getter for alpha
       double getAlpha();
       //getter for begin of vector belongToMe
       std::vector<unsigned int>::iterator getWordsBegin();
       //getter for end of vector belongToMe
       std::vector<unsigned int>::iterator getWordsEnd();
};

#endif
