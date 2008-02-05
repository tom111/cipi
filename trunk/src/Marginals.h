#ifndef MARGINALS_H
#define MARGINALS_H

#include <cstdlib>
#include <stdlib.h>
#include <math.h>
#include <iostream>
#include <fstream>
#include <string>
#include <vector>
#include <list>

#include "functions.h"
#include "Marginal.h"




class Marginals{

//list of all marginal
std::list<Marginal> marginals;
//iterator which will be set to the end of the list of marginals
std::list<Marginal>::iterator mEnd;

public:
       //create all marginals
       Marginals(std::vector<unsigned int> setsA,unsigned int N, std::vector<double>* pEmp, unsigned int order,unsigned int alphabetSize);
       //will do the one iterative scaling for the given distribution
       std::vector<double> doScaling(std::vector<double> pk);
};

#endif
