/*
 *  Marginal.h 
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
