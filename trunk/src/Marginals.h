/*
 *  Marginals.h 
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
