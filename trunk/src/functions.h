#ifndef FUNCTIONS_H
#define FUNCTIONS_H

#include <cstdio>
#include <stdlib.h>
#include <unistd.h>
#include <math.h>
#include <iostream>
#include <vector>
#include <set>
#include <map>
#include <string>
#include <fstream>
#include <sstream>
#include <cctype>


//calculate binary representation with size size of number as boolean array
//returns number of "1" in binary representation of number
unsigned int makeBinary(bool *binary,unsigned int length,unsigned int number);


//calculates the Kullback-Leiber Distance of two vectors<double>
//under the condition that word at pos i of vector 1 = w at pos i of vector 2
double KullbackLeiberDistance(std::vector<double> p, std::vector<double> q);

//read input from given file stream with sliding window of size windowSize
//returns vector of samles
std::vector<double> readInputViaSlidingWindow(std::ifstream* samplesfile,unsigned int windowsSize, unsigned int ele_nr, std::map<char,unsigned int>* alphabet);

//read input from given file stream line per line
std::vector<double> readInputFromLine(std::ifstream* samplesfile,unsigned int N, unsigned int ele_nr, std::map<char,unsigned int>* alphabet);

//read input from given file stream line per line and make them binary
//so there should be only one integer values per line
std::vector<double> readInputFromLineToBinaryString(std::ifstream* samplesfile, unsigned int ele_nr, std::map<char,unsigned int>* alphabet);

//calculate binary representation of the hypergraph given as string and add all sets containing just one element
std::vector< std::set<unsigned int> > makeHypergraphBinary(std::string set, unsigned int N);

//find max. elements of the hypergraph given bei the string
std::vector< std::vector<unsigned int> > findMaxElements(std::string hyp,unsigned  int N);

//write projection vector to file
void writeProjectionToFile(std::vector<double>* p,std::ofstream* pFile, std::map<unsigned int,char>* alphabet, unsigned int N);

//reverse the alphabet
std::map<unsigned int,char> reverseAlphabet(std::map<char,unsigned int>* alphabet);
#endif
