/*
 *  main.cpp 
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
#include <sstream>
#include <list>
#include <vector>
#include <map>
#include <ctime>
#include <cctype>
#include <boost/thread/thread.hpp>
#include <boost/bind.hpp>

#ifdef WIN32
       #include <windows.h>
#else 
      #include <unistd.h>
#endif

#include "ConfigFile.h"
#include "functions.h"
#include "Marginals.h"
#include "Marginal.h"



using namespace std;

//vector used for synchronize output
vector< double > results;
vector<bool> semaphore;



/*
 * PARAMS:
 * order - vector of unsigned int representing all "A's"
 * N - length of one "A" in binary representation
 * pEmp - empirical distribution vector 
 * order_i - size of A if representated as set
 * alphabetSize - size of th eused alphabet
 * p - pointer to the distribution vector which should be updated
 * IT - numbers of iteration which should be performed
 */
void calculateProjection(vector<unsigned int> order, int N, vector<double> pEmp, unsigned int order_i, unsigned int alphabetSize, vector<double>* p, int IT){
     //create Marginals
     Marginals marginals(order, N, &pEmp, order_i, alphabetSize);
     //do the itartiv scaling for p
     for(int j = 0; j < IT; j++){
         (*p) = (marginals).doScaling((*p));
     }//rof
     semaphore.pop_back();
}

/*
 * PARAMS:
 * n -  size of A
 * N - size of an A if represented as int
 * RETURN:
 * vector containing the all set A's with size n
 */
vector< vector<unsigned int> > createAn(unsigned int n, unsigned int N){
                               //initialize variables
                               vector< vector<unsigned int> > orders;
                               vector<unsigned int> dummy;
                               orders.reserve(N);
                               for(unsigned int i = 0; i <= N; i++)orders.push_back(dummy);
                               bool binary[N];
                               unsigned int order = 0;
                               unsigned int N2 = static_cast<unsigned int>(pow(2.0,static_cast<double>(N+1)))-1;
                               //find and save all A with n elements
                               for(unsigned int i = 1; i < N2 ; i++){
                                   order = makeBinary(binary, N, i);
                                   if(order == n)orders[order].push_back(i);
                               }//rof
                               return orders;
}


/*
 * PARAMS:
 * hypergraph - hypergraph used for scaling
 * orders - vector of unsigned int representing all "A's"
 * N - length of one "A" in binary representation
 * pEmp - empirical distribution vector 
 * order_i - size of A if representated as set
 * alphabetSize - size of th eused alphabet
 * p - pointer to the distribution vector which should be updated
 * IT - numbers of iteration which should be performed
 * res_i - position in result vector, used for synchronize output
 */
void calculateSpecialProjection(string hypergraph, vector< vector<unsigned int> >orders, int N, vector<double> pEmp, unsigned int alphabetSize, map<char,unsigned int>* alphabet, int IT, unsigned int res_i, string outputprefix){
     //initialize variables
     const unsigned int NX = static_cast<unsigned int>(pow(static_cast<double>(alphabetSize),static_cast<double>(N)));//intPower(alphabetComplete.size(),N);  //alphabet^N = number of possible different strings
     vector<double> current;
     for(unsigned int i = 0; i < NX; i++){
         current.push_back(1.0/NX);
     }//rof
     
     //write hypergraph into vector
     if(hypergraph[0] == 'A' || hypergraph[0] == 'a'){
        //given as A1 or A2 or a3 or ...
        string hypergraph2 = hypergraph;
        hypergraph2 = hypergraph2.erase(0,1);
        stringstream hypstream(hypergraph2);
        unsigned int n;
        hypstream >> n;
        orders = createAn(n, N);
     }//fi
     else{
          //given as (1),(1,2),...
          orders = findMaxElements(hypergraph,N);
     }//esle
     
     //create marginals
     list<Marginals> all;
     for(unsigned int i = 0; i < orders.size(); i++){
         if(orders[i].size()>0){
            Marginals marginals(orders[i], N, &pEmp, i, alphabetSize);
            all.push_back(marginals);
         }//fi
     }//rof
     
     
     //calculate projection
     for(int i = 0; i < IT; i++){
         list<Marginals>::iterator allend = all.end();
         for(list<Marginals>::iterator j = all.begin(); j != allend; ++j){
             current = (*j).doScaling(current);           
         }//rof
     }//rof
     
     //output projection to file if wanted
     if(outputprefix != ""){
        //write file for p
        //prepare filename
        string filename(outputprefix+"_"+hypergraph+"_p.txt");
        const char* fname = filename.c_str();
        //initalize file stream
        ofstream* out = new ofstream();
        (*out).open(fname);
        map<unsigned int,char> ralphabet = reverseAlphabet(alphabet);
        cerr << "write p into file \"" << filename << "\"" << endl;
        writeProjectionToFile(&current,out,&ralphabet,N);
        //write file fpr p empirical
        //prepare filename
        filename = outputprefix+"_"+hypergraph+"_p_empirical.txt";
        const char* fname2 = filename.c_str();
        //initalize file stream
        (*out).open(fname2);
        cerr << "write p empirical into file \"" << filename << "\"" << endl;
        writeProjectionToFile(&pEmp,out,&ralphabet,N);
     }//fi
     //write results into global result vector
     results[res_i]=KullbackLeiberDistance(pEmp, current);
     semaphore.pop_back();
}



int main(int argc, char *argv[]){
    if(argc<2){
       //print cl options and exit Cipi if there too few arguments
       cerr << "too few arguments" << endl;
       cerr << "use: Cipi [-?|<PARAMFILE> <INPUTFILE>]" << endl;
       return 1;
    }//fi
    string firstcloption(argv[1]);
    if(firstcloption == "-?" || firstcloption == "-h" || firstcloption == "--help"){
       //print cl options with explanation if first argument is "-?" and exit Cipi
       cerr << "use: Cipi [-?|<PARAMFILE> <INPUTFILE>]" << endl;
       cerr << "\t-h,--help,-?\t print this help" << endl;
       cerr << "\tPARAMFILE\t path to the parameter file" << endl;
       cerr << "\tINPUTFILE\t path to the input file"<< endl;        
       return 1;
    }//fi
    if(argc<3){
       //print cl options and exit Cipi if there too few arguments
       cerr << "too few arguments" << endl;
       cerr << "use: Cipi [-?|<PARAMFILE> <INPUTFILE>]" << endl;
       return 1;
    }//fi
    
    
    //get time of start
    time_t t1 = time(0);
    
    //initialize variables with value from cl and parameter file
    ifstream* sampleInFile = new ifstream(argv[2]);
    if(!(*sampleInFile).is_open()){
       cerr << "Could not open input file " << argv[2] << endl;
       return 1;
    }
    ifstream* confFile = new ifstream(argv[1]);
    if(!(*confFile).is_open()){
       cerr << "Could not open parameter file " << argv[1] << endl;
       return 1;
    }else{
          (*confFile).close();
    }
    ConfigFile config(argv[1]);
    const unsigned int N = config.read<unsigned int>("N");
    const int IT = config.read<int>("SetIterations");
    const string inputType = config.read<string>("InputType");
    const string alphabetComplete = config.read<string>("Alphabet");
    const string outputfileprefix = config.read<string>("OutputFilePrefix");
    string hypergraph = config.read<string>("Hypergraph");
    const unsigned int maxProcessNum = config.read<int>("NumberOfProcesses");
    map<char,unsigned int> alphabet;
    const unsigned int NX = static_cast<unsigned int>(pow(static_cast<double>(alphabetComplete.size()),static_cast<double>(N)));
    vector< vector<double> > p;
    vector<double> pEmp;
    
    //create encoding of the alphabet
    for(unsigned int i = 0; i< alphabetComplete.size(); i++){
        char letter = alphabetComplete[i];
        //convert upper case letters to lower case letters
        if(isalpha(letter)&& isupper(letter))letter = tolower(letter);
        alphabet[letter] = i;        
    }//rof
    
    //read input and calculate empirical distribution of the samples
    if(inputType == "CharacterSequence")pEmp = readInputViaSlidingWindow(sampleInFile,N,NX,&alphabet);
    else if(inputType == "Sample")pEmp = readInputFromLine(sampleInFile,N,NX,&alphabet);
    else if(inputType == "Integer")pEmp = readInputFromLineToBinaryString(sampleInFile,NX,&alphabet);
    else if(inputType == "Empirical")pEmp = readInputFromPEmp(sampleInFile,N,NX,&alphabet);
    else{
         cerr << "Wrong input typ! " << inputType << "is not an input type." << endl;
         return 1;
    }//esle
    //get time after read input
    time_t t3 = time(0);
    
    //initalize varibales
    vector< vector<unsigned int> > orders;
    vector<unsigned int> dummy;
    orders.reserve(N);
    for(unsigned int i = 0; i < N; i++){
        orders.push_back(dummy);
    }//rof
    vector<double> current;
    for(unsigned int i = 0; i < NX; i++){
        current.push_back(1.0/NX);
    }
    
     //choose right mode
    if(hypergraph == ""){
       //classical version: use all "A's"
       p.reserve(N);
       //create vector containg "A's"
       for(unsigned int i = 0; i < N; i++){
           p.push_back(current);        
       }//rof
       bool binary[N];
       unsigned int order = 0;
       unsigned int N2 = static_cast<unsigned int>(pow(2.0,static_cast<double>(N)))-1;//intPower(2,N)-1;
       for(unsigned int i = 1; i < N2 ; i++){
           order = makeBinary(binary, N, i);
           orders[order].push_back(i);
       }//rof
       boost::thread_group thrds;
       //calculate projection for every order(=subset containg all "A's" with same size)
       //use one thread per order
       for(unsigned int i = 1; i < orders.size(); i++){
           //while already max number of processes wait
           while(semaphore.size() >= maxProcessNum-1){
                 #ifdef WIN32
                 Sleep(10*1000); 
                 #else 
                 sleep(10); 
                 #endif
           }
           semaphore.push_back(1);
           thrds.create_thread(boost::bind(&calculateProjection,orders[i], N, pEmp, i, alphabet.size(), &(p[i]), IT));
           
       }//rof
       
       //wait until all threads finished there work
       thrds.join_all();
       if(outputfileprefix != ""){
          //reverse alphabet
          map<unsigned int,char> ralphabet = reverseAlphabet(&alphabet);
          //write p vectors to file
          for(unsigned int i = 0; i < p.size(); i++){
              //prepare filename
              stringstream i_str;
              i_str << i;
              string filename(outputfileprefix+"_p("+i_str.str()+").txt");
              const char* fname = filename.c_str();
               //initalize variables
              ofstream* out = new ofstream();
              (*out).open(fname);
              //write to file
              cerr << "writing p(" << i << ") to file " << filename << endl;
              writeProjectionToFile(&(p[i]),out,&ralphabet,N);
          }//rof
          //write p empirical to file 
          //prepare filename
          string filename(outputfileprefix +"_p_empirical.txt");
          const char* fname = filename.c_str();
          //initalize variables
          ofstream* out = new ofstream();
          (*out).open(fname);
          //write to file
          cerr << "writing p empirical to file " << filename << endl;
          writeProjectionToFile(&pEmp,out,&ralphabet,N);
       }//fi
       
       //print I-vector(=vector of Kullback-Leiber-Distances)
       cout << p.size() << " : " 
            << KullbackLeiberDistance(pEmp, p.at(N-1)) << endl;
       for(unsigned int i = p.size()-1; i > 0; i--){
           cout << i << " : " 
                << KullbackLeiberDistance(p.at(i), p.at(i-1)) << endl;
       }//rof
    }//fi
    else{
         //only same defined "A's", fullfilling the characteristic of max. elemtents
         boost::thread_group thrds2;
         vector<string> hyps;
         while(hypergraph != ""){
               //while still hypergraphs given execute Cipi for them
               
               //while already max number of processes wait
               while(semaphore.size() >= maxProcessNum-1){
                     #ifdef WIN32
                     Sleep(10*1000); 
                     #else 
                     sleep(10); 
                     #endif
               }
               semaphore.push_back(1);
               //initialize/reset same variables
               results.push_back(0.0);
               string hypergraphnew;
               unsigned int semicolon;
               semicolon = hypergraph.find_first_of(';');
               if(semicolon >= hypergraph.length())semicolon = hypergraph.length();
               hypergraphnew.assign(hypergraph,0,semicolon);
               hyps.push_back(hypergraphnew);
               hypergraph = hypergraph.erase(0,semicolon+1);
               //create new thread for hypergraph and start thread
               thrds2.create_thread(boost::bind(&calculateSpecialProjection,hypergraphnew,orders,N,pEmp,alphabet.size(),&alphabet,IT,hyps.size()-1,outputfileprefix));
         }//elihw
         
         //wait until all threads finished there work
         thrds2.join_all();
         
         //print all results
         for(unsigned int i = 0; i < results.size(); i++){
             cout << "result for hypergraph " << hyps[i] << endl;
             cout << endl << "D(p emp || p)=" 
                  << results[i] << endl;
         }//rof
    }//esle
    
    //get time after finishing all
    time_t t2 = time(0);
    //print time needed for calculation
    cerr << endl << endl;
    cerr << "run cipi in " << difftime(t2,t1) << "s" << endl;
    cerr << "read input in " << difftime(t3,t1) << "s" << endl;
    cerr << "calculate projection in " << difftime(t2,t3) << "s" << endl;
    return 0;
}
