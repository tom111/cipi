#include "functions.h"

using namespace std;

/*
 * PARAMS:
 * binary - array which will contasin binary representation
 * length - number of bit which should be used
 * number - number which should be make binary
 * RETURN:
 * number of "1" in binary
 */
unsigned int makeBinary(bool *binary,unsigned int length,unsigned int number){
             //number of "1" in binary representation of number
             static unsigned int size;
             size = 0;
             for(unsigned int i = 0; i < length ; i++ ){
                if (number & (1 << i)) {
                    binary[i] = 1; 
                    size++;
                }//fi
	               else binary[i] = 0;
             }//rof
             return size;
}

/*
 * PARAMS:
 * p - first vector
 * q - second vector
 * RETURN:
 * D( p || q ) - Kullback-Leiber-Distance of p and q 
 */
double KullbackLeiberDistance(vector<double> p, vector<double> q){
       double distance = 0.0;
       //sum up over all elements in p 
       for(unsigned int i = 0; i < p.size(); i++){
           //sum up if p not 0
           if (p[i] != 0) distance += ((log(p[i]/q[i])/0.6931471805599452862)*p[i]);
       }//rof
       return distance;
}

/*
 * PARAMS
 * samplesfile - input file stream for samples file
 * windowsSize - number of char's building one sample
 * ele_nr - number of different elements for distribution
 * alphabet - mapping from used alphabet to unsigned int's
 * RETURN:
 * samples distribution for all samples in file as vector of doubles
 */
vector<double> readInputViaSlidingWindow(ifstream* samplesfile, unsigned int windowsSize, unsigned int ele_nr, map<char,unsigned int>* alphabet){
               //initialize varibales
               //vector containg distribution of the samples 
               vector<double> samples;
               //count number of samples
               unsigned int samplesize = 0;
               //binary representation of one sample
               unsigned int oneSample = 0;
               //count letter currently read into sample
               unsigned int letterread = 0;
               //char currently read from file
               char nuc;
	           //error flag
               bool invalid = false;
               //needed to split first postion of a sample
               unsigned int mod = ele_nr/(*alphabet).size();
               //prepare samples vector
               for(unsigned int i = 0; i < ele_nr; i++)samples.push_back(0.0);
               
               //if file could be opened
               if((*samplesfile).is_open()){
                 //untill in file still left char's
                 while((*samplesfile) >> nuc ){
                       //needed to read a line
                       string line("");
                       //if line began with space skip line
                       if(nuc == ' '){
                          getline((*samplesfile),line);
                          continue;
                       }//fi
                       //if line is header
                       if(nuc == '>'){
                          //read complete line
                          getline((*samplesfile),line);
                          //new gene/gene part began
                          //reset variables
                          oneSample = 0;
                          letterread = 0;
                          continue;
                       }//fi
                       //converts all upper case letter to lower case letter
                       if(isalpha(nuc)&&isupper(nuc))nuc = tolower(nuc);
                       //push back "0"
                       oneSample *= (*alphabet).size();
                       //if nuc not in alphabet
                       if((*alphabet).find(nuc) == (*alphabet).end()){
		       	   if(!invalid){
				//set error flag 
				invalid = true;
				//inform user
				cerr << "Invalid sample! (All invalid samples will be ignored!)" << endl;
			    }//fi
                           //reset all and next char
                           oneSample = 0;
                           letterread = 0;
                           continue;
                       }//fi                
                       //add number of last read letter
                       oneSample += (*alphabet)[nuc];
                       //count letter
                       letterread++;
                       //if in current sample right number of letters
                       if(letterread == windowsSize){
                          //count sample and save it in samples vector
                          samples[oneSample]+=1.0;
                          samplesize++;
                          //split first position
                          oneSample %= mod;
                          //decrease count of letter in samples
                          letterread--;
                       }//fi
                 }//elihw      
               }//fi
               //close file
               (*samplesfile).close();
               //create ditribution of samples by dividing through numbers of samples
               for(unsigned int i = 0; i < ele_nr; i++)samples[i] /= samplesize;
               //notify user that no sample could be found
               if(samplesize == 0)cerr << "Could not find samples in input file!" << endl;
               return samples;
}

/*
 * PARAMS:
 * samplesfile - input file stream for file containing empirical p
 * N - length of one samples
 * ele_nr - number of different elements for distribution
 * alphabet - mapping from used alphabet to unsigned int's
 * RETURN:
 * samples given p empirical as vector of doubles
 */
vector<double> readInputFromPEmp(ifstream* samplesfile, unsigned int N, unsigned int ele_nr, map<char,unsigned int>* alphabet){
               //initialize variables
               //vector containing samples distribution
               vector<double> samples;
               //sample as string
               string sample("");
               //probability as string
               string prob_str("");
               //string stream for converting probability from string to double
               stringstream converter;
               //probability as double
               double probability = 0.0;
               //unsigned int coding sample
               unsigned int sampleCode = 0;
               //error flag
               bool error = false;
               //probability's sum 
               double sum = 0.0;
               //initialize samples vector
               for(unsigned int i = 0; i < ele_nr; i++)samples.push_back(0.0);
               
               //open file
               if((*samplesfile).is_open()){
                  //until reaching end of file
                  while(!(*samplesfile).eof()){
                       //initialize string for line
                       string line("");
                       //read line
                       getline((*samplesfile),line);
                       //if line is not in correct format go to next line
                       if(line == "" || line.find("\t") == line.size())continue;
                       //split into sample and probability
                       sample = line.substr(0,line.find("\t"));
                       prob_str = line.substr(line.find("\t") + 1);
                       //if the sample is too short or too large go to next line
                       if(sample.length() != N){
                          //if first error report it
                          if(!error){
                             cerr << "Invalid sample! (All invalid samples will be ignored!)" << endl;
                             error = true;
                          }
                          continue;
                       }
                       //if probability string is not correct got to next line
                       string numbers("0123456789e-");
                       if(prob_str.find_first_of("0") != 0 
                          || (prob_str.find_first_of(".") != 1 && prob_str.find_first_of(".")< prob_str.length())
                          || (prob_str.find_last_not_of(numbers) != 1 && prob_str.find_last_not_of(numbers) < prob_str.length())
                          || (prob_str.find_first_of("e") != prob_str.length() && prob_str.find_first_of("e") != prob_str.find_last_of("e"))
                          || (prob_str.find_first_of("-") != prob_str.length() && prob_str.find_first_of("-") != prob_str.find_last_of("-"))){
                          //if first error
                          if(!error){
                             cerr << "Invalid sample! (All invalid samples will be ignored!)" << endl;
                             error = true;
                          }
                          continue;
                       }
                       //reset sampleCode
                       sampleCode = 0;
                       //initalize baseExp
                       unsigned int baseExp = 1;
                       //unknow sign error
                       bool unknownSign = false;
                       for(unsigned int i = 0; i < N; i++){
                           //get char at current position
                           char letter = sample[sample.length() - 1 - i];
                           //convert to lower case
                           if(isalpha(letter) && isupper(letter))letter = tolower(letter);
                           //find letter
                           if((*alphabet).find(letter) == (*alphabet).end()){
                              //if first error report it
                              if(!error){
                                 cerr << "Invalid sample! (All invalid samples will be ignored!)" << endl;
                                 error = true;
                              }
                              unknownSign = true; 
                              break;
                           }//fi
                           //add to samplecode
                           sampleCode += ((*alphabet)[letter]*baseExp);
                           //update baseExp for next position
                           baseExp *= (*alphabet).size();
                       }//rof
                       //unknown letter in sample, go to next line
                       if(unknownSign)continue;
                       //convert probability
                       converter.clear();
                       converter.str(prob_str);
                       converter >> probability;
                       //sum up probability
                       sum += probability;
                       //save sample in vector
                       samples[sampleCode] = probability;
                  }//elihw
               }//fi
               //if sum of probability is not 1 than warn user
               if(sum < 1.0) {
		 cout << "Sum of probability < 1, normalizing !" << endl;
		 for(unsigned int i = 0; i < ele_nr; i++) samples[i] /= sum;
	       }
               if(sum > 1.0) {
		 cout << "Sum of probability > 1, normalizing !" << endl;
		 for(unsigned int i = 0; i < ele_nr; i++) samples[i] /= sum;
	       }
	       // normalize :
	       

               //return vector
               return samples;
}

/*
 * PARAMS:
 * samplesfile - input file stream for file containing samples
 * N - length of one samples
 * ele_nr - number of different elements for distribution
 * alphabet - mapping from used alphabet to unsigned int's
 * RETURN:
 * samples distribution for all samples in file as vector of doubles
 */
vector<double> readInputFromLine(ifstream* samplesfile,unsigned int N, unsigned int ele_nr, map<char,unsigned int>* alphabet){
               //initialize variables
               //vector containg samples distribution
               vector<double> samples;
               //number of samples
               unsigned int samplesize = 0;
               //one sample as string
               string onesample("");
	           //error flag 
	           bool invalid = false;
               //initialize vector for samples
               for(unsigned int i = 0; i < ele_nr; i++)samples.push_back(0.0);
               //if file could be opened
               if((*samplesfile).is_open()){
                  //while still lines left in file
                  while((*samplesfile) >> onesample){
                         //reverse string
                         reverse(onesample.begin(),onesample.end());
                        //intrepresentation of sample
                        unsigned int samplevalue = 0;
                        //current exponent
                        unsigned int baseExp = 1;
                        //for all postions in sample
                        for(unsigned int i = 0; i < onesample.size() ; i++){
                             //converts all upper case letter to lower case letter
                             if(isalpha(onesample[i])&&isupper(onesample[i]))onesample[i] = tolower(onesample[i]);
                             if((*alphabet).find(onesample[i]) == (*alphabet).end())break;
                            //add to int representation current postion
                            samplevalue += ((*alphabet)[onesample[i]]*baseExp);
                            //update exponent
                            baseExp *= (*alphabet).size();
                        }//rof
                        //if samplevalue is a valid sample
                        if(samplevalue >=0 && samplevalue < ele_nr && onesample.size() == N){
                           //count sample in distribution 
                           samples[samplevalue]++;
                           //increase samples counter
                           samplesize++;
                        }//fi
            else if(!invalid){
				//set error flag 
				invalid = true;
				//inform user
				cerr << "Invalid sample! (All invalid samples will be ignored!)" << endl;
			}//fi esle
                 }//elihw
             }//fi
             //close file
             (*samplesfile).close();
             //create samples distribution by dividing through samples counter
             for(unsigned int i = 0; i < ele_nr; i++)samples[i] /= samplesize;
             //notify user that no sample could be found
             if(samplesize == 0)cerr << "Could not find samples in input file!" << endl;  
             return samples;
}


/*
 * PARAMS:
 * samplesfile - input file stream for file containing samples
 * ele_nr - number of different elements for distribution
 * alphabet - mapping from used alphabet to unsigned int's
 * RETURN:
 * samples distribution for all samples in file as vector of doubles
 */
vector<double> readInputFromLineToBinaryString(ifstream* samplesfile, unsigned int ele_nr, map<char,unsigned int>* alphabet){
               //initialize variables
               //vector containing distribution of samples
               vector<double> samples;
               //counter for number aof samples
               unsigned int samplesize = 0;
               //one sample
               unsigned int onesample;
	       //error flag
	       bool invalid = false;
               //initalize distribution vector
               for(unsigned int i = 0; i < ele_nr; i++)samples.push_back(0.0);
               
               //if file could be opened
               if((*samplesfile).is_open()){
                  //while still samples left in file
                  while((*samplesfile) >> onesample){
                         if(onesample < samples.size() && onesample >= 0){
                            //count sample and increase samples counter
                            samples[onesample]++;
                            samplesize++;
                        }//fi
			else if(!invalid){
				//set error flag 
				invalid = true;
				//inform user
				cerr << "Invalid sample! (All invalid samples will be ignored!)" << endl;
			}//fi esle
                  }//elihw
               }//fi
               //close file
               (*samplesfile).close();
               //create distribution by dividing through number of samples
               for(unsigned int i = 0; i < ele_nr; i++)samples[i] /= samplesize;
               //notify user that no sample could be found
               if(samplesize == 0)cerr << "Could not find samples in input file!" << endl;
               return samples;
}

/*
 * PARAMS:
 * sets - hypergraph given as set of subsets of an alphabet
 * N - number of elements in the original alphabet
 * RETURN:
 * hypergraph representeted by unsigned int's given bei sets
 */
vector< set<unsigned int> > makeHypergraphBinary(string sets, unsigned int N){
       //prepare container for hypergraph
       vector< set<unsigned int> > returnMe;
       set<unsigned int> current;
       for(unsigned int i = 0; i <= N; i++)returnMe.push_back(current);
       
       //initialize variables
       unsigned int order = 0; 
       unsigned int binary = 0;
       //while there still subsets in hypergraph
       while(!sets.empty()){
           //extract subset from hypergraph
           unsigned int indexend = sets.find_first_of(")");
           string curElement = sets.substr(1,indexend-1);
           sets = sets.erase(0,indexend+1);
           unsigned int indexbegin = sets.find_first_of("(");
           sets = sets.erase(0,indexbegin);
           //while subset still have elements
           while(!curElement.empty()){
                      string curexpo;
                      //extrace element
                      if(curElement.find_first_of(",") < curElement.length()){
                          curexpo = curElement.substr(0,curElement.find_first_of(","));
                          curElement = curElement.erase(0,curElement.find_first_of(",")+1);
                      }//fi
                      else {
                              curexpo = curElement;
                              string dummy;
                              curElement = dummy;
                      }//esle
                      //increase counter for number of elements
                      order++;
                      //add element to binary representation of subset
                      istringstream expostream;
                      expostream.str(curexpo);
                      unsigned int exp;
                      expostream >> exp;
                      binary += static_cast<unsigned int>(pow(2.0,static_cast<double>(exp-1)));
           }//elihw
           //save binary subset in container for hypergraph
           returnMe[order].insert(binary);
           //reset
           order = 0;
           binary = 0;
       }//elihw
       //fill of with all one elements sets
       for(unsigned int i = 0; i < N; i++){
                    returnMe[1].insert( static_cast<unsigned int>(pow(2.0,static_cast<double>(i))));
       }//rof
       return returnMe;
}

/*
 * PARAMS:
 * hypergraph - hypergraph as string
 * N - number of different elements 
 * RETURN:
 * max. elements represented by unsigned ints of the given hypergraph
 */
vector< vector<unsigned int> > findMaxElements(string hypergraph, unsigned int N){
        //prepare container for max. elements of hypergraph
       vector<unsigned int> current;
       vector< vector<unsigned int> > returnMe;
       for(unsigned int i = 0; i <= N; i++)returnMe.push_back(current);
       vector < set <unsigned int > >temp;
       set<unsigned int> dummy;
       for(unsigned int i = 0; i <= N; i++)temp.push_back(dummy);
       vector< vector < set<unsigned int>::iterator> > deleteus;
       vector < set<unsigned int>::iterator> deldummy;
       for(unsigned int i = 0; i <= N; i++)deleteus.push_back(deldummy);
       //get binary representation of the hypergraph
       vector< set<unsigned int> > bin_hyp = makeHypergraphBinary(hypergraph,N);
       //for every possible number of elements
       for(unsigned int j = 0; j <= N; j++){
          set<unsigned int>::iterator orderend = bin_hyp[j].end();
          //for every subset with j elements
          for(set<unsigned int>::iterator i = bin_hyp[j].begin(); i != orderend; ++i){
             bool uncomparible = true;
             //for every element already max. element
             for(unsigned int l = 0; l < temp.size(); l++){
                 bool breakmax = false;
                 for(set<unsigned int>::iterator k = temp[l].begin(); k != temp[l].end(); ++k){
                     //a is part of b if a&~b==0
                     //if current subset part of current max elment
                     unsigned int test = ((*i)&~(*k));
                     if(test == 0){
                        //mark subset as comparible
                        uncomparible = false;
                        breakmax = true;
                        break;
                     }//fi
                     test = ((*k)&~(*i));
                     if(test == 0){
                        //if current max. element part of current subset
                        //mark subset as comparible, insert subset into max 
                        //elments and mark current max. element as deleted
                        uncomparible = false;
                        deleteus[l].push_back(k);
                        temp[j].insert((*i));
                     }//fi esle
                 }//rof
                 if(breakmax)break;
             }//rof
             //insert uncomparible subset to max. elements
             if(uncomparible)temp[j].insert((*i));
          }//rof
          //delete all subsets with j elemtents in bin_hyp
          //bin_hyp[j].clear();
       }//rof
       //delete deleteus in max elements
       for(unsigned int i = 0; i < deleteus.size(); i++){
           for(vector < set<unsigned int>::iterator>::iterator l = deleteus[i].begin(); l != deleteus[i].end(); ++l){
               if(temp[i].find((*(*l))) != temp[i].end())temp[i].erase(*l);
           }//rof
       }//rof
       //write into returnMe
       for(unsigned int i = 0; i < temp.size(); i++){
           for(set<unsigned int>::iterator l = temp[i].begin(); l != temp[i].end(); ++l){
               returnMe[i].push_back(*l);
           }//rof
       }//rof
       
       return returnMe;
}

/*
 * PARAMS:
 * p - vector which should be written to file
 * pFile - output file stream to file for output
 * alphabet - map assoziate integer with chars, contains reverse coding
 * N - length of one configuration
 */
void writeProjectionToFile(vector<double>* p,ofstream* pFile, map<unsigned int,char>* alphabet, unsigned int N){
                      //if output file isnt open print error and return
                      if(!(*pFile).is_open()){
                         cerr << "Could not open file for writing output!" << endl;
                         return;
                      }//fi
                      //initialize variables
                      unsigned int alSize = (*alphabet).size();
                      //for each position in i
                      for(unsigned int i = 0; i < (*p).size(); i++){
                          //initialize variables 
                          string conf("");
                          unsigned int cur = i;
                          //create configuration
                          for(unsigned int j = 0; j < N; j++){
                              unsigned int key = cur%alSize;
                              cur /= alSize;
                              conf = (*alphabet)[key] + conf;
                          }
                          stringstream val;
                          val << ((*p)[i]);
                          (*pFile) << (conf+"\t"+val.str()+"\n");
                      }//rof
                      (*pFile).close();
}


/*
 * PARAMS:
 * alphabet - alphabet which maps char on unsigend int
 * RETURN:
 * revered alphabet, maps unsigned int on char
 */
map<unsigned int,char> reverseAlphabet(map<char,unsigned int>* alphabet){
                       //map which will be returned
                       map<unsigned int,char> ralphabet;
                       //insert old keys als values and old values as keys
                       for(map<char,unsigned int>::iterator i = (*alphabet).begin(); i != (*alphabet).end(); ++i){
                           ralphabet[(*i).second] = (*i).first;
                       }//rof
                       //return result
                       return ralphabet;
}
