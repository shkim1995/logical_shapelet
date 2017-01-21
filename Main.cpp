
#include "Histogram.h"
#include <limits>
#include <cstddef>
#include <cstring>

#ifndef HEADERS

#define HEADERS

#include <iostream>
#include <fstream>
#include <string>
#include <cmath>
#include <vector>
#include <list>
#include <sstream>
#include "time.h"

#endif


using namespace std;


int MAXLEN = 60;
int MINLEN = 30;


//////////////////////////////////////////////////////////////////////////
///////////////////////////CLASS DEFINITION///////////////////////////////
//////////////////////////////////////////////////////////////////////////

class Data{

	public:
		int type;
		vector<double> data;
		int type2;

		void print(){
			cout<<type<<" | ";
			for(int i=0; i<data.size(); i++)
				cout<<" "<<data[i];
			cout<<endl;
			cout<<endl;
		}

		int size(){
			return data.size();
		}


};



//////////////////////////////////////////////////////////////////////////
///////////////////////////STRUCT DEFINITION//////////////////////////////
//////////////////////////////////////////////////////////////////////////


typedef struct _ShapeletInfo
{
	vector<double> best_S;
	double best_t;
	double maxGain;
	double maxGap;

}ShapeletInfo;


//////////////////////////////////////////////////////////////////////////
/////////////////////////////MAIN FUNCTIONS///////////////////////////////
//////////////////////////////////////////////////////////////////////////


//FOR DEBUGGINH
void printVector(vector<double> S){
	for(int i=0; i<S.size(); i++)
		cout<<S[i]<<" ";
	cout<<endl;
}


void printShapelet(ShapeletInfo info){
	cout<<"Shapelet : "<<endl;
	printVector(info.best_S);
	cout<<"Distance : "<<info.best_t<<endl;

}


//normalize dataVector[from] ~ dataVector[to-1]
vector<double> zNorm(vector<double>* dataVector, int from, int to){

    vector<double> normalizedVector;

	double sum = 0;
	int num = 0;
	double avg = 0;
	double var = 0;

	for(int i=from; i<to; i++){
		sum = sum + (*dataVector)[i];
		num++;
	}

	if(num!=0)	
		avg = sum/num;
	
	for(int i=from; i<to; i++){
		double x = (*dataVector)[i];
		var = var + (x-avg)*(x-avg);
	}
	var = var/(num-1);
	var = sqrt(var);

	for(int i=from; i<to; i++){
		normalizedVector.push_back(((*dataVector)[i]-avg)/var);
	}

	return normalizedVector;
	
}


double sdist(vector<double>* x, vector<double>* y){ // |x|<=|y|

	double minSum = 100000;
	double sum;

	vector<double> x_norm = zNorm(x, 0, x->size());
	int l = x->size();

	for(int i=0; i<y->size()-l; i++){
		sum = 0;
		vector<double> y_norm = zNorm(y, i, i+l);
		for(int k=0; k<l; k++){
			sum = sum+(y_norm[k]-x_norm[k])*(y_norm[k]-x_norm[k]);
		}

		if(minSum>sum)
			minSum = sum;
	}

	return sqrt(minSum/l);
}


//Entropy before division
// a, b : number of each elements 
double Entropy(int  a, int b){
	if(a==0) return 0;
	if(b==0) return 0;
	double pa = (a+0.0)/(a+b+0.0);
	double pb = (b+0.0)/(a+b+0.0);
	return -(pa*log(pa)+pb*log(pb));
}


//Entropy after division
double Entropy(int a1, int b1, int a2, int b2){
	return ((a1+b1)*Entropy(a1, b1)+(a2+b2)*Entropy(a2, b2))/(a1+a2+b1+b2);
}



bool BestIG(Histogram* hist, double* max_gain, double* max_gap){
	bool temp = false;

	list<HistData>::iterator it = hist->data.begin();
	
	int a1 = 0;
	int b1 = 0;
	int a2 = 0;
	int b2 = 0;

	double distSum_left = 0;
	double distSum_right = 0;
	
	while(it!=hist->data.end()){
		if(it->type==0)
			a2++;
		else
			b2++;
		distSum_right = distSum_right+it->value;
		it++;
	}

	cout<<"a "<<a2<<" b "<<b2<<" dist "<<distSum_right<<endl;
	
	it = hist->data.begin();

	double totalEntropy = Entropy(a2, b2);
	double entropy;
	double gain;

	while(it!=hist->data.end()){

		entropy = Entropy(a1, b1, a2, b2);
		gain = totalEntropy-entropy;

		if(gain>*max_gain){
			
		}
		it++;
	}

	return temp;
}

void Shapelet_Discovery(vector<Data>* dataVector){

	ShapeletInfo Shapelet;

	double max_gain = 0;
	double max_gap = 0;

	//for each data
	for(int k=0; k<dataVector->size(); k++){
		Data d = (*dataVector)[k];

		//for every subsequece
		for(int i=0; i<d.size(); i++){
			for(int j=i+MINLEN; j<fmin(i+MAXLEN, d.size()); j++){

				//make candidate
				vector<double> candidate;
				candidate.resize((j-i+1));
				memcpy(&(candidate[0]), &(d.data[i]), sizeof(double)*(j-i+1));

				//make Histogram
				Histogram hist;
				for(int t=0; t<dataVector->size(); t++){
					double distance = sdist(&candidate, &(*dataVector)[t].data);
					hist.insert((*dataVector)[t].type, distance);
				}

				BestIG(&hist, &max_gain, &max_gap);

			}
		}
	}
}


int main(){
	clock_t before;

	before = clock();

	vector <Data> dataVector;
 	ifstream inFile;
	inFile.open("Coffee_TRAIN", ios::in);
	string s;

	while(getline(inFile, s)){
		stringstream ss(s);
		double test;
		int t;
		Data temp;

		ss>>t;
		temp.type = t;
		
		while(ss>>test){
			temp.data.push_back(test);
		}
		dataVector.push_back(temp);
	}


	Shapelet_Discovery(&dataVector);

	// vector<double> test;
	// test.push_back(1);
	// test.push_back(2);
	// test.push_back(3);
	// test.push_back(2);
	// test.push_back(4);
	// test.push_back(1);
	// test.push_back(2);
	// test.push_back(3);
	// test.push_back(1);
	// test.push_back(4);



	// ShapeletInfo info;
	// info.best_S = test;
	// info.best_t = 10;
	// printShapelet(info);



	// // vector<double> test2 = zNorm(&test, 0, test.size());
	// // printVector(test2);
	

	// vector<double> test2;
	// test2.push_back(3);
	// test2.push_back(2);
	// test2.push_back(4);
	// test2.push_back(1);

	// cout<<sdist(&test2, &test);
}