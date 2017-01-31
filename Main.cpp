#include "Histogram.h"
#include <limits>
#include <cstddef>
#include <cstring>
 
#ifndef HEADERS
 
#define HEADERS
 
#include <iostream>
#include <iterator>
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
 
 
class Stats{
public:
 
	double** S; //array of SumArray
	double** S2; // array of SquareSumArray
	double**** M;
 
	Stats(vector<Data>* dataVector){
		
		//Inititate S, S2
		int dataLen = dataVector->size();
		S = new double* [dataLen];
		S2 = new double* [dataLen];
 
		int len;
		double sum;
		double squreSum;
		double temp;
 
		for(int i=0; i<dataLen; i++){
			len = (*dataVector)[i].size();
			sum = 0;
			squreSum = 0;
 
			S[i] = new double[len+1];
			S2[i] = new double[len+1];
			S[i][0] = 0;
			S2[i][0] = 0;
 
			for(int j=0; j<len; j++){
				temp = ((*dataVector)[i].data)[j];
				sum = sum + temp;
				squreSum = squreSum + temp*temp;
				S[i][j+1] = sum;
				S2[i][j+1] = squreSum;
			}
		}
 
		//Inititate M
		int x = 0;
		int y = 0;
 
		M = new double*** [dataLen];
		for(int i=0; i<dataLen; i++)
			M[i] = new double** [dataLen];
 
		for(int i=0; i<dataLen; i++){
			for(int j=0; j<dataLen; j++){
				
				x = (*dataVector)[i].size();
				y = (*dataVector)[j].size();
 
				//inititalization
				M[i][j] = new double* [x+1];
				for(int u=0; u<x+1; u++)
					M[i][j][u] = new double[y+1];
 
				//boundary condition
				for(int u=0; u<x+1; u++)
					M[i][j][u][0] = 0;
				for(int v=0; v<y+1; v++)
					M[i][j][0][v] = 0;
 
				for(int v=1; v<y+1; v++){
					for(int u=v; u<x+1; u++){
						M[i][j][u][v] = M[i][j][u-1][v-1]+((*dataVector)[i].data)[u-1]*((*dataVector)[j].data)[v-1];
						
					}
				}
 
				for(int u=1; u<x+1; u++){
					for(int v=u; v<y+1; v++){
						M[i][j][u][v] = M[i][j][u-1][v-1]+((*dataVector)[i].data)[u-1]*((*dataVector)[j].data)[v-1];
						
					}
				}
			}
		}
 
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
	var = var/(num);
	var = sqrt(var);
 
	for(int i=from; i<to; i++){
		normalizedVector.push_back(((*dataVector)[i]-avg)/var);
	}
 
	return normalizedVector;
	
}
 
//brute force
double sdist(vector<double>* x, vector<double>* y){ // |x|<=|y|
 
	double minSum = 100000;
	double sum;
 
	vector<double> x_norm = zNorm(x, 0, x->size());
	int l = x->size();
 
	for(int i=0; i<y->size()-l+1; i++){
		sum = 0;
		vector<double> y_norm = zNorm(y, i, i+l);
		for(int k=0; k<l; k++){
			sum = sum+(y_norm[k]-x_norm[k])*(y_norm[k]-x_norm[k]);
		}
		//cout<<sqrt(sum/l)<<" ";
		if(minSum>sum)
			minSum = sum;
	}
	//cout<<endl;
	return sqrt(minSum/l);
}
 
double sdist_new(vector<Data>* dataVector, Stats* stats, int i, int j, int u, int l){ 
	// x = Di, y = Dj
	// |x| <= |y|
 
	// cout<<stats->S2[i][0]<<" "<<stats->S2[i][1]<<" "<<stats->S2[i][2]<<" "<<stats->S2[i][3]<<" "<<stats->S2[i][4]<<" ";
	// cout<<stats->S2[i][5]<<" "<<stats->S2[i][6]<<" "<<stats->S2[i][7]<<" "<<stats->S2[i][8]<<" "<<endl;
	// cout<<endl;
 
	// cout<<stats->S2[j][0]<<" "<<stats->S2[j][1]<<" "<<stats->S2[j][2]<<" "<<stats->S2[j][3]<<" "<<stats->S2[j][4]<<" ";
	// cout<<stats->S2[j][5]<<" "<<stats->S2[j][6]<<" "<<stats->S2[j][7]<<" "<<stats->S2[j][8]<<" "<<endl;
	// cout<<endl;
 
	int x = ((*dataVector)[i]).data.size();
	int y = ((*dataVector)[j]).data.size();
 
	double avg_x = (stats->S[i][u+l]-stats->S[i][u])/l;
	double avg_y;
	// cout<<stats->S2[i][u+l]-stats->S2[i][u]<<endl;
	// cout<<endl;
 
	double stdev_x = (stats->S2[i][u+l]-stats->S2[i][u]+0.0)/l - avg_x*avg_x;
	stdev_x = sqrt(stdev_x);
	double stdev_y;
 
	// cout<<"avg_x :"<<avg_x<<endl;
	// cout<<"stdev_x : "<<sqrt(stdev_x)<<endl;
	// cout<<endl;
 
	double multiple;
	double Cs;
	double dist;
	double mindist = 100000;
 
	for(int v=0; v<=y-l; v++){
		avg_y = (stats->S[j][v+l]-stats->S[j][v])/l;
		stdev_y= (stats->S2[j][v+l]-stats->S2[j][v])/l - avg_y*avg_y;
		stdev_y = sqrt(stdev_y);
		// cout<<"avg_y :"<<avg_y<<endl;
		// cout<<"stdev_y : "<<sqrt(stdev_y)<<endl;
 
		multiple = (stats->M[i][j][u+l][v+l]-stats->M[i][j][u][v]);
 
		Cs = (multiple-l*avg_x*avg_y)/(l*stdev_x*stdev_y);
		if(abs(1-Cs)<0.001)
			dist = 0;
		else
			dist = sqrt(2*(1-Cs));
		//cout<<dist<<endl;
 
		//cout<<endl;
 
		if(dist<mindist)
			mindist = dist;
	}
	//cout<<endl;
 
	return mindist;
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
 
 
double Gap(double distSum_left, double distSum_right, int num_left, int num_right){
	if(num_left == 0)
		return distSum_right/num_right;
	if(num_right == 0)
		return 100000;
 
	return distSum_right/num_right-distSum_left/num_left;
}
 
 
bool BestIG(Histogram* hist, double* best_t, double* max_gain, double* max_gap){
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
 
	//cout<<"a "<<a2<<" b "<<b2<<" dist "<<distSum_right<<endl;
	
	it = hist->data.begin();
 
	double totalEntropy = Entropy(a2, b2);
	double entropy;
	double gain;
	double dist_before = -1;
	double dist = -1;
 
	while(it!=hist->data.end()){
 
		dist_before = dist;
		dist = it->value;
 
		entropy = Entropy(a1, b1, a2, b2);
		gain = totalEntropy-entropy;
		//cout<<a1<<" "<<b1<<" "<<a2<<" "<<b2<<endl;
		//cout<<dist_before<<" "<<dist<<" "<<gain<<" "<<Gap(distSum_left, distSum_right, a1+b1, a2+b2)<<endl;
		//cout<<endl;
		// if gain increases
		if(gain>*max_gain){
 
			temp = true;
			*max_gain = gain;
			*max_gap = Gap(distSum_left, distSum_right, a1+b1, a2+b2);
			*best_t = (dist_before+dist)/2;
		}
 
		// if gain ties
		else if(abs(gain-*max_gain)<0.0001){
 
			double gap = Gap(distSum_left, distSum_right, a1+b1, a2+b2);
			if(gap>*max_gap){

				temp = true;
				*max_gap = gap;
				*best_t = (dist_before+dist)/2;
			}
		}
 
		// change values
		if(it->type==0){
			a1++;
			a2--;
		}
		else{
			b1++;
			b2--;
		}
 
 
		distSum_right = distSum_right - dist;
		distSum_left = distSum_left + dist;
 
		it++;
 
	}
 
	return temp;
}
 
void Shapelet_Discovery(vector<Data>* dataVector){
 	
	Stats stats(dataVector);
 
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
					//double distance = sdist(&candidate, &(*dataVector)[t].data);
					double distance = sdist_new(dataVector, &stats, k, t, i, j-i+1);
					hist.insert((*dataVector)[t].type, distance);
				}
				if(BestIG(&hist, &(Shapelet.best_t), &max_gain, &max_gap)){
					// hist.print();
					// cout<<i<<" "<<j<<" :  gain = "<<max_gain<<", gap = "<<max_gap<<endl;
				}
			}
			//cout<<i<<endl;
		}
		cout<<"Round "<<k<<" : t = "<<(Shapelet.best_t)<<", gain = "<<max_gain<<", gap = "<<max_gap<<endl;
		
	}
}
 
 
int main(){
	clock_t before;
 
	before = clock();
 
	vector <Data> dataVector;
 	ifstream inFile;
	inFile.open("gun_train_2", ios::in);
	string s;
 
	// while(getline(inFile, s)){
	// 	stringstream ss(s);
	// 	double test;
	// 	int t;
	// 	Data temp;
 
	// 	ss>>t;
	// 	temp.type = t;
		
	// 	while(ss>>test){
	// 		temp.data.push_back(test);
	// 	}
	// 	dataVector.push_back(temp);
	// }
 
	//FOR DEBUGGING - DATA SET SIZE = 10
	for(int i=0; i<10; i++){
		getline(inFile, s);
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
 
	// for(int i=0; i<dataVector.size(); i++){
	// 	cout<<i<<endl;;
	// 	printVector(dataVector[i].data);
	// 	cout<<endl;
	// }
 
 
	Shapelet_Discovery(&dataVector);
 

 	//GET TIME
	double result = (double)(clock()-before)/CLOCKS_PER_SEC;
	cout<<"Time = "<<result<<endl;

	//////////////////////sdist_new and dist TESTING///////////////////////////
 
	// Data data1;
	// Data data2;
 
	// data1.data.push_back(1);
	// data1.data.push_back(2);
	// data1.data.push_back(3);
 
	// data1.data.push_back(4);
	// data1.data.push_back(5);
	// data1.data.push_back(4);
	// data1.data.push_back(5);
 
	// data1.data.push_back(6);
	// data1.data.push_back(7);
	// data1.data.push_back(8);
 
	// data2.data.push_back(9);
	// data2.data.push_back(8);
	// data2.data.push_back(6);
	// data2.data.push_back(7);
	// data2.data.push_back(4);
	// data2.data.push_back(3);
	// data2.data.push_back(5);
	// data2.data.push_back(2);
	// data2.data.push_back(0);
	// data2.data.push_back(1);
 
 
	// Data data3;
	// data3.data.push_back(4);
	// data3.data.push_back(5);
	// data3.data.push_back(4);
	// data3.data.push_back(5);
 
	// vector<Data> dataVector2;
	// dataVector2.push_back(data1);
	// dataVector2.push_back(data2);
 
	// Stats stats2(&dataVector2);
 
	// cout<<sdist(&data3.data, &data2.data)<<endl;
	// cout<<endl;
	// cout<<sdist_new(&dataVector2, &stats2, 0, 1, 3, 4);
 
	/////////////////////////////////////////////////////////////////////////////
 
 
 
	// Histogram hist_temp;
	// hist_temp.insert(1, 1.321);
	// hist_temp.insert(1, 2.321);
	// hist_temp.insert(1, 2.564);
	// hist_temp.insert(1, 3.124);
	// hist_temp.insert(0, 5.123);
	// hist_temp.insert(1, 6.123);
	// hist_temp.insert(0, 6.431);
	// hist_temp.insert(0, 7.123);
	// hist_temp.insert(0, 8.231);
	// hist_temp.insert(0, 9.231);
 
	// hist_temp.print();
	// double best_t = 0;
	// double max_gain = 0;
	// double max_gap = 0;
	// BestIG(&hist_temp, &best_t, &max_gain, &max_gap);
 
	// cout<<endl;
 
	// cout<<best_t<<" "<<max_gain<<" "<<max_gap<<endl;
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