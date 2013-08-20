//============================================================================
// Name        : SignificantHeight.cpp
// Author      : Konstantin Kuznetsov
// Version     :
// Copyright   : 
// Description : Anomalous waves, Ansi-style
//============================================================================

#include <iostream>
#include <ipp.h>
#include <fstream>
#include <cmath>
#include <time.h>
#include <stdio.h>
#include <stdlib.h>

using namespace std;


class Waves
{
private:
	Ipp32f* wavesHeight;
	unsigned int* wavesPeriod;
	unsigned int* wavesIndex;
	unsigned int wavesCount;
	unsigned int firstMinPeak;
	Ipp32f bdysh;
public:
	Waves(unsigned int waveCount = 1){
		wavesHeight = ippsMalloc_32f(waveCount);
		wavesIndex = new unsigned int[waveCount];
		wavesPeriod = new unsigned int[waveCount];
		wavesCount = waveCount;
	}
	~Waves(){
		ippsFree(wavesHeight);
		delete(wavesPeriod);
		delete(wavesIndex);
	}
	Waves(Ipp32f* x, unsigned int len, bool isCorrecting = false){
		IppStatus st;
		bdysh = 0;
		/**define average wave line*/
		Ipp32f tmp=0;
		st = ippsMean_32f(x,len,&tmp,ippAlgHintFast);
		st = ippsSubC_32f_I(tmp,x,len);
		/**define average wave line done*/
		/**define max and min peaks****/
		unsigned int maxPeaksCount = 0;
		unsigned int minPeaksCount = 0;
		unsigned int* maxPeaksId = new unsigned int[len/2];
		unsigned int* minPeaksId = new unsigned int[len/2];
		for(unsigned int i=0;i<len/2;i++)
		{
			maxPeaksId[i]=0;
			minPeaksId[i]=0;
		}
		//st = ippsSet_32u(0, maxPeaksId,len/2);
		//st = ippsSet_32u(0, minPeaksId,len/2);
		for(unsigned int i=0;i<len;i++)
		{
			/*if((x[i]>0)&&(maxPeaksId[minPeaksCount]==0))
				minPeaksId[minPeaksCount] = i;*/
			if((x[i]>=0)&&(maxPeaksId[maxPeaksCount]==0))
				maxPeaksId[maxPeaksCount] = i;
			else if((x[i]>=0)&&(x[i]>x[maxPeaksId[maxPeaksCount]]))
				maxPeaksId[maxPeaksCount] = i;
			else if((x[i]<0)&&(maxPeaksId[maxPeaksCount]!=0))
				maxPeaksCount++;
		}
		for(unsigned int i=0;i<len;i++)
		{
			//cout<<"x["<<i<<"]="<<x[i]<<"meanPeaksId["<<minPeaksCount<<"]="<<minPeaksId[minPeaksCount]<<endl;
			if((x[i]<0)&&(minPeaksId[minPeaksCount]==0))
				minPeaksId[minPeaksCount] = i;
			else if((x[i]<0)&&(x[i]<x[minPeaksId[minPeaksCount]]))
				minPeaksId[minPeaksCount] = i;
			else if((x[i]>=0)&&(minPeaksId[minPeaksCount]!=0))
				minPeaksCount++;
		}
		/**define max and min peaks done****/
		unsigned int waveCount=maxPeaksCount;
		/*if(maxPeaksCount<minPeaksCount)
			waveCount = minPeaksCount;*/
		//cout<<"maxPeaksCount="<<maxPeaksCount<<" minPeaksCount="<<minPeaksCount<<endl;
		wavesHeight = ippsMalloc_32f(waveCount);
		wavesIndex = new unsigned int[waveCount];
		wavesPeriod = new unsigned int[waveCount];
		wavesCount = waveCount;
		int wvPeriodTmp=0;

		if((minPeaksId[0]<maxPeaksId[0])&&(minPeaksId[0]!=0)){
			//cout<<"minPeaksId "<<waveCount<<" ";
			//cout<<minPeaksId[0]<<endl;
			for(unsigned int i=0;i<waveCount; i++){

				//wavesHeight[i] = max(abs(x[maxPeaksId[i]] - x[minPeaksId[i]]),abs(x[maxPeaksId[i]] - x[minPeaksId[i+1]]));
				wavesHeight[i] = x[maxPeaksId[i]] - x[minPeaksId[i]];
				//wavesHeight[i] = x[maxPeaksId[i]] - x[minPeaksId[i]];
				//wavesHeight[i] = 2*x[maxPeaksId[i]];
				//wavesHeight[i] = (abs(x[maxPeaksId[i]] - x[minPeaksId[i]])+abs(x[maxPeaksId[i]] - x[minPeaksId[i+1]]))/2;
				//wavesHeight[i] = abs(x[maxPeaksId[i]] - x[minPeaksId[i]]);
				wvPeriodTmp = (maxPeaksId[i+1] - minPeaksId[i]);
				if(wvPeriodTmp<-500)
					wvPeriodTmp=wavesPeriod[i-1];
				//wvPeriodTmp = 2*(maxPeaksId[i] - minPeaksId[i]);
				//if(wvPeriodTmp<0)	wvPeriodTmp=(-1)*wvPeriodTmp;
				wavesPeriod[i] = wvPeriodTmp;
				wavesIndex[i] = maxPeaksId[i];
				/*cout<<"x[maxPeaksId["<<i<<"]]="<<x[maxPeaksId[i]]<<"\t";
				cout<<"x[minPeaksId["<<i<<"]]="<<x[minPeaksId[i]]<<"\t";
				cout<<"maxPeaksId["<<i<<"]="<<maxPeaksId[i]<<"\t";
				cout<<"minPeaksId["<<i<<"]="<<minPeaksId[i]<<"\t";
				cout<<"waveHeight["<<i<<"]="<<wavesHeight[i]<<endl;*/
			}
		}
		/*else if((minPeaksId[0]<maxPeaksId[0])&&(minPeaksId[0]==0)){
			cout<<"minPeakId=0"<<waveCount<<endl;
			wavesHeight[0] = x[maxPeaksId[1]]-bdysh;//iz predydushego
			wavesPeriod[0] = maxPeaksId[1]*2;//primerno
			wavesIndex[0] = maxPeaksId[1];
			for(unsigned int i=1;i<waveCount; i++){
				//wavesHeight[i-1] = max(abs(x[maxPeaksId[i]] - x[minPeaksId[i-1]]),abs(x[maxPeaksId[i]] - x[minPeaksId[i]]));
				//wavesHeight[i-1] = x[maxPeaksId[i]] - x[minPeaksId[i-1]];
				wavesHeight[i] = x[maxPeaksId[i]] - x[minPeaksId[i-1]];
				//wavesHeight[i-1] = (abs()+abs(x[maxPeaksId[i]] - x[minPeaksId[i]]))/2;
				//wavesHeight[i-1] = 2*x[maxPeaksId[i]];
				//wvPeriodTmp = minPeaksId[i+1] - minPeaksId[i];
				wvPeriodTmp = 2*(maxPeaksId[i] - minPeaksId[i]);
				//if(wvPeriodTmp<0)	wvPeriodTmp=(-1)*wvPeriodTmp;
				wavesPeriod[i] = wvPeriodTmp;
				wavesIndex[i] = maxPeaksId[i];
				/*cout<<"x[maxPeaksId["<<i<<"]]="<<x[maxPeaksId[i]]<<"\t";
				cout<<"x[minPeaksId["<<i<<"]]="<<x[minPeaksId[i]]<<"\t";
				cout<<"maxPeaksId["<<i<<"]="<<maxPeaksId[i]<<"\t";
				cout<<"minPeaksId["<<i<<"]="<<minPeaksId[i]<<"\t";
				cout<<"waveHeight["<<i<<"]="<<wavesHeight[i]<<endl;
			}
		}*/
		else if((minPeaksId[0]>maxPeaksId[0]))//1 and 2
		{
			//cout<<"maxPeaksId "<<waveCount<<" ";
			//cout<<maxPeaksId[0]<<endl;
			wavesHeight[0] = 2*x[maxPeaksId[0]]-bdysh;//iz predydushego
			wavesPeriod[0] = maxPeaksId[0]*2;//primerno
			wavesIndex[0] = maxPeaksId[0];
			for(unsigned int i=1;i<waveCount; i++){
				//wavesHeight[i-1] = max(abs(x[maxPeaksId[i]] - x[minPeaksId[i-1]]),abs(x[maxPeaksId[i]] - x[minPeaksId[i]]));
				//wavesHeight[i-1] = x[maxPeaksId[i]] - x[minPeaksId[i-1]];
				wavesHeight[i] = x[maxPeaksId[i]] - x[minPeaksId[i-1]];
				//wavesHeight[i-1] = (abs()+abs(x[maxPeaksId[i]] - x[minPeaksId[i]]))/2;
				//wavesHeight[i-1] = 2*x[maxPeaksId[i]];
				wvPeriodTmp = maxPeaksId[i+1] - maxPeaksId[i];
				if(wvPeriodTmp<-500)
					wvPeriodTmp=wavesPeriod[i-1];
				//wvPeriodTmp = 2*(maxPeaksId[i] - minPeaksId[i-1]);
				//if(wvPeriodTmp<0)	wvPeriodTmp=(-1)*wvPeriodTmp;
				wavesPeriod[i] = wvPeriodTmp;
				wavesIndex[i] = maxPeaksId[i];
				/*cout<<"x[maxPeaksId["<<i<<"]]="<<x[maxPeaksId[i]]<<"\t";
				cout<<"x[minPeaksId["<<i<<"]]="<<x[minPeaksId[i]]<<"\t";
				cout<<"maxPeaksId["<<i<<"]="<<maxPeaksId[i]<<"\t";
				cout<<"minPeaksId["<<i<<"]="<<minPeaksId[i]<<"\t";
				cout<<"waveHeight["<<i<<"]="<<wavesHeight[i]<<endl;*/
			}
		}
		else{
			cout<<"Errrr"<<endl;
			for(unsigned int i=1;i<waveCount; i++){
				//wavesHeight[i-1] = max(abs(x[maxPeaksId[i]] - x[minPeaksId[i-1]]),abs(x[maxPeaksId[i]] - x[minPeaksId[i]]));
				//wavesHeight[i-1] = x[maxPeaksId[i]] - x[minPeaksId[i-1]];
				wavesHeight[i-1] = x[maxPeaksId[i]] - x[minPeaksId[i-1]];
				//wavesHeight[i-1] = (abs()+abs(x[maxPeaksId[i]] - x[minPeaksId[i]]))/2;
				//wavesHeight[i-1] = 2*x[maxPeaksId[i]];
				//wvPeriodTmp = minPeaksId[i+1] - minPeaksId[i];
				wvPeriodTmp = 2*(maxPeaksId[i] - minPeaksId[i]);
				//if(wvPeriodTmp<0)	wvPeriodTmp=(-1)*wvPeriodTmp;
				wavesPeriod[i-1] = wvPeriodTmp;
				wavesIndex[i-1] = maxPeaksId[i];
				/*cout<<"x[maxPeaksId["<<i<<"]]="<<x[maxPeaksId[i]]<<"\t";
				cout<<"x[minPeaksId["<<i<<"]]="<<x[minPeaksId[i]]<<"\t";
				cout<<"maxPeaksId["<<i<<"]="<<maxPeaksId[i]<<"\t";
				cout<<"minPeaksId["<<i<<"]="<<minPeaksId[i]<<"\t";
				cout<<"waveHeight["<<i<<"]="<<wavesHeight[i]<<endl;*/
			}
		}
		/*for(unsigned int i=0;i<waveCount; i++){
			wavesHeight[i] = max(abs(x[maxPeaksId[i]] - x[minPeaksId[i]]),abs(x[maxPeaksId[i]] - x[minPeaksId[i+1]]));

			wvPeriodTmp = maxPeaksId[i] - minPeaksId[i];
			//if(wvPeriodTmp<0)	wvPeriodTmp=(-1)*wvPeriodTmp;

			wavesPeriod[i] = wvPeriodTmp;
			wavesIndex[i] = maxPeaksId[i];
			/*cout<<"x[maxPeaksId["<<i<<"]]="<<x[maxPeaksId[i]]<<"\t";
			cout<<"x[minPeaksId["<<i<<"]]="<<x[minPeaksId[i]]<<"\t";
			cout<<"maxPeaksId["<<i<<"]="<<maxPeaksId[i]<<"\t";
			cout<<"minPeaksId["<<i<<"]="<<minPeaksId[i]<<"\t";
			cout<<"waveHeight["<<i<<"]="<<wavesHeight[i]<<endl;
		}*/
		delete(maxPeaksId);
		delete(minPeaksId);
		if(isCorrecting == true){

		}
	}
	unsigned int getLastMinPeak(){
		return 0;
	}
	double getSignHeight(){
		Ipp32f signHeight;
		Ipp32f* waveHeightSort = ippsMalloc_32f(wavesCount);
		ippsCopy_32f(wavesHeight, waveHeightSort, wavesCount);
		ippsSortDescend_32f_I(waveHeightSort,wavesCount);
		ippsMean_32f(waveHeightSort,wavesCount/3,&signHeight,ippAlgHintAccurate);
		ippsFree(waveHeightSort);
		return (double)signHeight;
	}
	double getMeanHeight(){
		Ipp32f meanHeight;
		ippsMean_32f(wavesHeight, wavesCount, &meanHeight,ippAlgHintAccurate);
		return (double)meanHeight;
	}
	double getMeanPeriod(){
		//Ipp32f meanPeriod=0;
		unsigned int sumT=0;
		//ippsMean_32f((Ipp32f*)wavesPeriod, wavesCount, &meanPeriod,ippAlgHintAccurate);
		for(unsigned int i=0;i<wavesCount;i++){
			sumT+=wavesPeriod[i];
		}
		return sumT/wavesCount;
	}
	void getIndexesOfAnomWaves(unsigned int* indexesAnom, double* coeffMass, unsigned int& lenAnomWaves, const double coeff=2){
		Ipp32f signHeight = getSignHeight();
		/*for(unsigned int i = 0; i<wavesCount; i++){
			if(wavesHeight[i]>(coeff*signHeight))
				lenAnomWaves++;
		}*/
		lenAnomWaves = 0;
		for(unsigned int i = 0; i<wavesCount; i++){
			if(wavesHeight[i]>(coeff*signHeight)){
				indexesAnom[lenAnomWaves] = wavesIndex[i];
				coeffMass[lenAnomWaves] = wavesHeight[i]/signHeight;
				lenAnomWaves++;
			}
		}
	}
	unsigned int getWaveNum(){
		return wavesCount;
	}
	void printHeights(){
		for(unsigned int i=0; i<wavesCount; i++){
			cout<<"wavesHeight["<<i<<"] = "<<wavesHeight[i]<<endl;
		}
	}
	void printIndexes(){
		for(unsigned int i=0; i<wavesCount; i++){
			cout<<"wavesIndexs["<<i<<"] = "<<wavesIndex[i]<<endl;
		}
	}
	void printPeriods(){
		for(unsigned int i=0; i<wavesCount; i++){
			cout<<"wavesPeriod["<<i<<"] = "<<wavesPeriod[i]<<endl;
		}
	}
	Ipp32f getHeight(unsigned int ind){
		return wavesHeight[ind];
	}
	unsigned int getPeriod(unsigned int ind){
		return wavesPeriod[ind];
	}
	unsigned int getIndex(unsigned int ind){
			return wavesIndex[ind];
	}
};

int main(int argc, char* argv[]) {
	if(argc!=4)
	{
		cout<<"Wrong input parameters!\n";
		cout<<"Usage: ./waveHeights infile.txt OutHeights SingHeights"<<endl;
		cout<<"file OutHeights consists of:"<<endl;
		cout<<"#1waveHeights #2waveIndex #3wavePeriods #4waveHeights/Hs #5waveHeights/Hmean #6waveHeights/Tmean"<<endl;
		return 1;
	}
	float datain;
	//timeMeter tm;//time checker, shows how mush time executes some code blocks
	//const char* filename="/home/konst/ProjCPP/SignificantHeight/sinTest.txt";
	//const char* filename="/media/KONST_500GB/anomWaves/Vzmorie#24_2007-07-14_1sec.txt";
	const char* filename=argv[1];
	const char* filenameOutHeights=argv[2];
	const char* filenameSingHeights=argv[3];
	//const char* filename="c:/Work/anomWaves/Vzmorie/modelRowRayl.txt";
	//const char* filenameOutHeights="c:/Work/anomWaves/Vzmorie/modelRowOutHeights.txt";
	//const char* filenameSingHeights="c:/Work/anomWaves/Vzmorie/modelRowSignHeights.txt";
	FILE*in=NULL;
	FILE*outHeight=NULL;
	FILE*outSignHeight=NULL;
	in=fopen(filename,"r");
	if(in==NULL)
	{
		printf("no input files\n");
		return 1;
	}
	unsigned int cnt=0;
	cout<<"scaning..."<<endl;
	while(fscanf(in,"%f",&datain)!=EOF)
	{
		cnt++;
	}
	fclose(in);
	unsigned int len=cnt; // length of input vector
	/****************************************/
	Ipp32f* ost=ippsMalloc_32f(len);//input vector
	/****loading file******/
	in=fopen(filename,"r");
	cnt=0;
	cout<<"reading..."<<endl;
	while(fscanf(in,"%f",&datain)!=EOF)
	{
		ost[cnt]=datain;
		cnt++;
	}
	fclose(in);
	cout<<len<<"  values read"<<endl;
	/****loading done******/
	unsigned int otrLenght = 20*60*1;
	unsigned int start, end = 0;
	Ipp32f* rowTmp = ippsMalloc_32f(otrLenght);
	/*unsigned int* anomWaves = new unsigned int[len/10];
	unsigned int lenAnomWaves = 0;*/
	outHeight=fopen(filenameOutHeights,"w");
	outSignHeight=fopen(filenameSingHeights,"w");
	if(outHeight==NULL)
	{
		printf("Could not open file %s for write\n",filenameOutHeights);
		return 1;
	}
	/*Waves wv(ost,len);
	*/
	//unsigned int* anomWavesTMP = new unsigned int[otrLenght/10];
	//unsigned int lenAnomWavesTMP = 0;
	unsigned int wavesNum = 0;
	double Hs=0;
	double Hmean=0;
	double Tmean=1;

	//double* coeffMass = new double[otrLenght/10];
	for(unsigned int i=0;i<len/otrLenght;i++){
		start = i*otrLenght; end = start+otrLenght;
		for(unsigned int j=start;j<end;j++){
			rowTmp[j-start] = ost[j];
		}
		Waves* wv = new Waves(rowTmp,otrLenght);
		//printf("Sizeof waves wv = %d\n",sizeof(&wv));
		//wv.printIndexes();
		//wv.printPeriods();
		//cout<<"wave num"<<wv.getWaveNum()<<endl;
		wavesNum = wv->getWaveNum();
		Hs = wv->getSignHeight();
		Hmean = wv->getMeanHeight();
		Tmean = wv->getMeanPeriod();
		fprintf(outSignHeight,"%8.12f\n",Hs);
		for(unsigned int i=0; i<wavesNum; i++){
			fprintf(outHeight,"%8.12f\t%d\t%d\t%8.12f\t%8.12f\t%8.12f\n",wv->getHeight(i),start+1+wv->getIndex(i),wv->getPeriod(i),wv->getHeight(i)/Hs,wv->getHeight(i)/Hmean,wv->getPeriod(i)/Tmean);
		}
		delete wv;
		//wv.getIndexesOfAnomWaves(anomWavesTMP, coeffMass, lenAnomWavesTMP,atof(argv[4]));
		//for(unsigned int k=0; k<lenAnomWavesTMP; k++){
		//	fprintf(outAnom,"%d\t%d\t%f\n",i,start+1+anomWavesTMP[k], coeffMass[k]);
		//}
	}
	fclose(outHeight);
	fclose(outSignHeight);
	ippsFree(ost);
	ippsFree(rowTmp);
	/*Waves wv(ost,len);
	unsigned int* anomWaves = new unsigned int[wv.getWaveNum()/5];
	unsigned int lenAnom=0;
	//signHeight(ost,swh,len);
	cout<<"Significat wave heights "<<wv.getSignHeight()<<endl;
	wv.getIndexesOfAnomWaves(anomWaves,lenAnom,1);
	for(unsigned int i=0; i<lenAnom; i++){
		cout<<"anomWaves "<<anomWaves[i]<<endl;
	}
	wv.printHeights();
	wv.printIndexes();
	wv.printPeriods();*/
	return 0;
}
