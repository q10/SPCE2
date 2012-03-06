#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <time.h>

float deltaf[20]={0};
float samples[20][550000]={0};
double weights[20][550000]={0};
int numberofsamples[20]={0};
int i,j,k,l,n=0;
float windowspacing, windowwidth, minimumh=0;
int nwindows=0;
int niterations=0;

void iteratedeltaf();
void readfiles();
void calculateweights();
double calcdenominator();
void createhistogram();
char histogramfile[128];
char datafile[128];

int main(int argc, char *argv[]){
	sscanf(argv[1], "%d", &nwindows);
	sscanf(argv[2], "%f", &windowspacing);
	sscanf(argv[3], "%f", &windowwidth);
	sscanf(argv[4], "%f", &minimumh);
	sscanf(argv[5], "%s", &histogramfile);
	sscanf(argv[6], "%d", &niterations);
	sscanf(argv[7], "%s", &datafile);
	readfiles();
	for(k=0; k<=niterations-1; k++){
		iteratedeltaf();
	}
	for(k=0; k<=nwindows-1; k++){
		printf("%f\n", deltaf[k]);

	}
	calculateweights();
	createhistogram();
	return 0;
	
}

void createhistogram(){
	FILE *histptr;
	double histogram[200]={0};
	float range, binspacing, average=0;
	int binindex=0;
	range = (nwindows-1)*windowspacing + windowwidth;
	binspacing = range/100;
	histptr=fopen(histogramfile, "w");
	for(i=0; i<=nwindows-1; i++){
		for(j=0; j<=numberofsamples[i]-1; j++){
			binindex= (samples[i][j]-minimumh)/binspacing;
			histogram[binindex]+=weights[i][j];
		}
	}
	for(i=0; i<=99; i++){
		fprintf(histptr, "%f %.12f %.12f\n", i*binspacing+minimumh, histogram[i], -log(histogram[i]));
	}
	fclose(histptr);
}

void readfiles(){
	char filename[128];
	char lineoffile[128];
	long int nlines=0;
	FILE *inputfile;
	for(i=0; i<=nwindows-1; i++){
		sprintf(filename, "%s%d.txt", &datafile, i+1);
		inputfile = fopen(filename, "r");
		nlines=0;
		while((fgets(lineoffile, sizeof(lineoffile), inputfile)) !=NULL){
			fscanf(inputfile, "%f\n", &samples[i][nlines]);
			numberofsamples[i]+=1;
			nlines+=1;
		}
		printf("%d\n", nlines);
	}
	
}


double calcdenominator(){
	double denominator=0;
	
	if(samples[j][n]< (minimumh + (j-1)*windowspacing + windowwidth) && j!=0){
		denominator = 1/(exp(deltaf[j])*numberofsamples[j]+exp(deltaf[j-1])*numberofsamples[j-1]);
	}
	else if(samples[j][n]> (minimumh + (j+1)*windowspacing)  && j!=nwindows-1){
		denominator = 1/(exp(deltaf[j])*numberofsamples[j] + exp(deltaf[j+1])*numberofsamples[j+1]);
	}
	else{
		denominator = 1/(exp(deltaf[j])*numberofsamples[j]);
	}
	return denominator;
	
}

void calculateweights(){
	for(j=0; j<=nwindows-1; j++){
			for(n=0; n<=numberofsamples[j]-1; n++){
				weights[j][n] = calcdenominator();
			}
		}
	
	
}

void iteratedeltaf(){
	double fenew=0;
	double sum=0;
	double deltafnew[20]={0};
	for(i=0; i<=nwindows-1; i++){
		sum=0;
		for(j=0; j<=nwindows-1; j++){
			for(n=0; n<=numberofsamples[j]-1; n++){
				if(samples[j][n]>(minimumh+i*windowspacing) && samples[j][n]<(minimumh+i*windowspacing+windowwidth)){
					sum += calcdenominator();
				}
			}
		}
		deltafnew[i] = -log(sum);
	}
	for(i=0; i<=nwindows-1; i++){
		deltaf[i] = deltafnew[i];
	}
	
}



