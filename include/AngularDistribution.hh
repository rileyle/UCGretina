#ifndef ANGULARDISTRIBUTION_HH
#define ANGULARDISTRIBUTION_HH
/* class for coulex angular distribution */
#include <stdio.h>
#include <math.h>
#include <stdio.h>
#include <time.h>
#include <stdlib.h>

class AngularDistribution {
public:
	// constructor
	AngularDistribution();
	AngularDistribution(double,double,double);
	~AngularDistribution();	
	double GetCoeff(int);
	void SetCoeff(int,double);
	void FillPDF();
	double PDF(double);	
	double* GetAngDis();
	double GetRandomAngle();
private:
	double ai[3];
	double angdis[1440];
	double angdismax;
};

#endif
