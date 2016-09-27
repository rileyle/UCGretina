/* 
*  class for coulex angular distributions a la Alder and Winther
*  w(th) = 1/4/pi*(a0 P0(cos(th)) + a2 P2(cos(th)) + a4 P4(cos(th))
*	 P(th) = W(th)*sin(th)
*/
#include "AngularDistribution.hh"
/***************************************************************/
AngularDistribution::AngularDistribution() {
	ai[0] = 1;
	ai[1] = 0;
	ai[2] = 0;
	//	FillPDF();
}
/***************************************************************/
AngularDistribution::AngularDistribution(double a0,double a2,double a4) {
	ai[0] = a0;
	ai[1] = a2;
	ai[2] = a4;
	//	FillPDF();
}
/***************************************************************/
AngularDistribution::~AngularDistribution() {
}
/***************************************************************/
double AngularDistribution::GetCoeff(int index) {
	return ai[index/2];
}
/***************************************************************/
void AngularDistribution::SetCoeff(int index,double value) {
  //	printf("inside SetCoeff(%i,%f)\n",index,value);
	ai[index/2] = value;
	//	FillPDF();
	return;
}
/***************************************************************/
void AngularDistribution::SetCoeffs(double a0, double a2, double a4) {
	ai[0] = a0;
	ai[1] = a2;
	ai[2] = a4;
	//	FillPDF();
	return;
}
/***************************************************************/
void AngularDistribution::FillPDF() {
	int i = 0;
	angdismax = 0;
	double stepsize = 3.14159/1440.0;
	for(i=0;i<1440;i++) {
		angdis[i] = PDF(i*stepsize);
		if(angdis[i]>angdismax) angdismax = angdis[i];
	//	printf("angdismax = %f\n",angdismax);
	}
	//	printf("built angular distribution with coefficients (a0,a2,a4) = (%f,%f,%f)\n",ai[0],ai[1],ai[2]);
	//printf("angdismax = %f\n",angdismax);
	return;
}
/***************************************************************/
double* AngularDistribution::GetAngDis() {
	return angdis;
}
/***************************************************************/
double AngularDistribution::PDF(double x) {
	return 1.0/4.0/3.14159*(ai[0] + 
		ai[1]*0.5*(3.0*pow(cos(x),2)-1) + 
		0.125*ai[2]*(35.0*pow(cos(x),4)-30.0*pow(cos(x),2)+3.0)) * sin(x);
}
/***************************************************************/
double AngularDistribution::GetRandomAngle() {
	//FILE* fp = fopen("internal_test.dat","w");
	//printf("in random angle\n");
	//	int i=0;
	//for(i=0;i<1440;i++)	fprintf(fp,"%i\t%f\n",i,angdis[i]);
	//fclose(fp);
	double r1   = 0;
	double rw   = 0;
	int    rbin = 0;
	double r2   = 0;
	do {
		r1          = rand();
		rw          = rand();
		r1          = (double)r1/RAND_MAX*1440;
		rbin        = (int)r1;
		rw          = (double)rw/RAND_MAX*angdismax;
		//printf("%i\t%f\n",rbin,rw);
	} while(angdis[rbin]<rw);
	r2 = rand();
	r2 = (double)r2/RAND_MAX;
	return (rbin+r2)*3.14159/1440; 
}
/***************************************************************/
void AngularDistribution::Report() {
  printf("built angular distribution with coefficients (a0,a2,a4) = (%f,%f,%f)\n",ai[0],ai[1],ai[2]);
}
/***************************************************************/
