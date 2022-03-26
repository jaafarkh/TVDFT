// TVDFT.cpp : Defines the entry point for the console application.
//demonstrate the application of Time Variant DFT with non-orthogonality compensation
//generate a simulated signal with close and cross orders and analyze it twice
//with non-orthog. compensation enabled and without it

#include "stdafx.h"
#include <iostream>
#include <chrono>
#include <cstdlib>
#include <iomanip>
#include <cmath>
#include <ctime>
#include <windows.h>
#include <sstream>
#define DbgMsg( s )            \
{                             \
std::wostringstream os_;    \
   os_ << s;                   \
   OutputDebugStringW( os_.str().c_str() );  \
}

const double pi = 3.1415926535897932384626433832795;
void Gauss(double **A, double B[], int m, double X[])
{
	//**A, B[] and X[] are 1-based
	double max = 0;
	double a1 = 0;
	int i = 0;
	int j = 0;
	int k = 0;
	int Imx = 0;
	//pivoting

	for (j = 1; j < m; j++)
	{
		//pivoting
		max = abs(A[j][j]);
		Imx = j;
		for (i = j; i < m; i++)
		{
			if (abs(A[i + 1][j]) > max)
			{
				max = abs(A[i + 1][j]);
				Imx = i + 1;
			}
		}
		if (Imx != j)
		{
			for (k = 1; k <= m; k++)
			{
				a1 = A[j][k];
				A[j][k] = A[Imx][k];
				A[Imx][k] = a1;
			}
			a1 = B[j];
			B[j] = B[Imx];
			B[Imx] = a1;
		}
		//elimination
		for (i = j + 1; i <= m; i++)
		{
			a1 = -A[i][j] / A[j][j];
			//If Math.Abs(A(j, j)) < 1.0E-20 Then Debug.Print("pivoting error")
			for (k = 1; k <= m; k++)
			{
				A[i][k] = A[i][k] + A[j][k] * a1;
			}
			B[i] = B[i] + B[j] * a1;
		}
		//Debug.Print(Str(A(m, m - 1)) + "  " + Str(A(m, m)))
	}

	//back substitution
	X[m] = B[m] / A[m][m];

	for (j = m - 1; j >= 1; j--)
	{
		for (i = j + 1; i <= m; i++)
		{
			B[j] = B[j] - X[i] * A[j][i];
		}
		X[j] = B[j] / A[j][j];
	}

}
//==============================================TVDFT=====================================================
void TVDFTOCM(float Ai[], int NP, float Ph[], int NOrder, float OA[], float OB[], bool OCMEnabled)
{
	double z;
	int k, i, m, i1, m1;
	int NO2 = NOrder * 2;
	double* g = new double[NOrder + 1];
	double** A = new double*[NO2 + 1];
	double* B = new double[NO2 + 1];
	double* AC = new double[NO2 + 1];
	double* Tf = new double[NO2 + 1];
	z = 1000.0;
	for (i = 1; i <= NO2; i++)
	{
		B[i] = 0;
		A[i] = new double[NO2 + 1];  ////NO2+1 columns: here we define the elements of each row
		for (k = 1; k <= NO2; k++)
		{
			A[i][k] = 0;
		}
	}

	for (k = 1; k <= NP; k++)
	{
		for (i = 0; i < NOrder; i++)
		{
			g[i] = Ph[i*(NP + 1) + k];  //convert to double for better accuracy!
		}

		for (i = 0; i < NOrder; i++)
		{
			Tf[i * 2 + 1] = cos(g[i]);
			Tf[i * 2 + 2] = sin(g[i]);
		}
		for (i = 0; i < NOrder; i++)
		{
			B[i * 2 + 1] = B[i * 2 + 1] + Ai[k] * z * cos(g[i]);
			B[i * 2 + 2] = B[i * 2 + 2] + Ai[k] * z * sin(g[i]);
		}



		for (i = 1; i <= NO2; i++)
		{
			if (i % 2 == 0)
			{
				i1 = i - 1;
			}
			else
			{
				i1 = i + 1;
			}
			for (m = i; m <= NO2; m++)
			{
				if (m % 2 == 0)
				{
					m1 = m - 1;
				}
				else
				{
					m1 = m + 1;
				}
				A[i][m] = A[i][m] + Tf[i] * z * Tf[m];

			}
		}
	}
	for (i = 2; i <= NO2; i++)
	{
		for (m = 1; m < i; m++)
		{
			A[i][m] = A[m][i];
		}
	}
	if (OCMEnabled)
	{
		Gauss(A, B, NO2, AC);
	}
	for (i = 0; i < NOrder; i++)
	{
		if (OCMEnabled)
		{
			OA[i] = (float) AC[2 * i + 1];
			OB[i] = (float) AC[2 * i + 2];
		}
		else
		{
			OA[i] = (float)(2.0 * B[2 * i + 1] / NP / z);
			OB[i] = (float)(2.0 * B[2 * i + 2] / NP / z);
		}

	}
	for (i = 1; i <= NO2; i++) delete[] A[i];
	delete[] A;
	delete[] B;
	delete[] g;
	delete[] AC;
	delete[] Tf;

}
void GetRealImag(float *Ai,float *Psi,float *OA, float *OB, int NP, float Speed, int StP, int NofTraces, 
		         float *TraceOrder, unsigned char *TraceType, bool OCMEnabled, float SFactor, float SamRate)
{
	float SumA, Sc, DT;
	int j, i,k;
	float *Ph = new float[NofTraces*(NP + 1)];	
	Sc = SFactor; //a Scaling factor
	DT = 1.0f / SamRate;
	Speed = Speed / 60.0f;	//convert to Hz
	Speed = 2.0f * (float)pi * Speed;	//convert to rad/sec
	
	for (j = 0; j < (NofTraces-1); j++) {
		k = TraceType[j];
		//if(j == 0) std::cout << "type, order: " << k <<","<<TraceOrder[j] <<  " \n";
		for (i = 1; i <= NP; i++) {
			switch (k) {
			case 0:				//order
				Ph[j*(NP + 1) + i] = TraceOrder[j] * Psi[i - 1];
				break;
			case 1:				// Hz
				Ph[j*(NP + 1) + i] = 2 * (float)pi * TraceOrder[j] * (StP + i - 1) * DT;	//fixed freq
				break;
			case 2:		//RPM
				Ph[j*(NP + 1) + i] = 2 * (float)pi * (TraceOrder[j] / 60) * (StP + i - 1) * DT;		// fixed freq
				break;
			}

		}

	}
	TVDFTOCM(Ai, NP, Ph, NofTraces - 1, OA, OB, OCMEnabled); //orders
	for (i = 0; i < NofTraces-1; i++) {
		OA[i] = Sc * OA[i];
		OB[i] = Sc * OB[i];
	}

	SumA = 0;
	for (j = 1; j <= NP; j++) {
		SumA = SumA + Ai[j] * Ai[j];
	}
	SumA = SumA / NP;
	SumA = sqrt(SumA);
	OA[NofTraces-1] = 1.414f * Sc * SumA; //overall vibration for the last Trace
	OB[NofTraces-1] = 0.0f;
	delete[] Ph;
}

float Phase(double R, double Im)
{
	//returns an angle from 0 to 2*pi

	return (float)atan2(Im, R);
}


void Analyse(float *x, float *Trig, int RecLen,int NofTraces,double **AA, double **Ph, int &NofReadings, float *TraceOrder, unsigned char *TraceType,
			 float SamRate, float PPR, unsigned char SpacingMode, float TimSpac, float RevSpac, bool OCMEnabled)
{
	/*
	x: input vibration signal
	Trig: input trigger signal
	RecLen: total number of points (length of x and Trig)
	NofTraces: number of traces including overall trace
	AA: amplitude data
	Ph: phase data
	NofReadings: total number of blocks
	TraceOrder: trace order
	TraceType: tracing type; 0: Order, 1: fixed freq in Hz, 2: fixed freq in RPM, 3: overall
	SamRate: sample rate
	PPR: pulse per revoultion (normally = 1)
	SpacingMode: spacing mode; 0: revoltions, 1:time
	TimSpac: time spacing in sec
	RevSpac: revolution spacing
	bool OCMEnabled: if true then enable Non-Orthogonality compensation

	*/
	float *OA = new float[NofTraces + 1]; //cosine coefficients
	float *OB = new float[NofTraces + 1]; //sine coefficients
	float *Ai = new float[10241]; //temporary processing array holding one block data
	float *Psi = new float[10241]; //radian angular displacement temporary array
	int i, Count1, k, kk, NAvgS,j;
	float  vibSam1, vibSam2;
	float spdSam1, spdSam2;
	bool FTrigDetected, DoAnalysis, FirstSpeedDetected;
	float Speed, LSpeed, HSpeed, Rev;
	float Psi2 = 0.0;
	float SFactor = 1.0f; //scale factor

	float DeltaT = 1 / SamRate;
	float t1, t2;
	float a1, AvSpd, ang1,w,w1,w2,tt;
	int CC, CPuls;
	CC = 0;
	LSpeed = 1; //lower speed
	HSpeed = 20000; //higher speed	
	Count1 = 0;
	Rev = 0;
	NAvgS = 0;
	j = 0;
	//CaptureSam = false;
	FTrigDetected = false;
	FirstSpeedDetected = false;

	NofReadings = 0;
	float Lev1 = 0.5; //trigger low state plus a threshold, assuming Trig is changing from 0 to 5
	float Lev2 = 4.5; //trigger high state minus a threshold
	vibSam1 = x[0];
	spdSam1 = Trig[0];
	AvSpd = 0;
	for (i = 1; i < (RecLen - 1); i++) {
		//
		//b1 = binReader.ReadByte : b2 = binReader.ReadByte
		
		 vibSam2 = x[i];
		 spdSam2 = Trig[i];

		if ((spdSam1 < Lev1) && (spdSam2 > Lev2)) { //assuming positive edge trigger
			if (!FTrigDetected) {
				//CaptureSam = true;
				FTrigDetected = true;
				t1 = 0;
				t2 = 0;
				Count1 = 0;
				CC = 0;
				j = 0;
				CC = 0;
				//NAV = 0;
				Psi[0] = 0;
				Psi2 = 0;
				CPuls = 0;
			}
			else {
				//FTrigDetected = False
				CPuls += 1;


				t2 = 1.0f * Count1;
				if (CPuls > 1) {		//linear speed variation between pulses
					w1 = (2.0f * (float)pi / PPR) * SamRate / t1; //assuming linear speed variation, calculate the slope
					w2 = (2.0f * (float)pi / PPR) * SamRate / t2;
					a1 = (w2 - w1) / t2; //dw / dt

				}
				else {	//constant speed between cycles
					w1 = (2.0f * (float)pi / PPR) * SamRate / t2;
					a1 = 0;
				}
				tt = 0;
				for (j = 1; j <= Count1; j++) {
					w = w1 + a1*j; //instantaneous speed
					tt = tt + w / SamRate; //to ensure complete 2pi per cycle
				}
				for (kk = 1; kk <= Count1; kk++) {
					w = w1 + a1 * kk; //instantaneous speed
					Psi[kk + CC - Count1] = Psi2 + (w / SamRate) * (2.0f * (float)pi/PPR) / tt;
					Psi2 = Psi[kk + CC - Count1];
				}
				//last angle
				t1 = t2;
				Rev += 1 / PPR;
				Speed = 60 * SamRate / Count1 / PPR;
				//qDebug()<<"current cycle Speed is:"<<Speed;
				NAvgS += 1;
				AvSpd = AvSpd + Speed;
				j = Count1;
				Count1 = 0;

				if ((Speed >= LSpeed) && (Speed <= HSpeed)) {
					if (!FirstSpeedDetected) {
						//FirstSpeed = Speed;
						FirstSpeedDetected = true;
					}
					DoAnalysis = false;
					switch (SpacingMode) {
					case 0:
						//rev spacing
						if (Rev >= RevSpac)
							DoAnalysis = true;
						break;
					case 1:
						//time spacing
						if ((CC * DeltaT) >= TimSpac)
							DoAnalysis = true;
						break;
					case 2: //both whichever true
						if (Rev >= RevSpac)
							DoAnalysis = true;
						if ((CC * DeltaT) >= TimSpac)
							DoAnalysis = true;
						break;

					}
					if (DoAnalysis) {						
						
						GetRealImag(Ai, Psi, OA, OB, CC, Speed, i - CC, NofTraces, TraceOrder, TraceType, OCMEnabled, SFactor, SamRate);
						for (k = 0; k < NofTraces; k++) {							
							AA[k][NofReadings] = sqrt(OA[k] * OA[k] + OB[k] * OB[k]);							
							ang1 = Phase(OA[k], OB[k]);
							if (ang1 > pi) ang1 = ang1 - 2 * (float) pi; // -Pi < ph < Pi
							Ph[k][NofReadings] = ang1;
						}
						AvSpd = AvSpd / NAvgS;
						NofReadings += 1;
						
						//CaptureSam = False
						FirstSpeedDetected = false;
						CC = 0;
						Rev = 0;
						NAvgS = 0;
						AvSpd = 0;
					}
				}
				else {
					j = 0;
					//CaptureSam = False
					FirstSpeedDetected = false;
					// Rev = 0
				}
			}
		}
		Count1 += 1;
		CC += 1;
		if (CC > 200000) {
			//msgbox("Too long processing length\n please reduce #Revolutions!", QMessageBox::Ok);
			//rv = false;
			break;
		}

		Ai[CC] = vibSam1;
		spdSam1 = spdSam2;
		vibSam1 = vibSam2;

	}
	delete[] OA;
	delete[] OB;
	delete[] Ai;
	delete[] Psi;

}
int main()
{
	int NofTraces = 5; //toral number of traces including overall
	const int maxReadings = 4000; //expected max number of blocks
	int RecLen = 8192; //signal length
	double **AA = new double* [NofTraces];//amp data
	double **Ph = new double*[NofTraces];//phase data
	for (int j = 0; j < NofTraces; j++) {
		AA[j] = new double[maxReadings];
		Ph[j] = new double[maxReadings];
	}
	
	
	float* x = new float[RecLen];
	float* Trig = new float[RecLen];

	float *TraceOrder = new float[NofTraces + 1];
	
	unsigned char *TraceType = new unsigned char[NofTraces + 1];
	TraceOrder[0] = 1.0f; //1X
	TraceOrder[1] = 1.2f; //1.2X
	TraceOrder[2] = 2.0f; //2X
	TraceOrder[3] = 30.0f; //30 Hz

	TraceType[0] = 0;  //order of shaft speed
	TraceType[1] = 0;  //order of shaft speed 
	TraceType[2] = 0;  //order of shaft speed
	TraceType[3] = 1; // fixed freq in Hz
	TraceType[4] = 3; //overall, should be the last one

	//generate the signal
	float f, SampleRate,dt,w,phi,ph1,ph2,ph3,ph4,a1,a2,a3,a4,phi0;
	int i;
	SampleRate = 1024;
	dt = 1 / SampleRate;
	phi = 0;//shaft angular displacement
	ph1 = 0;
	ph2 = 0;
	ph3 = 0;
	ph4 = 0;
	phi0 = 0;
	for (i = 0; i < RecLen; i++) {
		f = 10 + 4 * i*dt;
		w = f * 2.0f * (float)pi;
		phi = phi + w * dt;

		ph1 = ph1+ TraceOrder[0] * w * dt;
		ph2 = ph2 + TraceOrder[1] * w * dt;
		ph3 = ph3 + TraceOrder[2] * w * dt;
		ph4 = ph4 + TraceOrder[3] * 2 * (float)pi * dt;
		a1 = 3.0f + (3.0f * i / RecLen);
		a1 = a1 * sin(ph1);
		a2 = (2.5f - (1.0f * i / RecLen)) * sin(ph2);
		a3 = (1.5f + (0.5f * i / RecLen)) * sin(ph3);
		a4 = (2.6f) *sin(ph4);
		x[i] = a1 + a2 + a3 + a4;
	/*	if ((sin(phi0) <= 0) && (sin(phi) > 0)) {
			Trig[i] = 5.0f;
		}
		else {
			Trig[i] = 0.0f;
		}*/
		if (sin(phi)  > 0) {
			Trig[i] = 5.0f;
		}
		else {
			Trig[i] = 0.0f;
		}
		phi0 = phi;
	}
	
	int NReads;
	float RevSpc = 10.0f;
	float TimSpc = 0.3f; 
	unsigned char SpacMode = 0; //0: const number of Rev mode, 1: constant delt-time

	//now will analyze with Non-orthogonality compensation enabled
	Analyse(x, Trig, RecLen, NofTraces, AA, Ph, NReads, TraceOrder, TraceType, SampleRate, 1.0f, SpacMode, TimSpc, RevSpc,true);

	std::cout << "No. of Readings: " << NReads << "\n";
	std::cout << "OCM enabled \n";
	for (i = 0; i < NReads; i++) {
		std::cout << AA[0][i] << " "<< AA[1][i] <<" "<< AA[2][i] <<" "<< AA[3][i] <<" "<< AA[4][i] <<";\n";
	}

	//now will analyze with Non-orthogonality compensation turned-off
	Analyse(x, Trig, RecLen, NofTraces, AA, Ph, NReads, TraceOrder, TraceType, SampleRate, 1.0f, SpacMode, TimSpc, RevSpc, false);

	std::cout << "\n now with OCM disabled \n";
	for (i = 0; i < NReads; i++) {
		std::cout << AA[0][i] << " " << AA[1][i] << " " << AA[2][i] << " " << AA[3][i] << " " << AA[4][i] << ";\n";
	}
	//std::cout << "X: " << X[1] <<","<< X[2]<<"," << X[3] <<"\n";
	getchar();
	delete[] x; delete[] Trig; 
	delete[] TraceOrder; delete[] TraceType;
	for (int j = 0; j < NofTraces; j++) {
		delete[] AA[j];
		delete[] Ph[j];
	}
	delete[] AA; delete[] Ph;
    return 0;
}

