#include <stdio.h>
#include <stdlib.h>
#include <float.h>
#include <time.h>
#include <math.h>
#include <algorithm>
#include <random>

#define debug 0

int frame = 0;
int length = 64;
const int error_target = 200;
int times;
double stddev;
double Eb_N0_dB;
double code_rate = 1./2.;

using namespace::std;

double rand_normal(double mean, double stddev)
{//Box muller method


    static double n2 = 0.0;
    static int n2_cached = 0;
    if (!n2_cached)
    {
        double x, y, r;
        do
        {
            x = 2.0*rand()/RAND_MAX - 1;
            y = 2.0*rand()/RAND_MAX - 1;

            r = x*x + y*y;
        }
        while (r == 0.0 || r > 1.0);
        {
            double d = sqrt(-2.0*log(r)/r);
            double n1 = x*d;
            n2 = y*d;
            double result = n1*stddev + mean;
            n2_cached = 1;
            return result;
        }
    }
    else
    {
        n2_cached = 0;
        return n2*stddev + mean;
    }
}

double max_star(double a)
{
    return a;
}
double max_star(double a, double b)
{
	if (a>b)
		swap(a,b);

	if (b-a > 500)	return b;
	
	return a+log(1+exp(b-a));
}
double max_star(double a, double b,double c)
{
	if (a>b)
		swap(a,b);
	if (a>c)
		swap(a,c);
	if (b>c)
		swap(b,c);
	
	if (c-b > 500)	return c;
	
	if (c-a > 500)	return max_star(b,c);
	
    return a+log(1+exp(b-a)+exp(c-a));
}
double max_star(double a, double b,double c,double d)
{
    if (a>b)
		swap(a,b);
	if (a>c)
		swap(a,c);
	if (a>d)
		swap(a,d);
	if (b>c)
		swap(b,c);
	if (b>d)
		swap(b,d);
	if (c>d)
		swap(c,d);
	
	if (d-c > 500)	return d;
	
	if (d-b > 500)	return max_star(c,d); 
	
	if (d-a > 500)	return max_star(b,c,d);
	
    return a+log(1+exp(b-a)+exp(c-a)+exp(d-a));
}
int sign(double a)
{
    return a>0?1:0;
}


double max(double a)
{
    return a;
}
double max(double a,double b)
{
	return max_star(a,b);
}
double max(double a, double b,double c,double d)
{
    //return max(max(a,b),max(c,d));
    return max_star(a,b,c,d);
}

int main()
{
    srand(time(NULL));
    //std::default_random_engine generator;
    printf ("length = %d\n",length);
    
    FILE *f = fopen("BCJR.txt","w");
	fclose(f);

	for (Eb_N0_dB=-3;Eb_N0_dB<7.41;Eb_N0_dB+=0.1)
	{
		//std::normal_distribution<double> distribution(0,1./SNR);
		stddev = sqrt(pow(10,-Eb_N0_dB/10)/2/code_rate);
		times = 5;
		printf ("SNR = %lf, stddev = %lf\n",Eb_N0_dB,stddev);
		if (debug) printf ("\n");
	    while (times)
	    {
	    	int error_count = 0;
	    	frame = 0;
	    	while(1)
	    	{
	    		frame++;
		        int* message_array;
		        message_array = new int[length];
		        double* coded_array;
		        coded_array = new double[length*2];
		        double* noise_array;
		        noise_array = new double[length*2];
		        double* receive_array;
		        receive_array = new double[length*2];
		
		        int reg0=0;
		        int reg1=0;
		        int reg2=0;
		
		        for (int i=0;i<length;i++)      //construct initial message
		        {
		            if (i>length-3) message_array[i]=0;
		            else message_array[i] = rand()/0.5 > RAND_MAX ? 1 : 0;
		            if (debug) printf ("%d ",message_array[i]);
		            reg2 = reg1;
		            reg1 = reg0;
		            reg0 = message_array[i];
		            coded_array[2*i  ] = (reg0          )%2 ? 1 : -1;			//4
		            coded_array[2*i+1] = (reg0+reg1+reg2)%2 ? 1 : -1;			//7
		        }
		        if (debug) printf ("\n");
		        for (int i=0;i<2*length;i++)
		        {
		            noise_array[i] = rand_normal(0,stddev);
		            receive_array[i] = coded_array[i] + noise_array[i];
		            if (debug) printf ("%+.0lf ",receive_array[i]);
		        }
		        if (debug) printf ("\n");
		        delete [] coded_array;
		        delete [] noise_array;
		
				//BCJR
		        double** alpha;
		        alpha = new double* [4];
		        double** beta;
		        beta = new double* [4];
		        double** gamma;
		        gamma = new double* [8];
		        int u,p;
		        int route[4][4] = {{0,-1,1,-1},{2,-1,3,-1},{-1,4,-1,5},{-1,6,-1,7}};
		        for (int i=0;i<4;i++)
		        {
		            alpha[i] = new double [length];
		            beta[i] = new double [length];
		        }
		        for (int i=0;i<8;i++)
		            gamma[i] = new double [length];
		        //(i)
		        for (int i=0;i<length;i++)
		        {
		            for (int j=0;j<8;j++)
		            {
		            	u = j==0||j==2||j==4||j==6?-1.:1.;											//edit for change g
				        p = j==0||j==3||j==5||j==6?-1.:1.;											//edit for change g
		                gamma[j][i] = u*receive_array[2*i]/(stddev*stddev)+p*receive_array[2*i+1]/(stddev*stddev);				//check
		                if (debug) printf ("%+.6lf ",gamma[j][i]);
		            }
		            if (debug) printf ("\n");
		            for (int j=0;j<4;j++)
		            {
		                if (i==0)   alpha[j][i] = (j==0||j==2?max(0+gamma[route[0][j]][i]):-100.);
		                else
		                {
		                    int s1 = (j==0||j==2)?0:2;
		                    int s2 = (j==0||j==2)?1:3;
		                    int r1 = route[s1][j];
		                    int r2 = route[s2][j];
		                    alpha[j][i] = max(alpha[s1][i-1]+gamma[r1][i],alpha[s2][i-1]+gamma[r2][i]);
		                }
		                if (debug) printf ("%+.6lf ",alpha[j][i]);
		            }
		            if (debug) printf ("\n");
		        }
		        //(ii)
		        for (int i=length-1;i>0;i--)
		        {
		        	if (i==length-1)
		        	{
		        		for (int j=0;j<4;j++)
		        		{
		        			beta[j][i] = (j==0||j==1?0:-100.);
						}
					}
		            for (int j=0;j<4;j++)
		            {
		                    int s1 = (j==0||j==1)?0:1;
		                    int s2 = (j==0||j==1)?2:3;
		                    int r1 = route[j][s1];
		                    int r2 = route[j][s2];
		                    beta[j][i-1] = max(beta[s1][i]+gamma[r1][i],beta[s2][i]+gamma[r2][i]);
		            }
		        }
		        //(iii)
		        int* decoded_array;
		        decoded_array = new int[length];
		        for (int i=0;i<length;i++)
		        {
		            if (i==0)               decoded_array[i] = sign(max(0+gamma[route[0][2]][i]+beta[2][i])-max(0+gamma[route[0][0]][i]+beta[0][i]));
		            else if (i==length-1)   decoded_array[i] = sign(max(alpha[0][i-1]+(-100)+0)-max(alpha[0][i-1]+gamma[route[0][0]][i]+0,alpha[1][i-1]+gamma[route[1][0]][i]+0));
		            else                    decoded_array[i] = sign(max(alpha[0][i-1]+gamma[route[0][2]][i]+beta[2][i], alpha[1][i-1]+gamma[route[1][2]][i]+beta[2][i], alpha[2][i-1]+gamma[route[2][3]][i]+beta[3][i], alpha[3][i-1]+gamma[route[3][3]][i]+beta[3][i])
		                                           -max(alpha[0][i-1]+gamma[route[0][0]][i]+beta[0][i], alpha[1][i-1]+gamma[route[1][0]][i]+beta[0][i], alpha[2][i-1]+gamma[route[2][1]][i]+beta[1][i], alpha[3][i-1]+gamma[route[3][1]][i]+beta[1][i]));
		            if (debug) printf ("%d ",decoded_array[i]);
		        }
		        if (debug) printf ("\n");
		
				//Summary
		        for (int i=0;i<length;i++)
		            if (decoded_array[i]!=message_array[i])
					{
						error_count++;
						break;					//Frame error rate
					} 
		
		        delete [] message_array;
		        delete [] receive_array;
		        for (int i=0;i<4;i++)
		        {
		            delete [] alpha[i];
		            delete [] beta[i];
		        }
		        for (int i=0;i<8;i++)
		            delete [] gamma[i];
		        delete [] alpha;
		        delete [] beta;
		        delete [] gamma;
		        delete [] decoded_array;
		        //scanf ("%c",&temp);
		        if (error_count>=error_target)	break;   
			}
			printf ("frame count = %d\n",frame);	        
		    FILE *f = fopen("BCJR.txt","a");
		    //fprintf (f,"%.1lf %d %d\n",Eb_N0_dB,error_count,length*frame);		//bit error rate
		    fprintf (f,"%.1lf %d %d\n",Eb_N0_dB,error_count,frame);					//frame error rate
		    fclose(f);
		    
		    char temp;
			if (debug) scanf ("%c",&temp);
			
			times--;
	    }
	}

    return 0;
}

