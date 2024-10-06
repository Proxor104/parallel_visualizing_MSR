#include "headers.h"

int input_args (double & a, double & b, double & c, double & d, int & nx, int & ny, 
                int & mx, int & my, int & k, double & eps, int & mi, int & p, int argc, char * argv[])
{
	if(!((argc==13)&&
                (sscanf(argv[1],"%le",&a) == 1)&&
                (sscanf(argv[2],"%le",&b) == 1)&&
                (sscanf(argv[3],"%le",&c) == 1)&&
                (sscanf(argv[4],"%le",&d) == 1)&&
		(sscanf(argv[5],"%d",&nx) == 1)&&
		(sscanf(argv[6],"%d",&ny) == 1)&&
		(sscanf(argv[7],"%d",&mx) == 1)&&
		(sscanf(argv[8],"%d",&my) == 1)&&
		(sscanf(argv[9],"%d",&k) == 1)&&
                (sscanf(argv[10],"%lf",&eps) == 1)&&
		(sscanf(argv[11],"%d",&mi) == 1)&&
		(sscanf(argv[12],"%d",&p) == 1))||
                (k < 0 || k > 7)||
                (fabs(b - a) < 1.0e-06)||
                (fabs(d - c) < 1.0e-06))
		{
			if ((argc == 1))
				{
					a = -1;
					b = 1;
                                        c = -0.5;
                                        d = 0.5;
					nx = 10;
					ny = 10;
					mx = 2;
					my = 2;
					k = 2;
					eps = 1.0e-08;
					mi = 100;
					p = 3;
					return 0;
				}
			return -1;
		}
	return 0;
}

