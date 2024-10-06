#include "headers.h"


Trade::Trade (double * x_, double a_, double b_, 
							double c_, double d_, int nx_, int ny_, 
							int mx_, int my_, int k_, double eps_, 
							int mi_, int p_, int N_, int task_, char * arg_)
							: x (x_), a (a_), b (b_), c (c_), d (d_), nx (nx_),
								ny (ny_), mx (mx_), my (my_), k (k_), eps (eps_), mi (mi_),
								p (p_), N (N_), task (task_), arg (arg_), r1 (0), r2 (0), r3 (0), r4 (0), 
								t1 (0), t2 (0), it (0), pxy (0), fmax (0), first_step (1), state (0) {}

Data::Data (double * x_, int * I_, double * A_, double * B_, 
						double * r_, double * u_, double * v_, double * block_)
						: x (x_), I (I_), A (A_), B (B_), r (r_), u (u_), v (v_), block (block_){}
/*
Trade::~Trade()
	{
		if (x != nullptr)
			{
				delete [] x;
				x = nullptr;
			}
	}
	
Data::~Data ()
{
	if (x != nullptr)
			{
				delete [] x;
				x = nullptr;
			}
	if (I != nullptr)
			{
				delete [] I;
				I = nullptr;
			}
		if (A != nullptr)
			{
				delete [] A;
				A = nullptr;
			}
		if (B != nullptr)
			{
				delete [] B;
				B = nullptr;
			}
		if (r != nullptr)
			{
				delete [] r;
				r = nullptr;
			}
		if (u != nullptr)
			{
				delete [] u;
				u = nullptr;
			}
		if (v != nullptr)
			{
				delete [] v;
				v = nullptr;
			}
		if (block != nullptr)
			{
				delete [] block;
				block = nullptr;
			}
}


Args::~Args()
	{
		if (tr_gui != nullptr)
			{
				delete tr_gui;
				tr_gui = nullptr;
			}
		if (data != nullptr)
			{
				delete data;
				data = nullptr;
			}
	}
*/

void calc2gui (Trade * tr_calc, Trade * tr_gui)
{
	tr_gui->a = tr_calc->a; tr_gui->b = tr_calc->b; tr_gui->c = tr_calc->c; tr_gui->d = tr_calc->d; 
	tr_gui->nx = tr_calc->nx; tr_gui->ny = tr_calc->ny; 
	tr_gui->mx = tr_calc->mx; tr_gui->my = tr_calc->my;
	tr_gui->k = tr_calc->k; 
	tr_gui->N = tr_calc->N;
	tr_gui->r1 = tr_calc->r1; tr_gui->r2 = tr_calc->r2; tr_gui->r3 = tr_calc->r3; tr_gui->r4 = tr_calc->r4;
	tr_gui->t1 = tr_calc->t1; tr_gui->t2 = tr_calc->t2;
	tr_gui->it = tr_calc->it;
	tr_gui->pxy = tr_calc->pxy;
	tr_gui->fmax = tr_calc->fmax;
	//for (int i = 0; i < tr_calc->N; i++)
		//{
			//(tr_gui->x) [i] = (tr_calc->x) [i];
		//}
}

void calc2gui_mem (Trade * tr_calc, Trade * tr_gui)
{
	for (int i = 0; i < tr_calc->N; i++)
		{
			(tr_gui->x) [i] = (tr_calc->x) [i];
		}
}
