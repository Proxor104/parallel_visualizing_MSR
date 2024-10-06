#include "headers.h"
#define F(I,J) f(x0+(I)*hx,y0+(J)*hy)

static
double f0 (double /*x*/, double /*y*/)
{
	return 1;
}

static
double f1 (double x, double /*y*/)
{
	return x;
}

static
double f2 (double /*x*/, double y)
{
	return y;
}

static
double f3 (double x, double y)
{
	return x + y;
}

static
double f4 (double x, double y)
{
	return sqrt (x * x + y * y);
}

static
double f5 (double x, double y)
{
	return x * x + y * y;
}

static
double f6 (double x, double y)
{
	return exp (x * x - y * y);
}

static
double f7 (double x, double y)
{
	return 1 / (25 * (x * x + y * y) + 1);
}

static double * results = nullptr;


int init_reduce_sum (int p)
{
	//printf ("Hi1\n");
	results = new double [p];
	if (results == nullptr) return -1;
	return 0;
}


double reduce_max_det (int p, int k, double m)
{
	double max = 0; int l = 0;
	results [k] = m;
	reduce_sum (p);
	for ( ; l < p; l++)
		{
			 if (m > max) max = m;
		}
	reduce_sum (p);
	return max;
}


double reduce_sum_det (int p, int k, double s)
{
	double sum = 0; int l = 0;
	results [k] = s;
	reduce_sum (p);
	for ( ; l < p; l++)
		{
			sum += results [l];
		}
	reduce_sum (p);
	return sum;
}


int finilize_reduce_sum (int /*p*/)
{
	if(results)
		{
			//printf ("Hi aaaaa\n");
			delete [] results;
		}
	return 0;
}

void * thread_main(void * ptr)
{
	Args * Arg = (Args *) ptr;
	Trade * tr_gui = Arg->tr_gui;
	Data * data = Arg->data;
  static pthread_cond_t * cond_start = Arg->cond_start;
	static pthread_mutex_t * mut = Arg->mut;
  static pthread_cond_t * cond = Arg->cond;
	int & state = tr_gui->state;
	int & first_step = tr_gui->first_step;
	double a = 0, b = 0, c = 0, d = 0;
	int nx = 0, ny = 0;
	int k = 0;
	double eps = tr_gui->eps;
	int mi = tr_gui->mi;
	int p = tr_gui->p;
	int q = Arg->q;
	int main_q = 0;
	int err = 0;
	int len_msr = 0;
	int N = 0;
	double & r1 = tr_gui->r1;
	double & r2 = tr_gui->r2;
	double & r3 = tr_gui->r3;
	double & r4 = tr_gui->r4;
	double & t1 = tr_gui->t1;
	double & t2 = tr_gui->t2;
	int & it = tr_gui->it;
	int pxy = 0;
	double fmax = 0;
	
	double hx = 0;
	double hy = 0;
	
	double * x_tr = nullptr;
	double * x = nullptr;
	int * I = nullptr;
	double * A = nullptr;
	double * B = nullptr;
	double * u = nullptr;
	double * v = nullptr;
	double * r = nullptr;
	double * block = nullptr;
	double (*func) (double, double) = nullptr;
	
	cpu_set_t cpu;
	CPU_ZERO(&cpu);
	int n_cpus = get_nprocs();
	int cpu_id = n_cpus - 1 - (q % n_cpus);
	CPU_SET(cpu_id, &cpu);
	pthread_t tid = pthread_self();
	pthread_setaffinity_np(tid, sizeof(cpu), &cpu);
	
	int max_steps = 100;
	double time = 0;
	
	int i1 = 0, i2 = 0, i = 0;
	
	reduce_sum (p);		
	if (q == main_q)
		{
			init_reduce_sum (p);
		}
	reduce_sum (p);
	
	for (;state != -2;)
		{
			state = 1;
			a = tr_gui->a, b = tr_gui->b, c = tr_gui->c, d = tr_gui->d;
			nx = tr_gui->nx, ny = tr_gui->ny;
			hx = (b - a) / nx;
			hy = (d - c) / ny;
			k = tr_gui->k;
			len_msr = get_len_msr (nx, ny) + 1 + (nx + 1) * (ny + 1);
			N = (nx + 1) * (ny + 1);
			pxy = tr_gui->pxy;
			fmax = tr_gui->fmax;
			
			if (q == main_q)
				{

					data->x = new double [N];
					data->I = new int [len_msr];
					data->A = new double [len_msr];
					data->B = new double [N];
				 	data->u = new double [N]; 
				 	data->r = new double [N]; 
				 	data->v = new double [N]; 
					data->block = new double [N];
				}
			reduce_sum (p);
			
			x = data->x;
			x_tr = tr_gui->x;
			I = data->I;//
			A = data->A;//
			B = data->B;
			u = data->u;
			v = data->v;
			r = data->r;
			block = data->block;
			
			reduce_sum (p);
			
			thread_rows (N, p, q, i1, i2);
			for (i = i1; i < i2; i++)
				{
					x [i] = 0;
					B [i] = 0;
					u [i] = 0;
					r [i] = 0;
					u [i] = 0;
					v [i] = 0;
					block [i] = 0;
 				}
 			thread_rows (len_msr, p, q, i1, i2);
			for (i = i1; i < i2; i++)
				{
					A [i] = 0;
					I [i] = 0;
 				}
 			
			switch (k)
				{
					case 0:
						func = f0;
						break;
					case 1:
						func = f1;
						break;
					case 2:
						func = f2;
						break;
					case 3:
						func = f3;
						break;
					case 4:
						func = f4;
						break;
					case 5:
						func = f5;
						break;
					case 6:
						func = f6;
						break;
					case 7:
						func = f7;
						break;
				}
			
			reduce_sum (p);
			
			if (q == main_q)
				{
					if (fill_I (nx, ny, I) != len_msr) 
						{
							err = -1;
						}
				}
			reduce_sum (p, &err, 1);
			if (err < 0)
				{
					printf ("fill_I (nx, ny, I) != len_msr\n");
					state = -1;
					break;
				}
			err = 0;

			if (fill_IA (nx, ny, hx, hy, I, A, p, q) < 0) 
				{
					//printf ("Hi\n");
					err = -1;
				}
			reduce_sum (p, &err, 1);
			if (err < 0)
				{
					printf ("fill_IA (nx, ny, hx, hy, I, A, p, q) < 0\n");
					state = -1;
					break;
				}
			err = 0;
			reduce_sum (p);
			
			if (fill_B (nx, ny, hx, hy, a, c, func, B, p, q, pxy, fmax) < 0) 
				{
					err = -1;
				}
			reduce_sum (p, &err, 1);
			if (err < 0)
				{
					printf ("fill_B (nx, ny, hx, hy, a, c, func, B, p, q, pxy, fmax) < 0\n");
					state = -1;
					break;
				}
			err = 0;
			reduce_sum (p);

			
			time = get_full_time();
			it = minimal_residual_msr_matrix_full (N, A, I, B, x, r, u, v, block, eps, mi, max_steps, p, q);
			time = get_full_time() - time;
			t1 = time;

			reduce_sum (p);

			time = 0;
			time = get_full_time();
			get_value_r_1_2 (x, func, nx, ny, hx, hy, r1, r2, a, c, p, q);
			get_value_r_3_4 (x, func, nx, ny, hx, hy, r3, r4, a, c, p, q, pxy, fmax);
			time = get_full_time() - time;
			t2 = time;

			reduce_sum (p);
			
			thread_rows (N, p, q, i1, i2);
			for (i = i1; i < i2; i++)
				{
					 x_tr [i] = x [i];
 				}
			
			reduce_sum (p);
			
			if (first_step == 1)
				{
					reduce_sum (p);
					if (q == main_q)
						{
							first_step = 0;
							pthread_cond_broadcast(cond_start);
						}	
					reduce_sum (p);
				}
			
			if (q == main_q)
				{
					//if (data->x != nullptr)
					if (data->x != nullptr)
						{
							//printf ("Hi  20\n");
							delete [] data->x;
							data->x = nullptr;
							x = nullptr;
						}
					if (data->I != nullptr)
						{
							//printf ("Hi  20\n");
							delete [] data->I;
							data->I = nullptr;
							I = nullptr;
						}
					if (data->A != nullptr)
						{
							//printf ("Hi  20\n");
							delete [] data->A;
							data->A = nullptr;
							A = nullptr;
						}
					if (data->B != nullptr)
						{
							//printf ("Hi  20\n");
							delete [] data->B;
							data->B = nullptr;
							B = nullptr;
						}
					if (data->u != nullptr)
						{
							//printf ("Hi  20\n");
							delete [] data->u; 
						 	data->u = nullptr;
						 	u = nullptr;
						}
					if (data->r != nullptr )
						{
							//printf ("Hi  20\n");
							delete [] data->r; 
						 	data->r = nullptr;
						 	r = nullptr;
						}
					if (data->v != nullptr )
						{
							//printf ("Hi  20\n");
							delete [] data->v; 
						 	data->v = nullptr;
						 	v = nullptr;
						}
					if (data->block != nullptr )
						{
							//printf ("Hi  20\n");
							delete [] data->block;
							data->block = nullptr;
							block = nullptr;
						}							
				}
			reduce_sum (p);
			
			state = 0;
			//printf ("Hi6 %d\n", q);
			
			pthread_mutex_lock(mut);
			while (state == 0)
				{
					 pthread_cond_wait(cond, mut);
				}
			pthread_mutex_unlock(mut);
			/*
			if (first_step != 1)
				{
					state = -1;
					break;
				}
			*/
		}
	reduce_sum (p);
	if (q == main_q)
		{
			finilize_reduce_sum (p);	
		}
	reduce_sum (p);
	if (state == -1)
		{
			if (q == main_q)
				{
					printf ("We got an error!\n");
					
				}
			pthread_mutex_lock(mut);
			while (state != -2)
				{
					 pthread_cond_wait(cond, mut);
				}
			pthread_mutex_unlock(mut);
		}

	return  0;
}

double scalar_product (int n, double * x, double * y, int p, int k)
{
	int i = 0, i1 = 0, i2 = 0;
	double s = 0;
	thread_rows (n, p, k, i1, i2);
	for (i = i1; i < i2; i++)
		{
			//printf ();
			s += x [i] * y [i];
		}
	s = reduce_sum_det (p, k, s);
	//printf ("k = %d s = %10.3e\n", k, s);
	return s;
}


void l2ij_middle (int nx, int /* ny*/ , int l, int & i, int & j)
{
	j = l / nx;
	i = l - j * nx;
}


void thread_cords_r (int n, int p, int k, int & i1, int & i2)
{
	i1 = n * k; i1 /= p;
	i2 = n * (k + 1); i2 /= p;	
}


void get_value_r_1_2 (double * x, double (*f) (double, double), int nx, int ny, double hx, double hy, double & r1, double & r2, double a, double c, int p, int k)
{// нужен массив results
	int l = 0, l1 = 0, l2 = 0;
	int i = 0, j = 0;
	double sum = 0;
	double x1 = 0, y1 = 0, x2 = 0, y2 = 0;
	int num_points = nx * ny;
	double value1 = 0, value2 = 0;
	double max = 0;
	double res1 = 0, res2 = 0;
	thread_cords_r (num_points, p, k, l1, l2);
	for (l = l1; l < l2; l++)
		{
		
			l2ij_middle (nx, ny, l, i, j);
			x1 = (i * hx) + a + (2 * hx / 3); 
			y1 = (j * hy) + c + (hy / 3);
			
			l2ij_middle (nx, ny, l, i, j);
			x2 = (i * hx) + a + (hx / 3); 
			y2 = (j * hy) + c + (2 * hy / 3);
			
		  //printf ("%10.3e %10.3e %10.3e %10.3e\n", x1, y1, x2, y2);
			
			value1 = get_value_pf_center (x1, y1, nx, ny, hx, hy, a, c, x, 1);
			value2 = get_value_pf_center (x2, y2, nx, ny, hx, hy, a, c, x, -1);
			//printf ("k = %d %10.3e %10.3e\n", k, value1, value2);
			
			res1 = fabs (value1 - f (x1, y1));
			res2 = fabs (value2 - f (x2, y2));
			
			sum += res1 + res2;
			
			if (res1 > max) max = res1;
			if (res2 > max) max = res2;
		}
	reduce_sum_det (p, k, sum);
	r2 = sum * hx * hy / 2;
	reduce_max_det (p, k, max);
	r1 = max;
}



void get_value_r_3_4 (double * x, double (*f) (double, double), int nx, int ny, double hx, double hy, double & r3, double & r4, double a, double c, int p, int k, int pxy, double fmax)
{// нужен массив results
	int l = 0, l1 = 0, l2 = 0; int l3 = 0;
	int i = 0, j = 0;
	double sum = 0;
	double x1 = 0, y1 = 0;
	int N = (nx + 1) * (ny + 1);
	double value1 = 0;
	double max = 0;
	double res1 = 0;
	thread_rows (N, p, k, l1, l2);
	ij2l (nx, ny, nx / 2, ny / 2, l3);
	for (l = l1; l < l2; l++)
		{
			l2ij (nx, ny, i, j, l);
			x1 = (i * hx) + a; 
			y1 = (j * hy) + c;

			value1 = get_value_pf_corr (x1, y1, nx, ny, hx, hy, a, c, x);
			if (l == l3)
				{
					res1 = fabs (value1 - (f (x1, y1) + 0.1 * fmax * pxy));
				}
			else 
				{
					res1 = fabs (value1 - f (x1, y1));
				}

			//printf ("pf(%10.3e,%10.3e) = %10.3e | f(%10.3e,%10.3e) = %10.3e\n", x1, y1, value1, x1, y1, f (x1, y1));
			
			sum += res1;
			if (res1 > max) max = res1;
		}
	reduce_sum_det (p, k, sum);
	r4 = sum * hx * hy;
	reduce_max_det (p, k, max);
	r3 = max;
}

int minimal_residual_msr_matrix (int n, double * A, int * I, double * b, double * x, double * r, double * u, 
                  double * v, double * block, double eps, int max_it, int p, int k)
{
	double prec = 0, b_norm2 = 0, tau = 0, c1 = 0, c2 = 0;
	int it = 0;
	b_norm2 = scalar_product (n, b, b, p, k);
	prec = b_norm2 * eps * eps;
	matrix_mult_vector_msr (n, A, I, x, r, p, k);
	mult_sub_vector (n, r, b, 1, p, k); // r -= 1 * b
	for (it = 0; it < max_it; it++)
		{
			apply_precondition_msr_matrix (n, A, I, v, r, block, p, k); // Mv = r
			matrix_mult_vector_msr (n, A, I, v, u, p, k); // u = Av
			c1 = scalar_product (n, u, r, p, k); // c1 = (u,r)
			c2 = scalar_product (n, u, u, p, k); // c2 = (u,u)
			if (c1 < prec || c2 < prec) break;
			tau = c1 / c2;
			mult_sub_vector (n, x, v, tau, p, k); // x -= tau v
			mult_sub_vector (n, r, u, tau, p, k); // r -= tau u
		}
	if (it >= max_it) return -1;
	return it;
}


int minimal_residual_msr_matrix_full (int n, double * A, int * I, double * b, double * x, double * r, double * u, 
                                     double * v, double * block, double eps, int max_it, int max_steps, int p, int k)
{
	int step = 0, ret = 0, its = 0;
	for ( ; step < max_steps; step++)
		{
			ret = minimal_residual_msr_matrix (n, A, I, b, x, r, u, v, block, eps, max_it, p, k);
			if (ret >= 0) 
				{
					its += ret;
					break;
				}
			its += max_it; 
		}
	if (step >= max_steps) return -1;
	return its;
}

void reduce_sum (int p, double * a, int n)
{
  static pthread_mutex_t m = PTHREAD_MUTEX_INITIALIZER;
  static pthread_cond_t c_in = PTHREAD_COND_INITIALIZER;
  static pthread_cond_t c_out = PTHREAD_COND_INITIALIZER;
  static int thread_in = 0;
  static int thread_out = 0;
  static double * r = nullptr;
  //if(n == 0) a = nullptr;
  int i = 0;
  if (p <= 1) return ;

  pthread_mutex_lock(&m);
  if (r == nullptr)
    r = a;
  else
		for(i = 0; i < n; i++)
			r[i] += a[i];
  thread_in++;
  if (thread_in >= p)
    {
      thread_out = 0;
      pthread_cond_broadcast(&c_in);
    }
  else
    while (thread_in < p)
      pthread_cond_wait(&c_in, &m);

  if (r != a)
    for(i = 0; i < n; i++)
      {
        a[i] = r[i];
      }
  thread_out++;
  if(thread_out >= p)
    {
      thread_in = 0;
      r = nullptr;
      pthread_cond_broadcast(&c_out);
    }
  else
    while (thread_out < p)
      pthread_cond_wait(&c_out, &m);
  pthread_mutex_unlock(&m);
}


void reduce_sum (int p, int * a, int n)
{
  static pthread_mutex_t m = PTHREAD_MUTEX_INITIALIZER;
  static pthread_cond_t c_in = PTHREAD_COND_INITIALIZER;
  static pthread_cond_t c_out = PTHREAD_COND_INITIALIZER;
  static int thread_in = 0;
  static int thread_out = 0;
  static int * r = nullptr;
  //if(n == 0) a = nullptr;
  int i = 0;
  if (p <= 1) return ;

  pthread_mutex_lock(&m);
  if (r == nullptr)
   	r = a;
  else
		for(i = 0; i < n; i++)
			r[i] += a[i];
  thread_in++;
  if (thread_in >= p)
    {
      thread_out = 0;
      pthread_cond_broadcast(&c_in);
    }
  else
    while (thread_in < p)
      pthread_cond_wait(&c_in, &m);

  if (r != a)
    for(i = 0; i < n; i++)
      {
        a[i] = r[i];
      }
  thread_out++;
  if(thread_out >= p)
    {
      thread_in = 0;
      r = nullptr;
      pthread_cond_broadcast(&c_out);
    }
  else
    while (thread_out < p)
      pthread_cond_wait(&c_out, &m);
  pthread_mutex_unlock(&m);
}


int get_len_msr (int nx, int ny)
{
	return (nx - 1) * (ny - 1) * 6 + (2 * (nx - 1) + 2 * (ny - 1)) * 4 + 2 * 3 + 2 * 2;
}


int get_len_msr_off_diag (int nx, int ny, double hx, double hy)
{
	int i = 0, j = 0, res = 0;
	for ( ; i <= nx; i++)
		{
			for (j = 0; j <= ny; j++)
				{
					res += get_off_diag (nx, ny, hx, hy, i, j);
				}
		}
	return res;
}


void reduce_sum (int p)
{
	int i = 0;
	reduce_sum (p, &i, 0);
}


void thread_rows (int n, int p, int k, int & i1, int & i2)
{
	i1 = n * k; i1 /= p;
	i2 = n * (k + 1); i2 /= p;	
}


void mult_sub_vector (int n, double * x, double * y, double t, int p, int k)
{
	int i = 0, i1 = 0, i2 = 0;
	thread_rows (n, p, k, i1, i2);
	for (i = i1; i < i2; i++)
		{
			x [i] -= t * y [i];
		}
	reduce_sum (p);
}


void matrix_mult_vector_msr (int n, double * A, int * I, double * x, double * y, int p, int k)
{
	int i = 0, j = 0, l = 0, J = 0, i1 = 0, i2 = 0;
	double s = 0;
	thread_rows (n, p, k, i1, i2);
	for (i = i1; i < i2; i++)
		{
			s = A [i] * x [i];
			l = I [i + 1] - I [i];
			J = I [i];
			for (j = 0; j < l; j++)
				{
					s += A [J + j] * x [I [J + j]];
				}
			y [i] = s;
		}
	reduce_sum (p);
}

/*
void apply_precondition_msr_matrix (int n, double * A, int * I, double * v, double * r, double * block, int p, int k) // Mv = r
{
	int i = 0, j = 0, l = 0, J = 0, i1 = 0, i2 = 0;
	double s = 0;
	int main_k = 0;
	if (k == main_k)
		{
			for (i = 0; i < n; i++)
				{
					s = r [i];
					l = I [i + 1] - I [i];
					J = I [i];
					for (j = 0; j < l; j++)
						{
							if (I [J + j] < i)
								{
									s -= A [J + j] * block [I [J + j]];
								}
						}
					block [i] = s / A [i];
				}
		}
	reduce_sum (p);
	thread_rows (n, p, k, i1, i2);
	for (i = i1; i < i2; i++)
		{
			block [i] *= A [i] ;
		}
	reduce_sum (p);
	
	if (k == main_k)
		{
			for (i = n - 1; i >= 0; i--)
				{
					s = block [i];
					l = I [i + 1] - I [i];
					J = I [i];
					for (j = 0; j < l; j++)
						{
							if (I [J + j] > i)
								{
									s -= A [J + j] * v [I [J + j]];
								}
						}
					v [i] = s / A [i];
				}
		}
	reduce_sum (p);
}
*/



void apply_precondition_msr_matrix (int n, double * A, int * I, double * v, double * r, double * block, int p, int k) // Mv = r
{
	int i = 0, j = 0, l = 0, J = 0, i1 = 0, i2 = 0;
	double s = 0;
	//int main_k = 0;

	thread_rows (n, p, k, i1, i2);
	
	s = r [i1];
	l = I [i1 + 1] - I [i1];
	J = I [i1];
	block [i1] = s / A [i1];
	
	for (i = i1 + 1; i < i2; i++)
		{
			s = r [i];
			l = I [i + 1] - I [i];
			J = I [i];
			for (j = 0; j < l; j++)
				{
					if (I [J + j] >= i1 && I [J + j] < i)
						{
							s -= A [J + j] * block [I [J + j]];
						}
				}
			block [i] = s / A [i];
		}
	
	for (i = i1; i < i2; i++)
		{
			block [i] *= A [i];
		}

	s = block [i2 - 1];
	l = I [i2] - I [i2 - 1];
	J = I [i2 - 1];
	v [i2 - 1] = s / A [i2 - 1];

	for (i = i2 - 2; i >= i1; i--)
		{
			s = block [i];
			l = I [i + 1] - I [i];
			J = I [i];
			for (j = 0; j < l; j++)
				{
					if (I [J + j] > i && I [J + j] < i2)
						{
							s -= A [J + j] * v [I [J + j]];
						}
				}
			v [i] = s / A [i];
		}
	reduce_sum (p);
} 


/*
void apply_precondition_msr_matrix (int n, double * A, int * I, double * v, double * r, double * tmp, int p, int k) // Mv = r
{
        int i = 0, j = 0, l = 0, J = 0, i1 = 0, i2 = 0;
	//double s = 0;
	//int main_k = 0;
	thread_rows (n, p, k, i1, i2);
	for (i = i1; i < i2; i++)
		{
			v [i] = r [i] / A [i] ;
		}
	reduce_sum (p);
}
*/

double get_value_of_approx (double x, double y, int nx, int ny, double hx, double hy, double a, double c, double * vec)
{
	int i = (int) ((x - a) / hx);
	int j = (int) ((y - c) / hy);
	double is_i = (x - a) - i * hx;
	double is_j = (y - c) - j * hy;
	int l = 0;
	double x1 = a + i * hx;
	double y1 = c + j * hy;
	double x2 = a + (i + 1) * hx;
	double y2 = c + (j + 1) * hy;
	double var = (y2 * (x - x1) - y1 * (x - x2)) / (x2 - x1);
	//printf ("%d %d %10.3e %10.3e\n", i, j, y, var);
	double num1 = 0, num2 = 0, num3 = 0; 
	double res = 0;
	if (fabs (is_i) < 1.0e-15 && fabs (is_j) < 1.0e-15)
		{
			ij2l (nx, ny, i, j, l);
			return vec [l];
		}
	else if (var > y)
		{
			num1 = (x - x1) / hx; // i + 1 j + 1
			num3 = 1 - ((y - y1) / hy); // i j
			num2 = ((hx * (y - y1) - hy * (x - x1)) / (hx * hy)); // i j + 1
			ij2l (nx, ny, i, j, l);
			res += vec [l] * num3;
			ij2l (nx, ny, i + 1, j + 1, l);
			res += vec [l] * num1;
			ij2l (nx, ny, i, j + 1, l);
			res += vec [l] * num2;
		}
	else if (var < y)
		{
			num1 = 1 - ((x - x1) / hx); // i j
			num3 = (y - y1) / hy; // i + 1 j + 1
			num2 = ((hy * (x - x1) - hx * (y - y1)) / (hx * hy)); // i + 1 j
			ij2l (nx, ny, i, j, l);
			res += vec [l] * num1;
			ij2l (nx, ny, i + 1, j, l);
			res += vec [l] * num2;
			ij2l (nx, ny, i + 1, j + 1, l);
			res += vec [l] * num3;
		}
	else if (fabs (var - y) < 1.0e-15) 
		{
			num1 = 1 - ((x - x1) / hx); // i j
			num3 = (y - y1) / hy; // i + 1 j + 1
			num2 = ((hy * (x - x1) - hx * (y - y1)) / (hx * hy)); // i + 1 j
			ij2l (nx, ny, i, j, l);
			res += vec [l] * num1;
			ij2l (nx, ny, i + 1, j, l);
			res += vec [l] * num2;
			ij2l (nx, ny, i + 1, j + 1, l);
			res += vec [l] * num3;
		}
	return res;
}


double get_value_pf_corr (double x, double y, int nx, int ny, double hx, double hy, double a, double c, double * vec)
{
	int i = round ((x - a) / hx);
	int j = round ((y - c) / hy);
	int l = 0;
	ij2l (nx, ny, i, j, l);
	return vec [l];
}


double get_value_pf_center (double x, double y, int nx, int ny, double hx, double hy, double a, double c, double * vec, int up_down)
{
	int i = round ((x - a) / hx);
	int j = round ((y - c) / hy);
	int l = 0;
	double x1 = a + i * hx;
	double y1 = c + j * hy;
	double num1 = 0, num2 = 0, num3 = 0; 
	double res = 0;
	if (up_down > 0)
		{
			num1 = (-1) * (x - x1) / hx; // i - 1 j
			num3 = (y - y1) / hy; // i j + 1
			num2 = 1 + ((hy * (x - x1) - hx * (y - y1)) / (hx * hy)); // i j
			ij2l (nx, ny, i, j, l);
			res += vec [l] * num2;
			ij2l (nx, ny, i - 1, j, l);
			res += vec [l] * num1;
			ij2l (nx, ny, i, j + 1, l);
			res += vec [l] * num3;
		}
	else if (up_down < 0)
		{
			num1 = (x - x1) / hx; // i + 1 j
			num3 = (-1) * (y - y1) / hy; // i j - 1
			num2 = 1 + ((hx * (y - y1) - hy * (x - x1)) / (hx * hy)); // i j
			ij2l (nx, ny, i, j, l);
			res += vec [l] * num2;
			ij2l (nx, ny, i + 1, j, l);
			res += vec [l] * num1;
			ij2l (nx, ny, i, j - 1, l);
			res += vec [l] * num3;
		}
	return res;
}

void ij2l (int nx, int /* ny */, int i, int j, int & l)
{
	l = i + j * (nx + 1);
}


void l2ij (int nx, int /* ny */, int & i, int & j, int l)
{
	j = l / (nx + 1);
	i = l - j * (nx + 1);
}


void IA_ij (int nx, int ny, double hx, double hy, int i, int j, int is, int js, int s, int * I, double * A)
{
	int l = 0, ls = 0;
	ij2l (nx, ny, i, j, l);
	ij2l (nx, ny, is, js, ls);
	if (I) I [s] = ls;
	if (A)
		{
			if (l == ls)
				{
					A [s] = (
									 (i < nx && j > 0 ? 1 : 0) + (i > 0 && j > 0 ? 2 : 0) + 
					         (i > 0 && j < ny ? 1 : 0) + (i < nx && j < ny ? 2 : 0)
					        ) * hx * hy / 12;
				}
			else
				{
					A [s] = (
									(is == i + 1 && js == j ? ((j < ny ? 1 : 0) + (j > 0 ? 1 : 0)) : 0) + 
									(is == i && js == j - 1 ? ((i < nx ? 1 : 0) + (i > 0 ? 1 : 0)) : 0) +
									(is == i - 1 && js == j - 1 ? 2 : 0) +
									(is == i - 1 && js == j ? ((j > 0 ? 1 : 0) + (j < ny ? 1 : 0)) : 0) +
									(is == i && js == j + 1 ? ((i > 0 ? 1 : 0) + (i < nx ? 1 : 0)) : 0) +
									(is == i + 1 && js == j + 1 ? 2 : 0)
									) * hx * hy / 24;
				}
		}
}


int get_off_diag (int nx, int ny, double hx, double hy, int i, int j, int * I, double * A)
{
	int s = 0;
	if (i < nx) 
		{
			IA_ij (nx, ny, hx, hy, i, j, i + 1, j, s++, I, A);
		}
	if (j > 0) 
		{	
			IA_ij (nx, ny, hx, hy, i, j, i, j - 1, s++, I, A);
		}
	if (i > 0 && j > 0)
		{		
			IA_ij (nx, ny, hx, hy, i, j, i - 1, j - 1, s++, I, A);
		}
	if (i > 0) 
		{	
			IA_ij (nx, ny, hx, hy, i, j, i - 1, j, s++, I, A);
		}
	if (j < ny) 
		{	
			IA_ij (nx, ny, hx, hy, i, j, i, j + 1, s++, I, A);
		}
	if (i < nx && j < ny) 
		{	
			IA_ij (nx, ny, hx, hy, i, j, i + 1, j + 1, s++, I, A);
		}
	return s;
}


int fill_I (int nx, int ny, int * I)
{
	int i = 0, j = 0, l = 0, r = 0, s = 0;
	int N = (nx + 1) * (ny + 1);
	double hx = 0, hy = 0;
	r = N + 1;
	for ( ; l < N; l++)
		{
			l2ij (nx, ny, i, j, l);
			s = get_off_diag (nx, ny, hx, hy, i, j);
			I [l] = r;
			r += s;
		}
	I [l] = r;
	return r; // д.б. равно общей длине матрицы (N + 1 + get_len_msr ())
}


int get_diag (int nx, int ny, double hx, double hy, int i, int j, int * /* I */, double * A)
{
	int l = 0;
	ij2l (nx, ny, i, j, l);
	IA_ij (nx, ny, hx, hy, i, j, i, j, l /* s */, nullptr, A);
	return 0;
}


int fill_IA (int nx, int ny, double hx, double hy, int * I, double * A, int p, int k)
{ // fill_I уже вызвано
	int i = 0, j = 0, l = 0, l1 = 0, l2 = 0, N = (nx + 1) * (ny + 1), s = 0, r = 0, err = 0, len = 0, t = 0;
	l1 = N * k; l1 /= p; l2 = N * (k + 1); l2 /= p;
	for (l = l1; l < l2; l++)
		{
			l2ij (nx, ny, i, j, l);
			if (get_diag (nx, ny, hx, hy, i, j, I, A + 0) != 0)
				{
					err = -1;
					break;
				}
			r = I [l];
			s = I [l + 1] - I [l];
			t = get_off_diag (nx, ny, hx, hy, i, j, I + r, A + r);
			if (s != t)
				{
					//printf ("%d %d %d\n", k, s, t);
					err = -2;
					break;
				}
			len += s;
		}
	reduce_sum (p, &err, 1);
	if(err < 0)
		{
			//printf ("hi %d\n",err);
			return -1;
		}
	reduce_sum (p, &len, 1);
	if(I [N] != (N + 1) + len)
		{
			return -2;
		}
	return 0;
}


int fill_B (int nx, int ny, double hx, double hy, double x0, double y0, double (*f) (double, double), double * B, int p, int q, int pxy, double fmax)
{
	int l1 = 0, l2 = 0, l = 0, N = (nx + 1) * (ny + 1); int l3 = 0;
	l1 = N * q; l1 /= p; l2 = N * (q + 1); l2 /= p;
	ij2l (nx, ny, nx / 2, ny / 2, l3);
	for (l = l1; l < l2; l++)
		{
			if (l == l3)
				{
					B [l] = F_ij (nx, ny, hx, hy, x0, y0, f, l) + pxy * 0.1 * fmax;
				}
			else 
				{
					B [l] = F_ij (nx, ny, hx, hy, x0, y0, f, l);
				}
		}
	return 0;
}

/*
double F_ij (int nx, int ny, double hx, double hy, double x0, double y0, double (*f)(double, double), int l)
{
	int i = 0, j = 0;
	l2ij (nx, ny, i, j, l);
	return (
	       (i < nx && j > 0 ? (2 * F(i,j) + F(i+1,j) + F(i,j-1)) : 0) + 
	       (i > 0 && j > 0 ? (4 * F(i,j) + F(i,j-1) + 2 * F(i-1,j-1)+ F(i-1,j)) : 0) + 
	       (i > 0 && j < ny ? (2 * F(i,j) + F(i-1,j) + F(i,j+1)) : 0) +
	       (i < nx && j < ny ? (4 * F(i,j) + F(i,j+1) + 2 * F(i+1,j+1) + F(i+1,j)) : 0)
	       ) * hx * hy / 24;
}
*/
double F_ij (int nx, int ny, double hx, double hy, double x0, double y0, double (*f)(double, double), int l)
{
	int i = 0, j = 0;
	l2ij (nx, ny, i, j, l);
	return (
	       (i < nx && j > 0 ? (6 * F(i,j) + 10 * F(i+0.5,j) + F(i+1,j) + 4 * F(i+0.5,j-0.5) + F(i,j-1) + 10 * F(i,j-0.5)) : 0) + 
	       (i > 0 && j > 0 ? (12 * F(i,j) + 10 * F(i,j-0.5) + F(i,j-1) + 4 * F(i-0.5,j-1) + 2 * F(i-1,j-1) 
	                         + 20 * F(i-0.5,j-0.5) + 4 * F(i-1,j-0.5) + F(i-1,j) + 10 * F(i-0.5,j)) : 0) + 
	       (i > 0 && j < ny ? (6 * F(i,j) + 10 * F(i-0.5,j) + F(i-1,j) + 4 * F(i-0.5,j+0.5) + F(i,j+1) + 10 * F(i,j+0.5)) : 0) +
	       (i < nx && j < ny ? (12 * F(i,j) + 10 * F(i,j+0.5) + F(i,j+1) + 4 * F(i+0.5,j+1) + 2 * F(i+1,j+1) + 
	                           20 * F(i+0.5,j+0.5) + 4 * F(i+1,j+0.5) + F(i+1,j) + 10 * F(i+0.5,j)) : 0)
	       ) * hx * hy / 192;
}

