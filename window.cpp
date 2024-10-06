
#include <QPainter>
#include <stdio.h>

#include "headers.h"

#define DEFAULT_MODE 1
#define DEFAULT_TIMER 500
#define DEFAULT_IS_TIMER 0
#define DEFAULT_FMAX 0

#define L2G(X,Y) (l2g ((X), (Y)))

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

Window::Window (QWidget * parent)
  : QWidget (parent)
{
  mode = DEFAULT_MODE;
  is_timer = DEFAULT_IS_TIMER;
  fmax = DEFAULT_FMAX;
  s = 0;
  p = 0;
  par = parent;
}

void Window::set_param (Trade * tr_calc_, Trade * tr_gui_, pthread_mutex_t * mut_, pthread_cond_t * cond_)
{
	tr_calc = tr_calc_;
	tr_gui = tr_gui_;
  mut = mut_;
  cond = cond_;
	synhronize_param ();
	x = tr_gui->x;
  //func_id = k;
  arg = tr_gui->arg;
  
  //set_func ();
}

QSize Window::minimumSizeHint () const
{
  return QSize (300, 300);
}

QSize Window::sizeHint () const
{
  return QSize (1000, 1000);
}

void Window::set_func ()
{
	switch (func_id)
		{
			case 0:
				f_name = "f(x, y) = 1";
				f = f0;
				break;
			case 1:
				f_name = "f(x, y) = x";
				f = f1;
				break;
			case 2:
				f_name = "f(x, y) = y";
				f = f2;
				break;
			case 3:
				f_name = "f(x, y) = x + y";
				f = f3;
				break;
			case 4:
				f_name = "f(x, y) = sqrt (x * x + y * y)";
				f = f4;
				break;
			case 5:
				f_name = "f(x, y) = x * x + y * y";
				f = f5;
				break;
			case 6:
				f_name = "f(x, y) = exp (x * x - y * y)";
				f = f6;
				break;
			case 7:
				f_name = "f(x, y) = 1 / (25 * (x * x + y * y) + 1)";
				f = f7;
				break;
		}
}

void Window::synhronize_param ()
{
  a = tr_gui->a; b = tr_gui->b; c = tr_gui->c; d = tr_gui->d; 
	nx = tr_gui->nx; ny = tr_gui->ny;
	mx = tr_gui->mx; my = tr_gui->my;
	k = tr_gui->k;
	eps = tr_gui->eps;
	mi = tr_gui->mi;
	task = tr_gui->task;
	r1 = tr_gui->r1; r2 = tr_gui->r2; r3 = tr_gui->r3; r4 = tr_gui->r4;
	t1 = tr_gui->t1; t2 = tr_gui->t2;
	it = tr_gui->it;
	p = tr_gui->p;
	pxy = tr_gui->pxy;
	func_id = tr_gui->k;
	switch (k)
		{
			case 0:
				f_name = "f(x, y) = 1";
				f = f0;
				break;
			case 1:
				f_name = "f(x, y) = x";
				f = f1;
				break;
			case 2:
				f_name = "f(x, y) = y";
				f = f2;
				break;
			case 3:
				f_name = "f(x, y) = x + y";
				f = f3;
				break;
			case 4:
				f_name = "f(x, y) = sqrt (x * x + y * y)";
				f = f4;
				break;
			case 5:
				f_name = "f(x, y) = x * x + y * y";
				f = f5;
				break;
			case 6:
				f_name = "f(x, y) = exp (x * x - y * y)";
				f = f6;
				break;
			case 7:
				f_name = "f(x, y) = 1 / (25 * (x * x + y * y) + 1)";
				f = f7;
				break;
		}
}


void Window::synhronize_tr ()
{
  //tr_gui->a = a; tr_gui->b = b; tr_gui->c = c; tr_gui->d = d; 
	tr_gui->nx = nx; tr_gui->ny = ny;
	tr_gui->mx = mx; tr_gui->my = my;
	tr_gui->k = k;
	tr_gui->eps = eps;
	tr_gui->mi = mi;
	tr_gui->r1 = r1; tr_gui->r2 = r2; tr_gui->r3 = r3; tr_gui->r4 = r4;
	tr_gui->t1 = t1; tr_gui->t2 = t2;
	tr_gui->it = it;
	tr_gui->p = p;
	tr_gui->pxy = pxy;
	tr_gui->fmax = fmax;
}


void Window::zoom_out ()
{
	double left_g = 0;
	double right_g = 0;
	double left_v = 0, right_v = 0;
	if (is_timer == 0)
		{
			s--;
			left_g = (3 * a - b) / 2;
			right_g = (3 * b - a) / 2;
			tr_calc->a = left_g;
			tr_calc->b = right_g;
			left_v = (3 * c - d) / 2;
			right_v = (3 * d - c) / 2;
			tr_calc->c = left_v;
			tr_calc->d = right_v;
			tr_calc->state = 1;
			is_timer = 1;
			pthread_cond_broadcast(cond);
			QTimer::singleShot(DEFAULT_TIMER, this, &Window::timer);
		}
}


void Window::zoom_in ()
{
	double left_g = 0, right_g = 0;
	double left_v = 0, right_v = 0;
	if (is_timer == 0)
		{
			left_g = (3 * a + b) / 4;
			right_g = (3 * b + a) / 4;
			if (fabs (left_g - right_g) < 1.0e-06)
				return;
			
			left_v = (3 * c + d) / 4;
			right_v = (3 * d + c) / 4;
			if (fabs (left_v - right_v) < 1.0e-06)
				return;
			
			tr_calc->a = left_g;
			tr_calc->b = right_g;
			tr_calc->c = left_v;
			tr_calc->d = right_v;
			
			s++;
			//synhronize_tr ();
			tr_calc->state = 1;
			is_timer = 1;
			pthread_cond_broadcast(cond);
			QTimer::singleShot(DEFAULT_TIMER, this, &Window::timer);
		}
}


void Window::increase_nx_ny ()
{
	int N = 0;
	int r = nx, s = ny;
	if (is_timer == 0)
		{
			r *= 2;
			s *= 2;
			if (tr_calc->x)
				{
					//printf ("Hi  200\n");
					delete [] tr_calc->x;
					tr_calc->x = nullptr;
				}
			N = (r + 1) * (s + 1);
			tr_calc->N = N;
			tr_calc->nx = r;
			tr_calc->ny = s;
			//printf ("Hi  100\n");
			tr_calc->x = new double [N];
			memset (tr_calc->x, 0, N * sizeof (double));
			//synhronize_tr ();
			tr_calc->state = 1;
			is_timer = 2;
			pthread_cond_broadcast(cond);
			QTimer::singleShot(DEFAULT_TIMER, this, &Window::timer);
		}
}

void Window::reduce_nx_ny ()
{
	int N = 0;
	int r = nx, s = ny;
	if (is_timer == 0)
		{
			r /= 2;
			s /= 2;
			if (r < 2 || s < 2)
				{
					return;
				}
			if (tr_calc->x)
				{
					//printf ("Hi  200\n");
					delete [] tr_calc->x;
					tr_calc->x = nullptr;
				}
			N = (r + 1) * (s + 1);
			tr_calc->N = N;
			tr_calc->nx = r;
			tr_calc->ny = s;
			//printf ("Hi  100\n");
			tr_calc->x = new double [N];
			memset (tr_calc->x, 0, N * sizeof (double));
			//synhronize_tr ();
			tr_calc->state = 1;
			is_timer = 2;
			pthread_cond_broadcast(cond);
			QTimer::singleShot(DEFAULT_TIMER, this, &Window::timer);
			//printf ("Hi ========= %10.3e %d %d %d\n", tr_calc->x[0], tr_calc->nx, tr_calc->ny, tr_calc->N);
		}
}



void Window::change_func ()
{
  if (is_timer == 0)
  	{
  		//k = (k + 1) % 8;
  		tr_calc->k = (tr_calc->k + 1) % 8;
			//synhronize_tr ();
			tr_calc->state = 1;
			is_timer = 1;
			pthread_cond_broadcast(cond);
			QTimer::singleShot(DEFAULT_TIMER, this, &Window::timer);
  	}
}

void Window::increase_p ()
{
	if (is_timer == 0)
  	{
  		tr_calc->pxy++;
  		//if (pxy == 0)
  			//{
  		tr_calc->fmax = fmax;
  			//}
  		//synhronize_tr ();
  		tr_calc->state = 1;
			is_timer = 1;
			pthread_cond_broadcast(cond);
			QTimer::singleShot(DEFAULT_TIMER, this, &Window::timer);
  	}
}
void Window::reduce_p ()
{
	if (is_timer == 0)
  	{
  		tr_calc->pxy--;
  		//if (pxy == 0)
  			//{
  		tr_calc->fmax = fmax;
  			//}
  		//synhronize_tr ();
  		tr_calc->state = 1;
			is_timer = 1;
			pthread_cond_broadcast(cond);
			QTimer::singleShot(DEFAULT_TIMER, this, &Window::timer);
  	}
}


void Window::timer ()
{
	//printf ("Hi___ =========\n");
	if (tr_calc->state == 1)
		{
			QTimer::singleShot(DEFAULT_TIMER, this, &Window::timer);
		}
	else if (tr_calc->state == 0)
		{
			//printf ("Hi1 =========\n");
			calc2gui (tr_calc, tr_gui);
			//printf ("Hi2 =========\n");
			synhronize_param ();
			//printf ("Hi3 =========\n");
			if (is_timer == 2)
				{
					if (tr_gui->x)
						{
							//printf ("Hi  201\n");
							delete [] tr_gui->x;
							tr_gui->x = nullptr;
						}
					//printf ("Hi  101\n");
					tr_gui->x = new double [tr_gui->N];
					memset (tr_gui->x, 0, tr_gui->N * sizeof (double));
				}
			calc2gui_mem (tr_calc, tr_gui);
			x = tr_gui->x;
			is_timer = 0;
			update ();
		}
	else if (tr_calc->state == -1)
		{
			tr_calc->state = -2;
			pthread_cond_broadcast(cond);
			QApplication::quit();
		}
}


void Window::increase_mx_my ()
{
	mx *= 2;
	my *= 2;
        tr_gui->mx = mx;
        tr_gui->my = my;
        tr_calc->mx = mx;
        tr_calc->my = my;
	update ();
}

void Window::reduce_mx_my ()
{
	mx /= 2;
	my /= 2;
	if (mx < 2 || my < 2)
		{
			mx *= 2;
			my *= 2;
		}
        tr_gui->mx = mx;
        tr_gui->my = my;
        tr_calc->mx = mx;
        tr_calc->my = my;
	update ();
}

void Window::change_mode ()
{
	mode = (mode + 1) % 3;
	update ();
}


void Window::is_close ()
{
	if (is_timer == 0)
  	{
  		par->close();
  	}
}


void Window::delete_param ()
{
	if (tr_calc->x != nullptr)
		{
			//printf ("Hi  200\n");
			delete [] tr_calc->x;
			tr_calc->x = nullptr;
		}
	
	if (tr_gui->x != nullptr)
		{
			//printf ("Hi  201\n");
			delete [] tr_gui->x;
			tr_gui->x = nullptr;
		}
	tr_gui = nullptr;
	tr_calc = nullptr;
	data = nullptr;
	//x = nullptr;
	mut = nullptr;
	cond = nullptr;
	arg = nullptr;
}


QPointF Window::l2g (double x_loc, double y_loc)
{
  double x_gl = ((x_loc - a) / (b - a)) * (width () - 1);
  double y_gl = ((d - y_loc) / (d - c)) * (height () - 1);
  return QPointF (x_gl, y_gl);
}

void Window::draw_tr (QPainter & painter, QPointF qp1, QPointF qp2, QPointF qp3, QColor qc)
{
	QPolygonF polygon;
	polygon << qp1 << qp2 << qp3;
	painter.setBrush(qc);
	painter.drawPolygon(polygon);
} 

void Window::draw_graph (QPainter & painter, double min_z, double max_z)
{
	int M = mx * my;
	int l = 0; int i = 0; int j = 0;
	double x_loc1 = 0, y_loc1 = 0;
	double x_loc2 = 0, y_loc2 = 0;
	double x_loc3 = 0, y_loc3 = 0;
	double var;
	double hx = (b - a) / mx;
	double hy = (d - c) / my;
	double nhx = (b - a) / nx;
	double nhy = (d - c) / ny;
	double alpha = 0;
	int red = 50;
	int green = 50;
	int blue = 0;
	if (mode == 0)
		{
			for (l = 0; l < M; l++)
				{
					j = l / mx;
					i = l - j * mx;
					x_loc1 = (i * hx) + a + (hx / 3); 
					y_loc1 = d - (j * hy) - (hy / 3);
					var = f (x_loc1, y_loc1);
					f2alpha (var, alpha, min_z, max_z);
					x_loc1 = (i * hx) + a;
					y_loc1 = d - (j * hy);
					x_loc2 = (i * hx) + a + hx;
					y_loc2 = d - (j * hy);
					x_loc3 = (i * hx) + a;
					y_loc3 = d - (j * hy) - hy;
					blue = alpha * 180 + 50;
					draw_tr (painter, L2G(x_loc1, y_loc1), L2G(x_loc2, y_loc2), L2G(x_loc3, y_loc3), QColor(red, green, blue));
					
					x_loc1 = (i * hx) + a + (2 * hx / 3); 
					y_loc1 = d - (j * hy) - (2 * hy / 3);
					var = f (x_loc1, y_loc1);
					f2alpha (var, alpha, min_z, max_z);
					x_loc1 = (i * hx) + a;
					y_loc1 = d - (j * hy) - hy;
					x_loc2 = (i * hx) + a + hx;
					y_loc2 = d - (j * hy) - hy;
					x_loc3 = (i * hx) + a + hx;
					y_loc3 = d - (j * hy);
					blue = alpha * 180 + 50;
					draw_tr (painter, L2G(x_loc1, y_loc1), L2G(x_loc2, y_loc2), L2G(x_loc3, y_loc3), QColor(red, green, blue));
				}
		}
	else if (mode == 1)
		{
			for (l = 0; l < M; l++)
				{
					j = l / mx;
					i = l - j * mx;
					
					x_loc1 = (i * hx) + a + (hx / 3); 
					y_loc1 = d - (j * hy) - (hy / 3);
					var = get_value_of_approx (x_loc1, y_loc1, nx, ny, nhx, nhy, a, c, x);
					//printf ("%10.3e %10.3e     %10.3e\n", x_loc1, y_loc1, var);
					f2alpha (var, alpha, min_z, max_z);
					x_loc1 = (i * hx) + a;
					y_loc1 = d - (j * hy);
					x_loc2 = (i * hx) + a + hx;
					y_loc2 = d - (j * hy);
					x_loc3 = (i * hx) + a;
					y_loc3 = d - (j * hy) - hy;
					blue = alpha * 180 + 50;
					draw_tr (painter, L2G(x_loc1, y_loc1), L2G(x_loc2, y_loc2), L2G(x_loc3, y_loc3), QColor(red, green, blue));
					
					x_loc1 = (i * hx) + a + (2 * hx / 3); 
					y_loc1 = d - (j * hy) - (2 * hy / 3);
					var = get_value_of_approx (x_loc1, y_loc1, nx, ny, nhx, nhy, a, c, x);
					//printf ("%10.3e %10.3e     %10.3e\n", x_loc1, y_loc1, var);
					f2alpha (var, alpha, min_z, max_z);
					x_loc1 = (i * hx) + a;
					y_loc1 = d - (j * hy) - hy;
					x_loc2 = (i * hx) + a + hx;
					y_loc2 = d - (j * hy) - hy;
					x_loc3 = (i * hx) + a + hx;
					y_loc3 = d - (j * hy);
					blue = alpha * 180 + 50;
					draw_tr (painter, L2G(x_loc1, y_loc1), L2G(x_loc2, y_loc2), L2G(x_loc3, y_loc3), QColor(red, green, blue));
				}
		}
	else if (mode == 2)
		{
			for (l = 0; l < M; l++)
				{
					j = l / mx;
					i = l - j * mx;
					x_loc1 = (i * hx) + a + (hx / 3); 
					y_loc1 = d - (j * hy) - (hy / 3);
					var = fabs (f (x_loc1, y_loc1) - get_value_of_approx (x_loc1, y_loc1, nx, ny, nhx, nhy, a, c, x));
					//printf ("%10.3e %10.3e     %10.3e\n", x_loc1, y_loc1, var);
					f2alpha (var, alpha, min_z, max_z);
					x_loc1 = (i * hx) + a;
					y_loc1 = d - (j * hy);
					x_loc2 = (i * hx) + a + hx;
					y_loc2 = d - (j * hy);
					x_loc3 = (i * hx) + a;
					y_loc3 = d - (j * hy) - hy;
					blue = alpha * 180 + 50;
					draw_tr (painter, L2G(x_loc1, y_loc1), L2G(x_loc2, y_loc2), L2G(x_loc3, y_loc3), QColor(red, green, blue));
					
					x_loc1 = (i * hx) + a + (2 * hx / 3); 
					y_loc1 = d - (j * hy) - (2 * hy / 3);
					var = fabs (f (x_loc1, y_loc1) - get_value_of_approx (x_loc1, y_loc1, nx, ny, nhx, nhy, a, c, x));
					//printf ("%10.3e %10.3e     %10.3e\n", x_loc1, y_loc1, var);
					f2alpha (var, alpha, min_z, max_z);
					x_loc1 = (i * hx) + a;
					y_loc1 = d - (j * hy) - hy;
					x_loc2 = (i * hx) + a + hx;
					y_loc2 = d - (j * hy) - hy;
					x_loc3 = (i * hx) + a + hx;
					y_loc3 = d - (j * hy);
					blue = alpha * 180 + 50;
					draw_tr (painter, L2G(x_loc1, y_loc1), L2G(x_loc2, y_loc2), L2G(x_loc3, y_loc3), QColor(red, green, blue));
				}
		}
	
}

void Window::min_max(double & min_z, double & max_z)
{
	int M = mx * my;
	int l = 0; int i = 0; int j = 0;
	double x_loc = 0, y_loc = 0;
	double var;
	double hx = (b - a) / mx;
	double hy = (d - c) / my;
	double nhx = (b - a) / nx;
	double nhy = (d - c) / ny;
	min_z = 0; max_z = 0;
	//printf ("Mode == %d\n", mode);
	if (mode == 0)
		{
			for (l = 0; l < M; l++)
				{
					j = l / mx;
					i = l - j * mx;
					x_loc = (i * hx) + a + (hx / 3); 
					y_loc = d - (j * hy) - (hy / 3);
					var = f (x_loc, y_loc);
					if (var < min_z)
						min_z = var;
					if (var > max_z)
						max_z = var;

					x_loc = (i * hx) + a + (2 * hx / 3); 
					y_loc = d - (j * hy) - (2 * hy / 3);
					var = f (x_loc, y_loc);
					if (var < min_z)
						min_z = var;
					if (var > max_z)
						max_z = var;
				}
		}
	else if (mode == 1)
		{
			for (l = 0; l < M; l++)
				{
					j = l / mx;
					i = l - j * mx;
					x_loc = (i * hx) + a + (hx / 3); 
					y_loc = d - (j * hy) - (hy / 3);
					var = get_value_of_approx (x_loc, y_loc, nx, ny, nhx, nhy, a, c, x);
					//printf ("%10.3e %10.3e    %10.3e\n", x_loc, y_loc, var);
					if (var < min_z)
						min_z = var;
					if (var > max_z)
						max_z = var;

					x_loc = (i * hx) + a + (2 * hx / 3); 
					y_loc = d - (j * hy) - (2 * hy / 3);
					var = get_value_of_approx (x_loc, y_loc, nx, ny, nhx, nhy, a, c, x);
					//printf ("%10.3e %10.3e    %10.3e\n", x_loc, y_loc, var);
					if (var < min_z)
						min_z = var;
					if (var > max_z)
						max_z = var;
				}
		}
	else if (mode == 2)
		{
			for (l = 0; l < M; l++)
				{
					j = l / mx;
					i = l - j * mx;
					x_loc = (i * hx) + a + (hx / 3); 
					y_loc = d - (j * hy) - (hy / 3);
					var = fabs (f (x_loc, y_loc) - get_value_of_approx (x_loc, y_loc, nx, ny, nhx, nhy, a, c, x));
					if (var < min_z)
						min_z = var;
					if (var > max_z)
						max_z = var;

					x_loc = (i * hx) + a + (2 * hx / 3); 
					y_loc = d - (j * hy) - (2 * hy / 3);
					var = fabs (f (x_loc, y_loc) - get_value_of_approx (x_loc, y_loc, nx, ny, nhx, nhy, a, c, x));
					if (var < min_z)
						min_z = var;
					if (var > max_z)
						max_z = var;
				}
		}
	
	if (fabs (min_z) < fabs (max_z))
		{
			fmax = fabs (max_z);
		}
	else 
		{
			fmax = fabs (min_z);
		}
	
	if (fabs (max_z) < 1.0e-16 && fabs (min_z) < 1.0e-16)
		{
			max_z = 1;
			min_z = -1;
		}
}

void Window::set_text (QPainter & painter, int draw)
{
	QString qstr;
	std::stringstream strm;
	std::string str;
	if (draw == 0)
		{
			painter.drawText (0, 20, f_name);
		}
	else if (draw == 1)
		{
			str = "{|F_min|, |F_max|} = ";
			strm << str << fmax;
			str = strm.str();
			qstr = QString::fromStdString(str);
			painter.drawText (0, 40, qstr);
		}
	else if (draw == 2)
		{
			str = "s = ";
			strm << str << s;
			str = strm.str();
			qstr = QString::fromStdString(str);
			painter.drawText (0, 60, qstr);
		}
	else if (draw == 3)
		{
			str = "nx = ";
			strm << str << nx;
			str = "ny = ";
			strm << " " << str << ny;
			str = strm.str();
			qstr = QString::fromStdString(str);
			painter.drawText (0, 80, qstr);
		}
	else if (draw == 4)
		{
			str = "p = ";
			strm << str << pxy;
			str = strm.str();
			qstr = QString::fromStdString(str);
			painter.drawText (0, 100, qstr);
		}
	else if (draw == 5)
		{
			str = "mx = ";
			strm << str << mx;
			str = "my = ";
			strm << " " << str << my;
			str = strm.str();
			qstr = QString::fromStdString(str);
			painter.drawText (0, 120, qstr);
		}
	else if (draw == 6)
		{
			str = "a = ";
			strm << str << a;
			str = "b = ";
			strm << " " << str << b;
			str = "c = ";
			strm << " " << str << c;
			str = "d = ";
			strm << " " << str << d;
			str = strm.str();
			qstr = QString::fromStdString(str);
			painter.drawText (0, 140, qstr);
		}
	else if (draw == 7)
		{
			str = "mode = ";
			strm << str << mode;
			str = strm.str();
			qstr = QString::fromStdString(str);
			painter.drawText (0, 160, qstr);
		}
}

void Window::f2alpha (double x, double & alpha, double min_z, double max_z)
{
	alpha = (x - min_z) / (max_z - min_z);
}

/// render graph
void Window::paintEvent (QPaintEvent * /* event */)
{ 
	printf ("========================\n");
	printf ("mode = %d\n", mode);
	QPainter painter (this);
	double min_z = 0; double max_z = 0;
	QPen pen_black(Qt::black, 1, Qt::SolidLine);
	QPen pen_green(Qt::green, 3, Qt::SolidLine);
	painter.setPen (pen_black);

	painter.drawLine (L2G(a, c), L2G(b, c));
	painter.drawLine (L2G(b, c), L2G(b, d));
	painter.drawLine (L2G(b, d), L2G(a, d));
	painter.drawLine (L2G(a, d), L2G(a, c));

	min_max(min_z, max_z);
	//printf ("%10.3e %10.3e\n", min_z, max_z);
	draw_graph (painter, min_z, max_z);
	painter.setPen (pen_green);
	printf ("{|F_min|, |F_max|} = %10.3e\n", fmax);
	
	printf (
	"%s : Task = %d R1 = %10.3e R2 = %10.3e R3 = %10.3e R4 = %10.3e T1 = %.2f T2 = %.2f It = %d E = %10.3e K = %d Nx = %d Ny = %d P = %d\n",
	arg, task, r1, r2, r3, r4, t1, t2, it, eps, k, nx, ny, p);
	
	//printf ("a = %10.3e b = %10.3e c = %10.3e d = %10.3e\n", a, b, c, d);
	set_text (painter, 0);
	set_text (painter, 1);
	set_text (painter, 2);
	set_text (painter, 3);
	set_text (painter, 4);
	set_text (painter, 5);
	set_text (painter, 6);
	set_text (painter, 7);
}




