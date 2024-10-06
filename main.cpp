#include "headers.h"

int main (int argc, char *argv[])
{
	QApplication app (argc, argv);
	static pthread_mutex_t mut_start = PTHREAD_MUTEX_INITIALIZER;
  static pthread_cond_t cond_start = PTHREAD_COND_INITIALIZER;
  static pthread_mutex_t mut = PTHREAD_MUTEX_INITIALIZER;
  static pthread_cond_t cond = PTHREAD_COND_INITIALIZER;

	double * x_tr = nullptr;
	double * x_gui = nullptr;
	Args * Arg = nullptr;
	Trade * tr_calc = nullptr;
	Trade * tr_gui = nullptr;
	Data * data = nullptr;

	int task = 2; // номер метода
	double a = 0, b = 0, c = 0, d = 0; // границы прямоугольника
	int nx = 0, ny = 0; // количество точек разбиения интерполяции
	int mx = 0, my = 0; // количество точек разбиения на экране
	int k = 0; // номер формулы функции f
	double eps = 0; // точность решения системы линейных уравнений
	int mi = 0; // максимальное количество интераций итерационного метода
	int p = 0; // количество потоков

	double r1 = -1, r2 = -1, r3 = -1, r4 = -1; // невязки
	double t1 = -1, t2 = -1; // время
	int it = -1; // количество итераций
	
	if (input_args (a, b, c, d, nx, ny, mx, my, k, eps, mi, p, argc, argv) < 0)
		{
                        printf ("a = %10.3e b = %10.3e c = %10.3e d = %10.3e\n", a, b, c, d);
                        printf (
			"%s : Task = %d R1 = %e R2 = %e R3 = %e R4 = %e T1 = %.2f T2 = %.2f\
			It = %d E = %e K = %d Nx = %d Ny = %d P = %d\n",
			argv[0], task, r1, r2, r3, r4, t1, t2, it, eps, k, nx, ny, p);
			return -1;
		}

           //printf ("%10.3e %10.3e %lf %lf %d %d %d %d %d %10.3e %d %d\n", a, b, c, d, nx, ny, mx, my, k, eps, mi, p);
	//int len_msr = get_len_msr (nx, ny) + 1 + (nx + 1) * (ny + 1);
	int N = (nx + 1) * (ny + 1);
	int q = 0;

	//printf ("Hi  100\n");
	x_tr = new double [N];
	memset (x_tr, 0, N * sizeof (double));
	//printf ("Hi  101\n");
	x_gui = new double [N];
	memset (x_gui, 0, N * sizeof (double));

	//printf ("Hi  101\n");
	Arg = new Args[p];
	//printf ("Hi  102\n");
	tr_calc = new Trade (x_tr, a, b, c, d, nx, ny, mx, my, k, eps, mi, p, N, task, argv[0]);
	//printf ("Hi  103\n");
	tr_gui = new Trade (x_gui, a, b, c, d, nx, ny, mx, my, k, eps, mi, p, N, task, argv[0]);
	//printf ("Hi  104\n");
	data = new Data ();

	for(q = 0; q < p; q++)
		{
			//Arg[q].len_msr = len_msr;
			//Arg[q].N = N;
			Arg[q].q = q;
			Arg[q].tr_gui = tr_calc;
			Arg[q].cond_start = &cond_start;
			Arg[q].mut = &mut;
			Arg[q].cond = &cond;
			Arg[q].data = data;
 		}
	for(q = 0; q < p; q++)
    {
      if(pthread_create(&Arg[q].tid, 0, thread_main, Arg + q))
        {
          printf("Usage: can't create thread %d\n", q);
					printf (
					"%s : Task = %d R1 = %e R2 = %e R3 = %e R4 = %e T1 = %.2f T2 = %.2f\
					It = %d E = %e K = %d Nx = %d Ny = %d P = %d\n",
					argv[0], task, r1, r2, r3, r4, t1, t2, it, eps, k, nx, ny, p);
          abort();
        }
    }	
	
	pthread_mutex_lock(&mut_start);
	while (tr_calc->first_step == 1)
		{
			 pthread_cond_wait(&cond_start, &mut_start);
		}
	pthread_mutex_unlock(&mut_start);
	
	calc2gui (tr_calc, tr_gui);
	calc2gui_mem (tr_calc, tr_gui);
	
  //printf ("Hi  1\n");
  QMainWindow *window = new QMainWindow;
  //printf ("Hi  1\n");
  QMenuBar *tool_bar = new QMenuBar (window);
  //printf ("Hi  1\n");
  Window *graph_area = new Window (window);
  
  graph_area->set_param (tr_calc, tr_gui, &mut, &cond);
  QAction *action;

  action = tool_bar->addAction ("&Change function", graph_area, SLOT (change_func ()));
  action->setShortcut (QString ("0"));
  
  action = tool_bar->addAction ("&Change mode", graph_area, SLOT (change_mode ()));
  action->setShortcut (QString ("1")); 
  
  action = tool_bar->addAction ("&Zoom in", graph_area, SLOT (zoom_in ()));
  action->setShortcut (QString ("2")); 
  
  action = tool_bar->addAction ("&Zoom out", graph_area, SLOT (zoom_out ()));
  action->setShortcut (QString ("3"));
  
  action = tool_bar->addAction ("&Increase nx, ny", graph_area, SLOT (increase_nx_ny ()));
  action->setShortcut (QString ("4"));
  
  action = tool_bar->addAction ("&Reduce nx, ny", graph_area, SLOT (reduce_nx_ny ()));
  action->setShortcut (QString ("5"));
 
  action = tool_bar->addAction ("&p++", graph_area, SLOT (increase_p ()));
  action->setShortcut (QString ("6"));
  
  action = tool_bar->addAction ("&p--", graph_area, SLOT (reduce_p ()));
  action->setShortcut (QString ("7"));
  
  action = tool_bar->addAction ("&Increase mx, my", graph_area, SLOT (increase_mx_my ()));
  action->setShortcut (QString ("8")); 
  
  action = tool_bar->addAction ("&Reduce mx, my", graph_area, SLOT (reduce_mx_my ()));
  action->setShortcut (QString ("9"));
  
  action = tool_bar->addAction ("E&xit", graph_area, SLOT (is_close ()));
  action->setShortcut (QString ("Ctrl+X"));

  tool_bar->setMaximumHeight (30);

  window->setMenuBar (tool_bar);
  window->setCentralWidget (graph_area);
  window->setWindowTitle ("Graph");

  window->show ();
  app.exec ();
  graph_area->delete_param ();
 	//printf ("Hi  2\n");
  delete tool_bar;
  //printf ("Hi  3\n");
  delete graph_area;
  //printf ("Hi  4\n");
  delete window;

	tr_calc->state = -2;
	pthread_cond_broadcast(&cond);

//printf ("AAAAAAAAAAAAAAAAAAAAAAA\n");

  for(q = 0; q < p; q++)
    pthread_join(Arg[q].tid, 0); 

//printf ("AAAAAAAAAAAAAAAAAAAAAAA\n");

/**/
	if (tr_calc)
		{
			//printf ("Hi  202\n");
			delete tr_calc;
		}
	if (tr_gui)
		{
			//printf ("Hi  202\n");
			delete tr_gui;
		}
	if (data)
		{
			//printf ("Hi  203\n");
			delete data;
		}

/*
	if (x_tr)
		{
			printf ("Hi  200\n");
			delete [] x_tr;
			x_tr = nullptr;
		}
*/
	//printf ("Hi  201\n");
  delete [] Arg;
 	
  return 0;
}
