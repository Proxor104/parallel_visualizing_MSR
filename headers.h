#ifndef WINDOW_H
#define WINDOW_H

#include <stdio.h>
#include <pthread.h>
#include <iostream>
#include <cmath>
#include <fstream>
#include <string.h>
#include <sys/sysinfo.h>
#include <sys/time.h>
#include <string>
#include <sstream>
#include <memory>

#include <QtWidgets/QApplication>
#include <QtWidgets/QMainWindow>
#include <QtWidgets/QVBoxLayout>
#include <QtWidgets/QAction>
#include <QtWidgets/QMenuBar>
#include <QtWidgets/QMessageBox>
#include <QtWidgets/QtWidgets>
#include <QTimer>

class Data
{
public:
	double * x;
	int * I;
	double * A;
	double * B;
	double * r;
	double * u;
	double * v;
	double * block;
	Data (double * x_ = nullptr, int * I_ = nullptr, double * A_ = nullptr, double * B_ = nullptr, 
				double * r_ = nullptr, double * u_ = nullptr, double * v_ = nullptr, double * block_ = nullptr);
	~Data () = default;
};

class Trade
{
public:
	double * x = nullptr;
	double a = 0; double b = 0; double c = 0; double d = 0; 
	int nx = 0; int ny = 0; 
	int mx = 0; int my = 0;
	int k = 0; 
	double eps = 0; 
	int mi = 0;
	int p = 0;
	int N = 0;
	int task = 2;
	char * arg = nullptr;
	double r1 = 0; double r2 = 0; double r3 = 0; double r4 = 0;
	double t1 = 0; double t2 = 0;
	int it = 0;
	int pxy = 0;
	double fmax = 0;
	int first_step = 0; // первый ли запуск потоков?
	int state = 0; // состояние вычисления:
								 // 0 - вычисления не идут
								 // 1 - вычисления идут
								 // -1 - была ошибка в вычислениях
								 // -2 - надо завершать работу выч. потоков
	Trade (double * x_, double a_, double b_, double c_, double d_, int nx_, int ny_, 
				 int mx_, int my_, int k_, double eps_, int mi_, int p_, int N_, int task_, char * arg_);
	~Trade () = default;
};

class Args
	{
		public:

			pthread_t tid = -1; // id of thread  
  		pthread_cond_t * cond_start = nullptr;
			pthread_mutex_t * mut = nullptr;
  		pthread_cond_t * cond = nullptr;
			//int len_msr = 0; // длина шаблона
			//int N = 0; // длина правой части
			int q = 0; // номер потока
			Trade * tr_gui = nullptr;
			Data * data = nullptr;

			Args() = default;
			Args(const Args & r) = delete;
			Args & operator = (const Args & r) = delete;
			Args(Args && r) = default;
			Args & operator = (Args && r) = default;
			~Args() = default;
	};

// window.cpp
class Window : public QWidget
{
  Q_OBJECT

private:
	
	int mode = 0;	
	int is_timer = 0; // 0 - таймер свободен
										// 1 - таймер занят
	double fmax = 0;
	int s = 0;
	int pxy = 0;

	Trade * tr_gui = nullptr;
	Trade * tr_calc = nullptr;
	Data * data = nullptr;
	pthread_mutex_t * mut = nullptr;
  pthread_cond_t * cond = nullptr;
  
  char * arg = nullptr;
  double * x = nullptr;
	double a = 0; double b = 0; double c = 0; double d = 0;
	int nx = 0; int ny = 0;
	int mx = 0; int my = 0;
	int k = 0;
	int func_id = 0;
	double (*f) (double, double) = nullptr;
  const char *f_name = nullptr;
	double eps = 0;
	int mi = 0;
	int p = 0;
	int task = 2;
	double r1 = 0; double r2 = 0; double r3 = 0; double r4 = 0;
	double t1 = 0; double t2 = 0;
	int it = 0;
    QWidget * par = nullptr;

public:
  Window (QWidget * parent);
  void set_param (Trade * tr_calc_, Trade * tr_gui_, pthread_mutex_t * mut_, pthread_cond_t * cond_);
  void synhronize_tr ();
  void synhronize_param ();
	void timer ();


  QSize minimumSizeHint () const;
  QSize sizeHint () const;

	void draw_tr (QPainter & painter, QPointF qp1, QPointF qp2, QPointF qp3, QColor qc);

  QPointF l2g (double x_loc, double y_loc);
  
  void set_func ();
  
  void draw_graph (QPainter & painter, double min_z, double max_z);
	void min_max(double & min_z, double & max_z);
	void f2alpha (double x, double & alpha, double min_z, double max_z);

	void set_text (QPainter & painter, int draw);

	void delete_param ();

public slots:
  void change_func ();
	void increase_mx_my ();
	void reduce_mx_my ();
	void change_mode ();
	void zoom_in ();
	void zoom_out ();
	void increase_nx_ny ();
	void reduce_nx_ny ();
	void increase_p ();
	void reduce_p ();
        void is_close ();
protected:
  void paintEvent (QPaintEvent *event);
};

// in_out.cpp
int input_args (double & a, double & b, double & c, double & d, int & nx, int & ny, 
                int & mx, int & my, int & k, double & eps, int & mi, int & p, int argc, char * argv[]);

//proc_time.cpp
double get_full_time();
double get_cpu_time();

void calc2gui (Trade * tr_calc, Trade * tr_gui);
void calc2gui_mem (Trade * tr_calc, Trade * tr_gui);

// calculations.cpp
void l2ij (int nx, int /* ny */, int & i, int & j, int l);
void ij2l (int nx, int /* ny */, int i, int j, int & l);
void IA_ij (int nx, int ny, double hx, double hy, int i, int j, int is, int js, int s, int * I, double * A);
int get_off_diag (int nx, int ny, double hx, double hy, int i, int y, int * I = nullptr, double * A = nullptr);
int fill_I (int nx, int ny, int * I);
int get_diag (int nx, int ny, double hx, double hy, int i, int j, int * /* I */, double * A);
int fill_IA (int nx, int ny, double hx, double hy, int * I, double * A, int p, int k);
int fill_B (int nx, int ny, double hx, double hy, double x0, double y0, double (*f) (double, double), double * B, int p, int q, int pxy, double fmax);
double F_ij (int nx, int ny, double hx, double hy, double x0, double y0, double (*f)(double, double), int l);
int init_reduce_sum (int p);
double reduce_sum_det (int p, int k, double s);
int finilize_reduce_sum (int /*p*/);
void reduce_sum (int p, double * a, int n);
void reduce_sum (int p, int * a, int n);
void reduce_sum (int p);
int get_len_msr (int nx, int ny);
int get_len_msr_off_diag (int nx, int ny, double hx, double hy);
void thread_rows (int n, int p, int k, int & i1, int & i2);
double scalar_product (int n, double * x, double * y, int p, int k);
void mult_sub_vector (int n, double * x, double * y, double t, int p, int k);
void matrix_mult_vector_msr (int n, double * A, int * I, double * x, double * y, int p, int k);
void apply_precondition_msr_matrix (int n, double * A, int * I, double * v, double * r, double * block, int p, int k); // Mv = r
double get_value_pf_center (double x, double y, int nx, int ny, double hx, double hy, double a, double c, double * vec, int up_down);
double get_value_pf_corr (double x, double y, int nx, int ny, double hx, double hy, double a, double c, double * vec);
double get_value_of_approx (double x, double y, int nx, int ny, double hx, double hy, double a, double c, double * vec);
double reduce_max_det (int p, int k, double m);
void * thread_main(void * ptr);
int minimal_residual_msr_matrix (int n, double * A, int * I, double * b, double * x, double * r, double * u, 
                  double * v, double * block, double eps, int max_it, int p, int k);
int minimal_residual_msr_matrix_full (int n, double * A, int * I, double * b, double * x, double * r, double * u, 
                                     double * v, double * block, double eps, int max_it, int max_steps, int p, int k);
void l2ij_middle (int nx, int /* ny*/ , int l, int & i, int & j);
void thread_cords_r (int n, int p, int k, int & i1, int & i2);
void get_value_r_1_2 (double * x, double (*f) (double, double), int nx, int ny, double hx, double hy, double & r1, double & r2, double a, double c, int p, int k);
void get_value_r_3_4 (double * x, double (*f) (double, double), int nx, int ny, double hx, double hy, double & r3, double & r4, double a, double c, int p, int k, int pxy, double fmax);

#endif
