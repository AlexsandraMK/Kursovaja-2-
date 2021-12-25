#include "Header.h"

using namespace std;

int main()
{
	setlocale(LC_ALL, "Russian");
	initial_data* form = new initial_data;   // Объявили область
	read(form);       // Чтение исходных данных
	//calc_global_matrix(form);	// Вычисление глобальной матрицы
	double time = form->time_g->start;

	double* q0 = new double[form->num_nodes]{};
	for (int i = 0; i < form->num_nodes; i++)
	{
		q0[i] = calc_u(form->nodes[i], time);
		form->q[i] = q0[i];
	}
	write_result(*form, time);
	time += form->time_g->h0;

	double* q1 = new double[form->num_nodes]{};
	for (int i = 0; i < form->num_nodes; i++)
	{
		q1[i] = calc_u(form->nodes[i], time);
		form->q[i] = q1[i];
	}
	write_result(*form, time);
	time += form->time_g->h0 * form->time_g->q;

	double* q2 = new double[form->num_nodes]{};
	for (int i = 0; i < form->num_nodes; i++)
	{
		q2[i] = calc_u(form->nodes[i], time);
		form->q[i] = q2[i];
	}
	write_result(*form, time);
	calc_global_M(form);
	calc_global_G(form);

	double* q0_ = q0;
	double* q1_ = q1;
	double* q2_ = q2;

	double t0 = form->time_g->start;
	double t1 = form->time_g->start + form->time_g->h0;
	double t2 = form->time_g->h0 * form->time_g->q + t1;

	double h = form->time_g->h0 * form->time_g->q * form->time_g->q;
	
	// Поиск векторов весов на каждом шаге
	for (double time = t2 + h; time <= form->time_g->end;
		h *= form->time_g->q,
		t0 = t1, t1 = t2, t2 = time,
		q0_ = q1_,
		q1_ = q2_,
		q2_ = form->q,
		form->q = new double[form->num_nodes]{},
		time += h)
	{

		form->global_b = new double[form->num_nodes]{};
		form->global_vector = new double[form->num_nodes]{};
		calc_global_A(form, time, t2, t1, t0);
		calc_global_F(form, time);	// Вычисление глобального вектора правой части
		calc_global_d(form, time, t2, t1, t0, q0_, q1_, q2_);
		calc_first_boundary_conditions(form, time);	// Учет первых краевых условий
		LOC(form);	// Решение СЛАУ
		write_result(*form, time);
	}

	return 0;
}
