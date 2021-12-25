#pragma once
#include <iostream>
#include <fstream>
#include <cmath>

struct local_area //локальные области
{
    double h_x, h_y, h_z; // Длины сторон
    double lambda, sigma, hi; // lambda и gamma в области
    int mas[8]; // Номера узлов
};

struct cross // Координаты узлов
{
    double x_y_z[3]; // x, y, z соотв.                    
};

struct bound // Граница
{
    int mas[2]; // Номера узлов 
};
struct time_grid
{
    double start;	// Начальное время
    double end;	// Конечное время
    int nSteps;	// Количество шагов по времени
    double h0;
    double q;	// Множитель для подсчета следующих шагов
};

//time_grid* TIME_GRID;	// сетка по времени

struct initial_data // Исходные данные
{
    cross* nodes; // Массив координат узлов
    local_area* locals; // Массив локальных областей
    bound* bounds; // Массив границ
                            
    int num_nodes; // Количество узлов
    int num_locals; // Количество  локальных элементов
    int num_bounds; // Количество границ
    double** global_matrix; // Глобальная матрица
    double** global_M; // Глобальная M
    double** global_G; // Глобальная G
    double** global_A; // Глобальная A
    double* global_b;
    double* global_d;
    double* global_vector; // Глобальный вектор
    double* q; // Решение в узлах (вектор весов)
    time_grid* time_g;
};



void read(initial_data* form); // Чтение исходных данных
void calc_global_matrix(initial_data* form); // Вычисление глобальной матрицы

void calc_global_M(initial_data* form);
void calc_global_G(initial_data* form);
void calc_global_A(initial_data* form, double t_j, double t_j1, double t_j2, double t_j3);

void calc_global_F(initial_data* form, double time);	 // Вычисление глобального вектора правой части
void calc_global_d(initial_data* form, double t_j, double t_j1, double t_j2, double t_j3, double* q_3, double* q_2, double* q_1);

void calc_first_boundary_conditions(initial_data* form, double time); // Учет первых краевых условий
double calc_u(cross coord, double time); // Решение уравнения
double calc_f(cross coord, double time);
void write(initial_data form, double time); // Вывод
void write_result(initial_data form, double time);
void write_matrix(initial_data form); // Вывод матрицы
void write_vector(initial_data form); // Вывод вектора
void mult_matr_by_vect(int size, double** matrix, double* vector, double* result); // Умножение матрицы на вектор
double get_u(cross point, initial_data form); // Получение значения в произвольной точке
double scalar(int size, double* v1, double* v2); // Скалярное умножение векторов
void LOC(initial_data* form); // Локально-Оптимальная Схема