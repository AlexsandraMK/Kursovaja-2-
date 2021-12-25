#pragma once
#include <iostream>
#include <fstream>
#include <cmath>

struct local_area //��������� �������
{
    double h_x, h_y, h_z; // ����� ������
    double lambda, sigma, hi; // lambda � gamma � �������
    int mas[8]; // ������ �����
};

struct cross // ���������� �����
{
    double x_y_z[3]; // x, y, z �����.                    
};

struct bound // �������
{
    int mas[2]; // ������ ����� 
};
struct time_grid
{
    double start;	// ��������� �����
    double end;	// �������� �����
    int nSteps;	// ���������� ����� �� �������
    double h0;
    double q;	// ��������� ��� �������� ��������� �����
};

//time_grid* TIME_GRID;	// ����� �� �������

struct initial_data // �������� ������
{
    cross* nodes; // ������ ��������� �����
    local_area* locals; // ������ ��������� ��������
    bound* bounds; // ������ ������
                            
    int num_nodes; // ���������� �����
    int num_locals; // ����������  ��������� ���������
    int num_bounds; // ���������� ������
    double** global_matrix; // ���������� �������
    double** global_M; // ���������� M
    double** global_G; // ���������� G
    double** global_A; // ���������� A
    double* global_b;
    double* global_d;
    double* global_vector; // ���������� ������
    double* q; // ������� � ����� (������ �����)
    time_grid* time_g;
};



void read(initial_data* form); // ������ �������� ������
void calc_global_matrix(initial_data* form); // ���������� ���������� �������

void calc_global_M(initial_data* form);
void calc_global_G(initial_data* form);
void calc_global_A(initial_data* form, double t_j, double t_j1, double t_j2, double t_j3);

void calc_global_F(initial_data* form, double time);	 // ���������� ����������� ������� ������ �����
void calc_global_d(initial_data* form, double t_j, double t_j1, double t_j2, double t_j3, double* q_3, double* q_2, double* q_1);

void calc_first_boundary_conditions(initial_data* form, double time); // ���� ������ ������� �������
double calc_u(cross coord, double time); // ������� ���������
double calc_f(cross coord, double time);
void write(initial_data form, double time); // �����
void write_result(initial_data form, double time);
void write_matrix(initial_data form); // ����� �������
void write_vector(initial_data form); // ����� �������
void mult_matr_by_vect(int size, double** matrix, double* vector, double* result); // ��������� ������� �� ������
double get_u(cross point, initial_data form); // ��������� �������� � ������������ �����
double scalar(int size, double* v1, double* v2); // ��������� ��������� ��������
void LOC(initial_data* form); // ��������-����������� �����