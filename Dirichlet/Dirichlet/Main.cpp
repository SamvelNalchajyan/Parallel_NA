#include <iostream>
#include <iomanip>

class heat_task {
public:
	double X; // ширина пластины
	double Y; // длина пластины
	int n; // размер сетки по x
	int m; // размер сетки по y
	double left_condition(double y); // функция, задающая граничное условие при x = 0
	double right_condition(double y); // функция, задающая граничное условие при x = X
	double bottom_condition(double x); // функция, задающая граничное условие при y = 0
	double top_condition(double x); // функция, задающая граничное условие при y = Y
	double f(double x, double y); // функция, задающая внешнее воздействие
};

double heat_task::left_condition(double y) {
	return 1.0;
	//return 1.5 - 2 * (y - Y / 2.0) * (y - Y / 2.0);
}
double heat_task::right_condition(double y) {
	return 1.0;
	//return 1.5 - 2 * (y - Y / 2.0) * (y - Y / 2.0);
}
double heat_task::bottom_condition(double x) {
	return 1.0;
	//return 1.5 - 2 * (x - X / 2.0) * (x - X / 2.0);
}
double heat_task::top_condition(double x) {
	return 1.0;
	//return 1.5 - 2 * (x - X / 2.0) * (x - X / 2.0);
}
double heat_task::f(double x, double y) {
	return -(x - X / 2.0) * (x - X / 2.0) - (y - Y / 2.0) * (y - Y / 2.0) + 0.0;
	//return -10.0;
	//return 0.0;
}

// base
void heat_dirichlet_sor_(heat_task task, double* v) {
	int n = task.n;
	int m = task.m;
	double h_x = task.X / n;
	double h_y = task.Y / m;
	double* f = new double[(n + 1) * (m + 1)];
	for (int i = 0; i <= n; ++i) {
		for (int j = 0; j <= m; ++j) {
			f[i * (m + 1) + j] = task.f(i * h_x, j * h_y);
			v[i * (m + 1) + j] = 0.0;
			if (i == 0) {
				v[i * (m + 1) + j] = task.left_condition(j * h_y);
			}
			if (i == n) {
				v[i * (m + 1) + j] = task.right_condition(j * h_y);
			}
			if (j == 0) {
				v[i * (m + 1) + j] = task.bottom_condition(i * h_x);
			}
			if (j == m) {
				v[i * (m + 1) + j] = task.top_condition(i * h_x);
			}
		}
	}
	double D = 2.0 * (1.0 / (h_x * h_x) + 1.0 / (h_y * h_y));
	//double omega = 2.0 / (1.0 + 2.0 * sin(0.5 * 3.141592653589793 * h_x));
	//double omega = 2.0 / (1.0 + 2.0 * sin(0.25 * 3.141592653589793 * (h_x + h_y)));
	//double omega = 2.0 / (1.0 + 2.0 * sin(0.5 * 3.141592653589793 * sqrt(h_x * h_y)));
	double omega = 2.0 - 2.0 * sqrt(h_x * h_y);
	double h_x_2 = h_x * h_x;
	double h_y_2 = h_y * h_y;
	double r_max;
	do {
		r_max = 0.0;
		double v_prev, v_curr;
		double r_curr;
		for (int i = 1; i < n; ++i) {
			for (int j = 1; j < m; ++j) {
				v_prev = v[i * (m + 1) + j];

				v_curr = -omega * ((v[(i - 1) * (m + 1) + j] + v[(i + 1) * (m + 1) + j]) / h_x_2 +
					(v[i * (m + 1) + j - 1] + v[i * (m + 1) + j + 1]) / h_y_2);
				v_curr -= (1.0 - omega) * D * v[i * (m + 1) + j] - omega * f[i * (m + 1) + j];
				v_curr /= -D;

				v[i * (m + 1) + j] = v_curr;
				r_curr = fabs(v_curr - v_prev);
				if (r_curr > r_max) {
					r_max = r_curr;
				}
			}
		}
	} while (r_max > 0.00000001);
	delete[] f;
}

void heat_dirichlet_sor(heat_task task, double* v) {
	int n = task.n;
	int m = task.m;
	double h_x = task.X / n;
	double h_y = task.Y / m;
	double* f = new double[(n + 1) * (m + 1)];
	for (int i = 0; i <= n; ++i) {
		for (int j = 0; j <= m; ++j) {
			f[i * (m + 1) + j] = task.f(i * h_x, j * h_y);
			v[i * (m + 1) + j] = 0.0;
			if (i == 0) {
				v[i * (m + 1) + j] = task.left_condition(j * h_y);
			}
			if (i == n) {
				v[i * (m + 1) + j] = task.right_condition(j * h_y);
			}
			if (j == 0) {
				v[i * (m + 1) + j] = task.bottom_condition(i * h_x);
			}
			if (j == m) {
				v[i * (m + 1) + j] = task.top_condition(i * h_x);
			}
		}
	}
	double D = 2.0 * (1.0 / (h_x * h_x) + 1.0 / (h_y * h_y));
	double omega = 2.0 - 2.0 * sqrt(h_x * h_y);
	double h_x_2 = h_x * h_x;
	double h_y_2 = h_y * h_y;
	double r_max;
	int n_1 = n + 1;
	int m_1 = m + 1;
	do {
		r_max = 0.0;
		double v_prev, v_curr;
		double r_curr;
		if (n_1 <= m_1) {
//#pragma omp parallel for
			for (int k = 1; k <= n_1 - 2; ++k) {
//#pragma omp parallel for
				for (int p = 0; p < k; ++p) {
					int i = p + 1;
					int j = k - p;
					v_prev = v[i * (m_1)+j];
					v_curr = -omega * ((v[(i - 1) * (m_1)+j] + v[(i + 1) * (m_1)+j]) / h_x_2 +
						(v[i * (m_1)+j - 1] + v[i * (m_1)+j + 1]) / h_y_2);
					v_curr -= (1.0 - omega) * D * v[i * (m_1)+j] - omega * f[i * (m_1)+j];
					v_curr /= -D;
					v[i * (m_1)+j] = v_curr;
					r_curr = fabs(v_curr - v_prev);
					if (r_curr > r_max) {
						r_max = r_curr;
					}
				}
			}
//#pragma omp parallel for
			for (int k = 0; k < m - n; ++k) {
//#pragma omp parallel for
				for (int p = 0; p < n_1 - 2; ++p) {
					int i = p + 1;
					int j = n_1 - 1 + k - p;
					v_prev = v[i * (m_1)+j];
					v_curr = -omega * ((v[(i - 1) * (m_1)+j] + v[(i + 1) * (m_1)+j]) / h_x_2 +
						(v[i * (m_1)+j - 1] + v[i * (m_1)+j + 1]) / h_y_2);
					v_curr -= (1.0 - omega) * D * v[i * (m_1)+j] - omega * f[i * (m_1)+j];
					v_curr /= -D;
					v[i * (m_1)+j] = v_curr;
					r_curr = fabs(v_curr - v_prev);
					if (r_curr > r_max) {
						r_max = r_curr;
					}
				}
			}
//#pragma omp parallel for
			for (int k = n_1 - 3; k >= 1; --k) {
//#pragma omp parallel for
				for (int p = 0; p < k; ++p) {
					int i = n_1 - k + p - 1;
					int j = m_1 - 2 - p;
					v_prev = v[i * (m_1)+j];
					v_curr = -omega * ((v[(i - 1) * (m_1)+j] + v[(i + 1) * (m_1)+j]) / h_x_2 +
						(v[i * (m_1)+j - 1] + v[i * (m_1)+j + 1]) / h_y_2);
					v_curr -= (1.0 - omega) * D * v[i * (m_1)+j] - omega * f[i * (m_1)+j];
					v_curr /= -D;
					v[i * (m_1)+j] = v_curr;
					r_curr = fabs(v_curr - v_prev);
					if (r_curr > r_max) {
						r_max = r_curr;
					}
				}
			}
		}
		else
		{
//#pragma omp parallel for
			for (int k = 1; k <= m_1 - 2; ++k) {
//#pragma omp parallel for
				for (int p = 0; p < k; ++p) {
					int i = p + 1;
					int j = k - p;
					v_prev = v[i * (m_1)+j];
					v_curr = -omega * ((v[(i - 1) * (m_1)+j] + v[(i + 1) * (m_1)+j]) / h_x_2 +
						(v[i * (m_1)+j - 1] + v[i * (m_1)+j + 1]) / h_y_2);
					v_curr -= (1.0 - omega) * D * v[i * (m_1)+j] - omega * f[i * (m_1)+j];
					v_curr /= -D;
					v[i * (m_1)+j] = v_curr;
					r_curr = fabs(v_curr - v_prev);
					if (r_curr > r_max) {
						r_max = r_curr;
					}
				}
			}
//#pragma omp parallel for
			for (int k = 0; k < n - m; ++k) {
//#pragma omp parallel for
				for (int p = 0; p < m_1 - 2; ++p) {
					int i = p + 2 + k;
					int j = m_1 - 2 - p;
					v_prev = v[i * (m_1)+j];
					v_curr = -omega * ((v[(i - 1) * (m_1)+j] + v[(i + 1) * (m_1)+j]) / h_x_2 +
						(v[i * (m_1)+j - 1] + v[i * (m_1)+j + 1]) / h_y_2);
					v_curr -= (1.0 - omega) * D * v[i * (m_1)+j] - omega * f[i * (m_1)+j];
					v_curr /= -D;
					v[i * (m_1)+j] = v_curr;
					r_curr = fabs(v_curr - v_prev);
					if (r_curr > r_max) {
						r_max = r_curr;
					}
				}
			}
//#pragma omp parallel for
			for (int k = m_1 - 3; k >= 1; --k) {
//#pragma omp parallel for
				for (int p = 0; p < k; ++p) {
					int i = n_1 - k + p - 1;
					int j = m_1 - 2 - p;
					v_prev = v[i * (m_1)+j];
					v_curr = -omega * ((v[(i - 1) * (m_1)+j] + v[(i + 1) * (m_1)+j]) / h_x_2 +
						(v[i * (m_1)+j - 1] + v[i * (m_1)+j + 1]) / h_y_2);
					v_curr -= (1.0 - omega) * D * v[i * (m_1)+j] - omega * f[i * (m_1)+j];
					v_curr /= -D;
					v[i * (m_1)+j] = v_curr;
					r_curr = fabs(v_curr - v_prev);
					if (r_curr > r_max) {
						r_max = r_curr;
					}
				}
			}
		}
	} while (r_max > 0.00000001);
	delete[] f;
}

int main() {
	heat_task task;
	task.X = 5.0;
	task.Y = 5.0;
	task.n = 20;
	task.m = 20;
	double* v = new double[(task.n + 1) * (task.m + 1)];

	heat_dirichlet_sor(task, v);
	for (int i = 0; i <= task.n; ++i) {
		for (int j = 0; j <= task.m; ++j) {
			std::cout << std::setw(9) << v[i * (task.m + 1) + j] << " ";
		}
		std::cout << "\n";
	}

	delete[] v;
	return 0;
}
