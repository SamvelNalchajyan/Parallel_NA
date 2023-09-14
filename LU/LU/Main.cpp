#include <iostream>
#include <iomanip>

void LU_Decomposition(double* A, double* L, double* U, int n) {
	for (int i = 0; i < n; ++i) {
#pragma omp parallel for
		for (int j = i; j < n; ++j) {
			if (i != j) {
				L[i * n + j] = 0.0;
				U[j * n + i] = 0.0;
			}
			double sum = 0;
			for (int k = 0; k < i; ++k) {
				sum += L[i * n + k] * U[k * n + j];
			}
			U[i * n + j] = A[i * n + j] - sum;
		}
#pragma omp parallel for
		for (int j = i; j < n; ++j) {
			if (i == j) {
				L[i * n + j] = 1.0;
			}
			else
			{
				double sum = 0;
				for (int k = 0; k < i; ++k) {
					sum += L[j * n + k] * U[k * n + i];
				}
				L[j * n + i] = (A[j * n + i] - sum) / U[i * n + i];
			}
		}
	}
}

int main() {
	//srand(time(nullptr));
	int n = 10;
	double* A;
	double* L;
	double* U;
	double* mult;
	A = new double[n * n];
	L = new double[n * n];
	U = new double[n * n];
	mult = new double[n * n];

	for (int i = 0; i < n; ++i) {
		for (int j = 0; j < n; ++j) {
			A[i * n + j] = double(-100 + rand() % 201);
			L[i * n + j] = 0.0;
			U[i * n + j] = 0.0;
		}
		A[i * n + i]++;
	}
	// print A
	std::cout << "A = \n";
	for (int i = 0; i < n; ++i) {
		for (int j = 0; j < n; ++j) {
			std::cout << std::setw(12) << A[i * n + j] << " ";
		}
		std::cout << "\n";
	}
	std::cout << "\n\n";

	LU_Decomposition(A, L, U, n);

	std::cout << "L = \n";
	for (int i = 0; i < n; ++i) {
		for (int j = 0; j < n; ++j) {
			std::cout << std::setw(12) << L[i * n + j] << " ";
		}
		std::cout << "\n";
	}
	std::cout << "\n\n";
	std::cout << "U = \n";
	for (int i = 0; i < n; ++i) {
		for (int j = 0; j < n; ++j) {
			std::cout << std::setw(12) << U[i * n + j] << " ";
		}
		std::cout << "\n";
	}
	std::cout << "\n\n";

	// mult
	for (int i = 0; i < n; ++i) {
		for (int j = 0; j < n; ++j) {
			mult[i * n + j] = 0.0;
			for (int k = 0; k < n; ++k) {
				mult[i * n + j] += L[i * n + k] * U[k * n + j];
			}
		}
	}
	std::cout << "LU = \n";
	for (int i = 0; i < n; ++i) {
		for (int j = 0; j < n; ++j) {
			std::cout << std::setw(12) << mult[i * n + j] << " ";
		}
		std::cout << "\n";
	}
	std::cout << "\n\n";
	// ------------------------------------------------ //

	delete[] A;
	delete[] L;
	delete[] U;
	delete[] mult;
	return 0;
}
