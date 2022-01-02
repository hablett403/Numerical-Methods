#include <iostream>
#include <vector>
#include <cmath>


std::vector<double> gauss_jacobi_3d(std::vector<std::vector<double>> A, std::vector<double> b){
    int iteration = 0;
    double x = 2, y = 2, z = 2, x_new, y_new, z_new;
    while(iteration != 100){
        ++iteration;
        x_new = (b[0] - A[0][1] * y - A[0][2] * z)/A[0][0];
        y_new = (b[1] - A[1][0] * x - A[1][2] * z)/A[1][1];
        z_new = (b[2] - A[2][0] * x - A[2][1] * y)/A[2][2];
        x = x_new, y = y_new, z = z_new;
    }
    return {x,y,z};
}

std::vector<double> gauss_jacobi(std::vector<std::vector<double>> A, std::vector<double> b){
    int iteration = 0;
    int n = A[0].size();
    std::vector<double> x(n, 0);
    std::vector<double> x_new(n);
    while(iteration != 100){
        ++iteration;
        for(int i = 0; i < n; ++i){

            double s1 = 0, s2 =0;
            for(int j = 0; j<=i-1; ++j) s1 += A[i][j] * x[j];
            for(int j = i+1; j<n; ++j) s2 += A[i][j] * x[j];
            x_new[i] = 1/A[i][i] * (b[i] - s1 -s2);

        }
        x = x_new;
    }
    return x;
}

std::vector<double> gauss_seidel_3d(std::vector<std::vector<double>> A, std::vector<double> b){
    int iteration = 0;
    double x = 2, y = 2, z = 2, x_new, y_new, z_new;
    while(iteration != 100){
        ++iteration;
        x_new = (b[0] - A[0][1] * y - A[0][2] * z)/A[0][0];
        y_new = (b[1] - A[1][0] * x_new - A[1][2] * z)/A[1][1];
        z_new = (b[2] - A[2][0] * x_new - A[2][1] * y_new)/A[2][2];
        x = x_new, y = y_new, z = z_new;
    }
    return {x,y,z};
}

std::vector<double> gauss_seidel(std::vector<std::vector<double>> A, std::vector<double> b){
    int iteration = 0;
    int n = A[0].size();
    std::vector<double> x(n, 0);
    std::vector<double> x_new(n);
    while(iteration != 100){
        ++iteration;
        for(int i = 0; i < n; ++i){

            double s1 = 0, s2 =0;
            for(int j = 0; j<=i-1; ++j) s1 += A[i][j] * x_new[j];
            for(int j = i+1; j<n; ++j) s2 += A[i][j] * x[j];
            x_new[i] = 1/A[i][i] * (b[i] - s1 -s2);

        }
        x = x_new;
    }
    return x;
}

bool is_diagonal_dominant(std::vector<std::vector<double>> A){
    int n = A[0].size();
    for(int i = 0; i < n; ++i){
        //sum up all elements in the row
        double s = 0;
        for(int j = 0; j < n; ++j)s += std::abs(A[i][j]);
        // subtract diagonal element
        s -= std::abs(A[i][i]);

        if(s >= std::abs(A[i][i]))return false;
    }
    return true;
}
