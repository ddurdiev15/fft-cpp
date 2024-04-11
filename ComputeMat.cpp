#include <iostream>
#include <vector>
#include <cstdlib> // For the exit function
#include <cmath>

using namespace std;


class Compute {
    public:

    // Function to print a matrix
    static void PrintMatrix(const vector<vector<float>>& matrix) {

        int rows = matrix.size();
        int cols = matrix[0].size();

        for (int i = 0; i < rows; i++) {
            for (int j = 0; j < cols; j++) {
                cout << matrix[i][j] << " ";
            }
            cout << endl;
        }
    }

    // Return the shape of a matirx, needed for matrix inversion
    static int CheckShapeMatrix(const vector<vector<float>>& matrix){

        int rows = matrix.size();  // rows
        int cols = matrix[0].size();  // columns

        // Check if the matrix is square
        if (rows != cols) {
            cout << "Error: Matrix is NOT square" << endl;
            exit(1);
        }

        return rows;
    }

    // Identity matrix 
    static vector<vector<float>> IdentityMatrix(int N){

        // create a zero matrix
        vector<vector<float>> identityMatrix(N, vector<float>(N, 0.0));

        // create the identity matrix with shape N x N
        for (int i = 0; i < N; i++) {
            identityMatrix[i][i] = 1.0;
        }

        return identityMatrix;
    }

    // Determinant of a matrix  
    static float Determinant(const vector<vector<float>>& matrix){

        int N = CheckShapeMatrix(matrix); // length of the matrix

        // create a zero matrix - upper triangular matrix
        vector<vector<float>> UppMatrix(N, vector<float>(N, 0.0));
        for (int i = 0; i < N; i++){
            for (int j = 0; j < N; j++){
                UppMatrix[i][j] = matrix[i][j];
            }
        }

        // forward elimination 
        for (int k = 0; k < N-1; k++) {
            for (int i = k + 1; i < N; i++) {
                float m = UppMatrix[i][k] / UppMatrix[k][k];
                for (int j = 0; j < N; j++) {
                    UppMatrix[i][j] -= m * UppMatrix[k][j];
                }
            }
        }

        float det = 1;
        for (int i = 0; i < N; i++){
            det *= UppMatrix[i][i];
        }

        // cout << "det: " << det << endl;

        return det;
    }


    // Calculate the inverse of a matrix
    static vector<vector<float>> Inverse(const vector<vector<float>>& matrix){

        int N = CheckShapeMatrix(matrix);
        float detA = Determinant(matrix);
        // Check if detA is zero
        if (detA == 0.0 || isnan(detA)) {
            cout << "Error: Determinant is " << detA << endl;
            exit(1);
        }


        // create a zero matrix
        vector<vector<float>> I = IdentityMatrix(N);

        // create a zero matrix
        vector<vector<float>> UppMatrix(N, vector<float>(N, 0.0));
        for (int i = 0; i < N; i++){
            for (int j = 0; j < N; j++){
                UppMatrix[i][j] = matrix[i][j];
            }
        }

        // start the matrix inversion
        // forward elimination 
        for (int k = 0; k < N-1; k++) {
            for (int i = k + 1; i < N; i++) {
                float m = UppMatrix[i][k] / UppMatrix[k][k];
                for (int j = 0; j < N; j++) {
                    UppMatrix[i][j] -= m * UppMatrix[k][j];
                    I[i][j] -= m * I[k][j];
                }
            }
        }

        // backward elimination
        for (int i = N-1; i >= 0; i--) {
            for (int k = i + 1; k < N; k++) {
                float m = UppMatrix[i][k] / UppMatrix[k][k];
                for (int j = 0; j < N; j++) {
                    UppMatrix[i][j] -= m * UppMatrix[k][j];
                    I[i][j] -= m * I[k][j];
                }
            }
        }

        // create a zero matrix for the inverse
        vector<vector<float>> InvMatrix(N, vector<float>(N, 0.0));
        for (int i = 0; i < N; i++){
            for (int j = 0; j < N; j++){
                InvMatrix[i][j] = I[i][j]/UppMatrix[i][i];
            }
        }

        return InvMatrix;
    }

    // Add
    static vector<vector<float>> Add(const vector<vector<float>>& matrix1, const vector<vector<float>>& matrix2) {
        int rows = matrix1.size();
        int cols = matrix1[0].size();

        // Check if matrices have the same dimensions
        if (rows != matrix2.size() || cols != matrix2[0].size()) {
            cout << "Error: Matrices must have the same dimensions for addition." << endl;
            exit(1);
        }

        // Initialize the result matrix
        vector<vector<float>> result(rows, vector<float>(cols, 0.0));

        for (int i = 0; i < rows; i++) {
            for (int j = 0; j < cols; j++) {
                result[i][j] = matrix1[i][j] + matrix2[i][j];
            }
        }

        return result;
    }

    // Subtract
    static vector<vector<float>> Subtract(const vector<vector<float>>& matrix1, const vector<vector<float>>& matrix2) {
        int rows = matrix1.size();
        int cols = matrix1[0].size();

        // Check if matrices have the same dimensions
        if (rows != matrix2.size() || cols != matrix2[0].size()) {
            cout << "Error: Matrices must have the same dimensions for addition." << endl;
            exit(1);
        }

        // Initialize the result matrix
        vector<vector<float>> result(rows, vector<float>(cols, 0.0));

        for (int i = 0; i < rows; i++) {
            for (int j = 0; j < cols; j++) {
                result[i][j] = matrix1[i][j] - matrix2[i][j];
            }
        }

        return result;
    }

};

int main() {

    // Example matrices
    vector<vector<float>> A = {{1.0, 2.0, 3.0}, {4.0, 5.0, 6.0}, {7.0, 8.0, 9.0}};
    vector<vector<float>> B = {{9.0, 8.0, 7.0}, {6.0, 5.0, 4.0}, {3.0, 2.0, 1.0}};

    Compute ComputeMatrices;

    vector<vector<float>> C = ComputeMatrices.Subtract(A, B);
    ComputeMatrices.PrintMatrix(C);

    return 0;
}