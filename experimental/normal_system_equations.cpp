#include "normal_system_equations.h"
#include <math.h>

normal_system_equations::normal_system_equations()
{

}

QVector<QVector<double> > normal_system_equations::multiplication_matrix(QVector<QVector<double> > matrix1, QVector<QVector<double> > matrix2)
{
    QVector<QVector<double> > result(matrix2.size(),QVector<double> (matrix1[0].size(),0));
    for (int i = 0;i<matrix1[0].size();i++)//i - количество строк внутри столбца первой матрицы
        for (int j = 0;j<matrix2.size();j++)//j - количество столбцов второй матрицы
        {
            for (int k = 0;k<matrix1.size();k++)//k - количество столбцов первой матрицы
                 result[j][i]+=matrix1[k][i]*matrix2[j][k];
        }
    return result;
}

QVector<QVector<double> > normal_system_equations::transposition_matrix(QVector<QVector<double> > matrix)
{
    QVector<double> non;
    QVector<QVector<double> > result(matrix[0].size(),non);
    for (int i = 0;i<matrix.size();i++)
        for (int j = 0;j<matrix[i].size();j++)
            result[j].push_back(matrix[i][j]);
    return result;
}

QVector<QVector<double> > normal_system_equations::reverse_matrix_count(QVector<QVector<double> > matrix)
{
    QVector<QVector<double> > A11(matrix.size()/2,QVector<double>(matrix.size()/2,0)),A12(matrix.size()/2,QVector<double>(matrix.size()/2,0)),A21(matrix.size()/2,QVector<double>(matrix.size()/2,0)),A22(matrix.size()/2,QVector<double>(matrix.size()/2,0)),rA22(matrix.size()/2,QVector<double>(matrix.size()/2,0)),rA11(matrix.size()/2,QVector<double>(matrix.size()/2,0)),C11(matrix.size()/2,QVector<double>(matrix.size()/2,0)),C12(matrix.size()/2,QVector<double>(matrix.size()/2,0)),C21(matrix.size()/2,QVector<double>(matrix.size()/2,0)),C22 (matrix.size()/2,QVector<double>(matrix.size()/2,0));
    QVector<QVector<double> > nul_matrix (3,QVector<double> (3,0));
    for (int i = 0;i<matrix.size()/2;i++)
        for (int j = 0;j<matrix.size()/2;j++)
        {
            A11[i][j] = matrix[i][j];
            A12[i][j] = matrix[i+matrix.size()/2][j];
            A21[i][j] = matrix[i][j+matrix.size()/2];
            A22[i][j] = matrix[i+matrix.size()/2][j+matrix.size()/2];
        }
    rA11 = reverse_matrix_3x3_count(A11);
    rA22 = reverse_matrix_3x3_count(A22);
    C11 = reverse_matrix_3x3_count(
                difference__matrix(A11,multiplication_matrix(multiplication_matrix(A12,rA22),A21)));
    C21 = difference__matrix(
                nul_matrix,multiplication_matrix(
                    multiplication_matrix(rA22,A21),C11));
    C22 = reverse_matrix_3x3_count(
                difference__matrix(A22,multiplication_matrix(
                                       multiplication_matrix(A21,rA11),A12)));
    C12 = difference__matrix(
                nul_matrix,multiplication_matrix(
                    multiplication_matrix(rA11,A12),C22));
    for (int i = 0;i<matrix.size()/2;i++)
        for (int j = 0;j<matrix.size()/2;j++)
        {
            matrix[i][j] = C11[i][j];
            matrix[i+matrix.size()/2][j] = C12[i][j];
            matrix[i][j+matrix.size()/2] = C21[i][j];
            matrix[i+matrix.size()/2][j+matrix.size()/2] = C22[i][j];
        }
    return matrix;
}

QVector<QVector<double> > normal_system_equations::reverse_matrix_3x3_count(QVector<QVector<double> > matrix)
{
    QVector<QVector<double> > result(matrix.size(),QVector<double> (matrix[0].size(),1/determenant_matrix_3x3_count(matrix)));
    for (int i = 0;i<matrix.size();i++)
        for (int j = 0;j<matrix[0].size();j++)
            result[i][j] *= pow(-1,i+j)*determenant_matrix_2x2_count(i,j,matrix);
    return normal_system_equations::transposition_matrix(result);
}

double normal_system_equations::determenant_matrix_3x3_count(QVector<QVector<double> > matrix)
{
    return matrix[0][0]*matrix[1][1]*matrix[2][2]
            +matrix[1][0]*matrix[2][1]*matrix[0][2]
            +matrix[0][1]*matrix[1][2]*matrix[2][0]
            -matrix[2][0]*matrix[1][1]*matrix[0][2]
            -matrix[0][0]*matrix[1][2]*matrix[2][1]
            -matrix[0][1]*matrix[1][0]*matrix[2][2];
}

double normal_system_equations::determenant_matrix_2x2_count(int i,int j,QVector<QVector<double> > matrix_3x3)
{
    QVector<QVector<double> > result;
    for (int ii = 0;ii<matrix_3x3.size();ii++)
        if (ii!=i)
        {
            result.push_back(QVector<double> (0));
            for (int jj = 0;jj<matrix_3x3[0].size();jj++)
                if (jj!=j)
                {
                    result[result.size()-1].push_back(matrix_3x3[ii][jj]);
                }
        }
    return result[0][0]*result[1][1]-result[0][1]*result[1][0];
}

QVector<QVector<double> > normal_system_equations::difference__matrix(QVector<QVector<double> > matrix1, QVector<QVector<double> > matrix2)
{
    QVector<QVector<double> > result(matrix1.size(),QVector<double> (matrix1[0].size(),0));
    for (int i = 0;i<matrix1.size();i++)
        for (int j = 0;j<matrix1[0].size();j++)
            result[i][j]=matrix1[i][j]-matrix2[i][j];
    return result;
}

