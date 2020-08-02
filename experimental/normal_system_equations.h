#ifndef NORMAL_SYSTEM_EQUATIONS_H
#define NORMAL_SYSTEM_EQUATIONS_H
#include <QVector>

class normal_system_equations
{
public:
    normal_system_equations();
    static QVector<QVector<double> > multiplication_matrix(QVector<QVector<double> > matrix1, QVector<QVector<double> > matrix2);
    static QVector<QVector<double> > transposition_matrix(QVector<QVector<double> > matrix);
    static QVector<QVector<double> > reverse_matrix_count(QVector<QVector<double> > matrix);
    static QVector<QVector<double> > reverse_matrix_3x3_count(QVector<QVector<double> > matrix);
    static double determenant_matrix_3x3_count(QVector<QVector<double> > matrix);
    static double determenant_matrix_2x2_count(int i,int j,QVector<QVector<double> > matrix_3x3);
    static QVector<QVector<double> > difference__matrix(QVector<QVector<double> > matrix1, QVector<QVector<double> > matrix2);
};

#endif // NORMAL_SYSTEM_EQUATIONS_H
