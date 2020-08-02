#include "partial_derivative.h"
#include "model_movement.h"
partial_derivative::partial_derivative()
{

}

QVector <QVector <double> > partial_derivative::count_column_matrix_A(double deviation, QVector <QVector <QVector <double> > > plus_irritated_column, QVector <QVector <QVector <double> > > minus_irritated_column)
{
    int min,max;
    QVector <QVector <double> > result (2,QVector <double> ());
    if (plus_irritated_column[1][0][0]>minus_irritated_column[1][0][0])
        min = int(plus_irritated_column[1][0][0]);
    else
        min = int(minus_irritated_column[1][0][0]);
    if (plus_irritated_column[1][0][1]<minus_irritated_column[1][0][1])
        max = int(plus_irritated_column[1][0][1]);
    else
        max = int(minus_irritated_column[1][0][1]);
    result[1].push_back(min);
    result[1].push_back(max);
    for (int i = min;i<max;i++)
    {
        for (int j = 0;j<3;j++)
            result[0].push_back((plus_irritated_column[0][i-plus_irritated_column[1][0][0]][j]-minus_irritated_column[0][i-minus_irritated_column[1][0][0]][j])/2/deviation);
    }
    return result;
}

QVector<QVector<QVector<double> > > partial_derivative::count_matrix_A(int size_H, QVector<double> deviation, QVector<double> gsc, QVector<double> geo)
{
    QVector<QVector<QVector<double> > >derivative_matrix_intermediat;
    QVector <QVector <QVector <double> > > plus_irritated_column,minus_irritated_column;
    QVector<double> gsc_plus,gsc_minus;
    for (int i = 0;i<deviation.size();i++)
    {
        gsc_plus=gsc;
        gsc_minus=gsc;
        gsc_plus[i]=gsc[i]+deviation[i];
        gsc_minus[i]=gsc[i]-deviation[i];
        plus_irritated_column = model_movement::puling_movement_count(gsc_plus,geo);
        minus_irritated_column = model_movement::puling_movement_count(gsc_minus,geo);
        derivative_matrix_intermediat.push_back(count_column_matrix_A(deviation[i],plus_irritated_column,minus_irritated_column));
    }

    return derivative_matrix_intermediat;
}

QVector<double> partial_derivative::count_dH(QVector<QVector<double> > sample_real, QVector<QVector<double> > sample_count)
{

}

