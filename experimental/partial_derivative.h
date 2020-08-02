#ifndef PARTIAL_DERIVATIVE_H
#define PARTIAL_DERIVATIVE_H
#include <QVector>

class partial_derivative
{
public:
    partial_derivative();
    static QVector<QVector<QVector<double> > >  count_matrix_A (int size_H, QVector<double> deviation, QVector<double> gsc, QVector<double> geo);
    static QVector <double> count_dH (QVector <QVector <double> > sample_real,QVector <QVector <double> > sample_count);
private:
    static QVector <QVector <double> >  count_column_matrix_A (double deviation, QVector <QVector <QVector <double> > > plus_irritated_column, QVector <QVector <QVector <double> > > minus_irritated_column);

};

#endif // PARTIAL_DERIVATIVE_H
