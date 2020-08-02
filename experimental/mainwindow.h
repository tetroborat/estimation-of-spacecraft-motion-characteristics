#ifndef MAINWINDOW_H
#define MAINWINDOW_H

#include <QMainWindow>
#include <QTableView>
#include <qwt_plot.h>
#include "QStandardItemModel"

namespace Ui {
class MainWindow;
}

class MainWindow : public QMainWindow
{
    Q_OBJECT

public:
    explicit MainWindow(QWidget *parent = 0);
    ~MainWindow();

private slots:
    void on_commandLinkButton_clicked();

private:
    Ui::MainWindow *ui;
    int time_begin_H,time_end_H;
    QVector<double> deviation,gsc_real,gsc_count,gsc_clarify,geo,msd;
    double sig0 = 2;
    QVector<QVector<QVector<double> > > matrix_A_intermediat;
    QVector <QVector <double> > sample_real, sample_count, matrix_A, matrix_Ph, matrix_C, matrix_C_reverse, matrix_CxC_reverse, matrix_H, matrix_L, matrix_q, matrix_Kq;
    QStandardItemModel *model_H = new QStandardItemModel;
    void read_entry_conditions();
    void create_table_A();
    void create_table_Kq();
    void fill_table_A();
    void fill_table_C_etc();
    void create_table_H(QVector <QVector <QVector <double> > > sample);
    void create_table_C_and_L();
    void create_curves(QVector <QVector <double> > sample_real,QVector <QVector <double> > sample_count);
    void count_one_iteration();
    QwtPlot *d_plot = new QwtPlot();
};

#endif // MAINWINDOW_H
