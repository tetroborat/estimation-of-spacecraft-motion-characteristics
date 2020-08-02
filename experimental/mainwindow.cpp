#include "mainwindow.h"
#include "ui_mainwindow.h"
#include "model_movement.h"
#include "partial_derivative.h"
#include "convert.h"
#include <qdebug.h>
#include "QStandardItem"
#include "normal_system_equations.h"

MainWindow::MainWindow(QWidget *parent) :
    QMainWindow(parent),
    ui(new Ui::MainWindow)
{
    ui->setupUi(this);
}

MainWindow::~MainWindow()
{
    delete ui;
}

void MainWindow::read_entry_conditions()
{
    deviation.clear();
    deviation.push_back(1000*ui->dx->text().toDouble());
    deviation.push_back(1000*ui->dy->text().toDouble());
    deviation.push_back(1000*ui->dz->text().toDouble());
    deviation.push_back(1000*ui->dvx->text().toDouble());
    deviation.push_back(1000*ui->dvy->text().toDouble());
    deviation.push_back(1000*ui->dvz->text().toDouble());

    gsc_real.clear();
    gsc_real.push_back(1000*ui->xr->text().toDouble());
    gsc_real.push_back(1000*ui->yr->text().toDouble());
    gsc_real.push_back(1000*ui->zr->text().toDouble());
    gsc_real.push_back(1000*ui->vxr->text().toDouble());
    gsc_real.push_back(1000*ui->vyr->text().toDouble());
    gsc_real.push_back(1000*ui->vzr->text().toDouble());

    gsc_count.clear();
    gsc_count.push_back(1000*ui->xc->text().toDouble());
    gsc_count.push_back(1000*ui->yc->text().toDouble());
    gsc_count.push_back(1000*ui->zc->text().toDouble());
    gsc_count.push_back(1000*ui->vxc->text().toDouble());
    gsc_count.push_back(1000*ui->vyc->text().toDouble());
    gsc_count.push_back(1000*ui->vzc->text().toDouble());

    msd.clear();
    msd.push_back(ui->msdD->text().toDouble());
    msd.push_back(ui->msdaz->text().toDouble()*M_PI/180);
    msd.push_back(ui->msdym->text().toDouble()*M_PI/180);

    geo.clear();
    geo.push_back(ui->geo0->text().toDouble()*M_PI/180);
    geo.push_back(ui->geo1->text().toDouble()*M_PI/180);
}

void MainWindow::create_table_H(QVector <QVector <QVector <double> > > sample)
{
    sample_real = sample[0];
    sample_count = sample[1];
    time_begin_H = sample[2][0][0];
    time_end_H = sample[2][0][1];
    matrix_H.clear();
    matrix_H.push_back(QVector <double> (0,0));
    QStandardItem *item;

    if (ui->matrix_H->horizontalHeader()->count()==0)
    {//Заголовки столбцов
    QStringList horizontalHeader;
    horizontalHeader.append("δh");

    //Заголовки строк
    QStringList verticalHeader;
    model_H->setHorizontalHeaderLabels(horizontalHeader);

    for (int i = 0;i<sample_real.size() && i<sample_count.size();i++)
    {
        verticalHeader.append(QString::number(i));
        verticalHeader.append("");
        verticalHeader.append("");
        model_H->setVerticalHeaderLabels(verticalHeader);
        for (int j = 0;j<3;j++)
        {
            item = new QStandardItem(QString::number(sample_real[i][j]-sample_count[i][j]));
            model_H->setItem(3*i+j,0, item);
            matrix_H[0].append(sample_real[i][j]-sample_count[i][j]);
        }
    }}
    else
    {
        int size = ui->matrix_H->horizontalHeader()->count();
        for (int i = 0;i<sample_count.size() && i<sample_real.size();i++)
        {
            for (int j = 0;j<3;j++)
            {
                item = new QStandardItem(QString::number(sample_real[i][j]-sample_count[i][j]));
                model_H->setItem(3*i+j,size, item);
                matrix_H[0].append(sample_real[i][j]-sample_count[i][j]);
            }
        }
    }
    ui->matrix_H->setModel(model_H);

    ui->matrix_H->resizeRowsToContents();
    ui->matrix_H->resizeColumnsToContents();
}

void MainWindow::create_table_C_and_L()
{
    matrix_Ph.clear();
    for (int i = 0;i<matrix_A[0].size();i+=3)
    {
        matrix_Ph.push_back(QVector<double> (matrix_A[0].size(),0));
        matrix_Ph.push_back(QVector<double> (matrix_A[0].size(),0));
        matrix_Ph.push_back(QVector<double> (matrix_A[0].size(),0));
        matrix_Ph[i][i] = 4/msd[0]/msd[0];
        matrix_Ph[i+1][i+1] = 4/msd[1]/msd[1];
        matrix_Ph[i+2][i+2] = 4/msd[2]/msd[2];
    }
    QVector<QVector<double> > At_multi_Ph = normal_system_equations::multiplication_matrix(normal_system_equations::transposition_matrix(matrix_A),matrix_Ph);
    matrix_C = normal_system_equations::multiplication_matrix(At_multi_Ph,matrix_A);
    matrix_C_reverse = normal_system_equations::reverse_matrix_count(matrix_C);
    matrix_L = normal_system_equations::multiplication_matrix(At_multi_Ph,matrix_H);
    matrix_q = normal_system_equations::multiplication_matrix(matrix_C_reverse,matrix_L);
    matrix_CxC_reverse = normal_system_equations::multiplication_matrix(matrix_C,matrix_C_reverse);

}

void MainWindow::create_table_A()
{
    matrix_A.clear();
    matrix_A_intermediat = partial_derivative::count_matrix_A(matrix_H[0].size(), deviation,gsc_count,geo);
    QVector <QVector <double> > matrix_H_new (1,QVector <double> ());;
    int min,max;
    min = time_begin_H;
    max = time_end_H;
    for (int i = 0;i<matrix_A_intermediat.size();i++)
    {
        if (min<matrix_A_intermediat[i][1][0])
        {
            min = matrix_A_intermediat[i][1][0];
        }
        matrix_A.push_back(QVector <double> ());
    }
    for (int i = 0;i<matrix_A_intermediat.size();i++)
        if (max>matrix_A_intermediat[i][1][1])
        {
            max = matrix_A_intermediat[i][1][1];
        }

    for (int i = min;i<max;i++)
    {
        matrix_H_new[0].push_back(matrix_H[0][(i-time_begin_H)*3]);
        matrix_H_new[0].push_back(matrix_H[0][(i-time_begin_H)*3+1]);
        matrix_H_new[0].push_back(matrix_H[0][(i-time_begin_H)*3+2]);
        for (int j = 0;j<6;j++)
        {
            matrix_A[j].push_back(matrix_A_intermediat[j][0][(i-matrix_A_intermediat[j][1][0])*3]);
            matrix_A[j].push_back(matrix_A_intermediat[j][0][(i-matrix_A_intermediat[j][1][0])*3+1]);
            matrix_A[j].push_back(matrix_A_intermediat[j][0][(i-matrix_A_intermediat[j][1][0])*3+2]);
        }
    }

    matrix_H = matrix_H_new;
}

void MainWindow::create_table_Kq()
{
    matrix_Kq = normal_system_equations::difference__matrix(
                    normal_system_equations::multiplication_matrix(
                        normal_system_equations::multiplication_matrix(
                            normal_system_equations::transposition_matrix(matrix_H),matrix_Ph),matrix_H),
                                normal_system_equations::multiplication_matrix(
                                    normal_system_equations::transposition_matrix(matrix_L),matrix_q));
    double Kq = matrix_Kq[0][0]/(matrix_Ph.size()-6);
    for (int i = 0;i<6;i++)
    {
        matrix_Kq.push_back(QVector <double> ());
        for (int j = 0;j<6;j++)
            matrix_Kq[i].push_back(matrix_C_reverse[i][j]*Kq);
    }
}

void MainWindow::fill_table_A()
{
    QStandardItemModel *model = new QStandardItemModel;
    QStandardItem *item;

    //Заголовки столбцов
    QStringList horizontalHeader;
    horizontalHeader.append("∂h/∂X");
    horizontalHeader.append("∂h/∂Y");
    horizontalHeader.append("∂h/∂Z");
    horizontalHeader.append("∂h/∂Vx");
    horizontalHeader.append("∂h/∂Vy");
    horizontalHeader.append("∂h/∂Vz");
    model->setHorizontalHeaderLabels(horizontalHeader);

    //Заголовки строк
    QStringList verticalHeader;

    for (int i = 0;i<matrix_A[0].size();i++)
    {
        if (int(i%3) == 0)
            verticalHeader.append(QString::number(i/3));
        else
            verticalHeader.append("");
        model->setVerticalHeaderLabels(verticalHeader);
        for (int j = 0;j<6;j++)
        {
            item = new QStandardItem(QString::number(matrix_A[j][i]));
            model->setItem(i,j, item);
        }
    }
    ui->matrix_A->setModel(model);

    ui->matrix_A->resizeRowsToContents();
    ui->matrix_A->resizeColumnsToContents();
}

void MainWindow::fill_table_C_etc()
{
    QStandardItemModel *model_C = new QStandardItemModel;
    QStandardItemModel *model_C_reverse = new QStandardItemModel;
    QStandardItemModel *model_CxC_reverse = new QStandardItemModel;
    QStandardItemModel *model_Kq = new QStandardItemModel;
    QStandardItemModel *model_L = new QStandardItemModel;
    QStandardItemModel *model_q = new QStandardItemModel;
    QStandardItem *item;

    for (int i = 0;i<6;i++)
    {
        item = new QStandardItem(QString::number(matrix_L[0][i]));
        model_L->setItem(i,0, item);
    }
    ui->matrix_L->setModel(model_L);
    ui->matrix_L->resizeRowsToContents();
    ui->matrix_L->resizeColumnsToContents();

    for (int i = 0;i<6;i++)
    {
        item = new QStandardItem(QString::number(matrix_q[0][i]));
        model_q->setItem(i,0, item);
    }
    ui->matrix_q->setModel(model_q);
    ui->matrix_q->resizeRowsToContents();
    ui->matrix_q->resizeColumnsToContents();

    for (int i = 0;i<6;i++)
    {
        for (int j = 0;j<6;j++)
        {
            item = new QStandardItem(QString::number(matrix_C[j][i]));
            model_C->setItem(i,j, item);
        }
    }
    ui->matrix_C->setModel(model_C);
    ui->matrix_C->resizeRowsToContents();
    ui->matrix_C->resizeColumnsToContents();

    for (int i = 0;i<6;i++)
    {
        for (int j = 0;j<6;j++)
        {
            item = new QStandardItem(QString::number(matrix_C_reverse[j][i]));
            model_C_reverse->setItem(i,j, item);
        }
    }
    ui->matrix_C_reverse->setModel(model_C_reverse);
    ui->matrix_C_reverse->resizeRowsToContents();
    ui->matrix_C_reverse->resizeColumnsToContents();

    for (int i = 0;i<6;i++)
    {
        for (int j = 0;j<6;j++)
        {
            item = new QStandardItem(QString::number(matrix_CxC_reverse[j][i]));
            model_CxC_reverse->setItem(i,j, item);
        }
    }
    ui->matrix_CxC_reverse->setModel(model_CxC_reverse);
    ui->matrix_CxC_reverse->resizeRowsToContents();
    ui->matrix_CxC_reverse->resizeColumnsToContents();

    for (int i = 0;i<6;i++)
    {
        for (int j = 0;j<6;j++)
        {
            item = new QStandardItem(QString::number(matrix_Kq[j][i]));
            model_Kq->setItem(i,j, item);
        }
    }
    ui->matrix_Kq->setModel(model_Kq);
    ui->matrix_Kq->resizeRowsToContents();
    ui->matrix_Kq->resizeColumnsToContents();
}
#include <qwt_plot_curve.h>
void MainWindow::create_curves(QVector <QVector <double> > sample_real,QVector <QVector <double> > sample_count)
{
    d_plot->setAutoReplot(true);
    // Создать поле со шкалами для отображения графика
    ui->horizontalLayout->addWidget(d_plot);// привязать поле к границам окна
    d_plot->setTitle( "Трасса полета" ); // заголовок

    d_plot->setAxisTitle(QwtPlot::yLeft, "Широта, град");
    d_plot->setAxisTitle(QwtPlot::xBottom, "Долгота, град");



    ////////////////////////////////////////////////////////////////////////////
    QwtPlotCurve * curve = new QwtPlotCurve();
    QwtPlotCurve * curve1 = new QwtPlotCurve();
    QwtPlotCurve * dot = new QwtPlotCurve();
    QPolygonF points,points1,dot1;
    dot1<<QPointF(geo[0]*180/M_PI,geo[1]*180/M_PI);
    dot->setPen(QPen(Qt::red,10,Qt::SolidLine));
    dot->setStyle(QwtPlotCurve::Dots);
    //dot->setSamples(dot1);
    curve->setPen(QPen(Qt::red,2,Qt::SolidLine));
    curve->attach(d_plot);
    curve1->attach(d_plot);
    dot->attach(d_plot);
    double x,y,x1,y1;
    for (int zz =0;zz<sample_real.size();zz++)
    {
        //x = sample_real[zz][3];
        //y = sample_real[zz][4];
        x = zz;
        y = sample_real[zz][0];
        points<<QPointF(x,y);

//        x1 = sample_count[zz][3];
//        y1 = sample_count[zz][4];
        x1 = zz;
        y1 = sample_count[zz][2]*180/M_PI;
        //points1<<QPointF(x1,y1);
    }
    curve->setSamples(points);
    curve1->setSamples(points1);
        ///////////////////////////////////////////////////////////////////////////
}

void MainWindow::count_one_iteration()
{
    create_table_H(model_movement::puling_movement_discrepancy(gsc_real,gsc_clarify,geo,msd));
                    qApp->processEvents();
    create_table_A();
                    qApp->processEvents();
    create_table_C_and_L();
                    qApp->processEvents();
    for (int i = 0;i<gsc_count.size();i++)
        gsc_clarify[i] += matrix_q[0][i];
    qDebug()<<gsc_clarify<<matrix_q;
}

void MainWindow::on_commandLinkButton_clicked()
{
    model_H->clear();
    read_entry_conditions();
    gsc_clarify = gsc_count;

    while (sqrt(pow(gsc_real[0]-gsc_clarify[0],2)+pow(gsc_real[1]-gsc_clarify[1],2)+pow(gsc_real[2]-gsc_clarify[2],2))>3000
           ||
           sqrt(pow(gsc_real[3]-gsc_clarify[3],2)+pow(gsc_real[4]-gsc_clarify[4],2)+pow(gsc_real[5]-gsc_clarify[5],2))>3)
        count_one_iteration();
    create_table_Kq();
    fill_table_A();
    fill_table_C_etc();
    ui->dr->setText(QString::number(sqrt(pow(gsc_real[0]-gsc_clarify[0],2)+pow(gsc_real[1]-gsc_clarify[1],2)+pow(gsc_real[2]-gsc_clarify[2],2))));
    ui->dv->setText(QString::number(sqrt(pow(gsc_real[3]-gsc_clarify[3],2)+pow(gsc_real[4]-gsc_clarify[4],2)+pow(gsc_real[5]-gsc_clarify[5],2))));

    create_curves(sample_real,sample_count);
    ui->xc1->setText(QString::number(gsc_clarify[0]/1000));
    ui->yc1->setText(QString::number(gsc_clarify[1]/1000));
    ui->zc1->setText(QString::number(gsc_clarify[2]/1000));
    ui->vxc1->setText(QString::number(gsc_clarify[3]/1000));
    ui->vyc1->setText(QString::number(gsc_clarify[4]/1000));
    ui->vzc1->setText(QString::number(gsc_clarify[5]/1000));
}
