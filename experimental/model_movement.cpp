#include "model_movement.h"
#include <math.h>
#include <qdebug.h>
#include "convert.h"

model_movement::model_movement()
{

}


QVector<double> model_movement::rung4(QVector<double> gsc)
{
    QVector <double> k1,k2,k3,k4,nu;
    nu = gsc;
    k1 = integration_movement(nu);
    for (int i =0;i<6;i++)
        nu[i] = gsc[i]+k1[i]/2;
    k2 = integration_movement(nu);
    for (int i =0;i<6;i++)
        nu[i] = gsc[i]+k2[i]/2;
    k3 = integration_movement(nu);
    for (int i =0;i<6;i++)
        nu[i] = gsc[i]+k3[i];
    k4 = integration_movement(nu);
    for (int i =0;i<6;i++)
        nu[i] = gsc[i] + (k1[i]+2*k2[i]+2*k3[i]+k4[i])/6;
    return nu;
}

QVector<double> model_movement::integration_movement(QVector<double> gsc)
{
    double count = 1;
    QVector <double> result;
    for (int i = 3;i<6;i++)
        result.push_back(count*gsc[i]);
    double a,b,c,d,ro,v,r,h;

    v = R(gsc[3],gsc[4],gsc[5]);
    r = R(gsc[0],gsc[1],gsc[2]);
    h = r - aKrasov*(1-0.2*alzg*d);
    b = r0/r/r/r;
    c = 1.5*al20*r0*r0/r/r;
    d = 5*gsc[2]*gsc[2]/r/r;
    a = b*(al00+c*(d-1));
    h = r - aKrasov*(1-0.2*alzg*d);
    ro = Ro(h);
    if (((w0*w0-a)*gsc[0]+2*w0*gsc[4]-sb*ro*v*gsc[3])*count==NAN)
        qDebug()<<(w0*w0-a)*gsc[0]<<2*w0*gsc[4]<<sb*ro*v*gsc[3];
    result.push_back(((w0*w0-a)*gsc[0]+2*w0*gsc[4]-sb*ro*v*gsc[3])*count);
    result.push_back(((w0*w0-a)*gsc[1]-2*w0*gsc[3]-sb*ro*v*gsc[4])*count);
    result.push_back(((2*b*c-a)*gsc[2]-sb*ro*v*gsc[5])*count);
    return result;
}

double model_movement::R(double x, double y, double z)
{
    return sqrt(x*x+y*y+z*z);
}

double model_movement::Ro(double h)
{
    double H[8] = {0,0.2,0.6,1,1.5,3,6,9};
    double kk0[8] = {1.225,0.891e-2,2.578e-4,4.061e-7,2.13e-9,4.764e-11,8.726e-12,6.367e-13};
    double kk1[8] = {-0.2639e-8,0.4407e-9,-0.256e-8,0.14688e-8,0.8004e-10,0.7111e-11,0.1831e-11,0};
    double kk2[8] = {0.7825e-4,0.16375e-3,0.5905e-4,0.1787e-3,0.37336e-3,0.15467e-4,0.9275e-5,0.954e-5};
    int i=0;
    while (h>H[i]*100000)
        i+=1;
    i-=1;
    return kk0[i]*exp(kk1[i]*pow((h-H[i]),2)-kk2[i]*(h-H[i]));
}

QVector <QVector <QVector <double> > > model_movement::puling_movement_discrepancy(QVector <double> gsc_real, QVector <double> gsc_count, QVector<double> geo, QVector <double> msd)
{
    QVector <QVector <QVector <double> > > selection_measurements (3,QVector <QVector <double> >());
    selection_measurements[2].push_back(QVector <double> ());
    QVector <double> DAzYm_real,DAzYm_count;
    bool check = false;
    bool check1 = false;
    int time = 0;
    while (!check1)
        {
            DAzYm_real = convert::isk2DAzYm(convert::gsk2isk(gsc_real,geo));
            DAzYm_count = convert::isk2DAzYm(convert::gsk2isk(gsc_count,geo));
            if (check)
            {
                for (int i = 0;i<DAzYm_real.size();i++)
                    DAzYm_real[i] += norm(msd[i]);
                /*долгота*/
                if (gsc_real[1]>=0)
                    DAzYm_real.push_back(acos(gsc_real[0]/pow(gsc_real[0]*gsc_real[0]+gsc_real[1]*gsc_real[1],0.5))/M_PI*180.);
                else
                    DAzYm_real.push_back(360.-acos(gsc_real[0]/pow(gsc_real[0]*gsc_real[0]+gsc_real[1]*gsc_real[1],0.5))/M_PI*180.);
                /*широта*/
                DAzYm_real.push_back(asin(gsc_real[2]/pow((gsc_real[0]*gsc_real[0]+gsc_real[1]*gsc_real[1]+gsc_real[2]*gsc_real[2]),0.5))/M_PI*180.);

                /*долгота*/
                if (gsc_count[1]>=0)
                    DAzYm_count.push_back(acos(gsc_count[0]/pow(gsc_count[0]*gsc_count[0]+gsc_count[1]*gsc_count[1],0.5))/M_PI*180.);
                else
                    DAzYm_count.push_back(360.-acos(gsc_count[0]/pow(gsc_count[0]*gsc_count[0]+gsc_count[1]*gsc_count[1],0.5))/M_PI*180.);
                /*широта*/
                DAzYm_count.push_back(asin(gsc_count[2]/pow((gsc_count[0]*gsc_count[0]+gsc_count[1]*gsc_count[1]+gsc_count[2]*gsc_count[2]),0.5))/M_PI*180.);

                selection_measurements[0].push_back(DAzYm_real);
                selection_measurements[1].push_back(DAzYm_count);
            }
            gsc_real = rung4(gsc_real);
            gsc_count = rung4(gsc_count);
            time ++;
            if ((!check_zone_overview(DAzYm_real) || !check_zone_overview(DAzYm_count)) && check)
            {
                check1 = true;
                selection_measurements[2][0].push_back(time-1);
            }
            if (check_zone_overview(DAzYm_real) && check_zone_overview(DAzYm_count) && selection_measurements[2][0].size() == 0)
            {
                check = true;
                selection_measurements[2][0].push_back(time-1);
            }
        }
    return selection_measurements;
}
QVector <QVector <QVector <double> > > model_movement::puling_movement_count(QVector<double> gsc, QVector<double> geo)
{
    QVector <QVector <QVector <double> > > selection_measurements (2,QVector <QVector <double> > ());
    selection_measurements[1].push_back(QVector<double> ());
    QVector <double> DAzYm;
    bool check = false;
    bool check1 = false;
    int time = 0;
    while (!check1)
        {
            DAzYm = convert::isk2DAzYm(convert::gsk2isk(gsc,geo));
            if (check)
            {
               selection_measurements[0].push_back(DAzYm);
            }
            gsc = rung4(gsc);
            time ++;
            if (!check_zone_overview(DAzYm) && check)
            {
                check1 = true;
                selection_measurements[1][0].push_back(time-1);
            }
            if (check_zone_overview(DAzYm) && selection_measurements[1][0].size() == 0)
            {
                check = true;
                selection_measurements[1][0].push_back(time-1);
            }
        }
    return selection_measurements;
}


bool model_movement::check_zone_overview(QVector<double> DAzYm)
{
    if (DAzYm[2]>(7/180*M_PI)) return true;
    else return false;
}

double model_movement::norm(double sig)
{
    return (rand()%10001-5000)*0.0002*sig;
}

