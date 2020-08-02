#include "convert.h"
#include <qdebug.h>
convert::convert()
{

}

QVector<double> convert::gsk2isk(QVector<double> gsc, QVector<double> geo)
{
    QVector < double > ty;
    ty.push_back(-(gsc[0]*cos(geo[0])+gsc[1]*sin(geo[0]))*sin(geo[1])+gsc[2]*cos(geo[1]));
    ty.push_back((gsc[0]*cos(geo[0])+gsc[1]*sin(geo[0]))*cos(geo[1])+gsc[2]*sin(geo[1])-pow(Rze,2)*(1-alzg)/sqrt(pow(Rze*(1-alzg)*cos(geo[1]),2)+pow(Rze*sin(geo[1]),2)));
    ty.push_back(-gsc[0]*sin(geo[0])+gsc[1]*cos(geo[0]));
    return ty;
}

QVector<double> convert::isk2DAzYm(QVector<double> isc)
{
    QVector <double> result;
    result.push_back(sqrt(isc[0]*isc[0]+isc[1]*isc[1]+isc[2]*isc[2]));
    if (isc[0]>0 && isc[2]<0) result.push_back(fmod(5*M_PI_2-atan2(isc[2],isc[0]),2*M_PI));
    else result.push_back(fmod(M_PI_2-atan2(isc[2],isc[0]),2*M_PI));
    result.push_back(asin(isc[1]/result[0]));

    return result;
}

