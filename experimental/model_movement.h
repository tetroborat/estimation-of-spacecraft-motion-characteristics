#ifndef MODEL_MOVEMENT_H
#define MODEL_MOVEMENT_H
#include <QVector>
#include "constants.h"

class model_movement
{
public:
    model_movement();
    static QVector <QVector <QVector <double> > > puling_movement_count (QVector <double> gsc, QVector <double> geo);
    static QVector <QVector <QVector <double> > > puling_movement_discrepancy (QVector <double> gsc_real, QVector <double> gsc_count, QVector <double> geo,QVector <double> msd);
private:
    static QVector <double> rung4 (QVector <double> gsc);
    static QVector <double> integration_movement (QVector <double> gsc);
    static double R (double x, double y, double z);
    static double Ro (double h);
    static bool check_zone_overview (QVector <double> DAzYm);
    static double norm (double sig);
};

#endif // MODEL_MOVEMENT_H
