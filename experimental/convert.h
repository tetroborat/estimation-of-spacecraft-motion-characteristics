#ifndef CONVERT_H
#define CONVERT_H
#include <qvector.h>
#include "constants.h"
#include "math.h"

class convert
{
public:
    convert();
    static QVector <double> gsk2isk (QVector <double> gsc, QVector <double> geo);
    static QVector <double> isk2DAzYm (QVector <double> isc);
};

#endif // CONVERT_H
