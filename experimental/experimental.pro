#-------------------------------------------------
#
# Project created by QtCreator 2019-04-20T08:38:39
#
#-------------------------------------------------

QT       += core gui

INCLUDEPATH += D:/Qt/qwt-6.1.4/src
LIBS += -LD:/Qt/qwt-6.1.4/lib -lqwtd


greaterThan(QT_MAJOR_VERSION, 4): QT += widgets

TARGET = experimental
TEMPLATE = app


SOURCES += main.cpp\
        mainwindow.cpp \
    model_movement.cpp \
    convert.cpp \
    partial_derivative.cpp \
    normal_system_equations.cpp

HEADERS  += mainwindow.h \
    model_movement.h \
    convert.h \
    constants.h \
    partial_derivative.h \
    normal_system_equations.h

FORMS    += mainwindow.ui
