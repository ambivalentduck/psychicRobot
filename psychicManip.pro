######################################################################
# Automatically generated by qmake (2.01a) Sat Nov 22 17:22:58 2008
######################################################################

TEMPLATE = app
TARGET = manipulandumDisplay
DEPENDPATH += . 
INCLUDEPATH += .
QT += opengl network
CONFIG += debug
LIBS += -lm -lGL -lGLU -lgsl -lgslcblas
#CUDA_DIR = $$system(which nvcc | sed 's,/bin/nvcc$,,')
#INCLUDEPATH += $$CUDA_DIR/include
#QMAKE_LIBDIR += $$CUDA_DIR/lib
 
# Input
HEADERS += displaywidget.h timestuff.h point.h randb.h controlwidget.h arm.h mat2.h armsolver.h
SOURCES += main.cpp displaywidget.cpp timestuff.cpp randb.cpp controlwidget.cpp arm.cpp mat2.cpp armsolver.cpp
