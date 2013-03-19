#ifndef MAT2WIDGET_H
#define MAT2WIDGET_H

#include <QWidget>
#include <QtGui>
#include "mat2.h"

class Mat2Widget : public QWidget
{
	Q_OBJECT
	
public:
	Mat2Widget(QString str, mat2 m,double min=-50, double max=50, int dec=3);
	mat2 mat;
signals:
	void valueChanged(mat2 m);
public slots:
	void setm0(double d) {mat.m[0]=d; emit(valueChanged(mat))};
	void setm1(double d) {mat.m[1]=d; emit(valueChanged(mat))};
	void setm2(double d) {mat.m[2]=d; emit(valueChanged(mat))};
	void setm3(double d) {mat.m[3]=d; emit(valueChanged(mat))};
private:
	QDoubleSpinBox * boxes[4];
	QGridLayout * layout;
	QLabel label;
};

#endif
