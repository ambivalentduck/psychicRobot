#include "mat2widget.h"

Mat2Widget::Mat2Widget(QString str, mat2 m, double min, double max, int dec)
{
	QLabel label(str);
	layout=new QGridLayout(this);
	layout->addWidget(&label,0,0);
	for(int k=0;k<4;k++)
	{
		boxes[k]=new QDoubleSpinBox();
		boxes[k]->setValue(m.m[k]);
		boxes[k]->setMaximum(max);
		boxes[k]->setMinimum(min);
		layout->addWidget(boxes[k],1+floor(k/2),k%2);
	}
	connect(boxes[0],SIGNAL(valueChanged(double)),this,SLOT(setm0(double)));
	connect(boxes[1],SIGNAL(valueChanged(double)),this,SLOT(setm1(double)));
	connect(boxes[2],SIGNAL(valueChanged(double)),this,SLOT(setm2(double)));
	connect(boxes[3],SIGNAL(valueChanged(double)),this,SLOT(setm3(double)));
}
