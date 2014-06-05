#ifndef DISPLAYWIDGET_H
#define DISPLAYWIDGET_H

#include <QGLWidget>
#include <QRect>
#include <QBasicTimer>
#include <QMutex>
#include <QGLPixelBuffer>
#include <vector>
#include <deque>
#include "point.h"


#define LEFT .32l
#define RIGHT -.37l
#define TOP .17l
#define BOTTOM .68l
#define PROJECTORX .025l
#define PROJECTORY .66l
#define PROJECTORZ 1.35l

#define LOWERBAR .50l
#define UPPERBAR .65l

class DisplayWidget : public QGLWidget
{
	Q_OBJECT 
	
public:
	struct Sphere
	{
		point position;
		double radius;
		point color;
	};
	
	enum Shapes {TRIANGLE=0, SQUARE=1, CIRCLE=2, INFSIGN=3};
	
	DisplayWidget(QWidget *parent=0, bool FullScreen=false);
	~DisplayWidget();
	void initializeGL();
	void paintGL();
	void resizeGL(int w, int h);
	void timerEvent(QTimerEvent * event) {updateGL();}
	
	void setDeepBGColor(point color) {dataMutex.lock(); deepBackgroundColor=color; dataMutex.unlock();}
	void setBGColor(point color) {dataMutex.lock(); backgroundColor=color; dataMutex.unlock();}
	void setSpheres(std::vector<Sphere> s) {dataMutex.lock(); spheres=s; dataMutex.unlock();}
	void setBars(std::deque<double> t) {dataMutex.lock(); times=t; dataMutex.unlock();}
	void setText(QString t) {dataMutex.lock(); text=t; dataMutex.unlock();}
	void setShape(Shapes s, bool on) {drawShapes[s]=on;}
	void setShapes(bool triangle, bool square, bool circle, bool infsign) {drawShapes[0]=triangle; drawShapes[1]=square; drawShapes[2]=circle; drawShapes[3]=infsign;}
	
	
private:
	GLuint sphereList, dyntexture, shapeList[4];
	bool drawShapes[4];
	int W, H;
	QBasicTimer timer;
	std::vector<Sphere> spheres;
	std::deque<double> times;
	point backgroundColor,deepBackgroundColor;
	QMutex dataMutex;
	double min;
	QString text;
	QGLPixelBuffer * pbuffer;
};

#endif
