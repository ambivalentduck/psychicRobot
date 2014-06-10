#include "displaywidget.h"
#include <GL/glu.h>
#include <cmath>
#include <iostream>

#define BUFFSIZE 1920
#define LINEWIDTH 5

DisplayWidget::DisplayWidget(QWidget *parent,bool FullScreen)
:QGLWidget(QGLFormat(QGL::DoubleBuffer|QGL::AlphaChannel|QGL::SampleBuffers|QGL::AccumBuffer), parent, 0, FullScreen?Qt::X11BypassWindowManagerHint:Qt::Widget)
{
	//Take care of window and input initialization.
	timer.start(16, this); //Draw again shortly after constructor finishes
	if(FullScreen)
	{
		setWindowState(Qt::WindowFullScreen); 
		setCursor(QCursor(Qt::BlankCursor)); //Hide the cursor
		raise(); //Make sure it's the top window
	}
	setPalette(QPalette(QColor(0, 0, 0))); //IF the background draws, draw it black.
	setAutoFillBackground(false); //Try to let glClear work...
	setAutoBufferSwap(false); //Don't let QT swap automatically, we want to control timing.
	backgroundColor=point(0,0,0);
	deepBackgroundColor=point(0,0,0);

	for(int k=0;k<4;k++) drawShapes[k]=false;
	calibrationMode=true;
	
	//Set up a "calibration" field. Should be a 1/4 circle in each corner
	Sphere sphere;
	spheres.clear();
	sphere.color=point(1,0,0);
	sphere.position=point(0,0,HANDLEDEPTH);
	sphere.radius=.018;
	spheres.push_back(sphere);
	sphere.color=point(0,1,0);
	sphere.position=point(LEFTPROBE,0,HANDLEDEPTH);
	spheres.push_back(sphere);
	sphere.color=point(0,0,1);
	sphere.position=point(0, UPPROBE,HANDLEDEPTH);
	spheres.push_back(sphere);
}

DisplayWidget::~DisplayWidget()
{
	makeCurrent();
	glDeleteLists(sphereList,1);
	for(int k=0;k<4;k++) glDeleteLists(shapeList[k],1);
}

void DisplayWidget::initializeGL()
{  
	glClearColor(0,0,0,1);
	glClear(GL_COLOR_BUFFER_BIT);
	glShadeModel(GL_FLAT);
	sphereList = glGenLists(1);
	GLUquadricObj *qobj=gluNewQuadric();
	gluQuadricDrawStyle(qobj, GLU_FILL);
	glNewList(sphereList, GL_COMPILE);
	gluSphere(qobj, 1, 100, 100); //Arbitrary defaults "grid" size: 100 x 100		
	glEndList();
	
	shapeList[TRIANGLE] = glGenLists(1);
	glNewList(shapeList[TRIANGLE],GL_COMPILE); //Triangle
	glLineWidth(LINEWIDTH);
	glBegin(GL_LINE_LOOP);
		glVertex2f(1,-1);
		glVertex2f(1,1);
		glVertex2f(-sqrt(2.0)/2.0,0);
	glEnd();
	glEndList();
	
	shapeList[SQUARE] = glGenLists(1);
	glNewList(shapeList[SQUARE],GL_COMPILE); //Square
	glLineWidth(LINEWIDTH);
	glBegin(GL_LINE_LOOP);
		glVertex2f(-1,-1);
		glVertex2f(-1,1);
		glVertex2f(1,1);
		glVertex2f(1,-1);
	glEnd();
	glEndList();
	
	double t;
	shapeList[CIRCLE] = glGenLists(1);
	glNewList(shapeList[CIRCLE],GL_COMPILE); //Circle
	glLineWidth(LINEWIDTH);
	glBegin(GL_LINE_LOOP);
		for(int k=0;k<100;k++) 
		{
			t=6.283185307*double(k)/100.0;
			glVertex2f(cos(t),sin(t));
		}
	glEnd();
	glEndList();
	
	shapeList[INFSIGN] = glGenLists(1);
	glNewList(shapeList[INFSIGN],GL_COMPILE); //Infinity sign
	glLineWidth(LINEWIDTH);
	glBegin(GL_LINE_LOOP);
		for(int k=0;k<100;k++) 
		{
			t=6.283185307*double(k)/100.0;
			glVertex2f(sin(t),.5*sin(2.0*t));
		}
	glEnd();
	glEndList();
	
	glEnable(GL_POINT_SMOOTH);
	glPointSize(1);
}

void DisplayWidget::calibrate(point Center, point probe1, point probe2)
{
	center=Center;
	double rot=std::atan2(probe1.Y()-center.Y(),probe1.X()-center.X())-std::atan2(0,LEFTPROBE); //In theory, any one point is enough
	screenRotation=(180l/M_PI)*rot; //OpenGL likes degrees.
	point one=(probe1-center).rotateZero(-rot);
	point two=(probe2-center).rotateZero(-rot);
	std::cout << rot << " " << one.X() << " " << one.Y() << " " << two.X() << " " << two.Y() << std::endl; //Why not?
	
	width=SCREENWIDTH/LEFTPROBE*one.X();
	height=SCREENHEIGHT/UPPROBE*two.Y();
	
	calibrationMode=false;
	std::cout << width << " " << height << std::endl; //Why not?
	
	glMatrixMode(GL_PROJECTION);
	glLoadIdentity();
	glOrtho(-width/2l+center.X(),width/2l+center.X(),-height/2l+center.Y(),height/2l+center.Y(),MINDEPTH,MAXDEPTH);
	glTranslated(center.X(),center.Y(),0);
	glRotated(fmod(screenRotation,180),0,0,1);
	glTranslated(-center.X(),-center.Y(),0);
}

void DisplayWidget::resizeGL(int w, int h)
{
	makeCurrent();
	W=w; H=h;
	glViewport(0,0, w, h);
	glMatrixMode(GL_PROJECTION);
	glLoadIdentity();
	if(calibrationMode)
	{
		//Naively assume the screen is the center of the universe.
		glOrtho(-SCREENWIDTH/2l,SCREENWIDTH/2l,-SCREENHEIGHT/2l,SCREENHEIGHT/2l,MINDEPTH,MAXDEPTH);
	}
	else
	{
		//Post-calibration, assume you know where the center of the screen really is.
		glOrtho(-width/2l+center.X(),width/2l+center.X(),-height/2l+center.Y(),height/2l+center.Y(),MINDEPTH,MAXDEPTH);
		glTranslated(-center.X(),-center.Y(),0);
		glRotated(-screenRotation,0,0,1);
		glTranslated(center.X(),center.Y(),0);
	}
	update();
}

void DisplayWidget::paintGL()
{
	timer.stop();
	dataMutex.lock();
	
	makeCurrent();
	glMatrixMode(GL_MODELVIEW);
	glLoadIdentity();
	if(calibrationMode) gluLookAt(0,0,0,0,0,HANDLEDEPTH,0,-1,0); //Camera has opposite default orientation
	//else gluLookAt(center.X(),center.Y(),0,center.X(),center.Y(),HANDLEDEPTH,0,1,0); //Camera has opposite default orientation
	else gluLookAt(0,0,0,0,0,HANDLEDEPTH,0,1,0); //Camera has opposite default orientation
	
	//Unused area is unlit by default
	if(calibrationMode) glClearColor(deepBackgroundColor.X(), deepBackgroundColor.Y(), deepBackgroundColor.Z(),1);
	else glClearColor(backgroundColor.X(), backgroundColor.Y(), backgroundColor.Z(),1);
	glClear(GL_COLOR_BUFFER_BIT);
	
	glShadeModel(GL_FLAT);
	glEnable(GL_LINE_SMOOTH);
	glEnable(GL_POLYGON_SMOOTH);
	glEnable(GL_POINT_SMOOTH);
	glEnable(GL_BLEND);
	glBlendFunc(GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA);
	glDisable(GL_DEPTH_TEST);
		
	glPushMatrix();
	glTranslated(center.X(),center.Y()+.05,0);
	glScaled(min*(1.0/6.0),min*(1.0/6.0),1.0);
	glColor3d(.5,.5,.5); //Grey because...why not?
	for(int k=0;k<4;k++)
	{
		if(drawShapes[k])
		{
			glCallList(shapeList[k]);
		}
	}
	glPopMatrix();
	
	for(std::vector<Sphere>::iterator it=spheres.begin();it!=spheres.end();++it)
	{
		if(it->radius<=0) continue;
		glColor3dv(it->color);
		glPushMatrix();
		if(!calibrationMode) glTranslated(0,0,HANDLEDEPTH);
		glTranslated(it->position.X(),it->position.Y(),it->position.Z());
		glScaled(it->radius,it->radius,it->radius);
		glCallList(sphereList);
		glPopMatrix();
	}

	renderText(textLocation.X(),textLocation.Y(),textLocation.Z(),text);
	
	dataMutex.unlock();
	swapBuffers();
	glFinish();  //Get precise timing by recording time after this, blocks until swap succeeds.  Swap happens during refresh.
	timer.start(15, this); //60 Hz = 16.6 ms, guarantee a paint in each refresh and almost immediately before refresh to minimize lag.
}

