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
	sphere.color=point(0,.5,0);
	sphere.position=point(0,0,HANDLEDEPTH);
	sphere.radius=.1;
	spheres.push_back(sphere);
	sphere.position=(LEFTPROBE,0,HANDLEDEPTH);
	spheres.push_back(sphere);
	sphere.position=(0, UPPROBE,HANDLEDEPTH);
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
	double rot=std::atan2(probe1.Y()-center.Y(),probe1.X()-center.X()); //In theory, any one point is enough
	screenRotation=(180l/M_PI)*rot; //OpenGL likes degrees.
	point one=(probe1-center).rotateZero(-rot);
	point two=(probe2-center).rotateZero(-rot);
	std::cout << (180l/M_PI)*std::atan2(one.Y(),one.X()) << " " << (180l/M_PI)*std::atan2(two.Y(),two.X()) << std::endl; //Why not?
	
	if(two.mag()>one.mag()) // make sure one > two
	{
		point temp=two;
		two=one;
		one=temp;
	}
	
	if(std::abs(LEFTPROBE)>std::abs(UPPROBE))
	{
		width=SCREENWIDTH/LEFTPROBE*one.X();
		height=SCREENHEIGHT/UPPROBE*two.Y();
	}
	else
	{
		width=SCREENWIDTH/LEFTPROBE*two.X();
		height=SCREENHEIGHT/UPPROBE*one.Y();
	}
	calibrationMode=false;
	
	pbuffer = new QGLPixelBuffer(QSize(BUFFSIZE,BUFFSIZE),QGLFormat(QGL::DoubleBuffer|QGL::AlphaChannel),this);
	pbuffer->makeCurrent();
	glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MAG_FILTER,GL_NEAREST);
	glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MIN_FILTER,GL_NEAREST);
	//glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MAG_FILTER,GL_LINEAR);
	//glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MIN_FILTER,GL_LINEAR);
	glTexEnvf(GL_TEXTURE_ENV, GL_TEXTURE_ENV_MODE, GL_REPLACE);
	
	glClearColor(0,0,0,1);
	glClear(GL_COLOR_BUFFER_BIT);
	glShadeModel(GL_FLAT);
	glViewport(0,0,BUFFSIZE,BUFFSIZE);
	glMatrixMode(GL_PROJECTION);
	glLoadIdentity();
	glOrtho(-width/2l,width/2l,height/2l,height/2l,MINDEPTH,MAXDEPTH);

	dyntexture=pbuffer->generateDynamicTexture();
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
		glFrustum(-SCREENWIDTH/2l,SCREENWIDTH/2l,-SCREENHEIGHT/2l,SCREENHEIGHT/2l,MINDEPTH,MAXDEPTH);
	}
	else
	{
		//Render from projector's perspective, projector must be 0,0,0, looking down -Z
		glFrustum(-width/2l-center.X(),width/2l-center.X(),-height/2l-center.Y(),height/2l-center.Y(),MINDEPTH,MAXDEPTH);
		glRotated(-screenRotation,0,0,1);
	}
	update();
}

void DisplayWidget::paintGL()
{
	timer.stop();
	dataMutex.lock();
	
	if(calibrationMode)
	{
		makeCurrent();
		glMatrixMode(GL_MODELVIEW);
		glLoadIdentity();
		//gluLookAt(0,0,0,0,0,HANDLEDEPTH,0,-1,0);
	}
	else
	{
		pbuffer->makeCurrent();
		glClearColor(backgroundColor.X(), backgroundColor.Y(), backgroundColor.Z(),1);
		glClear(GL_COLOR_BUFFER_BIT);
		glShadeModel(GL_FLAT);
		glMatrixMode(GL_MODELVIEW);
		glLoadIdentity(); 
		
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
	}
	
	for(std::vector<Sphere>::iterator it=spheres.begin();it!=spheres.end();++it)
	{
		if(it->radius<=0) continue;
		glColor3dv(it->color);
		glPushMatrix();
		glTranslated(it->position.X(),it->position.Y(),it->position.Z());
		glScaled(it->radius,it->radius,it->radius);
		glCallList(sphereList);
		glPopMatrix();
	}

	if(!calibrationMode)
	{
		pbuffer->updateDynamicTexture(dyntexture);

		makeCurrent();
		glMatrixMode(GL_MODELVIEW);
		glLoadIdentity();
		//glScaled(1,-1,1);
		glTranslated(-center.X(),-center.Y(),HANDLEDEPTH); //Projector now ignored
			
		glClearColor(deepBackgroundColor.X(), deepBackgroundColor.Y(), deepBackgroundColor.Z(),1);  //Unused area is unlit by default
		glClear(GL_COLOR_BUFFER_BIT);
		glEnable(GL_TEXTURE_2D);
		glTexEnvf(GL_TEXTURE_ENV, GL_TEXTURE_ENV_MODE, GL_REPLACE);
		glBindTexture(GL_TEXTURE_2D,dyntexture);
		glBegin(GL_POLYGON);
			glTexCoord2f(0,0); glVertex2f(-width/2l,-height/2l);
			glTexCoord2f(0,1); glVertex2f(-width/2l,height/2l);
			glTexCoord2f(1,1); glVertex2f(width/2l,height/2l);
			glTexCoord2f(1,0); glVertex2f(width/2l,height/2l);
		glEnd();
		glDisable(GL_TEXTURE_2D);
		
		renderText(textLocation.X(),textLocation.Y(),textLocation.Z(),text);
	}
	dataMutex.unlock();
	swapBuffers();
	glFinish();  //Get precise timing by recording time after this, blocks until swap succeeds.  Swap happens during refresh.
	timer.start(15, this); //60 Hz = 16.6 ms, guarantee a paint in each refresh and almost immediately before refresh to minimize lag.
}

