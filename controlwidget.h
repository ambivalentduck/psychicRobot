#ifndef CONTROLWIDGET_H
#define CONTROLWIDGET_H

#include <QtGui>
#include <QByteArray>
#include <QUdpSocket>
#include <QFile>
#include <QTextStream>
#include <QString>
#include <vector>
#include <deque>
#include <iostream>
#include "displaywidget.h"
#include "timestuff.h"
#include <cmath>
#include "armsolver.h"
#include "arm.h"


class ControlWidget : public QWidget
{
	Q_OBJECT

public:
	ControlWidget(QDesktopWidget * qdw);
	
private:
	QSpinBox *trialNumBox, *subjectBox;
	QDoubleSpinBox *virtualMassBox, *earlyPulseGainBox, *latePulseGainBox, *blwnGainBox;
	QDoubleSpinBox *eaGainBox, *cursorFadeBox, *extractionFadeBox, *rawFadeBox;
	QDoubleSpinBox *heightBox, *weightBox, *armL1Box, *armL2Box;
	QDoubleSpinBox *x0xBox, *x0yBox;
	QPushButton *startButton; 	
	QComboBox *stimulusBox;
	QFormLayout * layout;
		
	DisplayWidget * userWidget;
	ArmSolver * armsolver;

	void closeEvent(QCloseEvent *event);
	void goGray() {for(std::vector<QWidget*>::iterator it=grayList.begin();it!=grayList.end();++it) (*it)->setEnabled(false); }
	void unGray() {for(std::vector<QWidget*>::iterator it=grayList.begin();it!=grayList.end();++it) (*it)->setEnabled(true); }
	point loadTrial(int T);
	
	QByteArray in,out;
	int inSize, outSize;
	QUdpSocket * us;
	QFile contFile, trialFile;
	QTextStream outStream;
	
	double min, perturbGain, weight, eaGain, xpcTime, virtualMass, cursorFadeTime, extractionFadeTime, rawFadeTime, blwnGain, earlyPulseGain, latePulseGain;
	enum GameState {acquireTarget=0, inTarget=1, hold=2} state;
	std::vector<QWidget*> grayList;
	std::vector<DisplayWidget::Sphere> sphereVec;
	std::deque<timespec> times;
	std::deque<QByteArray> data;
	std::deque<point> extractedD, handleD, cursorD;
	DisplayWidget::Sphere sphere;
	
	twoLinkArm::ArmParams params;
	timespec zero, now, trialStart, targetAcquired, holdStart, lastFade;
	bool ExperimentRunning, inputReady, outputReady, ignoreInput, leftOrigin, firstpush, leftSide;
	int trial, subject, pulls;
	point x0, origin, cursor, desposition, position, velocity, accel, target, force, center;
	
signals:
	void endApp();
	
public slots:
	void readPending();
	void startClicked();
	void setTrialNum(int i) {trial=i;}
	void setSubject(int i) {subject=i;}
	double evalSigmoid(double t, double risetime) {double a=10l/risetime; t-=risetime/2l; return (a*t/sqrt(1l+pow(a*t,2))+1l)/2l;}
	void setEAGain(double g) {eaGain=g;}
	void setCursorFade(double f) {cursorFadeTime=f;}
	void setRawFade(double f) {rawFadeTime=f;}
	void setExtractionFade(double f) {extractionFadeTime=f;}
	void setl1(double l) {params=twoLinkArm::calcParams(weight,l,params.l2,x0);}
	void setl2(double l) {params=twoLinkArm::calcParams(weight,params.l1,l,x0);}
	void setHeight(double h) {params=twoLinkArm::calcParams(weight,h,x0);}
	void setWeight(double w) {weight=w; params=twoLinkArm::calcParams(w,params.l1,params.l2,x0);}
	void setVirtualMass(double m) {virtualMass=m;}
	void setEarlyPulseGain(double g) {earlyPulseGain=g;}
	void setLatePulseGain(double g) {latePulseGain=g;}
	void setBLWNGain(double g) {blwnGain=g;}
	void setX0x(double p) {x0.X()=p;}
	void setX0y(double p) {x0.Y()=p;}
};

#endif
