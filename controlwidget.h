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
	QPushButton *startButton, *resetTGButton; 	
	QComboBox *stimulusBox;
	QFormLayout * layout;
	QLabel * probe0Label, *probe1Label, *probe2Label;
	QPushButton *probe0Button, *probe1Button, *probe2Button;
		
	DisplayWidget * userWidget;
	ArmSolver * armsolver;

	void closeEvent(QCloseEvent *event);
	void goGray() {for(std::vector<QWidget*>::iterator it=grayList.begin();it!=grayList.end();++it) (*it)->setEnabled(false); }
	void unGray() {for(std::vector<QWidget*>::iterator it=grayList.begin();it!=grayList.end();++it) (*it)->setEnabled(true); }
	void loadTrial(int T);
	QString makeProbeText(point P, int N);
	void writeCalib2File();
	
	QByteArray in,out;
	int inSize, outSize;
	QUdpSocket * us;
	QFile contFile, trialFile;
	QTextStream outStream;
	
	double resetTG, min, perturbGain, weight, eaGain, xpcTime, virtualMass, cursorFadeTime, extractionFadeTime, rawFadeTime, blwnGain, earlyPulseGain, latePulseGain;
	enum GameState {acquireTarget=0, inTarget=1, hold=2} state;
	std::vector<QWidget*> grayList;
	std::vector<DisplayWidget::Sphere> sphereVec;
	std::deque<timespec> times;
	std::deque<QByteArray> data;
	std::deque<point> extractedD, handleD, cursorD;
	DisplayWidget::Sphere sphere;
	
	twoLinkArm::ArmParams params;
	timespec zero, now, trialStart, targetAcquired, holdStart, lastFade;
	bool ExperimentRunning, inputReady, outputReady, ignoreInput, leftOrigin, firstpush, leftSide, hideCursor;
	int trial, subject, pulls, acquisitionsNeeded;
	point x0, origin, cursor, desposition, position, velocity, accel, claimedTarget, target, force, center;
	point probe0,probe1,probe2;
	
signals:
	void endApp();
	
public slots:
	void readPending();
	void startClicked();
	void resetTGClicked() {resetTG=1;}
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
	void setX0x(double p) {x0.X()=p; params=twoLinkArm::calcParams(weight,params.l1,params.l2,x0);}
	void setX0y(double p) {x0.Y()=p; params=twoLinkArm::calcParams(weight,params.l1,params.l2,x0);}
	void acquireProbe0() {probe0=position; writeCalib2File(); probe0Label->setText(makeProbeText(probe0,0)); center=probe0;}
	void acquireProbe1() {probe1=position; writeCalib2File(); probe1Label->setText(makeProbeText(probe1,1));}
	void acquireProbe2() {probe2=position; writeCalib2File(); probe2Label->setText(makeProbeText(probe2,2));}
};

#endif
