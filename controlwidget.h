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
	QDoubleSpinBox *gainBox, *eaGainBox, *l1Box, *l2Box, *massBox;
	QPushButton *startButton; 	
	QComboBox *stimulusBox, *acidBox;
	QFormLayout * layout;
		
	DisplayWidget * userWidget;
	ArmSolver * armsolver;

	void closeEvent(QCloseEvent *event);
	void goGray() {for(std::vector<QWidget*>::iterator it=grayList.begin();it!=grayList.end();++it) (*it)->setEnabled(false); }
	void unGray() {for(std::vector<QWidget*>::iterator it=grayList.begin();it!=grayList.end();++it) (*it)->setEnabled(true); }
	point loadTrial(int T);
	void noConsecutive(bool * array, int n);
	
	QByteArray in,out;
	int inSize, outSize;
	QUdpSocket * us;
	QFile contFile;
	QTextStream outStream;
	
	double sigGain, gain, min, smalls, bigs, eaGain, xpcTime, mass;
	enum stimuli {UNSTIMULATED=0, STIMULATED=1} stimulus;
	enum GameState {acquireTarget=0, inTarget=1, hold=2} state;
	enum AcidTrails {NEITHER=0, EXTRACTED=1, BOTH=2} trails;
	std::vector<QWidget*> grayList;
	std::vector<DisplayWidget::Sphere> sphereVec;
	std::deque<timespec> times;
	std::deque<QByteArray> data;
	std::deque<point> extracted, handle;
	DisplayWidget::Sphere sphere;
	
	twoLinkArm::ArmParams params;
	timespec zero, now, trialStart, targetAcquired, holdStart;
	bool ExperimentRunning, inputReady, outputReady, ignoreInput, leftOrigin, leftSide, firstpush;
	int trial, subject,lastStim, pulls;
	point origin, cursor, desposition, position, velocity, accel, target, force, center;
	
signals:
	void endApp();
	
public slots:
	void readPending();
	void startClicked();
	void setTrialNum(int i) {trial=i;}
	void setSubject(int i) {subject=i;}
	double evalSigmoid(double t, double risetime) {double a=10l/risetime; t-=risetime/2l; return (a*t/sqrt(1l+pow(a*t,2))+1l)/2l;}
	void setStimulus(int i) {stimulus=stimuli(i);}
	void setGain(double g) {sigGain=g;}
	void setEAGain(double g) {eaGain=g;}
	void setl1(double l) {params.l1=l;}
	void setl2(double l) {params.l2=l;}
	void setMass(double m) {mass=m;}
	void setAcid(int i) {trails=AcidTrails(i);}
};

#endif
