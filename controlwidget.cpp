#include "controlwidget.h"
#include "randb.h"
#include <QDesktopWidget>
#include <cstdio>
#include <sstream>

#define targetDuration .5
#define HOLDTIME .5
#define oRadius .01
#define cRadius .01
#define tRadius .01
#define TAB << "\t" <<

ControlWidget::ControlWidget(QDesktopWidget * qdw) : QWidget(qdw->screen(qdw->primaryScreen()))
{
	pulls=0;
	//Take care of window and input initialization.
	setFocus(); //Foreground window that gets all X input
	
	min=.49;
	cursor=point(0,.5);
	origin=point(0,.5);
	//state=acquireTarget;

	target=point(5,5);
	
	//Snag a UDP Socket and call a function (readPending) every time there's a new packet.
	ignoreInput=true;
	us=new QUdpSocket(this);
	us->bind(QHostAddress("192.168.1.1"),25000,QUdpSocket::DontShareAddress); //Bind the IP and socket you expect packets to be received from XPC on.
	connect(us, SIGNAL(readyRead()),this, SLOT(readPending()));
	
	//Develop GUI window
	layout = new QFormLayout(this);

	layout->addRow(startButton=new QPushButton("Start"));
	connect(startButton, SIGNAL(clicked()), this, SLOT(startClicked()));
	
	layout->addRow(tr("Subject Number:"), subjectBox=new QSpinBox(this));
	subjectBox->setValue(0);
	subjectBox->setMaximum(1000);
	subjectBox->setMinimum(0);
	grayList.push_back(subjectBox);
	subject=0;
	connect(subjectBox, SIGNAL(valueChanged(int)), this, SLOT(setSubject(int)));
	
	layout->addRow(tr("Trial Number:"), trialNumBox=new QSpinBox(this));
	trialNumBox->setValue(1);
	trialNumBox->setMaximum(10000);
	trialNumBox->setMinimum(0);
	grayList.push_back(trialNumBox);
	trial=1;
	connect(trialNumBox, SIGNAL(valueChanged(int)), this, SLOT(setTrialNum(int)));
	
	layout->addRow("EA Gain:",eaGainBox=new QDoubleSpinBox(this));
	eaGainBox->setValue(1);
	eaGainBox->setMaximum(5);
	eaGainBox->setMinimum(0);
	eaGainBox->setDecimals(3);
	eaGain=1;
	connect(eaGainBox, SIGNAL(valueChanged(double)), this, SLOT(setEAGain(double)));
	
	layout->addRow(tr("Cursor Fade Time (s):"), cursorFadeBox=new QDoubleSpinBox(this));
	cursorFadeBox->setValue(0);
	cursorFadeBox->setMaximum(5);
	cursorFadeBox->setMinimum(0);
	cursorFadeBox->setDecimals(3);
	cursorFadeTime=0;
	connect(cursorFadeBox, SIGNAL(valueChanged(double)), this, SLOT(setCursorFade(double)));
	
	layout->addRow(tr("Extraction Fade Time (s):"), extractionFadeBox=new QDoubleSpinBox(this));
	extractionFadeBox->setValue(0);
	extractionFadeBox->setMaximum(5);
	extractionFadeBox->setMinimum(0);
	extractionFadeBox->setDecimals(3);
	extractionFadeTime=0;
	connect(extractionFadeBox, SIGNAL(valueChanged(double)), this, SLOT(setExtractionFade(double)));
	
	layout->addRow(tr("Handle Fade Time (s):"), rawFadeBox=new QDoubleSpinBox(this));
	rawFadeBox->setValue(0);
	rawFadeBox->setMaximum(5);
	rawFadeBox->setMinimum(0);
	rawFadeBox->setDecimals(3);
	rawFadeTime=0;
	connect(rawFadeBox, SIGNAL(valueChanged(double)), this, SLOT(setRawFade(double)));
	
	layout->addRow(tr("Forearm Length (meters):"), armL2Box=new QDoubleSpinBox(this));
	armL2Box->setValue(.34);
	armL2Box->setMaximum(2);
	armL2Box->setMinimum(0);
	armL2Box->setDecimals(4);
	connect(armL2Box, SIGNAL(valueChanged(double)), this, SLOT(setl2(double)));
	
	layout->addRow(tr("Upper Arm Length (meters):"), armL1Box=new QDoubleSpinBox(this));
	armL1Box->setValue(.33);
	armL1Box->setMaximum(2);
	armL1Box->setMinimum(0);
	armL1Box->setDecimals(4);
	connect(armL1Box, SIGNAL(valueChanged(double)), this, SLOT(setl1(double)));
	
	layout->addRow(tr("Subject Weight (lbs):"), weightBox=new QDoubleSpinBox(this));
	weightBox->setMaximum(400);
	weightBox->setMinimum(0);
	weightBox->setDecimals(4);
	weightBox->setValue(160);
	weight=160;
	connect(weightBox, SIGNAL(valueChanged(double)), this, SLOT(setWeight(double)));
	
	layout->addRow(tr("Subject Shoulder Pos X (m):"), x0xBox=new QDoubleSpinBox(this));
	x0xBox->setValue(0);
	x0xBox->setMaximum(2);
	x0xBox->setMinimum(-2);
	x0xBox->setDecimals(4);
	x0.X()=0;
	connect(x0xBox, SIGNAL(valueChanged(double)), this, SLOT(setX0x(double)));
	
	layout->addRow(tr("Subject Shoulder Pos Y (m):"), x0yBox=new QDoubleSpinBox(this));
	x0yBox->setValue(.9);
	x0yBox->setMaximum(2);
	x0yBox->setMinimum(0);
	x0yBox->setDecimals(4);
	x0.Y()=.9;
	connect(x0yBox, SIGNAL(valueChanged(double)), this, SLOT(setX0y(double)));
	params=twoLinkArm::calcParams(160,.33,.34,x0);
	
	layout->addRow(tr("Virtual Mass (kg):"), virtualMassBox=new QDoubleSpinBox(this));
	virtualMassBox->setValue(5);
	virtualMassBox->setMaximum(10);
	virtualMassBox->setMinimum(-1);
	virtualMassBox->setDecimals(3);
	virtualMass=5;
	connect(virtualMassBox, SIGNAL(valueChanged(double)), this, SLOT(setVirtualMass(double)));
	                                                  
	layout->addRow(tr("Early Pulse Gain (N):"), earlyPulseGainBox=new QDoubleSpinBox(this));
	earlyPulseGainBox->setValue(0);
	earlyPulseGainBox->setMaximum(40);
	earlyPulseGainBox->setMinimum(-40);
	earlyPulseGainBox->setDecimals(3);
	earlyPulseGain=0;
	connect(earlyPulseGainBox, SIGNAL(valueChanged(double)), this, SLOT(setEarlyPulseGain(double)));
	
	layout->addRow(tr("Late Pulse Gain (N):"), latePulseGainBox=new QDoubleSpinBox(this));
	latePulseGainBox->setValue(0);
	latePulseGainBox->setMaximum(40);
	latePulseGainBox->setMinimum(-40);
	latePulseGainBox->setDecimals(3);
	latePulseGain=0;
	connect(latePulseGainBox, SIGNAL(valueChanged(double)), this, SLOT(setLatePulseGain(double)));
	
	layout->addRow(tr("Band-Limited White Stdev (N):"), blwnGainBox=new QDoubleSpinBox(this));
	blwnGainBox->setValue(0);
	blwnGainBox->setMaximum(20);
	blwnGainBox->setMinimum(0);
	blwnGainBox->setDecimals(3);
	blwnGain=0;
	connect(blwnGainBox, SIGNAL(valueChanged(double)), this, SLOT(setBLWNGain(double)));
	
	layout->addRow(resetTGButton=new QPushButton("Reset Target"));
	connect(resetTGButton, SIGNAL(clicked()), this, SLOT(resetTGClicked()));
	resetTG=1;
	
	//Need to populate with calibration data or make it
	QFile calibFile("calibration.dat");
	QTextStream calibStream;
	bool success=false;
	bool firstLine, secondLine, thirdLine;
	if(calibFile.exists())
	{
		calibFile.open(QIODevice::ReadOnly);
		char line[41];
		double tempx, tempy;
		calibFile.readLine(line,40);
		firstLine=sscanf(line, "%lf\t%lf",&tempx,&tempy);
		probe0=point(tempx,tempy);
		calibFile.readLine(line,40);
		secondLine=sscanf(line, "%lf\t%lf",&tempx,&tempy);
		probe1=point(tempx,tempy);
		calibFile.readLine(line,40);
		thirdLine=sscanf(line, "%lf\t%lf",&tempx,&tempy);
		probe2=point(tempx,tempy);
		if(firstLine&&secondLine&&thirdLine) success=true;
	}
	
	if(!success) //The file doesn't exist or has formatting issues
	{
		calibFile.open(QIODevice::WriteOnly);
		calibStream.setDevice(&calibFile);
		double zero=0;
		calibStream << zero TAB zero << endl;
		calibStream << double(LEFTPROBE) TAB zero << endl;
		calibStream << zero TAB double(UPPROBE) << endl;
		probe0=point(0,0);
		probe1=point(LEFTPROBE,0);
		probe2=point(0,UPPROBE);
	}
	center=probe0;
	
	//With data in hand, fill text labels with it, make buttons, and connect them
	layout->addRow(probe0Label=new QLabel(makeProbeText(probe0,0),this), probe0Button=new QPushButton("Acquire"));
	connect(probe0Button, SIGNAL(clicked()), this, SLOT(acquireProbe0()));
	layout->addRow(probe1Label=new QLabel(makeProbeText(probe1,1),this), probe1Button=new QPushButton("Acquire"));
	connect(probe1Button, SIGNAL(clicked()), this, SLOT(acquireProbe1()));
	layout->addRow(probe2Label=new QLabel(makeProbeText(probe2,2),this), probe2Button=new QPushButton("Acquire"));
	connect(probe2Button, SIGNAL(clicked()), this, SLOT(acquireProbe2()));		
	setLayout(layout);
	
	//Plop window in a sane place on the primary screen	
	QRect geo=qdw->screenGeometry();
	geo.setWidth(1*geo.width()/3);
	geo.setHeight(4*geo.height()/5);
	geo.translate(80,80);
	setGeometry(geo);
	
	//Take care of plopping subject's interface onto the secondary screen.
	int notprimary=qdw->primaryScreen()==0?1:0;
	
	userWidget=new DisplayWidget(qdw->screen(notprimary), true);
	userWidget->setGeometry(qdw->screenGeometry(notprimary));
	userWidget->show();
	
	userWidget->setDeepBGColor(point(0,0,0));
	
	inSize=0;
	perturbGain=1;
	
	//Initialize everything UPD-related to values that prevent problems
	ExperimentRunning=false;
}

QString ControlWidget::makeProbeText(point P, int N)
{
	QString buildme;
	std::ostringstream stringStream;
	switch(N)
	{
		case 0:
			stringStream << "Red";
			break;
		case 1:
			stringStream << "Green";
			break;
		case 2:
			stringStream << "Blue";
			break;
	}
	stringStream << " Target. X="<<P.X()<<", Y=" << P.Y();

	return buildme.fromStdString(stringStream.str());	
}

void ControlWidget::writeCalib2File()
{
	QFile calibFile("calibration.dat");
	QTextStream calibStream;
	calibFile.open(QIODevice::WriteOnly);
	calibStream.setDevice(&calibFile);
	calibStream << probe0.X() TAB probe0.Y() << endl;
	calibStream << probe1.X() TAB probe1.Y() << endl;
	calibStream << probe2.X() TAB probe2.Y() << endl;
}


void ControlWidget::readPending()
{
	now=getTime();
	
	int s=us->pendingDatagramSize();
	if(inSize != s) in.resize(s);
	us->readDatagram(in.data(), in.size());
	
	//Make sure pva is seeded.
	xpcTime=*reinterpret_cast<double*>(in.data());
	position.X()=*reinterpret_cast<double*>(in.data()+sizeof(double));
	position.Y()=*reinterpret_cast<double*>(in.data()+2*sizeof(double));
	velocity.X()=*reinterpret_cast<double*>(in.data()+3*sizeof(double));
	velocity.Y()=*reinterpret_cast<double*>(in.data()+4*sizeof(double));
	accel.X()=*reinterpret_cast<double*>(in.data()+5*sizeof(double));
	accel.Y()=*reinterpret_cast<double*>(in.data()+6*sizeof(double));
	force.X()=*reinterpret_cast<double*>(in.data()+7*sizeof(double));
	force.Y()=*reinterpret_cast<double*>(in.data()+8*sizeof(double));
	cursor=position;
	
	double white=blwnGain;
	double originTargetLine[4];
	originTargetLine[0]=0;
	originTargetLine[1]=0;
	originTargetLine[2]=0;
	originTargetLine[3]=-1; //Initialized pointing into the wall AND unreachable
	
	if(ignoreInput) //Send something back out so that XPC doesn't choke/stall/worse
	{
		out=QByteArray(in.data(),sizeof(double));//Copy the timestamp from the input
		out.append(reinterpret_cast<char*>(&virtualMass),sizeof(double));
		out.append(reinterpret_cast<char*>(&resetTG),sizeof(double));
		out.append(reinterpret_cast<char*>(&white),sizeof(double));
		out.append(reinterpret_cast<char*>(&earlyPulseGain),sizeof(double));
		out.append(reinterpret_cast<char*>(&latePulseGain),sizeof(double));
		out.append(reinterpret_cast<char*>(originTargetLine),4*sizeof(double));
		us->writeDatagram(out.data(),out.size(),QHostAddress("192.168.1.2"),25000);
		
		if(resetTG>0) resetTG=0;
		return;
	}
	
	armsolver->push(xpcTime, position, velocity, accel, accel*-virtualMass-force);
	//armsolver->push(xpcTime, position, velocity, accel, -force);
	armsolver->solve();

	if (!leftOrigin) trialStart=now;
	if ((!leftOrigin)&&(cursor.dist(origin)>(oRadius+cRadius))) leftOrigin=true;
	
	sphereVec.clear();
	//Target
	if(state!=hold)
	{
		if((state!=inTarget)||(!leftOrigin)) sphere.color=point(1,0,0); //Red
		else //Yellow -> too slow, White -> too fast
		{
			/* double acquire_time=targetAcquired-trialStart;
			if (acquire_time<.4) sphere.color=point(1,1,0);
			else if (acquire_time>.6) sphere.color=point(1,1,1);
			else sphere.color=point(0,1,0); */
			sphere.color=point(0,1,0); //Switched to "natural" pace.
		}
		sphere.position=target;
		sphere.radius=tRadius;
		sphereVec.push_back(sphere);
	} //In the hold case, the target is not drawn.
	
	//Cursor
	bool cursorDodgy;
	while(armsolver->pull(desposition, cursorDodgy, 0)) pulls++;
	//cursor=desposition*(1l-eaGain)+position*eaGain;
	double effectiveEA=1l+(eaGain-1l)*(1l/(1l+exp(-294.4439*(cursor.dist(target)-.02l))));
	cursor=desposition*(1l-effectiveEA)+position*effectiveEA;
	
	if((eaGain!=1)&&(cursorDodgy)) sphere.color=point(.5,.5,1); //Dark blue
	else sphere.color=point(0,0,1); //Blue
	sphere.position=cursor;
	sphere.radius=cRadius;
	if((!hideCursor)||(cursor.dist(target)<(tRadius+cRadius)))
		sphereVec.push_back(sphere);
	
	
	double fade;
	double fade2;
	double fadeL;
	if((now-lastFade)>=(1l/60l)) //Ie. Resolution no GREATER than seconds since we don't sample faster.
	{
		lastFade=now;
		handleD.push_back(position);
		fadeL=floor(60l*rawFadeTime);
		while(handleD.size()>fadeL) handleD.pop_front();
		extractedD.push_back(desposition);
		fadeL=floor(60l*extractionFadeTime);
		while(extractedD.size()>fadeL) extractedD.pop_front();
		cursorD.push_back(cursor);
		fadeL=floor(60l*cursorFadeTime);
		while(cursorD.size()>fadeL) cursorD.pop_front();
	}
	
	std::deque<point>::iterator it=handleD.begin();
	fade=0;
	fadeL=floor(60l*rawFadeTime);
	while(it!=handleD.end())
	{
		fade2=.5*(1.0+fade/fadeL);
		sphere.color=point(1,1,1)*fade2;
		sphere.position=*it;
		sphere.radius=cRadius*fade2;
		sphereVec.push_back(sphere);
		fade++;
		it++;
	}

	it=extractedD.begin();
	fade=0;
	fadeL=floor(60l*extractionFadeTime);
	while(it!=extractedD.end())
	{
		fade2=.5*(1.0+fade/fadeL);
		sphere.color=point(0,1,1)*fade2;
		sphere.position=*it;
		sphere.radius=cRadius*fade2;
		sphereVec.push_back(sphere);
		fade++;
		it++;
	}
	
	it=cursorD.begin();
	fade=0;
	fadeL=floor(60l*cursorFadeTime);
	while(it!=cursorD.end())
	{
		fade2=.5*(1.0+fade/fadeL);
		sphere.color=point(0,1,1)*fade2;
		sphere.position=*it;
		sphere.radius=cRadius*fade2;
		sphereVec.push_back(sphere);
		fade++;
		it++;
	}
	
	userWidget->setSpheres(sphereVec);
	
	QString showText;
		
	switch(state)
	{
	case acquireTarget:
		perturbGain=1;
		if ((cursor.dist(target)<(tRadius+cRadius))&&leftOrigin)
		{
			state=inTarget;
			targetAcquired=now;
		}
		break;
	case inTarget:
		perturbGain=0;
		if (acquisitionsNeeded<=1)
		{
			showText.clear();
			userWidget->setText(showText,point()); //Offscreen empty string
			if (cursor.dist(target)<(tRadius+cRadius))
			{
				if((now-targetAcquired)>=targetDuration)
				{
					origin=target;
					loadTrial(trial+1);
					leftOrigin=false;
				}
			}
			else state=acquireTarget;
		}
		else if (cursor.dist(target)>(tRadius+cRadius)) //Decrement on leaving and get state machine ready for next encounter with target.
		{
			acquisitionsNeeded--;
			showText.setNum(acquisitionsNeeded);
			userWidget->setText(showText,point(center.X()+min/6.0+.03,center.Y()+.05));
			state=acquireTarget;
		}
		break;
	case hold:
		if((now-holdStart)>HOLDTIME) state=acquireTarget;
		perturbGain=0; //evalSigmoid((now-holdStart),.4);
		break;
	}
	
	white=perturbGain*blwnGain;

	originTargetLine[0]=origin.X();
	originTargetLine[1]=origin.Y();
	originTargetLine[2]=claimedTarget.X();
	originTargetLine[3]=claimedTarget.Y();
	
	out=QByteArray(in.data(),sizeof(double));//Copy the timestamp from the input
	out.append(reinterpret_cast<char*>(&virtualMass),sizeof(double));
	out.append(reinterpret_cast<char*>(&resetTG),sizeof(double));
	out.append(reinterpret_cast<char*>(&white),sizeof(double));
	out.append(reinterpret_cast<char*>(&earlyPulseGain),sizeof(double));
	out.append(reinterpret_cast<char*>(&latePulseGain),sizeof(double));
	out.append(reinterpret_cast<char*>(&originTargetLine),4*sizeof(double));
	//This will require additional appends for other stimuli
	us->writeDatagram(out.data(),out.size(),QHostAddress("192.168.1.2"),25000);
	
	if (resetTG>0) resetTG=0;
		
	outStream << trial TAB now-zero TAB position.X() TAB position.Y() TAB velocity.X() TAB velocity.Y() TAB accel.X() TAB accel.Y() TAB force.X() TAB force.Y() TAB desposition.X() TAB desposition.Y() TAB xpcTime << endl;
}

void ControlWidget::startClicked()
{
	//Get a file open and recording...or not.
	if(subject>0)
	{
		char fname[200];
		std::sprintf(fname, "./Data/output%i.dat",subject);
		contFile.setFileName(fname);
		contFile.open(QIODevice::Append);
		outStream.setDevice(&contFile);
		std::sprintf(fname, "./Data/input.dat");
		trialFile.setFileName(fname);
		if(trialFile.exists()) trialFile.open(QIODevice::ReadOnly);
		else 
		{
			QMessageBox::critical(this, "File Not Found!", "File not found, please select a different file.");
			return; //Should not start experiment with bad input file.
		}
	}
	else
	{
		contFile.setFileName("/dev/null");
		contFile.open(QIODevice::Append);
		outStream.setDevice(&contFile);
	}
	userWidget->calibrate(probe0,probe1,probe2);
	//Get the armsolver class initialized with default blah.
	//armsolver=new ArmSolver(params,ArmSolver::CONSTIMP);
	armsolver=new ArmSolver(params,ArmSolver::CONSTIMP);
	
	//Make UI Changes
	userWidget->setDeepBGColor(point(0,0,0));
	userWidget->setShapes(false,false,false,false);
	goGray();
	startButton->setText("Experiment running...");
	
	loadTrial(trial);
	ExperimentRunning=true;
	ignoreInput=false;
	zero=getTime(); //Get first time point.
	lastFade=zero;
}

void ControlWidget::closeEvent(QCloseEvent *event)
{
	if(ExperimentRunning)
		if(QMessageBox::question(this, tr("For realz?"), tr("Do you really want to shutdown the experiment?"), QMessageBox::Yes| QMessageBox::Cancel)!=QMessageBox::Yes)
		{
			event->ignore();
			return;
		}
	emit(endApp());
	event->accept();
}

void ControlWidget::loadTrial(int T)
{
	state=hold;
	holdStart=now;
	
	if (subject>0)
	{
		trial=T;
		if(trialFile.atEnd()) {emit(endApp()); return;} //Hold begins so that app can end gracefully
	
		char line[201];
		std::string qline;
		int tempShowCursor,tempshape,temptrial,tempAcquisitionsNeeded;
		double tempx, tempy, tempEarly, tempLate, tempWhite, tempeaGain;
		std::cout << "Loading Trial " << T << std::endl;
		do
		{
			trialFile.readLine(line,200);
			std::cout << line << std::endl;
			if(sscanf(line, "%d\t%lf\t%lf\t%lf\t%lf\t%lf\t%d\t%d\t%d\t%lf",&temptrial,&tempx,&tempy,&tempEarly,&tempLate,&tempWhite,&tempshape,&tempShowCursor, &tempAcquisitionsNeeded,&tempeaGain));
			else
			{
				std::cout << "Complete failure to read line: " << line << std::endl; emit(endApp()); return;
			}
		} while ((temptrial < T)&&(!trialFile.atEnd()));
		acquisitionsNeeded=tempAcquisitionsNeeded;
		
		QString showText;
		if(acquisitionsNeeded>1)
		{
			showText.setNum(acquisitionsNeeded);
			userWidget->setText(showText,point(center.X()+min/6.0+.03,center.Y()+.05));
		}
		else
		{
			showText.clear();
			userWidget->setText(showText,point());
		}
				
		origin=target;
		target=point(tempx,tempy);
		earlyPulseGain=tempEarly;
		earlyPulseGainBox->setValue(tempEarly);
		latePulseGain=tempLate;
		latePulseGainBox->setValue(tempLate);
		blwnGain=tempWhite;
		blwnGainBox->setValue(tempWhite);
		eaGain=tempeaGain;
		eaGainBox->setValue(eaGain);
		userWidget->setShapes(tempshape==0,tempshape==1,tempshape==2,tempshape==3);
		if(tempshape>=0) {target=point(center.X()+min/6.0,center.Y()+.05); claimedTarget=center;}
		else claimedTarget=target;
		if(tempShowCursor>0) hideCursor=false;
		else hideCursor=true;
		
		trialNumBox->setValue(T);
		std::cout << "Finished Loading Trial " << temptrial << " " << acquisitionsNeeded << std::endl;
		
		armsolver->setParams(params); //Any weirdness with the shoulder should be gone.
	}
	else
	{
		userWidget->setShapes(false,false,false,false);
		if(trial==1)
		{
			leftSide=true;
		}
		
		double mean=.1;
		if(leftSide) {mean*=-1; leftSide=false;}
		else leftSide=true;
		origin=target;
		target=point(mean+randb(-.05,.05),center.Y());
		acquisitionsNeeded=1;
		trial=T;
		hideCursor=false;
		claimedTarget=target;
		armsolver->setParams(params); //Any weirdness with the shoulder should be gone.
	}
}



