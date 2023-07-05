#include "teg-peak.hpp"

//-----------------------------------------------------------------------------

void TegPeak::PrintFirstLineGraf()
{
	int i = (itsHeight-1)*2*itsBase-2;  // vertical arcs
		i += (itsBase-1)*2*itsHeight-2;	// horiziontal arcs
    cout<<"f 5 "<<itsDimension<<" "<<i<<" \n";
}

//-----------------------------------------------------------------------------

void TegPeak::CheckInput()
{
    bool ok=true;

    if (itsFS_s_Dep<0 || itsFS_s_Dep>3) ok=false;
	if (itsHArcDep==1 && itsSym) ok=false;
    if (itsWaitDep!=1 && itsWaitDep!=0) ok=false;
	if (itsFSCmin>itsFSCmax) ok=false;
	if (itsTmin>itsTmax) ok=false;
	if (itsTType<0 || itsTType>2) ok=false;
	if (itsTType==1 && itsSym) {itsSym=false; cout << "Can not use symmetric travel times with flag_T = 1. Set flag_sym = 0.";}
    if (itsStatic<0 || itsStatic>1) ok=false;
	if (itsCmin>itsCmax) ok=false;

    itsPeakLength = 2*itsPeakTrans+itsPeakPure;
    itsNonPeakLength = (itsCycleLength - itsPeakStart - itsPeaks*itsPeakLength)/(itsPeaks);     // peaks start after itsPeakStart and distributed evenly afterwards (modified v1.6)
    if (itsNonPeakLength<0) {cout << "Peaks to long to fit in the cycle!\n"; ok=false;}

    if (!ok)
    {
        cout<< "Error - Wrong input parameters specified\n"<<
               "See the manual for further details.\n";
        exit(1);
    }
}

//-----------------------------------------------------------------------------

void TegPeak::GenPenalties()
{
    int rc,t;

    if (itsVerbose)
    {
        cout << "\nGenerate penalty costs for node 1 ...\n";
        cout << "------------------------------------------------------------\n";
    }

    if (itsXml) {
        rc = xmlTextWriterStartElement(pXmlWriter, BAD_CAST "node"); PrintError(rc);
        rc = xmlTextWriterWriteAttribute(pXmlWriter, BAD_CAST "number",BAD_CAST xmlXPathCastNumberToString(1)); PrintError(rc);
    }

	for(t=0; t<itsT; t++)
    {
        if (itsFSCmax!=0)
        {
            if (itsFS_s_Dep==0)
		    {
			    WPair(itsFSCmin,itsFSCmax);
		    }
            if (itsFS_s_Dep==1 || itsFS_s_Dep==2 || itsFS_s_Dep==3)
		    {
			    PenaltyW(t);
		    }
        }
        else itsW1=itsW2=0;

        if (itsVerbose) cout << "Arrival at time "<<t<<" - costs ("<<itsW1<<","<<itsW2<<")\n";
        if (itsXml) {
            rc = xmlTextWriterStartElement(pXmlWriter, BAD_CAST "penalty"); PrintError(rc);
            rc = xmlTextWriterWriteAttribute(pXmlWriter, BAD_CAST "t",BAD_CAST xmlXPathCastNumberToString(t)); PrintError(rc);
            rc = xmlTextWriterWriteAttribute(pXmlWriter, BAD_CAST "c1",BAD_CAST xmlXPathCastNumberToString(itsW1)); PrintError(rc);
            rc = xmlTextWriterWriteAttribute(pXmlWriter, BAD_CAST "c2",BAD_CAST xmlXPathCastNumberToString(itsW2)); PrintError(rc);
            rc = xmlTextWriterEndElement(pXmlWriter); PrintError(rc);
        }
        if (itsF5) {
            f5File << " " << HNumber(1,t) << " -1 0 " << itsW1 << " " << itsW2 << " \n";
            itsStat.ma++;
        }
    }
    if (itsXml) rc = xmlTextWriterEndElement(pXmlWriter); PrintError(rc);
    if (itsVerbose) cout << "... done\n";
}

//-----------------------------------------------------------------------------

void TegPeak::PenaltyW(int t)
{
  int size;
  float T2, i1, i2;

  t=itsT-1-t;   // since in old version use reverse time interval

  size = itsFSCmax - itsFSCmin + 1; // size of the cost interval
  T2=(float)itsT-1; T2=T2/2; // T2 may be fractional

  if (itsFS_s_Dep==3)
  {
	  /* treat special cases separately */
	  if (t==0) { itsW1=itsFSCmax; itsW2=itsFSCmin; return; };
	  if (t==itsT-1) { itsW1=itsFSCmax; itsW2=itsFSCmin; return; };
	  if (t==T2) { itsW1=itsFSCmin; itsW2=itsFSCmax; return; };
	  if ( (t>T2) && (t-1<T2) )
		then { /* T odd, i in the middle */
				itsW1=itsFSCmin; itsW2=itsFSCmax; return;
			 };
	  /* now treat standard cases */
	  if (t > T2)
	  then { /* early arrival */ i1 = t-T2; }
	  else { /* late arrival  */ i1 = T2-t; }
	  i2 = T2 - i1; /* note i1+i2==T2 and i1,i2 < T2 */
	  itsW1 = itsFSCmin + (size * i1)/T2;
	  itsW2 = itsFSCmin + (size * i2)/T2;
  }
  if (itsFS_s_Dep==1)
  {
	  /* treat special cases separately */
	  if (t==0) { itsW1=itsW2=itsFSCmin; return; };
	  if (t==itsT-1) { itsW1=itsW2=itsFSCmin; return; };
	  if (t==T2) { itsW1=itsW2=itsFSCmax; return; };
	  if ( (t>T2) && (t-1<T2) )
		then { /* T odd, i in the middle */
				itsW1=itsW2=itsFSCmax; return;
			 };
	  /* now treat standard cases */
	  if (t > T2)
	  then { /* early arrival */ i1 = t-T2; }
	  else { /* late arrival  */ i1 = T2-t; }
	  i2 = T2 - i1; /* note i1+i2==T2 and i1,i2 < T2 */
	  itsW1 = itsW2 = itsFSCmin + (size * i2)/T2;
  }
  if (itsFS_s_Dep==2)
  {
	  /* treat special cases separately */
	  if (t==0) { itsW2=itsW1=itsFSCmax; return; };
	  if (t==itsT-1) { itsW2=itsW1=itsFSCmax; return; };
	  if (t==T2) { itsW2=itsW1=itsFSCmin; return; };
	  if ( (t>T2) && (t-1<T2) )
		then { /* T odd, i in the middle */
				itsW2=itsW1=itsFSCmin; return;
			 };
	  /* now treat standard cases */
	  if (t > T2)
	  then { /* early arrival */ i1 = t-T2; }
	  else { /* late arrival  */ i1 = T2-t; }
	  i2 = T2 - i1; /* note i1+i2==T2 and i1,i2 < T2 */
	  itsW2=itsW1 = itsFSCmin + (size * i1)/T2;
  }
}

//-----------------------------------------------------------------------------

int TegPeak::Peak(int t)
{
	int i,ub;

	ub = 0;
	t = CycleTime(t);
	for (i=1; i<=itsPeaks+1;i++)	// number of nonpeaks
	{
		// add nonpeak length
		if (i==1) ub = itsFirstNonPeakLength;
		if (i==itsPeaks+1) ub += itsLastNonPeakLength;
		if (i!=itsPeaks+1 && i!=1) ub += itsNonPeakLength;
		if (t<ub) return 0;
		// add first transient peak length
		ub += itsPeakTrans;
		if (t<ub) return 1;
		// add pure peak length
		ub += itsPeakPure;
		if (t<ub) return 2;
		// add last transient peak length
		ub += itsPeakTrans;
		if (t<ub) return 3;

	}
	cout << "Error in TegPeak::Peak\n";
	exit(1);
}

//-----------------------------------------------------------------------------

double TegPeak::FindMean(int t,int meanOffPeak)
{
	int i,ub;
	double mean,ps,f;

	t = CycleTime(t);	// time in first cycle

	ub = 0;
	mean = (double)meanOffPeak;
	for (i=1; i<=itsPeaks+1;i++)	// number of nonpeaks
	{
		// add nonpeak length
		if (i==1) ub = itsFirstNonPeakLength;
		if (i==itsPeaks+1) ub += itsLastNonPeakLength;
		if (i!=itsPeaks+1 && i!=1) ub += itsNonPeakLength;
		if (t<=ub) return mean;// if t in nonpeak
		// add first transient peak length
		ps = (double)ub;
		ub += itsPeakTrans;
		if (t<=ub)
		{
			f = (double)itsMeanIncrease*((double)t - ps)/(double)itsPeakTrans;
			return mean*(1+f/100);
		}
		// add pure peak part
		ub += itsPeakPure;
		if (t<=ub)
		{
			f = (double)itsMeanIncrease;
			return mean*(1+f/100);
		}
		// add last transient peak part
		ub += itsPeakTrans;
		if (t<=ub)
		{
			f = (double)itsMeanIncrease*( ((double)ub) - (double)t )/(double)itsPeakTrans;
			return mean*(1+f/100);
		}
	}
	cout << "Error in FindCenter\n";
	exit(1);
}

//-----------------------------------------------------------------------------

double TegPeak::FindCostPlus(int t,int costOffPeak)
{
	int i,ub;		// cycle number
	double cost_d,ps,f;

	ub = 0;
	cost_d = (double)costOffPeak;
	t=CycleTime(t);
	for (i=1; i<=itsPeaks+1;i++)	// number of nonpeaks
	{
		// add nonpeak length
		if (i==1) ub = itsFirstNonPeakLength;
		if (i==itsPeaks+1) ub += itsLastNonPeakLength;
		if (i!=itsPeaks+1 && i!=1) ub += itsNonPeakLength;
		if (t<ub) // if t in nonpeak
		{
			return cost_d;
		}
		ps = (double)ub;
		// add first transient peak length
		ub += itsPeakTrans;
		if (t<ub)
		{
			f = (double)itsMeanIncrease/100*(((double)t - ps)/(double)itsPeakTrans);
			return cost_d*(1+f);
		}
		// add pure peak part
		ub += itsPeakPure;
		if (t<ub)
		{
			return cost_d*(1+(double)itsMeanIncrease/100);
		}
		// add last transient peak part
		ub += itsPeakTrans;
		if (t<ub)
		{
			f = (double)itsMeanIncrease/100*(( ((double)ub) - (double)t )/(double)itsPeakTrans);
			return cost_d*(1+f);
		}
	}
	cout << "Error in FindCostPlus\n";
	exit(1);
}

//-----------------------------------------------------------------------------

double TegPeak::FindCostMinus(int t,int costOffPeak)
{
	int i,ub;
	double cost_d,ps,f;

	ub = 0;
	cost_d = (double)costOffPeak;
	t=CycleTime(t);
	for (i=1; i<=itsPeaks+1;i++)	// number of nonpeaks
	{
		// add nonpeak length
		if (i==1) ub = itsFirstNonPeakLength;
		if (i==itsPeaks+1) ub += itsLastNonPeakLength;
		if (i!=itsPeaks+1 && i!=1) ub += itsNonPeakLength;
		if (t<ub) // if t in nonpeak
		{
			return cost_d;
		}
		ps = (double)ub;
		// add first transient peak length
		ub += itsPeakTrans;
		if (t<ub)
		{
			f = (double)itsMeanIncrease/100*(((double)t - ps)/(double)itsPeakTrans);
			return cost_d*(1-f);
		}
		// add pure peak part
		ub += itsPeakPure;
		if (t<ub)
		{
			return cost_d*(1-(double)itsMeanIncrease/100);
		}
		// add last transient peak part
		ub += itsPeakTrans;
		if (t<ub)
		{
			f = (double)itsMeanIncrease/100*(( ((double)ub) - (double)t )/(double)itsPeakTrans);
			return cost_d*(1-f);
		}
	}
	cout << "Error in FindCostMinus\n";
	exit(1);
}

//-----------------------------------------------------------------------------

int TegPeak::GenTravelArc(int gTail,int gHead,int mean,bool out)
{
	int i,j,t,lb,ub,tmp1,tmp2,top,dir,rc,ctr;
	double meanNow,varNow,ubd,lbd;
	bool skip;

    ctr=rc=0;
	dir = Direction(gHead,gTail);
	if (itsVerbose && out)
	{
		cout << "\nGrid arc: (" << gTail << "," << gHead << ") - direction: " << dir;
		cout << " - mean off-peak travel time: " << mean << "\n";
	}

	for (t=0;t<itsT;t++)
	{
		tmp1=tmp2=0;
		skip = false;
		meanNow = FindMean(t,mean);
		if (itsTType==1 && (dir==1 || dir==2)) meanNow = mean;  // no peak effect on south and east arcs
		if (itsTType==2 && (dir==0 || dir==2)) meanNow = mean;  // no peak effect on north and south arcs
		varNow = (double)itsVarMean/100*meanNow;
		lbd = meanNow - varNow;
		if (lbd<1) lbd = 1;
		ubd = meanNow + varNow;
		lb = floor(lbd);
		ub = ceil(ubd);
		if (itsStatic==1 && Peak(t)==0) ub=lb=meanNow;    // static travel time in non-peak
		top = ub-lb+1;  // number of elements
		for (i=0;i<top;i++)	// generate travel times
		{
			itsHArcArray[i] = lb+i;
        }
        if (t+ub>=itsT) continue;  // if above timehorizon don't output/generate anything
        /* Stat */ itsStat.avePdfUb+=ub;
        ctr++;
        if (!out) continue;     // if preprocess then continue
        /* Stat */ itsStat.avePdfSize+=top; itsStat.minT=MIN(itsStat.minT,lb); itsStat.maxT=MAX(itsStat.maxT,ub);

        // generate prob
        if (ub-lb==0) itsMultArray[0]=1; // an arc
        else
        {
		    for (i=0;i<top;i++)
			   itsMultArray[i] = (int)1000000*NumberGen.BinomPdf(ub-lb,(meanNow-lb)/(ub-lb),i);
        }

		// output weights
		if (itsHArcDep==0)	// don't use random weights
		{
			itsW1 = itsGrid[gTail].c1[t][dir];
			itsW2 = itsGrid[gTail].c2[t][dir];
		}
		if (itsHArcDep==1 || itsHArcDep==2)
		{
			itsW1 = Rand(itsGrid[gTail].c1[1][dir]);
			itsW2 = Rand(itsGrid[gTail].c2[1][dir]);
		}
		if (itsHArcDep==5)
		{
			if (Peak(t)>0)
			{
				itsW1 = Rand(itsGrid[gTail].c1[2][dir]);
				itsW2 = Rand(itsGrid[gTail].c2[2][dir]);
			}
			else
			{
				itsW1 = Rand(itsGrid[gTail].c1[1][dir]);
				itsW2 = Rand(itsGrid[gTail].c2[1][dir]);
			}
		}
		if (itsHArcDep==3 || itsHArcDep==4 || itsHArcDep==6)
		{
			itsW1 = Rand(itsGrid[gTail].c1[t][dir]);
			itsW2 = Rand(itsGrid[gTail].c2[t][dir]);
		}
		if (itsStatic==1 && Peak(t)==0) {   // if static use the same off-peak cost no rand effect
			itsW1 = itsGrid[gTail].c1[1][dir];
			itsW2 = itsGrid[gTail].c2[1][dir];
		}
		/* Stat */ itsStat.minC=MIN(itsStat.minC,itsW1);itsStat.minC=MIN(itsStat.minC,itsW2);itsStat.maxC=MAX(itsStat.maxC,itsW1);itsStat.maxC=MAX(itsStat.maxC,itsW2);

        // output
        if (itsVerbose && out) cout<<"   t="<<t<<" mean="<<meanNow<<" [lb,ub]=["<<lb<<","<<ub<<"] peakIdx="<<Peak(t)<<" costs ("<<itsW1<<","<<itsW2<<")\n      ";

        if (out && itsXml) {
            rc = xmlTextWriterStartElement(pXmlWriter, BAD_CAST "leavingTime"); PrintError(rc);
            rc = xmlTextWriterWriteAttribute(pXmlWriter, BAD_CAST "t",BAD_CAST xmlXPathCastNumberToString(t)); PrintError(rc);
            rc = xmlTextWriterWriteAttribute(pXmlWriter, BAD_CAST "c1",BAD_CAST xmlXPathCastNumberToString(itsW1)); PrintError(rc);
            rc = xmlTextWriterWriteAttribute(pXmlWriter, BAD_CAST "c2",BAD_CAST xmlXPathCastNumberToString(itsW2)); PrintError(rc);
        }
        if (out & itsF5) f5File << " "<< HNumber(gTail,t) << " ";
		for (j=0;j<top;j++)   // removed (j=0;j<top && top!=1;j++) 2013-06-03, i.e. can print arcs now
		{
		    if (out && itsXml) {
                rc = xmlTextWriterStartElement(pXmlWriter, BAD_CAST "travelTime"); PrintError(rc);
                rc = xmlTextWriterWriteAttribute(pXmlWriter, BAD_CAST "t",BAD_CAST xmlXPathCastNumberToString(itsHArcArray[j])); PrintError(rc);
                rc = xmlTextWriterWriteAttribute(pXmlWriter, BAD_CAST "prob",BAD_CAST xmlXPathCastNumberToString(itsMultArray[j])); PrintError(rc);
                rc = xmlTextWriterEndElement(pXmlWriter); PrintError(rc);
		    }
            if (itsVerbose && out) cout <<"("<<itsHArcArray[j]<<","<<itsMultArray[j]<< ") ";
            if (out && itsF5) f5File << -HNumber(gHead,t+itsHArcArray[j]) << " ";
		}
		if (out && itsF5) {
            f5File << "0 " << itsW1 << " " << itsW2 << " ";
            for (j=0;j<top && top!=1;j++) f5File << itsMultArray[j] << " ";
            f5File << "\n";
            if (top>1) {itsStat.mh++; itsStat.hsize+=1+top;} else itsStat.ma++;
		}
        if (out && itsXml) rc = xmlTextWriterEndElement(pXmlWriter); PrintError(rc);
		if (itsVerbose && out) cout << "\n";
	}
    return ctr;
}

//-----------------------------------------------------------------------------

void TegPeak::GenWaiting()
{
	int i,t,rc;

	if (itsWaitMax<0) return; // generate no waiting arcs

    if (itsVerbose)
    {
        cout << "\nGenerate waiting costs ...\n";
        cout << "------------------------------------------------------------";
    }

    // generate waiting arcs, except at source and destination
    for(i=2; i < itsDimension; i++) // scan grid node numbers
    {
        if (itsXml) {
            rc = xmlTextWriterStartElement(pXmlWriter, BAD_CAST "node"); PrintError(rc);
            rc = xmlTextWriterWriteAttribute(pXmlWriter, BAD_CAST "number",BAD_CAST xmlXPathCastNumberToString(i));
        }

		if (itsVerbose)
        {
            cout<<"\nNode "<<i;
            if (itsWaitDep==0) cout<<" (off-peak cost = ("<<itsGrid[i].cWait1
                <<","<<itsGrid[i].cWait2<<"))\n";
            else cout<<endl;
        }

		for( t=0; t<itsT-1; t++)
		{
			// generate arc from (i,t) to (i,t+1): currently, k==number(i,t)
			if (itsWaitDep==1)
				WaitPair(itsWaitMin,itsWaitMax);	// W1 and W2 are now calculated
			else
			{
				itsW1 = RandWait(itsGrid[i].cWait1);
				itsW2 = RandWait(itsGrid[i].cWait2);
			}

		    if (itsVerbose) cout<<"   Wait from "<<t<<" to "<<t+1<< " - cost ("
                <<itsW1<<","<<itsW2<<")\n";

            if (itsXml) {
                rc = xmlTextWriterStartElement(pXmlWriter, BAD_CAST "wait"); PrintError(rc);
                rc = xmlTextWriterWriteAttribute(pXmlWriter, BAD_CAST "t",BAD_CAST xmlXPathCastNumberToString(t)); PrintError(rc);
                rc = xmlTextWriterWriteAttribute(pXmlWriter, BAD_CAST "time",BAD_CAST "1"); PrintError(rc);
                rc = xmlTextWriterWriteAttribute(pXmlWriter, BAD_CAST "c1",BAD_CAST xmlXPathCastNumberToString(itsW1)); PrintError(rc);
                rc = xmlTextWriterWriteAttribute(pXmlWriter, BAD_CAST "c2",BAD_CAST xmlXPathCastNumberToString(itsW2)); PrintError(rc);
                rc = xmlTextWriterEndElement(pXmlWriter); PrintError(rc);
            }
            if (itsF5) {
                f5File << " " << HNumber(i,t) << " " << -HNumber(i,t+1) << " 0 " << itsW1 << " " << itsW2 << "\n";
                itsStat.ma++;
            }
		}
        if (itsXml) rc = xmlTextWriterEndElement(pXmlWriter); PrintError(rc);
    }
    if (itsVerbose) cout << "... done\n";
}

//-----------------------------------------------------------------------------

TegPeak::TegPeak(string filename,bool verbose,bool f5,bool xml)
{
    //if (xml) itsFilename=filename+".xml";
    //if (f5) itsFilename=filename+".f5";
    itsFilename=filename;
    itsVerbose=verbose;
    itsF5=f5;
    itsXml=xml;

    cout<<"TEGP generator\n";
    cout<<"==============\n";
    cout<<"Reading input data ... ";

	// read input data
	cin >>
		itsBase >> itsHeight >>
        itsCycleLength >>
		itsPeaks >> itsPeakTrans >> itsPeakPure >> itsPeakStart >> itsT >>
		itsMeanIncrease >> itsVarMean >>
		itsFSCmin >> itsFSCmax >> itsFS_s_Dep >>
		itsTmin >> itsTmax >> itsTType >> itsStatic >>
		itsWaitMin >> itsWaitMax >> itsWaitDep >>
		itsCmin >> itsCmax >> itsHArcDep >>
		itsSym >>
		itsCorType >>
		itsRand >> itsSeed;

    cout<<"done. ";
    if (itsVerbose) PrintInput();
    CheckInput();
    Initialize();
}

//-----------------------------------------------------------------------------

void TegPeak::Run()
{
	double avePdfUb,mean,ub;
	int tmp;

    // first set the time horizon so all paths can be travelled
    mean = (double)itsTmax*(1+(double)itsMeanIncrease/100);
    ub = mean + ((double)itsVarMean*mean)/100;
    itsPeriods = MAX(1,((double)(itsBase+itsHeight-2)*ub)/(double)itsCycleLength);    // use this number as starting point
    tmp=itsT;   // store input value in tmp
    itsT = itsPeriods*itsCycleLength;
    /*if (tmp>itsT) {
        cout<<"To high time horizon used as input can be at most "<<itsT<<" due preallocation of memory.\n"
            <<"You may increase it in the source code (Run method).\n";
        exit(1);
    }*/

    AllocMemGrid();
    GenOffPeakMeanTraveltimes();
    GenWeights();
    avePdfUb=Preprocess();     // preprocess to find a lower time-horizon
    //cout<<"Ave pdf ub = "<<avePdfUb<<endl;

    itsT=tmp;
	// Now set the right number of time instances
	if (itsT<=0) {
        itsPeriods = ((double)(itsBase+itsHeight)*avePdfUb)/(double)itsCycleLength;
        itsT = (int)(itsPeriods*itsCycleLength);   // if time horizon not input
        if (itsVerbose) cout << "\nSet the time horizon to "<<itsT<<"\n";
	}
	else {
        itsPeriods = ((double)itsT)/(double)itsCycleLength;
        if (itsVerbose) cout << "\nUse time horizon "<<itsT<<"\n";
	}

	// regenerate weights, since itsT may be higher than the one used in the preprocessing step
    FreeMemGrid();
    AllocMemGrid();
    GenWeights();


	/*// calculate number of nodes, arcs, harcs and hsize
	itsN = (itsBase * itsHeight * itsT) + 1;	// number of nodes
	itsNArcs += itsT;		// itsT arcs in FS(s)
	if (itsWaitMax>=0) itsNArcs += (itsT-1)*(itsBase*itsHeight-2);	// waiting arcs
	CalcHsize();	// calc for harcs*/

    //itsStat.Reset();
    if (itsXml) xmlInit();
    if (itsF5) f5Init();
	//PrintData();			// print comments
	//PrintFirstLine();
	GenPenalties();
	GenWaiting();
	GenTravel(true);
    if (itsXml) {
        WriteStat();
        WriteTegpInfo();
        xmlFree();
    }
    if (itsF5) f5Free();

	cout<<"STDN generated.\n";
}

//-----------------------------------------------------------------------------


TegPeak::~TegPeak()
{
	delete [] itsGrid;
	delete [] itsHArcArray;
	delete [] itsMultArray;
}

//-----------------------------------------------------------------------------

void TegPeak::GenOffPeakMeanTraveltimes()
{
	int i;

    if (itsVerbose) cout << "Generate off-peak mean travel times ... ";
	// generate mean traveltimes
	for(i=itsDimension;i>0;i--)
	{
		if (itsSym)
		{
			if (North(i)) itsGrid[i].traveltime[0] = NumberGen.Int_length(itsTmin,itsTmax);
			if (West(i)) itsGrid[i].traveltime[3] = NumberGen.Int_length(itsTmin,itsTmax);
			NumberGen.Int_length(itsTmin,itsTmax);	// just so have same traveltimes under sym and not sym
			NumberGen.Int_length(itsTmin,itsTmax);
			if (East(i)) itsGrid[i].traveltime[1] = itsGrid[Head(i,1)].traveltime[3];
			if (South(i)) itsGrid[i].traveltime[2] = itsGrid[Head(i,2)].traveltime[0];
		}
		else
		{
			if (North(i)) itsGrid[i].traveltime[0] = NumberGen.Int_length(itsTmin,itsTmax);
			if (West(i)) itsGrid[i].traveltime[3] = NumberGen.Int_length(itsTmin,itsTmax);
			if (East(i)) itsGrid[i].traveltime[1] = NumberGen.Int_length(itsTmin,itsTmax);
			if (South(i)) itsGrid[i].traveltime[2] = NumberGen.Int_length(itsTmin,itsTmax);
		}
	}
    if (itsVerbose) cout << "done\n";
}


//-----------------------------------------------------------------------------

void TegPeak::GenWeights()
{
	int i,t,cUb1,cUb2,wTmp1,wTmp2;

    if (itsVerbose) cout << "\nGenerate weights (before rand applied) ... ";

	// generate and store waiting costs if time-independent
	for(i=itsDimension;i>0;i--)
	{
		WaitPair(itsWaitMin,itsWaitMax);
		itsGrid[i].cWait1 = itsW1;
		itsGrid[i].cWait2 = itsW2;
	}

	// generate off-peak costs
	for(i=itsDimension;i>0;i--)
	{
		if (itsHArcDep==0)
		{
			for(t=0;t<itsT;t++)
			{
				if (itsSym)
				{
					if (North(i))
					{
						WPair(itsCmin,itsCmax);
						itsGrid[i].c1[t][0] = itsW1; // cost of leaving node i at time t in dir 0
						itsGrid[i].c2[t][0] = itsW2;
					}
					if (West(i))
					{
						WPair(itsCmin,itsCmax);
						itsGrid[i].c1[t][3] = itsW1;
						itsGrid[i].c2[t][3] = itsW2;
					}
					if (East(i))
					{
						WPair(itsCmin,itsCmax);		// so get the same cost in the sym and nonsym case
						itsGrid[i].c1[t][1] = itsGrid[Head(i,1)].c1[t][3];
						itsGrid[i].c2[t][1] = itsGrid[Head(i,1)].c2[t][3];
					}
					if (South(i))
					{
						WPair(itsCmin,itsCmax);
						itsGrid[i].c1[t][2] = itsGrid[Head(i,2)].c1[t][0];
						itsGrid[i].c2[t][2] = itsGrid[Head(i,2)].c2[t][0];
					}
				}
				else
				{
					if (North(i))
					{
						WPair(itsCmin,itsCmax);
						itsGrid[i].c1[t][0] = itsW1; // cost of leaving node i at time t in dir 0
						itsGrid[i].c2[t][0] = itsW2;
					}
					if (West(i))
					{
						WPair(itsCmin,itsCmax);
						itsGrid[i].c1[t][3] = itsW1;
						itsGrid[i].c2[t][3] = itsW2;
					}
					if (East(i))
					{
						WPair(itsCmin,itsCmax);
						itsGrid[i].c1[t][1] = itsW1;
						itsGrid[i].c2[t][1] = itsW2;
					}
					if (South(i))
					{
						WPair(itsCmin,itsCmax);
						itsGrid[i].c1[t][2] = itsW1;
						itsGrid[i].c2[t][2] = itsW2;
					}
				}
			}
			continue;
		}

		if (itsHArcDep==1)
		{
			if (!itsSym)
			{
				t = 1;	// store in t=1

				// north arc
				WPair(itsCmin,itsCmax);
				itsGrid[i].c1[t][0] = itsW1;
				itsGrid[i].c2[t][0] = itsW2;
				cUb1 = itsW1;
				cUb2 = itsW2;

				// west arc
				WPair(itsCmin,itsCmax);
				itsGrid[i].c1[t][3] = itsW1;
				itsGrid[i].c2[t][3] = itsW2;
				cUb1 = MIN(itsW1,cUb1);
				cUb2 = MIN(itsW2,cUb2);

				// east arc
				itsGrid[i].c1[t][1] = NumberGen.Int_length(itsCmin,cUb1);
				itsGrid[i].c2[t][1] = NumberGen.Int_length(itsCmin,cUb2);

				// south arc
				itsGrid[i].c1[t][2] = NumberGen.Int_length(itsCmin,cUb1);
				itsGrid[i].c2[t][2] = NumberGen.Int_length(itsCmin,cUb2);
			}
			continue;;
		}

		if (itsHArcDep==2)
		{
			t=1;	// store in t=1
			if (itsSym)
			{
				if (North(i))
				{
					WPair(itsCmin,itsCmax);
					itsGrid[i].c1[t][0] = itsW1; // cost of leaving node i at time t in dir 0
					itsGrid[i].c2[t][0] = itsW2;
				}
				if (West(i))
				{
					WPair(itsCmin,itsCmax);
					itsGrid[i].c1[t][3] = itsW1;
					itsGrid[i].c2[t][3] = itsW2;
				}
				if (East(i))
				{
					WPair(itsCmin,itsCmax);
					itsGrid[i].c1[t][1] = itsGrid[Head(i,1)].c1[t][3];
					itsGrid[i].c2[t][1] = itsGrid[Head(i,1)].c2[t][3];
				}
				if (South(i))
				{
					WPair(itsCmin,itsCmax);
					itsGrid[i].c1[t][2] = itsGrid[Head(i,2)].c1[t][0];
					itsGrid[i].c2[t][2] = itsGrid[Head(i,2)].c2[t][0];
				}
			}
			else
			{
				if (North(i))
				{
					WPair(itsCmin,itsCmax);
					itsGrid[i].c1[t][0] = itsW1; // cost of leaving node i at time t in dir 0
					itsGrid[i].c2[t][0] = itsW2;
				}
				if (West(i))
				{
					WPair(itsCmin,itsCmax);
					itsGrid[i].c1[t][3] = itsW1;
					itsGrid[i].c2[t][3] = itsW2;
				}
				if (East(i))
				{
					WPair(itsCmin,itsCmax);
					itsGrid[i].c1[t][1] = itsW1;
					itsGrid[i].c2[t][1] = itsW2;
				}
				if (South(i))
				{
					WPair(itsCmin,itsCmax);
					itsGrid[i].c1[t][2] = itsW1;
					itsGrid[i].c2[t][2] = itsW2;
				}
			}
			continue;;
		}

		if (itsHArcDep==5)
		{
			if (itsSym)
			{
				if (North(i))
				{
					// store nonpeak cost in t=1 and peak cost in t=2
					WPair(itsCmin,itsCmax);
					itsGrid[i].c1[1][0] = itsW1; // cost nonpeak of leaving node i at time t in dir 0
					itsGrid[i].c2[1][0] = itsW2;
					wTmp1 = itsW1;
					wTmp2 = itsW2;
					if (wTmp1<wTmp2)
					{
						itsW1 = NumberGen.Int_length(wTmp1,wTmp2);
						itsW2 = NumberGen.Int_length(wTmp2,itsCmax);
					}
					else
					{
						itsW2 = NumberGen.Int_length(wTmp2,wTmp1);
						itsW1 = NumberGen.Int_length(wTmp1,itsCmax);
					}
					itsGrid[i].c1[2][0] = itsW1; // cost peak
					itsGrid[i].c2[2][0] = itsW2;
				}
				if (West(i))
				{
					WPair(itsCmin,itsCmax);
					itsGrid[i].c1[1][3] = itsW1;
					itsGrid[i].c2[1][3] = itsW2;
					wTmp1 = itsW1;
					wTmp2 = itsW2;
					if (wTmp1<wTmp2)
					{
						itsW1 = NumberGen.Int_length(wTmp1,wTmp2);
						itsW2 = NumberGen.Int_length(wTmp2,itsCmax);
					}
					else
					{
						itsW2 = NumberGen.Int_length(wTmp2,wTmp1);
						itsW1 = NumberGen.Int_length(wTmp1,itsCmax);
					}
					itsGrid[i].c1[2][3] = itsW1;
					itsGrid[i].c2[2][3] = itsW2;
					//cout << "i:" << i << " w1:" << itsW1 << " w2:" << itsW2 << "\n";
				}
				if (East(i))
				{
					WPair(itsCmin,itsCmax);
					wTmp1 = itsW1;
					wTmp2 = itsW2;
					if (wTmp1<wTmp2)
					{
						itsW1 = NumberGen.Int_length(wTmp1,wTmp2);
						itsW2 = NumberGen.Int_length(wTmp2,itsCmax);
					}
					else
					{
						itsW2 = NumberGen.Int_length(wTmp2,wTmp1);
						itsW1 = NumberGen.Int_length(wTmp1,itsCmax);
					}
					itsGrid[i].c1[1][1] = itsGrid[Head(i,1)].c1[1][3];
					itsGrid[i].c2[1][1] = itsGrid[Head(i,1)].c2[1][3];
					itsGrid[i].c1[2][1] = itsGrid[Head(i,1)].c1[2][3];
					itsGrid[i].c2[2][1] = itsGrid[Head(i,1)].c2[2][3];
				}
				if (South(i))
				{
					WPair(itsCmin,itsCmax);
					wTmp1 = itsW1;
					wTmp2 = itsW2;
					if (wTmp1<wTmp2)
					{
						itsW1 = NumberGen.Int_length(wTmp1,wTmp2);
						itsW2 = NumberGen.Int_length(wTmp2,itsCmax);
					}
					else
					{
						itsW2 = NumberGen.Int_length(wTmp2,wTmp1);
						itsW1 = NumberGen.Int_length(wTmp1,itsCmax);
					}
					itsGrid[i].c1[1][2] = itsGrid[Head(i,2)].c1[1][0];
					itsGrid[i].c2[1][2] = itsGrid[Head(i,2)].c2[1][0];
					itsGrid[i].c1[2][2] = itsGrid[Head(i,2)].c1[2][0];
					itsGrid[i].c2[2][2] = itsGrid[Head(i,2)].c2[2][0];
				}
			}
			else
			{
				if (North(i))
				{
					WPair(itsCmin,itsCmax);
					itsGrid[i].c1[1][0] = itsW1; // cost nonpeak of leaving node i at time t in dir 0
					itsGrid[i].c2[1][0] = itsW2;
					wTmp1 = itsW1;
					wTmp2 = itsW2;
					if (wTmp1<wTmp2)
					{
						itsW1 = NumberGen.Int_length(wTmp1,wTmp2);
						itsW2 = NumberGen.Int_length(wTmp2,itsCmax);
					}
					else
					{
						itsW2 = NumberGen.Int_length(wTmp2,wTmp1);
						itsW1 = NumberGen.Int_length(wTmp1,itsCmax);
					}
					itsGrid[i].c1[2][0] = itsW1; // cost peak
					itsGrid[i].c2[2][0] = itsW2;
				}
				if (West(i))
				{
					WPair(itsCmin,itsCmax);
					itsGrid[i].c1[1][3] = itsW1;
					itsGrid[i].c2[1][3] = itsW2;
					wTmp1 = itsW1;
					wTmp2 = itsW2;
					if (wTmp1<wTmp2)
					{
						itsW1 = NumberGen.Int_length(wTmp1,wTmp2);
						itsW2 = NumberGen.Int_length(wTmp2,itsCmax);
					}
					else
					{
						itsW2 = NumberGen.Int_length(wTmp2,wTmp1);
						itsW1 = NumberGen.Int_length(wTmp1,itsCmax);
					}
					itsGrid[i].c1[2][3] = itsW1;
					itsGrid[i].c2[2][3] = itsW2;
					//cout << "i:" << i << " w1:" << itsW1 << " w2:" << itsW2 << "\n";
				}
				if (East(i))
				{
					WPair(itsCmin,itsCmax);
					itsGrid[i].c1[1][1] = itsW1;
					itsGrid[i].c2[1][1] = itsW2;
					wTmp1 = itsW1;
					wTmp2 = itsW2;
					if (wTmp1<wTmp2)
					{
						itsW1 = NumberGen.Int_length(wTmp1,wTmp2);
						itsW2 = NumberGen.Int_length(wTmp2,itsCmax);
					}
					else
					{
						itsW2 = NumberGen.Int_length(wTmp2,wTmp1);
						itsW1 = NumberGen.Int_length(wTmp1,itsCmax);
					}
					itsGrid[i].c1[2][1] = itsW1;
					itsGrid[i].c2[2][1] = itsW2;
				}
				if (South(i))
				{
					WPair(itsCmin,itsCmax);
					itsGrid[i].c1[1][2] = itsW1;
					itsGrid[i].c2[1][2] = itsW2;
					wTmp1 = itsW1;
					wTmp2 = itsW2;
					if (wTmp1<wTmp2)
					{
						itsW1 = NumberGen.Int_length(wTmp1,wTmp2);
						itsW2 = NumberGen.Int_length(wTmp2,itsCmax);
					}
					else
					{
						itsW2 = NumberGen.Int_length(wTmp2,wTmp1);
						itsW1 = NumberGen.Int_length(wTmp1,itsCmax);
					}
					itsGrid[i].c1[2][2] = itsW1;
					itsGrid[i].c2[2][2] = itsW2;
				}
			}
			continue;
		}

		if (itsHArcDep==3)
		{
			if (itsSym)
			{
				if (North(i))
				{
					WPair(itsCmin,itsCmax);	// static costs nonpeak
					wTmp1 = itsW1;
					wTmp2 = itsW2;
					for(t=0;t<itsT;t++)
					{
						itsW1 = Round(FindCostPlus(t,wTmp1));
						itsW2 = Round(FindCostPlus(t,wTmp2));
						itsGrid[i].c1[t][0] = itsW1;
						itsGrid[i].c2[t][0] = itsW2;
					}
				}
				if (West(i))
				{
					WPair(itsCmin,itsCmax);	// static costs nonpeak
					wTmp1 = itsW1;
					wTmp2 = itsW2;
					for(t=0;t<itsT;t++)
					{
						itsW1 = Round(FindCostPlus(t,wTmp1));
						itsW2 = Round(FindCostPlus(t,wTmp2));
						itsGrid[i].c1[t][3] = itsW1;
						itsGrid[i].c2[t][3] = itsW2;
					}
				}
				if (East(i))
				{
					WPair(itsCmin,itsCmax);	// static costs nonpeak
					wTmp1 = itsW1;
					wTmp2 = itsW2;
					for(t=0;t<itsT;t++)
					{
						itsGrid[i].c1[t][1] = itsGrid[Head(i,1)].c1[t][3];
						itsGrid[i].c2[t][1] = itsGrid[Head(i,1)].c2[t][3];
					}
				}
				if (South(i))
				{
					WPair(itsCmin,itsCmax);	// static costs nonpeak
					wTmp1 = itsW1;
					wTmp2 = itsW2;
					for(t=0;t<itsT;t++)
					{
						itsGrid[i].c1[t][2] = itsGrid[Head(i,2)].c1[t][0];
						itsGrid[i].c2[t][2] = itsGrid[Head(i,2)].c2[t][0];
					}
				}
			}
			else
			{
				if (North(i))
				{
					WPair(itsCmin,itsCmax);	// static costs nonpeak
					wTmp1 = itsW1;
					wTmp2 = itsW2;
					for(t=0;t<itsT;t++)
					{
						itsW1 = Round(FindCostPlus(t,wTmp1));
						itsW2 = Round(FindCostPlus(t,wTmp2));
						itsGrid[i].c1[t][0] = itsW1;
						itsGrid[i].c2[t][0] = itsW2;
					}
				}
				if (West(i))
				{
					WPair(itsCmin,itsCmax);	// static costs nonpeak
					wTmp1 = itsW1;
					wTmp2 = itsW2;
					for(t=0;t<itsT;t++)
					{
						itsW1 = Round(FindCostPlus(t,wTmp1));
						itsW2 = Round(FindCostPlus(t,wTmp2));
						itsGrid[i].c1[t][3] = itsW1;
						itsGrid[i].c2[t][3] = itsW2;
					}
				}
				if (East(i))
				{
					WPair(itsCmin,itsCmax);	// static costs nonpeak
					wTmp1 = itsW1;
					wTmp2 = itsW2;
					for(t=0;t<itsT;t++)
					{
						itsW1 = Round(FindCostPlus(t,wTmp1));
						itsW2 = Round(FindCostPlus(t,wTmp2));
						itsGrid[i].c1[t][1] = itsW1;
						itsGrid[i].c2[t][1] = itsW2;
					}
				}
				if (South(i))
				{
					WPair(itsCmin,itsCmax);	// static costs nonpeak
					wTmp1 = itsW1;
					wTmp2 = itsW2;
					for(t=0;t<itsT;t++)
					{
						itsW1 = Round(FindCostPlus(t,wTmp1));
						itsW2 = Round(FindCostPlus(t,wTmp2));
						itsGrid[i].c1[t][2] = itsW1;
						itsGrid[i].c2[t][2] = itsW2;
					}
				}
			}
			continue;
		}

		if (itsHArcDep==4)
		{
			if (itsSym)
			{
				if (North(i))
				{
					WPair(itsCmin,itsCmax);	// static costs nonpeak
					wTmp1 = itsW1;
					wTmp2 = itsW2;
					for(t=0;t<itsT;t++)
					{
						itsW1 = Round(FindCostPlus(t,wTmp1));
						itsW2 = Round(FindCostMinus(t,wTmp2));
						itsGrid[i].c1[t][0] = itsW1;
						itsGrid[i].c2[t][0] = itsW2;
					}
				}
				if (West(i))
				{
					WPair(itsCmin,itsCmax);	// static costs nonpeak
					wTmp1 = itsW1;
					wTmp2 = itsW2;
					for(t=0;t<itsT;t++)
					{
						itsW1 = Round(FindCostPlus(t,wTmp1));
						itsW2 = Round(FindCostMinus(t,wTmp2));
						itsGrid[i].c1[t][3] = itsW1;
						itsGrid[i].c2[t][3] = itsW2;
					}
				}
				if (East(i))
				{
					WPair(itsCmin,itsCmax);	// static costs nonpeak
					wTmp1 = itsW1;
					wTmp2 = itsW2;
					for(t=0;t<itsT;t++)
					{
						itsGrid[i].c1[t][1] = itsGrid[Head(i,1)].c1[t][3];
						itsGrid[i].c2[t][1] = itsGrid[Head(i,1)].c2[t][3];
					}
				}
				if (South(i))
				{
					WPair(itsCmin,itsCmax);	// static costs nonpeak
					wTmp1 = itsW1;
					wTmp2 = itsW2;
					for(t=0;t<itsT;t++)
					{
						itsGrid[i].c1[t][2] = itsGrid[Head(i,2)].c1[t][0];
						itsGrid[i].c2[t][2] = itsGrid[Head(i,2)].c2[t][0];
					}
				}
			}
			else
			{
				if (North(i))
				{
					WPair(itsCmin,itsCmax);	// static costs nonpeak
					wTmp1 = itsW1;
					wTmp2 = itsW2;
					for(t=0;t<itsT;t++)
					{
						itsW1 = Round(FindCostPlus(t,wTmp1));
						itsW2 = Round(FindCostMinus(t,wTmp2));
						itsGrid[i].c1[t][0] = itsW1;
						itsGrid[i].c2[t][0] = itsW2;
					}
				}
				if (West(i))
				{
					WPair(itsCmin,itsCmax);	// static costs nonpeak
					wTmp1 = itsW1;
					wTmp2 = itsW2;
					for(t=0;t<itsT;t++)
					{
						itsW1 = Round(FindCostPlus(t,wTmp1));
						itsW2 = Round(FindCostMinus(t,wTmp2));
						itsGrid[i].c1[t][3] = itsW1;
						itsGrid[i].c2[t][3] = itsW2;
					}
				}
				if (East(i))
				{
					WPair(itsCmin,itsCmax);	// static costs nonpeak
					wTmp1 = itsW1;
					wTmp2 = itsW2;
					for(t=0;t<itsT;t++)
					{
						itsW1 = Round(FindCostPlus(t,wTmp1));
						itsW2 = Round(FindCostMinus(t,wTmp2));
						itsGrid[i].c1[t][1] = itsW1;
						itsGrid[i].c2[t][1] = itsW2;
					}
				}
				if (South(i))
				{
					WPair(itsCmin,itsCmax);	// static costs nonpeak
					wTmp1 = itsW1;
					wTmp2 = itsW2;
					for(t=0;t<itsT;t++)
					{
						itsW1 = Round(FindCostPlus(t,wTmp1));
						itsW2 = Round(FindCostMinus(t,wTmp2));
						itsGrid[i].c1[t][2] = itsW1;
						itsGrid[i].c2[t][2] = itsW2;
					}
				}
			}
			continue;
		}

		if (itsHArcDep==6)
		{
			if (itsSym)
			{
				if (North(i))
				{
					WPair(itsCmin,itsCmax);	// static costs nonpeak
					wTmp1 = itsW1;
					wTmp2 = itsW2;
					for(t=0;t<itsT;t++)
					{
						itsGrid[i].c1[t][0] = wTmp1;
						itsGrid[i].c2[t][0] = wTmp2;
					}
				}
				if (West(i))
				{
					WPair(itsCmin,itsCmax);	// static costs nonpeak
					wTmp1 = itsW1;
					wTmp2 = itsW2;
					for(t=0;t<itsT;t++)
					{
						itsW1 = Round(FindCostPlus(t,wTmp1));
						itsW2 = Round(FindCostPlus(t,wTmp2));
						itsGrid[i].c1[t][3] = itsW1;
						itsGrid[i].c2[t][3] = itsW2;
					}
				}
				if (East(i))
				{
					WPair(itsCmin,itsCmax);	// static costs nonpeak
					wTmp1 = itsW1;
					wTmp2 = itsW2;
					for(t=0;t<itsT;t++)
					{
						itsGrid[i].c1[t][1] = itsGrid[Head(i,1)].c1[t][3];
						itsGrid[i].c2[t][1] = itsGrid[Head(i,1)].c2[t][3];
					}
				}
				if (South(i))
				{
					WPair(itsCmin,itsCmax);	// static costs nonpeak
					wTmp1 = itsW1;
					wTmp2 = itsW2;
					for(t=0;t<itsT;t++)
					{
						itsGrid[i].c1[t][2] = itsGrid[Head(i,2)].c1[t][0];
						itsGrid[i].c2[t][2] = itsGrid[Head(i,2)].c2[t][0];
					}
				}
			}
			else
			{
				if (North(i))
				{
					WPair(itsCmin,itsCmax);	// static costs nonpeak
					wTmp1 = itsW1;
					wTmp2 = itsW2;
					for(t=0;t<itsT;t++)
					{
						itsGrid[i].c1[t][0] = wTmp1;
						itsGrid[i].c2[t][0] = wTmp2;
					}
				}
				if (West(i))
				{
					WPair(itsCmin,itsCmax);	// static costs nonpeak
					wTmp1 = itsW1;
					wTmp2 = itsW2;
					for(t=0;t<itsT;t++)
					{
						itsW1 = Round(FindCostPlus(t,wTmp1));
						itsW2 = Round(FindCostPlus(t,wTmp2));
						itsGrid[i].c1[t][3] = itsW1;
						itsGrid[i].c2[t][3] = itsW2;
					}
				}
				if (East(i))
				{
					WPair(itsCmin,itsCmax);	// static costs nonpeak
					wTmp1 = itsW1;
					wTmp2 = itsW2;
					for(t=0;t<itsT;t++)
					{
						itsW1 = Round(FindCostPlus(t,wTmp1));
						itsW2 = Round(FindCostPlus(t,wTmp2));
						itsGrid[i].c1[t][1] = itsW1;
						itsGrid[i].c2[t][1] = itsW2;
					}
				}
				if (South(i))
				{
					WPair(itsCmin,itsCmax);	// static costs nonpeak
					wTmp1 = itsW1;
					wTmp2 = itsW2;
					for(t=0;t<itsT;t++)
					{
						itsGrid[i].c1[t][2] = wTmp1;
						itsGrid[i].c2[t][2] = wTmp2;
					}
				}
			}
			continue;
		}

		cout << "ERROR: Non-existing harcDep type.\n";
		exit(1);
	}
    if (itsVerbose) cout << "done\n";
}

//-----------------------------------------------------------------------------

void TegPeak::WPair(int lb,int ub)
{
  int ub1, lb1, size;     // intervals for NETMAKER costs
  int SWAP;

  if (ub==0)
  {
    //cout << "\n\n\nub=0 used\n\n\n";
	itsW1 = itsW2 = 0;
	return;
  }

  switch (itsCorType)
   {
     case 0:  itsW1 = NumberGen.Int_weight(lb,ub);
              itsW2 = NumberGen.Int_weight(lb,ub);
              break;
     case 1:  itsW1 = NumberGen.Int_weight(lb,ub);
              itsW2 = ub + lb - itsW1;
              break;
     case 2:  itsW1 = NumberGen.Int_weight(lb,ub);
              itsW2 = ub + lb - itsW1;      /* complement of W1 */
              if (itsW2 < itsW1) then itsW2 = NumberGen.Int_weight(lb,itsW2);
              if (itsW2 > itsW1) then itsW2 = NumberGen.Int_weight(itsW2,ub);
              /* note that may be W1=W2=(ub+lb)2 if (ub+lb) is even */
              break;
     case 3:  if (ub == lb) lb1=ub1=lb;
              if (ub == lb+1) {ub1=lb; lb1=ub;};
              if (ub == lb+2) {ub1=lb; lb1=ub;};
              if (ub == lb+3) {ub1=lb+1; lb1=ub-1;};
              if (ub == lb+4) {ub1=lb+1; lb1=ub-1;};
              if (ub >= lb+5)
              then { size = (ub - lb +1)/3; /* size of the interval >= 2 */
                     ub1 = lb + (size-1);   /* size = ub1 - lb + 1 */
                     lb1 = ub - (size-1);   /* size = ub - lb1 + 1 */
                   };
              itsW1 = NumberGen.Int_weight(lb,ub1);
              itsW2 = NumberGen.Int_weight(lb1,ub);
              /* switch the two values with probability 1/2 */
              if (NumberGen.Sign()) then {SWAP = itsW1; itsW1=itsW2; itsW2=SWAP; };
              break;
	 default: cout << "corType wrong!!\n";
			  exit(1);
   }
   /*if (itsW1>itsCmax) {cout<<"Error in setting weights (above ub_T)!\n"; exit(1);}
   if (itsW2>itsCmax) {cout<<"Error in setting weights (above ub_T)!\n"; exit(1);}
   if (itsW1<itsCmin) {cout<<"Error in setting weights (above ub_T)!\n"; exit(1);}
   if (itsW2<itsCmin) {cout<<"Error in setting weights (above ub_T)!\n"; exit(1);}*/
}

//-----------------------------------------------------------------------------

void TegPeak::WaitPair(int lb,int ub)
{
  int ub1, lb1, size;     // intervals for NETMAKER costs
  int SWAP;

  if (ub==0)
  {
	itsW1 = itsW2 = 0;
	return;
  }

  switch (itsCorType)
   {
     case 0:  itsW1 = NumberGen.Int_number(lb,ub);
              itsW2 = NumberGen.Int_number(lb,ub);
              break;
     case 1:  itsW1 = NumberGen.Int_number(lb,ub);
              itsW2 = ub + lb - itsW1;
              break;
     case 2:  itsW1 = NumberGen.Int_number(lb,ub);
              itsW2 = ub + lb - itsW1;      /* complement of W1 */
              if (itsW2 < itsW1) then itsW2 = NumberGen.Int_number(lb,itsW2);
              if (itsW2 > itsW1) then itsW2 = NumberGen.Int_number(itsW2,ub);
              /* note that may be W1=W2=(ub+lb)2 if (ub+lb) is even */
              break;
     case 3:  if (ub == lb) lb1=ub1=lb;
              if (ub == lb+1) {ub1=lb; lb1=ub;};
              if (ub == lb+2) {ub1=lb; lb1=ub;};
              if (ub == lb+3) {ub1=lb+1; lb1=ub-1;};
              if (ub == lb+4) {ub1=lb+1; lb1=ub-1;};
              if (ub >= lb+5)
              then { size = (ub - lb +1)/3; /* size of the interval >= 2 */
                     ub1 = lb + (size-1);   /* size = ub1 - lb + 1 */
                     lb1 = ub - (size-1);   /* size = ub - lb1 + 1 */
                   };
              itsW1 = NumberGen.Int_number(lb,ub1);
              itsW2 = NumberGen.Int_number(lb1,ub);
              /* switch the two values with probability 1/2 */
              if (NumberGen.Sign()) then {SWAP = itsW1; itsW1=itsW2; itsW2=SWAP; };
              break;
	 default: cout << "itsCorType wrong!!\n";
			  exit(1);
   }
}

//-----------------------------------------------------------------------------

void TegPeak::WriteStat()
{
    int rc;
    string tmp;
    char tmp1[20];

    rc = xmlTextWriterStartElement(pXmlWriter, BAD_CAST "statistics"); PrintError(rc);
    rc = xmlTextWriterStartElement(pXmlWriter, BAD_CAST "pdfSize"); PrintError(rc);
    rc = xmlTextWriterWriteAttribute(pXmlWriter, BAD_CAST "min",BAD_CAST xmlXPathCastNumberToString(itsStat.minT));
    rc = xmlTextWriterWriteAttribute(pXmlWriter, BAD_CAST "ave",BAD_CAST xmlXPathCastNumberToString(itsStat.avePdfUb));
    rc = xmlTextWriterWriteAttribute(pXmlWriter, BAD_CAST "max",BAD_CAST xmlXPathCastNumberToString(itsStat.maxT));
    rc = xmlTextWriterEndElement(pXmlWriter); PrintError(rc);

    rc = xmlTextWriterStartElement(pXmlWriter, BAD_CAST "travelTimes"); PrintError(rc);
    rc = xmlTextWriterWriteAttribute(pXmlWriter, BAD_CAST "min",BAD_CAST xmlXPathCastNumberToString(itsStat.minT));
    rc = xmlTextWriterWriteAttribute(pXmlWriter, BAD_CAST "max",BAD_CAST xmlXPathCastNumberToString(itsStat.maxT));
    sprintf(tmp1,"%d",itsStat.minI);
    tmp = tmp1;
    tmp += "-";
    sprintf(tmp1,"%d",itsStat.maxI);
    tmp += tmp1;
    rc = xmlTextWriterWriteAttribute(pXmlWriter, BAD_CAST "maxInterval",BAD_CAST tmp.c_str());
    rc = xmlTextWriterEndElement(pXmlWriter); PrintError(rc);
    rc = xmlTextWriterEndElement(pXmlWriter); PrintError(rc);
}

//-----------------------------------------------------------------------------

void TegPeak::WriteTegpInfo()
{
    int rc;
    string tmp;

    rc = xmlTextWriterStartElement(pXmlWriter, BAD_CAST "tegpInfo"); PrintError(rc);
        rc = xmlTextWriterStartElement(pXmlWriter, BAD_CAST "gridInfo"); PrintError(rc);
        rc = xmlTextWriterWriteAttribute(pXmlWriter, BAD_CAST "base",BAD_CAST xmlXPathCastNumberToString(itsBase));
        rc = xmlTextWriterWriteAttribute(pXmlWriter, BAD_CAST "height",BAD_CAST xmlXPathCastNumberToString(itsHeight));
        rc = xmlTextWriterEndElement(pXmlWriter); PrintError(rc);

        rc = xmlTextWriterStartElement(pXmlWriter, BAD_CAST "peakInfo"); PrintError(rc);
        rc = xmlTextWriterWriteAttribute(pXmlWriter, BAD_CAST "peaks",BAD_CAST xmlXPathCastNumberToString(itsPeaks));
        rc = xmlTextWriterWriteAttribute(pXmlWriter, BAD_CAST "tTransistent",BAD_CAST xmlXPathCastNumberToString(itsPeakTrans));
        rc = xmlTextWriterWriteAttribute(pXmlWriter, BAD_CAST "tPure",BAD_CAST xmlXPathCastNumberToString(itsPeakPure));
        rc = xmlTextWriterWriteAttribute(pXmlWriter, BAD_CAST "tFirst",BAD_CAST xmlXPathCastNumberToString(itsPeakStart));
        rc = xmlTextWriterEndElement(pXmlWriter); PrintError(rc);

        rc = xmlTextWriterStartElement(pXmlWriter, BAD_CAST "offPeakInfo"); PrintError(rc);
        rc = xmlTextWriterWriteAttribute(pXmlWriter, BAD_CAST "static",BAD_CAST xmlXPathCastNumberToString(itsStatic));
        rc = xmlTextWriterEndElement(pXmlWriter); PrintError(rc);

        rc = xmlTextWriterStartElement(pXmlWriter, BAD_CAST "cycleInfo"); PrintError(rc);
        rc = xmlTextWriterWriteAttribute(pXmlWriter, BAD_CAST "cycleLength",BAD_CAST xmlXPathCastNumberToString(itsCycleLength));
        rc = xmlTextWriterWriteAttribute(pXmlWriter, BAD_CAST "cycles",BAD_CAST xmlXPathCastNumberToString(itsPeriods));
        rc = xmlTextWriterEndElement(pXmlWriter); PrintError(rc);

    rc = xmlTextWriterEndElement(pXmlWriter); PrintError(rc);

}

//-----------------------------------------------------------------------------

void TegPeak::PrintInput()
{
  cout << "\nInput data for the TEGP generator\n";
  cout << "------------------------------------------------------------\n\n";
  cout << "Grid base and height           : " << itsBase << " " << itsHeight << "\n";
  cout << "Time instances in a cycle      : " << itsCycleLength << "\n";
  cout << "Peaks in a cycle               : " << itsPeaks<< "\n";
  cout << "Transient period length        : " << itsPeakTrans<<"\n";
  cout << "Pure peak length               : " << itsPeakPure<< "\n";
  cout << "First Peak starts at time      : " << itsPeakStart<< "\n";
  if (itsT>0) cout << "Use fixed time horizon of      : " << itsT<< "\n";
  cout << "Mean increase in peaks         : " << itsMeanIncrease << "\n";
  cout << "Std deviation/mean ratio       : " << itsVarMean << "\n";
  cout << "Penalty cost interval          : [" << itsFSCmin << "," << itsFSCmax << "]\n";
  cout << "Penalty flag                   : " << itsFS_s_Dep<< "\n";
  cout << "Travel time interval           : [" << itsTmin  << "," << itsTmax << "]\n";
  cout << "Static flag                    : " << itsStatic << "\n";
  cout << "Waitingarcs cost interval      : [" << itsWaitMin << "," << itsWaitMax << "]\n";
  cout << "Waiting flag                   : " << itsWaitDep << "\n";
  cout << "Travel cost interval           : [" << itsCmin << "," << itsCmax << "]\n";
  cout << "Travel cost flag               : " << itsHArcDep<< "\n";
  cout << "Symmetric flag                 : " << itsSym << "\n";
  cout << "Correlation flag               : " << itsCorType<< "\n";
  cout << "Random element (promille)      : "<<itsRand<<"\n";
  cout << "Seed                           : "<<itsSeed<<"\n";
  cout << "------------------------------------------------------------\n\n";
}

//-----------------------------------------------------------------------------

double TegPeak::GenTravel(bool out)
{
    int x,y,gTail,gHead,mean,rc,ctr;

    if (itsVerbose && out)
    {
        cout << "\nGenerate costs and travel time distributions ...\n";
        cout << "------------------------------------------------------------";
    }
    ctr=rc=0;
    /* Stat */ itsStat.avePdfUb=0;
    // generate other arcs and hyperarcs
    for(x=1; x <= itsBase; x++) // column counter
    {
        for(y=1; y <= itsHeight; y++)   // row counter
        {
            gTail = GridNumber(x,y); // gTail is the number of grid point (x,y)
            if (y > 1)
            {	// north arc in the grid,
                gHead = gTail - 1;
                mean = itsGrid[gTail].traveltime[0];
                //cout << "mean: " << mean << "\n";
                if (out && itsXml) startXmlArc(gHead,gTail);
                ctr+=GenTravelArc(gTail,gHead,mean,out);
                if (out && itsXml) endXmlArc();
            }
            if (x > 1)
            {	// west arc in the grid,
                gHead = gTail - itsHeight;
                mean = itsGrid[gTail].traveltime[3];
                if (out && itsXml) startXmlArc(gHead,gTail);
                ctr+=GenTravelArc(gTail,gHead,mean,out);
                if (out && itsXml) endXmlArc();
            }
            if ((y < itsHeight) && (gTail != itsDimension-1) && (gTail != 1) )
            {	// south arc in the grid,
                // NOTE: no south arc into source (gsize) and from destination (1)
                gHead = gTail + 1;
                if (itsSym) mean = itsGrid[gHead].traveltime[0]; // consider only north traveltimes
                else mean = itsGrid[gTail].traveltime[2];
                if (out && itsXml) startXmlArc(gHead,gTail);
                ctr+=GenTravelArc(gTail,gHead,mean,out);
                if (out && itsXml) endXmlArc();
            }
            if ((x < itsBase) && (gTail != itsDimension-itsHeight) && (gTail != 1) )
            {	// east arc in the grid
                // NOTE: no east arc into source (gsize) and from destination (1)
                gHead = gTail + itsHeight;
                if (itsSym) mean = itsGrid[gHead].traveltime[3]; // consider only west traveltimes
                else mean = itsGrid[gTail].traveltime[1];
                //cout << "mean: " << mean << "\n";
                if (out && itsXml) startXmlArc(gHead,gTail);
                ctr+=GenTravelArc(gTail,gHead,mean,out);
                if (out && itsXml) endXmlArc();
            }
        }
    }
    itsStat.avePdfUb=itsStat.avePdfUb/(double)ctr;
    itsStat.avePdfSize=itsStat.avePdfSize/(double)ctr;
    if (itsVerbose && out) cout << "done\n";;
    return itsStat.avePdfUb;
}

//-----------------------------------------------------------------------------

void TegPeak::GenerateGraf()
{
	int x,y,gTail,gHead;

    // generate arcs for graf representation
    for(x=1; x <= itsBase; x++) // column counter
	{
		for(y=1; y <= itsHeight; y++)   // row counter
		{
			gTail = GridNumber(x,y); // gTail is the number of grid point (x,y)

			if (y > 1)
			{	// north arc in the grid,
				gHead = gTail - 1;
				cout << " " << gHead << " " << gTail << "\n";
			}
			if (x > 1)
			{	// west arc in the grid,
				gHead = gTail - itsHeight;
				cout << " " << gHead << " " << gTail << "\n";
			}
			if ((y < itsHeight) && (gTail != itsDimension-1) && (gTail != 1) )
			{	// south arc in the grid,
				// NOTE: no south arc into source (gsize) and from destination (1)
				gHead = gTail + 1;
				cout << " " << gHead << " " << gTail << "\n";
			}
			if ( (x < itsBase) && (gTail != itsDimension-itsHeight) && (gTail != 1) )
			{	// east arc in the grid
				// NOTE: no east arc into source (gsize) and from destination (1)
				gHead = gTail + itsHeight;
				cout << " " << gHead << " " << gTail << "\n";
			}
       }
	}
}

//----------------------------------------------------------------------------

void TegPeak::xmlInit()
{
    int arcs;
    const char * MY_ENCODING="ISO-8859-1";
    int rc;
    string filename = itsFilename + ".xml";

    LIBXML_TEST_VERSION

    pXmlWriter = xmlNewTextWriterFilename(filename.c_str(), 0);
    if (pXmlWriter == NULL) {
        cout<<"Error creating the xml writer for output.\n";
        return;
    }
    xmlTextWriterSetIndent(pXmlWriter, TRUE);
    const xmlChar MY_INDENT_STR[5] = "    ";	// Your indent String is 5 blanks
    xmlTextWriterSetIndentString(pXmlWriter, MY_INDENT_STR);

    /* Start the document with the xml default for the version,
     * encoding ISO 8859-1 and the default for the standalone
     * declaration. */
    rc = xmlTextWriterStartDocument(pXmlWriter, NULL, MY_ENCODING, NULL); PrintError(rc);

    rc = xmlTextWriterStartElement(pXmlWriter, BAD_CAST "stdn"); PrintError(rc);
    rc = xmlTextWriterWriteAttribute(pXmlWriter, BAD_CAST "nodes",BAD_CAST xmlXPathCastNumberToString(itsDimension)); PrintError(rc);
    arcs=2*itsHeight*(itsBase-1)+2*itsBase*(itsHeight-1)-4;
    rc = xmlTextWriterWriteAttribute(pXmlWriter, BAD_CAST "arcs",BAD_CAST xmlXPathCastNumberToString(arcs)); PrintError(rc);
    rc = xmlTextWriterWriteAttribute(pXmlWriter, BAD_CAST "timeHorizon",BAD_CAST xmlXPathCastNumberToString(itsT)); PrintError(rc);

    //char file[50];
    //if (strrchr(filename.c_str(),'/')!=NULL) strcpy(file,(strrchr(filename.c_str(),'/')+1));
    rc = xmlTextWriterWriteAttribute(pXmlWriter, BAD_CAST "name",BAD_CAST filename.c_str()); PrintError(rc);
}

//----------------------------------------------------------------------------

void TegPeak::xmlFree()
{
    int rc;
    rc = xmlTextWriterEndDocument(pXmlWriter); PrintError(rc);
    xmlFreeTextWriter(pXmlWriter);
}

//----------------------------------------------------------------------------

void TegPeak::startXmlArc(int gHead,int gTail) {
    int rc;
    rc = xmlTextWriterStartElement(pXmlWriter, BAD_CAST "arc"); PrintError(rc);
    rc = xmlTextWriterWriteAttribute(pXmlWriter, BAD_CAST "head",BAD_CAST xmlXPathCastNumberToString(gHead)); PrintError(rc);
    rc = xmlTextWriterWriteAttribute(pXmlWriter, BAD_CAST "tail",BAD_CAST xmlXPathCastNumberToString(gTail)); PrintError(rc);
}

//----------------------------------------------------------------------------

void TegPeak::endXmlArc() {
    int rc;
    rc = xmlTextWriterEndElement(pXmlWriter); PrintError(rc);
}

//----------------------------------------------------------------------------

void TegPeak::f5Init()
{
    f5File.open("stdn.tmp");
}

//----------------------------------------------------------------------------

void TegPeak::f5Free()
{
    f5File.close();
    ifstream source("stdn.tmp");
    string dName = itsFilename+".f5";
    ofstream dest(dName.c_str());
    dest << "c Time Expanded Generator with Peaks (TEGP)             \n";
    dest << "c                                                       \n";
    if (itsSym) dest << "c      Use symmetric weights\n";
    else dest << "c      Use non symmetric weights\n";
    if (itsWaitMax<0) dest << "c      No waiting arcs allowed\n";
    dest << "c\n";
    dest << "c      Grid base and height           : " << itsBase << " " << itsHeight << "\n";
    dest << "c      Time instances in a cycle      : " << itsCycleLength << "\n";
    dest << "c      Peaks in a cycle               : " << itsPeaks<< "\n";
    dest << "c      Time instances in a Peak       : " << itsPeakLength<< "\n";
    dest << "c      Transient period in the peak   : " << itsPeakTrans<< "\n";
    dest << "c      Pure peak period               : " << itsPeakPure<< "\n";
    dest << "c      First Peak starts at time      : " << itsFirstNonPeakLength+1<< "\n";
    dest << "c      Time between Peaks             : " << itsNonPeakLength<< "\n";
    dest << "c      Number of cycles               : " << itsPeriods << "\n";
    dest << "c      Total number of time instances : " << itsT << "\n";
    dest << "c\n";
 	dest << "f 5 " << itsDimension*itsT+1 << " " <<  itsStat.ma << " " << itsStat.mh <<
		" " << itsStat.maxPdfSize << " " << itsStat.hsize << " " << itsBase << " " << itsHeight <<
		" " << itsCycleLength << " " << itsPeaks << " " << itsPeakTrans <<
		" " << itsPeakPure << " " << itsPeakStart << " " << itsT <<
		" " << (itsWaitMax>=0) << " \n";
    char line[500];
    if (source.is_open()){
        while (!source.eof()){
                source.getline(line,500,'\n');
            dest << line << endl;
        }
    }
    else{
        cout << "File couldn't be opened." << endl;
    }
    source.close();
    dest.close();
}

// -----------------------------------------------------------------------------

void TegPeak::Initialize()
{
    double mean,var;
    int lb,ub;

    if (itsSeed == 0) then itsSeed = NumberGen.Clock_seed();
    itsDimension = itsBase * itsHeight;
    itsGrid = new GridNode[itsBase*itsHeight+1];    // don't use entry 0
    itsPeakLength = 2*itsPeakTrans+itsPeakPure;
    itsNonPeakLength = (itsCycleLength - itsPeakStart - itsPeaks*itsPeakLength)/(itsPeaks);     // peaks start after itsPeakStart and distributed evenly afterwards (modified v1.6)
    itsFirstNonPeakLength = itsPeakStart;
    itsLastNonPeakLength = (itsCycleLength - itsPeakStart - itsPeaks*itsPeakLength) - (itsPeaks-1)*itsNonPeakLength;    // in v1.5 = itsNonPeakLength - itsFirstNonPeakLength;

    // calc max pdf size
    mean = (double)itsTmax*(1+(double)itsMeanIncrease/100);
    var = itsVarMean;
    lb = floor(mean - ((double)var*mean)/100);
    if (lb<1) lb = 1;
    ub = ceil(mean + ((double)var*mean)/100);
    itsHArcArray = new int[ub-lb+2];
    itsMultArray = new int[ub-lb+2];
    /* Stat */ itsStat.minI=lb; itsStat.maxI=ub; itsStat.maxPdfSize=ub-lb+1;

    NumberGen.Init_len(itsSeed);   // used to generate center traveltimes
    NumberGen.Init_num(itsSeed);   // used to generate waiting arc costs
    NumberGen.Init_sign(itsSeed);  // used for if cor = 3
    NumberGen.Init_w(itsSeed);     // used to generate weights (both for harcs and FS(s) arcs)
}

//-----------------------------------------------------------------------------









//-----------------------------------------------------------------------------

/** Print commandline error. */
void PrintCmdErr()
{
    cout << "Commandline parameters are not specified proberly!\n";
    cout << "See the manual for further details\n";
    exit(1);
}

// -----------------------------------------------------------------------------

/** Main function for the program. */
int main(int argc,char **argv)
{
	int i=1;
    bool verbose=false;
    bool f5=false;
    bool xml=false;
    string filename="";    // file of output

	while (i<argc)
	{
		if (!strcmp(argv[i],"-verbose")) verbose = true;
		if (!strcmp(argv[i],"-f5")) f5=true;
		if (!strcmp(argv[i],"-xml")) xml=true;
		if (!strcmp(argv[i],"-out")) filename=argv[++i];
		i++;
	}

    if (filename=="") PrintCmdErr();

	TegPeak TEGP(filename,verbose,f5,xml);
    TEGP.Run();
	return 0;
}




