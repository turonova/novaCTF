//============================================================================
// Name        : NovaCTF
// Author      : Beata Turonova
// Version     : 1.0
//============================================================================

#include <iostream>
#include <sstream>
#include <string>
#include <vector>
#include <ctime>
#include "parameterSetup.h"
#include "mrcStack.h"
#include "exception.h"
#include "volumeIO.h"

#include "fftRoutines.h"

#include "defocus.h"
#include "filterProjections.h"
#include "ctf3d.h"
#include "ctfCorrection.h"

using namespace std;

string dateAndTimeStamp(time_t startTime)
{
    time_t endTime;
    //struct tm * timeinfo;
    //char buffer[80];

    time(&endTime);

    double seconds_total = difftime(endTime,startTime);
    unsigned int hours = ((unsigned int)seconds_total)/3600;
    unsigned int minutes = (((int)seconds_total)%3600)/60;
    unsigned int seconds = (((int)seconds_total)%3600)%60;
    //timeinfo = localtime(&rawtime);
    stringstream timeDifference;


    if(hours>0)
        timeDifference << hours << "h ";

    if(minutes > 0 || hours > 0)
        timeDifference << minutes << "m ";

    timeDifference << seconds << "s";

    //strftime(buffer,80,"%d-%m-%Y %I:%M:%S",timeinfo);
    //std::string stamp(buffer);

    return timeDifference.str();
}

int main (int argc, const char* argv[]) {

    time_t startTime;
    time(&startTime);

    vector<std::string> argList;

    cout << "Reading parameters... " <<  endl;

    for(int i=1;i<argc;i++)
    {
        argList.push_back(argv[i]);
    }

    ParameterSetup params(argList);

    cout << "...done" <<  endl << endl;

    if(params.Algorithm()=="3dctf")
    {
        cout << "Starting WBP reconstruction with 3D-CTF correction." << endl << endl;
        CTF3d* newRec = new CTF3d(params);
        newRec->run();
        cout << "Reconstruction finished successfully! \n Required time: " << dateAndTimeStamp(startTime) << endl;
        delete newRec;
    }
    else if(params.Algorithm()=="defocus")
    {
        cout << "Starting computation of defocus." << endl << endl;
        Defocus* newRec = new Defocus(params);
        newRec->run();
        cout << "Finished successfully! \n Required time: " << dateAndTimeStamp(startTime) << endl;
        delete newRec;
    }
    else if(params.Algorithm()=="filterProjections")
    {
        cout << "Starting filtering projections." << endl << endl;
        FilterProjections* newRec = new FilterProjections(params);
        newRec->run();
        cout << "Filtering of projections finished successfully! \n Required time: " << dateAndTimeStamp(startTime) << endl;
        delete newRec;
    }
    else if(params.Algorithm()=="ctfCorrection")
    {
        cout << "Starting CTF correction." << endl << endl;
        CTFCorrection* newRec = new CTFCorrection(params);
        newRec->run();
        cout << "CTF correction finished successfully! \n Required time: " << dateAndTimeStamp(startTime) << endl;
        delete newRec;
    }
    else
    {
        cout << "An algorithm was not specified! Use -Algorithm defocus/ctfCorrection/filterProjections/3dctf option to specify the algorithm you want to use." << endl << endl;
    }

    return 0;
}
