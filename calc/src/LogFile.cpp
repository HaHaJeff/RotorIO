/*---------------------------------------------------------------------------*\


\*---------------------------------------------------------------------------*/


#include "LogFile.h"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

void LogInfo::Log(std::string loginfo) 
{  
    std::ofstream ofs;  
    time_t t = time(0);  
    char tmp[64];  
    strftime(tmp, sizeof(tmp), "\t[%Y.%m.%d %X %A]", localtime(&t));  
    ofs.open("LogInfo.log", std::ofstream::app);  
    ofs.write(loginfo.c_str(), loginfo.size());  
    ofs << tmp << '\n';  
    ofs.close();  
}  
// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //



