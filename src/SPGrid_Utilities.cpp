//#####################################################################
// Utility classes/functions
//#####################################################################
#include <iostream>
#include <stdexcept>
#include <stdio.h>
#include <string.h>
#include <stdint.h>
#include <errno.h>

#include <SPGrid_Utilities.h>

namespace SPGrid{

//#####################################################################
// Functions Bit_Spread/Bit_Pack
//#####################################################################

int Bit_Pack(const uint64_t data, const uint64_t mask)
{
    union{ int64_t slresult; uint64_t ulresult; };    

    uint64_t uldata=data; int count=0; ulresult=0;
    for(uint64_t bit=1;bit;bit<<=1)
        if(bit & mask) ulresult |= (uldata & bit)>>count; else count++;
    return (int)slresult;
}

//#####################################################################
// Function Fatal_Error
//#####################################################################
void Fatal_Error(const char* function,const char* file,unsigned int line)
{
    Fatal_Error(function,file,line,"Fatal error");
}
void Fatal_Error(const char* function,const char* file,unsigned int line,const char* message)
{
    Fatal_Error(function,file,line,std::string(message));
}
void Fatal_Error(const char* function,const char* file,unsigned int line,const std::string& message)
{
    static char buffer[2048];
    sprintf(buffer,"%s:%s:%d: %s",file,function,line,message.c_str());
    std::string error(buffer);
    std::cout<<std::flush;std::cerr<<"\n";
    std::cerr<<"\n*** ERROR: "<<error<<'\n'<<std::endl;
    throw std::runtime_error(error);
}

//#####################################################################
// Function Check_Compliance
//#####################################################################
void Check_Compliance()
{
    //if(sysconf(_SC_PAGESIZE)!=4096) FATAL_ERROR("Page size different than 4KB detected");
    if(sizeof(unsigned long)!=8) FATAL_ERROR("unsigned long is not 64-bit integer");
    if(sizeof(size_t)!=8) FATAL_ERROR("size_t is not 64-bit long");
    if(sizeof(void*)!=8) FATAL_ERROR("void* is not 64-bit long");
    typedef enum {dummy=0xffffffffffffffffUL} Long_Enum;
    if(sizeof(Long_Enum)!=8) FATAL_ERROR("Missing support for 64-bit enums");
}

//#####################################################################
}

