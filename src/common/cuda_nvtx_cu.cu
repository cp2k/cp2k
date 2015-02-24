/*****************************************************************************
 *  CP2K: A general program to perform molecular dynamics simulations        *
 *  Copyright (C) 2000 - 2015  CP2K developers group                         *
 *****************************************************************************/

#if defined( __CUDA_PROFILING)

#include <nvToolsExt.h>
#include <stdio.h>
#include <pthread.h>

const uint32_t colormap[] = { 0xFFFFFF00,  // Yellow
                              0xFFFF00FF,  // Fuchsia
                              0xFFFF0000,  // Red
                              0xFFC0C0C0,  // Silver
                              0xFF808080,  // Gray
                              0xFF808000,  // Olive
                              0xFF800080,  // Purple
                              0xFF800000,  // Maroon
                              0xFF00FFFF,  // Aqua
                              0xFF00FF00,  // Lime
                              0xFF008080,  // Teal
                              0xFF008000,  // Green
                              0xFF0000FF,  // Blue
                              0xFF000080}; // Navy

//==============================================================================
extern "C" int cuda_nvtx_range_push_cu(const char* message) {

    //assembling event attribute
    nvtxEventAttributes_t eventAttrib = {0};
    eventAttrib.version = NVTX_VERSION;
    eventAttrib.size = NVTX_EVENT_ATTRIB_STRUCT_SIZE;
    eventAttrib.messageType = NVTX_MESSAGE_TYPE_ASCII;
    eventAttrib.message.ascii = message;

    // colors are picked based on a (very simple) hash value of the message
    int hash=0;
    for (int i=0; i < strlen(message); i++)
        hash += i*message[i]*message[i];
    eventAttrib.colorType = NVTX_COLOR_ARGB;
    eventAttrib.color = colormap[hash%14];

    //these field could be fild with useful stuff
    eventAttrib.payloadType = NVTX_PAYLOAD_TYPE_INT64;
    eventAttrib.payload.llValue = 123;
    eventAttrib.category = 42;

    int level = nvtxRangePushEx(&eventAttrib);
    return(level);
}

//==============================================================================
extern "C" int cuda_nvtx_range_pop_cu() {
    int level = nvtxRangePop();
    return(level);
}

//==============================================================================
extern "C" void cuda_nvtx_name_osthread_cu(char* name){
    nvtxNameOsThread(pthread_self(), name);
}

#endif
