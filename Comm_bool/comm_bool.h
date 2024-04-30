#ifndef COMM_BOOL_H
#define COMM_BOOL_H

#include "mpi.h"
using namespace std;

uint64_t send_bytes = 0;
uint64_t recv_bytes = 0;
void init(int argc, char** argv);

void close();


template<typename T>
void boolSendToParty(T* data, int len, int PartyID, int shareTag);

template<typename T>
void boolRecvFromParty(T* data, int len, int PartyID, int shareTag);

uint64_t getsendbytes(){
    return send_bytes;
}

uint64_t getrecvbytes(){
    return recv_bytes;
}


#endif
