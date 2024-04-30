#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include "comm_bool.h"




#define PRIVATE static

#define NUM_PARTIES 3
#define XCHANGE_MSG_TAG 7
#define OPEN_MSG_TAG 777



PRIVATE void check_init(const char* f);

// initialize MPI, rank, num_parties
/*void init(int argc, char** argv) {
  MPI_Init(&argc, &argv);
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);
  MPI_Comm_size(MPI_COMM_WORLD, &num_parties);

  // this protocol works with 3 parties only
  if (rank == 0 && num_parties != NUM_PARTIES) {
    fprintf(stderr, "ERROR: The number of MPI processes must be %d for %s\n", NUM_PARTIES, argv[0]);
    MPI_Abort(MPI_COMM_WORLD, 1);
  }
  initialized = 1;
}*/

template<typename T>
void boolSendToParty(T* data, int len, int PartyID, int shareTag){
    if((std::is_same<T, uint64_t>::value)){
      MPI_Send(data, len, MPI_UNSIGNED_LONG_LONG, PartyID, shareTag, MPI_COMM_WORLD);
      send_bytes += len*sizeof(uint64_t);
    }else if((std::is_same<T, u_int32_t>::value)){
      MPI_Send(data, len, MPI_UINT32_T, PartyID, shareTag, MPI_COMM_WORLD);
      send_bytes += len*sizeof(uint32_t);
    }
      
    
}

template<typename T>
void boolRecvFromParty(T* data, int len, int PartyID, int shareTag){
    if((std::is_same<T, uint64_t>::value)){
      MPI_Recv(data, len, MPI_UNSIGNED_LONG_LONG, PartyID, shareTag, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
      recv_bytes += len*sizeof(uint64_t);
    }else if((std::is_same<T, u_int32_t>::value)){
      MPI_Recv(data, len, MPI_UINT32_T, PartyID, shareTag, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
      recv_bytes += len*sizeof(uint32_t);
    }
      
    
}
