#ifndef BoolSH2_H
#define BoolSH2_H
#include "../ABY/src/abycore/circuit/booleancircuits.h"

#include "../ABY/src/abycore/circuit/arithmeticcircuits.h"
#include "../ABY/src/abycore/circuit/circuit.h"
#include "../ABY/src/abycore/aby/abyparty.h"

#include<vector>
#include "../Comm_bool/comm_bool.cpp"
#include <sodium.h>

typedef uint64_t Data;
//#include "../Socket/Socket.h"

//#define p 10000


using namespace std;

class boolSh2 {

public:
    Data* s2x;

    
    int len;

    boolSh2(){};
    boolSh2(int length);
    
    ~boolSh2(){};

    void init(int length);
    template<typename T>
    void share(int party, int dest_role, int data_owner, T* data);
    void restruct(int party, Data* r_data, int dest_role);

    void add();
    void multiply();
};

#endif

