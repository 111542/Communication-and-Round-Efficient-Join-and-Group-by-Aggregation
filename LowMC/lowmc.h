#ifndef LOWMC_H
#define LOWMC_H

#include "../ABY/src/abycore/circuit/booleancircuits.h"
#include "../ABY/src/abycore/aby/abyparty.h"
#include <ENCRYPTO_utils/cbitvector.h>
#include <ENCRYPTO_utils/typedefs.h>
#include <cassert>
#include <vector>
#include "../ABY/src/abycore/sharing/sharing.h"
#include <ENCRYPTO_utils/crypto/crypto.h>
#include "../BoolSh2/boolSh2.cpp"



using namespace std;
static const BYTE mpccseed[16] = { 0x00, 0x11, 0x22, 0x33, 0x44, 0x55, 0x66, 0x77, 0x88, 0x99, 0xAA, 0xBB, 0xCC, 0xDD, 0xEE, 0xFF };

/* versions of the lowmc cipher: short term security, short term wide, long term security, long term wide */
enum LowMCVersion {
	STLowMC = 0, LTLowMC = 1
};

uint32_t blocksize = 64;
    uint32_t sboxes = 13, keysize = 128, rounds = 12;

struct LowMCParams {
	uint32_t nsboxes;
	uint32_t keysize;
	uint32_t blocksize;
	uint32_t data;
	uint32_t nrounds;
};

struct matmul {
	UGATE_T** matrix;
	uint32_t column;
};

//parameters: sboxes (m), key-length (k), statesize (n), data (d), rounds (r)
static const LowMCParams stp = { 49, 80, 256, 64, 12 };
static const LowMCParams ltp = { 63, 128, 256, 128, 14 };

static const LowMCParams lowmcparamlookup[] = { stp, ltp};
 
static CBitVector m_vRandomBits ;
static CBitVector key;

share* lowmc_exe(int role, uint32_t bitlen, int nvals, boolSh2& s, CBitVector& key ,LowMCParams* param, crypto* crypt, ABYParty* party);
share* Buildlowmccir(int role, share* val, share* key, BooleanCircuit* bc, LowMCParams* param,  uint32_t zerogate, crypto* crypt);
void LowMCXORMultipliedKey(std::vector<uint32_t>& state, std::vector<uint32_t> key, uint32_t lowmcstatesize, uint32_t round, BooleanCircuit* circ);
void LowMCXORConstants(std::vector<uint32_t>& state, uint32_t lowmcstatesize, BooleanCircuit* circ);
void FourRussiansMatrixMult(std::vector<uint32_t>& state, uint32_t lowmcstatesize, BooleanCircuit* circ);
void LowMCAddRoundKey(std::vector<uint32_t>& val, std::vector<uint32_t> key, uint32_t lowmcstatesize, uint32_t round, BooleanCircuit* circ);
void LowMCPutSBoxLayer(std::vector<uint32_t>& input, uint32_t nsboxes, BooleanCircuit* circ);
void LowMCPutSBox(uint32_t& o1, uint32_t& o2, uint32_t& o3, BooleanCircuit* circ);
uint32_t* BuildGrayCode(uint32_t length);
uint32_t* BuildGrayCodeIncrement(uint32_t length);
void LowMCMultiplyState(std::vector<uint32_t>& state, uint32_t lowmcstatesize, BooleanCircuit* circ);

#endif
