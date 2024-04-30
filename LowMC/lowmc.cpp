#include "lowmc.h"
#include "../ABY/src/abycore/sharing/sharing.h"
#include <ENCRYPTO_utils/crypto/crypto.h>
//#include "../BoolSh2/boolSh2.cpp"




static uint32_t m_nRndCtr;
static uint32_t* m_tGrayCode;
static uint32_t* m_tGrayCodeIncrement;
static uint32_t m_nZeroGate;

//x

share* lowmc_exe(int role, uint32_t bitlen, int nvals, boolSh2& s, CBitVector& key , LowMCParams* param, crypto* crypt, ABYParty* party){
    uint32_t ctr = 0, exp_key_bitlen = param->blocksize * (param->nrounds+1), zero_gate;
    
	int bytenum = param->blocksize/8;
	

	uint8_t* output;

    share *s_in, *s_key, *s_ciphertext;
//输入转换成电路输入形式
//cout<<"======================circuit========================"<<endl;
	vector<Sharing*>& sharings = party->GetSharings();

    //e_sharing sharing = S_BOOL;
    Circuit* circ = sharings[S_BOOL]->GetCircuitBuildRoutine();

 	//Circuit* ac = sharings[S_ARITH]->GetCircuitBuildRoutine();
	//Circuit* yc = sharings[S_YAO]->GetCircuitBuildRoutine();

	//share* s_ina;
	s_in = circ->PutSharedSIMDINGate(nvals, s.s2x, param->blocksize);
	        //s_ina = ac->PutB2AGate(s_in);
        //s_in = circ->PutA2BGate(s_ina, yc);
//circ->PutPrintValueGate(s_ina, "input");
	//s_in = circ->PutA2BGate(s_ina, yc);

	//s_in = circ->PutSharedSIMDINGate(nvals, in, param->blocksize);
//circ->PutPrintValueGate(s_in, "input_bool");

	s_key = circ->PutINGate(key.GetArr(), exp_key_bitlen, SERVER);
	s_key = circ->PutRepeaterGate(nvals, s_key);

	zero_gate = circ->PutConstantGate(0, nvals);
//circ->PutPrintValueGate(s_key, "key");
//cout<<"build bf"<<endl;
    s_ciphertext = Buildlowmccir(role, s_in, s_key, (BooleanCircuit*) circ, param, zero_gate, crypt);

//circ->PutPrintValueGate(s_ciphertext, "output");
	s_ciphertext = circ->PutOUTGate(s_ciphertext, ALL);//怎么给P0
//	 s_ciphertext = circ->PutOUTGate(s_in, ALL);
//cout<<"over"<<endl;
	
	return s_ciphertext;

	
}

share* Buildlowmccir(int role, share* val, share* key, BooleanCircuit* bc, LowMCParams* param,  uint32_t zerogate, crypto* crypt){


    uint32_t round, byte, i, j, k;
	m_nRndCtr = 0;
	uint32_t nsboxes = param->nsboxes;
	uint32_t statesize = param->blocksize;
	uint32_t nrounds = param->nrounds;

	std::vector<uint32_t> state(statesize);
	
	
	//m_vRandomBits.Create(2 * statesize * statesize * nrounds + nrounds * statesize, crypt);

	//cout<<"m_vRandom.size()"<<m_vRandomBits.GetSize()<<endl;
	

	m_nZeroGate = zerogate;

	//Build the GrayCode for the optimal window-size
	m_tGrayCode = BuildGrayCode(statesize);
	m_tGrayCodeIncrement = BuildGrayCodeIncrement(statesize);

	//copy the input to the current state
	for (i = 0; i < statesize; i++)
		state[i] = val->get_wire_id(i);

    LowMCAddRoundKey(state, key->get_wires(), statesize, 0, bc); //ARK

	//boolshare* bs = new boolshare(state, bc);
//bc->PutPrintValueGate(bs, "bs");
	for (round = 0; round < nrounds; round++) {

		//substitution via 3-bit SBoxes
		LowMCPutSBoxLayer(state, nsboxes, bc);

//boolshare* rs = new boolshare(state, bc);
//bc->PutPrintValueGate(rs, "rs");

		//multiply state with GF2Matrix
		LowMCMultiplyState(state, statesize, bc);//Naive version of the state multiplication
		//FourRussiansMatrixMult(state, statesize, bc);//4 Russians version of the state multiplication
		//LowMCMultiplyStateCallback(state, statesize, circ); //use callbacks to perform the multiplication in plaintext

		//XOR constants
		LowMCXORConstants(state, statesize, bc);

		//XOR with multiplied key
		LowMCXORMultipliedKey(state, key->get_wires(), statesize, round, bc);



	}



	free(m_tGrayCode);
	free(m_tGrayCodeIncrement);

    return new boolshare(state, bc);

    
}

//Multiply the key with a 192x192 matrix and XOR the result on the state.
void LowMCXORMultipliedKey(std::vector<uint32_t>& state, std::vector<uint32_t> key, uint32_t lowmcstatesize, uint32_t round, BooleanCircuit* circ) {
	uint32_t tmp;
	/*for(uint32_t i = 0; i < MPCC_STATE_SIZE; i++) {
	 tmp = 0;
	 for(uint32_t j = 0; j < MPCC_STATE_SIZE; j++, m_nRndCtr++) {
	 if(m_vRandomBits.GetBit(m_nRndCtr)) {
	 tmp = PutXORGate(tmp, key[j]);
	 }
	 }
	 state[i] = PutXORGate(state[i], tmp);
	 }*/
	//Assume outsourced key-schedule
	for (uint32_t i = 0; i < lowmcstatesize; i++) {
		state[i] = circ->PutXORGate(state[i], key[i+(1+round) * lowmcstatesize]);
	}

}

//XOR constants on the state
void LowMCXORConstants(std::vector<uint32_t>& state, uint32_t lowmcstatesize, BooleanCircuit* circ) {
	for (uint32_t i = 0; i < lowmcstatesize; i++, m_nRndCtr++) {
		if (m_vRandomBits.GetBit(m_nRndCtr)) {
			state[i] = circ->PutINVGate(state[i]);
		}
	}

}


void FourRussiansMatrixMult(std::vector<uint32_t>& state, uint32_t lowmcstatesize, BooleanCircuit* circ) {
	//round to nearest square for optimal window size
	uint32_t wsize = floor_log2(lowmcstatesize) - 2;

	//will only work if the statesize is a multiple of the window size
	uint32_t* lutptr;
	uint32_t* lut = (uint32_t*) malloc(sizeof(uint32_t) * (1 << wsize));
	uint32_t i, j, bitctr, tmp = 0;

	lut[0] = m_nZeroGate;	//circ->PutConstantGate(0, 1);

	std::vector<uint32_t> tmpstate(ceil_divide(lowmcstatesize, wsize) * wsize, lut[0]);
	//pad the state to a multiple of the window size and fill with zeros
	std::vector<uint32_t> state_pad(ceil_divide(lowmcstatesize, wsize) * wsize, lut[0]);
	for (i = 0; i < lowmcstatesize; i++)
		state_pad[i] = state[i];

	for (i = 0, bitctr = 0; i < ceil_divide(lowmcstatesize, wsize); i++) { //for each column-window
		for (j = 1; j < (1 << wsize); j++) {
			lut[m_tGrayCode[j]] = circ->PutXORGate(lut[m_tGrayCode[j - 1]], state_pad[i * wsize + m_tGrayCodeIncrement[j - 1]]);
		}

		for (j = 0; j < lowmcstatesize; j++, bitctr += wsize) {
			m_vRandomBits.GetBits((BYTE*) &tmp, bitctr, wsize);
			tmpstate[i] = circ->PutXORGate(tmpstate[j], lut[tmp]);
		}
	}

	for (i = 0; i < lowmcstatesize; i++)
		state[i] = tmpstate[i];

	free(lut);
}


void LowMCAddRoundKey(std::vector<uint32_t>& val, std::vector<uint32_t> key, uint32_t lowmcstatesize, uint32_t round, BooleanCircuit* circ) {
	for (uint32_t i = 0; i < lowmcstatesize; i++) {
		val[i] = circ->PutXORGate(val[i], key[i+(1+round) * lowmcstatesize]);
	}
}
//Put a layer of 3-bit LowMC SBoxes
void LowMCPutSBoxLayer(std::vector<uint32_t>& input, uint32_t nsboxes, BooleanCircuit* circ) {
	for (uint32_t i = 0; i < nsboxes * 3; i += 3) {
		LowMCPutSBox(input[i], input[i + 1], input[i + 2], circ);
	}
}
//Put a 3-bit LowMC SBoxes
void LowMCPutSBox(uint32_t& o1, uint32_t& o2, uint32_t& o3, BooleanCircuit* circ) {
	uint32_t i1 = o1;
	uint32_t i2 = o2;
	uint32_t i3 = o3;

	uint32_t ni1 = circ->PutINVGate(i1);
	uint32_t ni2 = circ->PutINVGate(i2);
	uint32_t ni3 = circ->PutINVGate(i3);

	//C = B * C + A
	o1 = circ->PutXORGate(circ->PutANDGate(i2, i3), i1);

	//E = A * (NOT C) + B
	o2 = circ->PutXORGate(circ->PutANDGate(i1, ni3), i2);

	//F = (NOT ((NOT B) * (NOT A))) + C
	o3 = circ->PutXORGate(circ->PutINVGate(circ->PutANDGate(ni2, ni1)), i3);
}


uint32_t* BuildGrayCode(uint32_t length) {
	uint32_t* gray_code = (uint32_t*) malloc(sizeof(uint32_t) * length);
	for(uint32_t i = 0; i < length; ++i) {
		gray_code[i] = i ^ (i >> 1);
	}
	return gray_code;
}

uint32_t* BuildGrayCodeIncrement(uint32_t length) {
	uint32_t* gray_code_increment = (uint32_t*) malloc(sizeof(uint32_t) * length);
	for(uint32_t i = 0; i < length; ++i) {
		gray_code_increment[i] = 0;
	}
	uint32_t length_inc = 2;
	while(length_inc < length) {
		uint32_t length_count = length_inc - 1;
		while(length_count <= length) {
			(gray_code_increment[length_count])++;
			length_count += length_inc;
		}
		length_inc <<= 1; 
	}
	return gray_code_increment;
}


void LowMCMultiplyState(std::vector<uint32_t>& state, uint32_t lowmcstatesize, BooleanCircuit* circ) {
	std::vector<uint32_t> tmpstate(lowmcstatesize);
	for (uint32_t i = 0; i < lowmcstatesize; i++) {
		tmpstate[i] = 0;
		for (uint32_t j = 0; j < lowmcstatesize; j++, m_nRndCtr++) {
			if (m_vRandomBits.GetBit(m_nRndCtr)) {
				tmpstate[i] = circ->PutXORGate(tmpstate[i], state[j]);
			}
		}
	}
}
