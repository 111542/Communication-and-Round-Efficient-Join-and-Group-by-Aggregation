#ifndef TABLE_H
#define TABLE_H

#include "../BoolSh2/boolSh2.h"
#include "../ABY/src/abycore/aby/abyparty.h"
#include "../ABY/src/abycore/sharing/sharing.h"
#include <ENCRYPTO_utils/crypto/crypto.h>
#include "../LowMC/lowmc.cpp"
#include <random>
#include<algorithm>
#include <utility>
#include <bitset>
#include <math.h>
#include <sys/time.h>
#include "../Sh2/Sh2.cpp"
#define p 10000
#define XSHARETAG 7


using namespace std;


//A*B: rowA*colB
void Multi_matrix(bool** A, int rowA, int colA, bool** B, int rowB, int colB, bool** res){
int count = 0;

    for(int i=0; i<rowA; i++){
        for(int j=0; j<colB; j++){
            res[i][j] = 0;
            for(int k1=0; k1<colA; k1++){

                res[i][j] |= (A[i][k1] & B[k1][j]);

            }
            //cout<<res[i][j];
        }
       // cout<<endl;
    }
    
}

bool isPowerOfTwo(int n) {
    return n > 0 && std::floor(std::log2(n)) == std::ceil(std::log2(n));
}


class STable{

public:
    boolSh2 s_K;
    //vector<Sh2> s_S();
    boolSh2 s_S, s_S1;


    int len;

    STable(int length){s_K.init(length);s_S1.init(length);s_S.init(length); len = length;};
    STable(){len = 0;s_K.init(len);s_S.init(len);s_S1.init(len);};
    //~STable(){};
};

void reRandomize(int role, boolSh2 s, int PartyID){
    Data* r1 = new Data[s.len];
    Data* r2 = new Data[s.len];
    

    for(int i=0; i<s.len; i++){
        randombytes_buf(&r1[i], sizeof(Data));
    }
    if(role > PartyID){
        
        
        boolSendToParty(r1, s.len, PartyID, XSHARETAG);
        boolRecvFromParty(r2, s.len, PartyID, XSHARETAG);

    }else{
        boolRecvFromParty(r2, s.len, PartyID, XSHARETAG);
        boolSendToParty(r1, s.len, PartyID, XSHARETAG);
        
    }
        


    for(int i=0; i<s.len; i++){
        s.s2x[i] = s.s2x[i] ^ r1[i] ^ r2[i];
    }

}

void reRandomize(int role, Sh2 s, int PartyID){
    uint64_t* r1 = new uint64_t[s.len];
    uint64_t* r2 = new uint64_t[s.len];
    

    for(int i=0; i<s.len; i++){
        r1[i] = rand()%1000;
    }
    if(role > PartyID){
        
        
        boolSendToParty(r1, s.len, PartyID, XSHARETAG);
        boolRecvFromParty(r2, s.len, PartyID, XSHARETAG);

    }else{
        boolRecvFromParty(r2, s.len, PartyID, XSHARETAG);
        boolSendToParty(r1, s.len, PartyID, XSHARETAG);
        
    }
        


    for(int i=0; i<s.len; i++){
        s.s2x[i] = s.s2x[i] - r1[i] + r2[i];
    }

}

//
void mpc_agg(int role, int operation, uint32_t bitlen, Data* a, Data* b1, Data* a2, Data* b2, Data** Mb, int rows, int j, crypto* crypt){
    
    //uint32_t bitlen = 64;
    string address = "127.0.0.1";
	e_mt_gen_alg mt_alg = MT_OT;
    uint32_t nthreads = 1;
	seclvl seclvl = crypt->get_seclvl();
	ABYParty* party;
    party = new ABYParty((e_role)(role-1), address, 7766, seclvl, 32, nthreads, mt_alg);

    vector<Sharing*>& sharings = party->GetSharings();

    Circuit* bc = sharings[S_BOOL]->GetCircuitBuildRoutine();

    share* s_a, *s_b1, *s_a2, *s_b2, *s_b, *s_tmp, *s_sel;
    s_a = bc->PutSharedSIMDINGate(rows, &a[1], bitlen);
    
    s_b1 = bc->PutSharedSIMDINGate(rows, &b1[1], bitlen);
    
    s_a2 = bc->PutSharedSIMDINGate(rows, &a2[1], bitlen); 
    s_b2 = bc->PutSharedSIMDINGate(rows, &b2[1], bitlen);

//ac->PutPrintValueGate(s_b1, "input_b1");
//ac->PutPrintValueGate(s_b2, "input_b2");

    //s_a = bc->PutMULGate(s_a1, s_a2);

    s_tmp = bc->PutMULGate(s_a2, s_b1);
    
    
    //cout<<"111"<<endl;
//A转B吧，没有A的GTGate
    //s_tmp = bc->PutA2BGate(s_tmp, yc);
    //s_a = ac->PutA2BGate(s_a, yc);
    //s_b2 = bc->PutA2BGate(s_b2, yc);
    

    s_sel = bc->PutGTGate(s_tmp, s_b2);
    
    if(operation == 1){
        s_b = bc->PutMUXGate(s_tmp, s_b2, s_sel);
    }else{
        s_b = bc->PutMUXGate(s_b2, s_tmp, s_sel);
    }
 //bc->PutPrintValueGate(s_b, "mggb");
    //s_b = ac->PutB2AGate(s_b);
    //share* a, *b;
    //a = circ->PutSharedOUTGate(s_a);
    //b = circ->PutSharedOUTGate(s_b);
    Data* r = new Data[rows];
    for(int i=0; i<rows; i++){
        randombytes_buf(&r[i], 0);
    }
    
    share* s_r = bc->PutSharedSIMDINGate(rows, r, bitlen);
    //s_a = bc->PutXORGate(s_a, s_r);
    s_b = bc->PutXORGate(s_b, s_r);
    //s_a = bc->PutOUTGate(s_a, ALL);
    s_b = bc->PutOUTGate(s_b, ALL);

    party->ExecCircuit();

    //Data* res_a = new Data[rows];
    Data *res_b = new Data[rows];
    
    //s_a->get_clear_value_vec(&res_a, &bitlen, &nvals);
    s_b->get_clear_value_vec(&res_b, &bitlen, (uint32_t*)&rows);
    std::cout <<"j: "<<j<<" "<< party->GetTiming(P_SETUP) << "\t" << party->GetTiming(P_ONLINE) << "\t" << party->GetTiming(P_TOTAL) << std::endl;


    delete party;
    boolSh2 r_a(rows), r_b(rows);
    //r_a.share(role, role*2%3, 1, res_a);
    r_b.share(role, role*2%3, 1, res_b);

    

    for(int i=1; i<rows+1; i++){
        //Ma[i][j] = r_a.s2x[i-1] ^ r[i-1];
        Mb[i][j] = r_b.s2x[i-1] ^ r[i-1];        
    }

}











void F_shuffle(int role, int bitlen, vector<uint32_t>& perm, boolSh2& V, boolSh2& R, uint32_t dsize){  //dsize的取值是n-m的m？

//判断n和m的大小，决定tsize的值
    uint32_t tsize = V.len;

    vector<uint32_t> rPerm1(tsize);
    boolSh2 temp(tsize), tmp1(tsize), tmp2(tsize);
    vector<uint32_t> reverse_rPerm1(tsize);
    vector<uint32_t> rPerm2(dsize);
    //P0
    if(role == 0){

        boolRecvFromParty(temp.s2x, temp.len, (role+1)%3, XSHARETAG);

        for(int i=0; i<tsize; i++){
            rPerm1[i] = i;
        }
        random_shuffle(rPerm1.begin(), rPerm1.end());

        boolSendToParty(rPerm1.data(), rPerm1.size(), (role-1+3)%3, XSHARETAG);
        
        for(int i=0; i<rPerm1.size(); i++){
            tmp1.s2x[i] = temp.s2x[rPerm1[i]];

            //cout<<temp.s2x[i]<<" "<<tmp1.s2x[i]<<endl;
        }
        //重随机化
        //reRandomize(tmp1);
        reRandomize(role, tmp1, 2);


        for(int i=0; i<tsize; i++){
            reverse_rPerm1[rPerm1[i]] = i;

        }

        
        for(int i=0; i<dsize; i++){
            rPerm2[i] = reverse_rPerm1[perm[i]];
           
        }

        boolSendToParty(rPerm2.data(), rPerm2.size(), (role+1)%3, XSHARETAG);

        for(int i=0; i<rPerm2.size(); i++){
            tmp2.s2x[i] = tmp1.s2x[rPerm2[i]];
        }

        //reRandomize(tmp2);
        reRandomize(role, tmp2, 1);

        boolSendToParty(tmp2.s2x, tmp2.len, (role+2)%3, XSHARETAG);
       
    }else if(role == 1){
    //P1
        boolSendToParty(V.s2x, V.len, (role+2)%3, XSHARETAG);
        boolRecvFromParty(tmp1.s2x, tmp1.len, (role+1)%3, XSHARETAG);

        boolRecvFromParty(rPerm2.data(), rPerm2.size(), (role+2)%3, XSHARETAG);
        
        for(int i=0; i<rPerm2.size(); i++){
            R.s2x[i] = tmp1.s2x[rPerm2[i]];
            //cout<<tmp1.s2x[i]<<" "<<tmp2.s2x[i]<<endl;
        }
        //reRandomize(R);

        reRandomize(role, R, 0);
    }else{        
    //P2
        boolRecvFromParty(rPerm1.data(), rPerm1.size(), (role+1)%3, XSHARETAG);
        for(int i=0; i<rPerm1.size(); i++){
            tmp1.s2x[i] = V.s2x[rPerm1[i]];
            
        }
        //reRandomize(tmp1);
        reRandomize(role, tmp1, 0);

        boolSendToParty(tmp1.s2x, tmp1.len, (role+2)%3, XSHARETAG);

        boolRecvFromParty(R.s2x, R.len, (role+1)%3, XSHARETAG);

    }

}
//
void F_encode(int role, uint32_t bitlen, CBitVector& key, LowMCParams* param , crypto* crypt, boolSh2& s, uint32_t len, Data* E){
//cout<<"111"<<endl;
    uint32_t blocksize = param->blocksize;
    int nvals = s.len;
    //uint8_t** E = new uint8_t*[nvals];

    string address = "127.0.0.1";
	e_mt_gen_alg mt_alg = MT_OT;
    uint32_t nthreads = 1;
	seclvl seclvl = crypt->get_seclvl();
	
	
	

    e_role er;
	ABYParty* party;
	if(role == 1){
		er = (e_role)0;
	}else if(role == 2){
		er = (e_role)1;
	}
	if(role != 0)
		party = new ABYParty(er, address, 7766, seclvl, 32, nthreads, mt_alg);


//P0
    if(role == 0){


        //uint8_t* Ex = new uint8_t[nvals*blocksize/8];
        //接收E并做E+r
        //comm_1.Receive(&E, sizeof(E));
        //if(recv(comm_1, Ex, nvals*blocksize/8, 0)==-1)  
            //cout<<"E"<<endl;
        
        //MPI_Recv(Ex, nvals*blocksize/8, MPI_UINT8_T, (role+1)%3, XSHARETAG, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
        boolRecvFromParty(E, len, 1, XSHARETAG);
        /*for(int i=0; i<E.size(); i++){
            E[i]+=rand()%p;         
        }*/
        

        //cout<<"=========encode: ========="<<endl;
        /*for(int i=0; i<nvals; i++){
            E[i] = new uint8_t[blocksize/8];
            for(int j=0; j<blocksize/8; j++){
                E[i][j] = Ex[i*blocksize/8+j];
                //cout<<(int)E[i][j];
            }
            //cout<<endl;
        }*/
        //cout<<"=======end========="<<endl;
        //return E;
    }else{
//P0，P1
        //生成随机矩阵I，l*lamda，l是输入的比特长，lamda是lowmc的块大小:P1生成之后发给P2
        //l < lambda的情况可能更多
   
    share* s_ciphertext = lowmc_exe(role, bitlen, s.len, s, key, param, crypt, party); //发送E给P0

  
    party->ExecCircuit();
cout<<"111111111111111111111111111111111111111111111111111111111111111111111111"<<endl;
    //Data* output;
    s_ciphertext->get_clear_value_vec(&E, &bitlen, &len);
    std::cout <<"encode: "<< party->GetTiming(P_SETUP) << "\t" << party->GetTiming(P_ONLINE) << "\t" << party->GetTiming(P_TOTAL) << std::endl;
    delete party;

    /*for(int i=0; i<nvals; i++){
        for(int j=0; j<param->blocksize/8; j++)
            cout<<(int)output[i*param->blocksize/8+j];
        cout<<endl;
    }*/

        if(role == 1){
            //comm_0.Send(&E, sizeof(E));
            //send(comm_0, output, nvals*blocksize/8, 0);
            
            //MPI_Send(output, nvals*blocksize/8, MPI_UINT8_T, (role+2)%3, XSHARETAG, MPI_COMM_WORLD);
            boolSendToParty(E, len, 0, XSHARETAG);
            
        }
        //return E;
    }
}

typedef pair<uint32_t, Data> index_value;

bool compare(index_value v1, index_value v2){
    if(v1.second<v2.second)
        return 1;
    else
        return 0;
}

bool compare_less(uint8_t** Ex, uint8_t** Ey, int inx, int iny){
    int nbytes = blocksize/8;
    for(int i=0; i<nbytes; i++){
        if(Ex[inx][i] >= Ey[iny][i])
            return 0;
    }
    return 1;
}
bool cmp(Data x, Data y){
    return x<y;
}
void join(int role, int bitlen, int operation, crypto* crypt, LowMCParams& param, STable& X, STable& Y, STable& res){
//operation = 0: inner 
//operation = 1: semi-join， 

    
    int nx = X.len;
    int ny = Y.len;


    //uint8_t** Ex;
    //uint8_t** Ey;
    Data* Ex, *Ey;
    //vector<uint32_t[2]> index(nx*ny);
    vector<uint32_t> ind1, ind2;
    boolSh2 Vx(nx), Rx(nx), Vy(ny), Ry(ny);
    struct timeval databegin, dataend;
    long datasec;
    //cout<<"P"<<role<<"   nx:"<<nx<<", ny:"<<ny<<endl;
    uint32_t k = 0;
    //P0
    if(role == 0){
            //随机编码 
        struct timeval p0begin, p0end;
        //cout<<"before: P0 nx "<<nx<<endl;
        //cout<<"P0 recv data"<<endl;
        //gettimeofday(&databegin, 0);
        
        boolRecvFromParty(&nx, 1, 1, XSHARETAG);
        boolRecvFromParty(&ny, 1, 1, XSHARETAG+1);
        
        //gettimeofday(&dataend, 0);

        //datasec = dataend.tv_sec-databegin.tv_sec;

        //cout<<"P0 recv end time: "<<datasec<<endl;
        //cout<<"P0     nx: "<<nx<<endl;
        //cout<<"P0     ny: "<<ny<<endl;
        ind1.resize(nx);
        ind2.resize(ny);

        Ex = new Data[nx];
        Ey = new Data[ny];
        struct timeval exbegin, exend, eybegin, eyend;
        //
        //cout<<"==================encode begin=================="<<endl;
        F_encode(role, bitlen, key, &param, crypt, X.s_K, nx, Ex);   //return Ex
        //
        F_encode(role, bitlen, key, &param, crypt, Y.s_K, ny, Ey);
        //
        //cout<<"==============encode end======================="<<endl;
        
        //排序
        
        for(int i=0; i<nx; i++){
            ind1[i]= i;
        }
        for(int i=0; i<ny; i++){
            ind2[i] = i;
        }
        
        //
        //int flag[Ex.size()] = {0};
        
        Data* E = new Data[nx+ny];
        for(int i=0; i<nx; i++){
            E[i] = Ex[i];
        }
        for(int i=0; i<ny; i++){
            E[i+nx] = Ey[i];
        }
        gettimeofday(&p0begin, 0);

        sort(&E[0], &E[nx+ny-1], cmp);
        gettimeofday(&p0end, 0);
        cout<<"P0 kkkkkk time: "<<p0end.tv_sec-p0begin.tv_sec+(p0end.tv_usec-p0begin.tv_usec)*1e-6<<endl;

        //int k = 0;
        for(int i=0; i<nx+ny-2; i++){
            if(E[i]==E[i+1])
                k++;
        }
        /*for(int i=0; i<nx; i++){
            for(int j=0; j<ny; j++){
                //cout<<i<<"-"<<j<<endl;
                if(Ex[i] == Ey[j]){
                    ind1[k] = i;
                    ind2[k++] = j;
                    //cout<<j<<"-"<<i<<endl;
    //！！！！！！！就是这！！！测试之后去掉            
                    break;
                }
                
            }
        }*/
        
        for(int i=k; i<nx; i++){
            ind1[i] = 0;
        }
        for(int i=k; i<ny; i++){
            ind2[i] = 0;
        }

        cout<<"P0 K: "<<k<<endl;
        res.len = k;
        res.s_K.init(k);
        res.s_S.init(k);
        res.s_S1.init(k);

        
        boolSendToParty(&k, 1, 2, XSHARETAG);
        
        //SendToParty(ind1.data(), ind1.size(), (role+2)%3, XSHARETAG);
        //SendToParty(ind2.data(), ind2.size(), (role+2)%3, XSHARETAG);
        //F_shuffle(role, bitlen, ind1, Vx, Rx, nx);
        //int role, int bitlen, vector<uint32_t>& perm, Sh2& V, Sh2& R, uint32_t dsize, int comm_0, int comm_1
        boolSendToParty(&k, 1, 1, XSHARETAG);
        //SendToParty(ind1.data(), ind1.size(), (role+1)%3, XSHARETAG);
        //SendToParty(ind2.data(), ind2.size(), (role+1)%3, XSHARETAG);
        //F_shuffle(role, bitlen, ind2, Vy, Ry, ny);

        struct timeval shufflebegin, shuffleoneend,shuffleend;
        long shufflesec, shuffleonesec;

       // 
        F_shuffle(role, bitlen, ind1, X.s_K, Vx, nx);
        //
        F_shuffle(role, bitlen, ind1, X.s_S, Rx, nx);
        F_shuffle(role, bitlen, ind2, Y.s_K, Vy, ny);
        F_shuffle(role, bitlen, ind2, Y.s_S, Ry, ny);
        //

        //
        


    }else{
        if(role == 1){
            struct timeval sendbegin, sendend;
            gettimeofday(&sendbegin, 0);
            boolSendToParty(&nx, 1, 0, XSHARETAG);
            boolSendToParty(&ny, 1, 0, XSHARETAG+1);      
            gettimeofday(&sendend, 0);
            cout<<"P1 send : "<<sendend.tv_sec-sendbegin.tv_sec+(sendend.tv_usec-sendbegin.tv_usec)*1e-6<<endl;
            //cout<<"P1 nx: "<<nx<<"            send end"<<endl; 
            //cout<<"P1 send"<<endl;
        }

        
        //P1，P2        
        struct timeval exbegin, exend, eybegin, eyend, shufflebegin, shuffleoneend, shuffleend, kbegin, kend;
        gettimeofday(&exbegin, 0);
        F_encode(role, bitlen, key, &param, crypt, X.s_K, nx, Ex);
        gettimeofday(&exend, 0);
        cout<<"P"<<role<<" ex time: "<<exend.tv_sec-exbegin.tv_sec+(exend.tv_usec-exbegin.tv_usec)*1e-6<<endl;
        F_encode(role, bitlen, key, &param, crypt, Y.s_K, ny, Ey);
        gettimeofday(&eyend, 0);
        cout<<"P"<<role<<" ey time: "<<eyend.tv_sec-exend.tv_sec+(eyend.tv_usec-exend.tv_usec)*1e-6<<endl;
        
        boolRecvFromParty(&k, 1, 0, XSHARETAG);




        cout<<"kkkkkkkkkkkk: "<<k<<endl;
        gettimeofday(&kend, 0);
        cout<<"P"<<role<<" k time: "<<kend.tv_sec-eyend.tv_sec+(kend.tv_usec-eyend.tv_usec)*1e-6<<endl;
        //cout<<"P1/2 k: "<<k<<endl;
        gettimeofday(&shufflebegin, 0);
        F_shuffle(role, bitlen, ind1, X.s_K, Vx, nx);
        gettimeofday(&shuffleoneend, 0);
        F_shuffle(role, bitlen, ind1, X.s_S, Rx, nx);
        F_shuffle(role, bitlen, ind2, Y.s_K, Vy, ny);
        F_shuffle(role, bitlen, ind2, Y.s_S, Ry, ny);
        gettimeofday(&shuffleend, 0);
        cout<<"P"<<role<<" shuffle 1 time: "<<shuffleoneend.tv_sec-shufflebegin.tv_sec+(shuffleoneend.tv_usec-shufflebegin.tv_usec)*1e-6<<endl;
        cout<<"P"<<role<<" shuffle all time: "<<shuffleend.tv_sec-shufflebegin.tv_sec+(shuffleend.tv_usec-shufflebegin.tv_usec)*1e-6<<endl;

        //对K，S执行pi，合并K，S
    //要对res重新初始化  res.init(ind1.size());
        struct timeval resbegin, resend;
        gettimeofday(&resbegin, 0);
        res.len = k;
        res.s_K.init(k);
        res.s_S.init(k);
        res.s_S1.init(k);

        
        for(int i=0; i<k; i++){
            
            res.s_K.s2x[i] = Vx.s2x[i];
            res.s_S.s2x[i] = Rx.s2x[i]; 
            
            res.s_S1.s2x[i] = Ry.s2x[i];   

            //cout<<Vx.s2x[i]<<endl;
            //cout<<i<<": "<<res.s_K.s2x[i]<<" "<<res.s_S.s2x[i]<<" "<<res.s_S1.s2x[i]<<endl;
        }
        gettimeofday(&resend, 0);
        cout<<"P"<<role<<" res time: "<<resend.tv_sec-resbegin.tv_sec+(resend.tv_usec-resbegin.tv_usec)*1e-6<<endl;

    }
    
}

uint32_t group_by(int role, int bitlen,vector<uint32_t>& Ind, LowMCParams* param, crypto* crypt, STable& T, STable& R){   //operation是aggregation的函数

    
    int len = T.len;
    vector<uint32_t> index(len);
    vector<index_value> iv(len);
    //iv.index.resize(len);
    uint32_t m = 1;
    Data* ET = new Data[len];
    F_encode(role, bitlen, key, param, crypt, T.s_K, len, ET);
    //cout<<"P"<<role<<" encode end"<<endl;
    //P0
    if(role == 0){

     

        //int flag[len];
        for(int i=0; i<len; i++){
            iv[i].first = i;
            iv[i].second = ET[i];
        }
        int tag = 0;
        //index[tag++] = 0;
        //flag[0] = 1;
        struct timeval begin, end;
        gettimeofday(&begin, 0);
        /*for(int i=0; i<len; i++){
            if(flag[i] == 1)
                continue;
            else{
                index[tag++] = i;
                flag[i] = 1;
            }
                
            for(int j=i+1; j<len; j++){
                if(flag[j]!=1 && ET[i] == ET[j]){
                    index[tag++] = j;
                    flag[j] = 1;
                }
                    
            }
            
        }*/
        sort(iv.begin(), iv.end(), compare);
        

        gettimeofday(&end, 0);
        //sort(&ET[0], &ET[len-1], cmp);
        cout<<"time: "<<end.tv_sec-begin.tv_sec+(end.tv_usec-begin.tv_usec)*1e-6<<endl;

        for(int i=0; i<index.size(); i++){
            index[i] = iv[i].first;
        }
        
        boolSh2 R(len);
        //cout<<"K shuffle"<<endl;
        
        F_shuffle(role, bitlen, index, T.s_K, R, index.size());
        //cout<<"S shuffle"<<endl;
        F_shuffle(role, bitlen, index, T.s_S, R, index.size());



        //m
        Ind[0] = 1;
        for(int i=1; i<len; i++){
            //if(!compare(ET, ET, index[i], index[i-1])){  
            if(!(ET[index[i]]==ET[index[i-1]])){              //!!!
                Ind[i] = 0;
                m++;
            }else{
                Ind[i] = 1;
            }
        }

        //send(comm_0, &m, sizeof(int), 0);
        //send(comm_1, &m, sizeof(int), 0);

        //cout<<"P0 m: "<<m<<endl;
        boolSendToParty(&m, 1, 2, XSHARETAG);
        boolSendToParty(&m, 1, 1, XSHARETAG);
        //cout<<m<<endl;
        return m;
    }else{
        //P1,P2
        //uint32_t m;
        //cout<<"K shuffle"<<endl;
        //F_encode(role, bitlen, key, param, crypt, T.s_K);
        F_shuffle(role, bitlen, index, T.s_K, R.s_K, len);
        //cout<<"S shuffle"<<endl;
        F_shuffle(role, bitlen, index, T.s_S, R.s_S, len);

        boolRecvFromParty(&m, 1, 0, XSHARETAG);

        //cout<<"P"<<role<<" m: "<<m<<endl;

        return m;


    }
}

void F_shuffle(int role, int bitlen, vector<uint32_t>& perm, Sh2& V, Sh2& R, uint32_t dsize){  //dsize的取值是n-m的m？

//判断n和m的大小，决定tsize的值
    uint32_t tsize = V.len;

    vector<uint32_t> rPerm1(tsize);
    Sh2 temp(tsize), tmp1(tsize), tmp2(tsize);
    vector<uint32_t> reverse_rPerm1(tsize);
    vector<uint32_t> rPerm2(dsize);
    //P0
    if(role == 0){

        boolRecvFromParty(temp.s2x, temp.len, (role+1)%3, XSHARETAG);

        for(int i=0; i<tsize; i++){
            rPerm1[i] = i;
        }
        random_shuffle(rPerm1.begin(), rPerm1.end());

        boolSendToParty(rPerm1.data(), rPerm1.size(), (role-1+3)%3, XSHARETAG);
        
        for(int i=0; i<rPerm1.size(); i++){
            tmp1.s2x[i] = temp.s2x[rPerm1[i]];

            //cout<<temp.s2x[i]<<" "<<tmp1.s2x[i]<<endl;
        }
        //重随机化
        //reRandomize(tmp1);
        reRandomize(role, tmp1, 2);


        for(int i=0; i<tsize; i++){
            reverse_rPerm1[rPerm1[i]] = i;

        }

        
        for(int i=0; i<dsize; i++){
            rPerm2[i] = reverse_rPerm1[perm[i]];
           
        }

        boolSendToParty(rPerm2.data(), rPerm2.size(), (role+1)%3, XSHARETAG);

        for(int i=0; i<rPerm2.size(); i++){
            tmp2.s2x[i] = tmp1.s2x[rPerm2[i]];
        }

        //reRandomize(tmp2);
        reRandomize(role, tmp2, 1);

        boolSendToParty(tmp2.s2x, tmp2.len, (role+2)%3, XSHARETAG);
       
    }else if(role == 1){
    //P1
        boolSendToParty(V.s2x, V.len, (role+2)%3, XSHARETAG);
        boolRecvFromParty(tmp1.s2x, tmp1.len, (role+1)%3, XSHARETAG);

        boolRecvFromParty(rPerm2.data(), rPerm2.size(), (role+2)%3, XSHARETAG);
        
        for(int i=0; i<rPerm2.size(); i++){
            R.s2x[i] = tmp1.s2x[rPerm2[i]];
            //cout<<tmp1.s2x[i]<<" "<<tmp2.s2x[i]<<endl;
        }
        //reRandomize(R);

        reRandomize(role, R, 0);
    }else{        
    //P2
        boolRecvFromParty(rPerm1.data(), rPerm1.size(), (role+1)%3, XSHARETAG);
        for(int i=0; i<rPerm1.size(); i++){
            tmp1.s2x[i] = V.s2x[rPerm1[i]];
            
        }
        //reRandomize(tmp1);
        reRandomize(role, tmp1, 0);

        boolSendToParty(tmp1.s2x, tmp1.len, (role+2)%3, XSHARETAG);

        boolRecvFromParty(R.s2x, R.len, (role+1)%3, XSHARETAG);

    }

}


// sum
void aggre_sum(int role, uint32_t bitlen, vector<uint32_t>& Ind, STable& T, STable& R, int m, crypto* crypt){
    bitlen = 64;
    string address = "127.0.0.1";
	e_mt_gen_alg mt_alg = MT_OT;
    uint32_t nthreads = 1;
	seclvl seclvl = crypt->get_seclvl();
    vector<uint32_t> Idx(Ind.size());
    Sh2 u(T.len);
    Sh2 s_res(T.len);
    boolSh2 k_res(T.len);

    //P1，P2
    if(role != 0){
        ABYParty* party = new ABYParty((e_role)(role-1), address, 7766, seclvl, 32, nthreads, mt_alg);
        vector<Sharing*>& sharings = party->GetSharings();
        Circuit* circ = sharings[S_BOOL]->GetCircuitBuildRoutine();
        Circuit* ac = sharings[S_ARITH]->GetCircuitBuildRoutine();
        Circuit* yc = sharings[S_YAO]->GetCircuitBuildRoutine();

        Data* r = new Data[T.len];
        for(int i=0; i<T.len; i++){
            r[i] = rand()%1000;
        }
        share* val, *rl;

        
        val = circ->PutSharedSIMDINGate(T.len, T.s_S.s2x, bitlen);
        rl = ac->PutSharedSIMDINGate(T.len, r, bitlen);
        val = ac->PutB2AGate(val);
        val = ac->PutADDGate(val, rl);
        //val = circ->PutA2BGate(val, yc);
        share* s_out = ac->PutOUTGate(val, ALL);
        cout<<"111"<<endl;
        party->ExecCircuit();
        std::cout << party->GetTiming(P_SETUP) << "\t" << party->GetTiming(P_ONLINE) << "\t" << party->GetTiming(P_TOTAL) << std::endl;

        Data* sval = new Data[T.len];
        s_out->get_clear_value_vec(&sval, &bitlen, (uint32_t*)&T.len);
        cout<<"s_out"<<endl;
        //Sh2 tmp(T.len);
        u.share(role, 2*role%3, 1, sval);


        for(int i=0; i<u.len; i++){
            u.s2x[i] -= r[i];
        }
       
        delete party;


        //u.s2x[0] = T.s_S.s2x[0];
        //for(int i=1; i<T.len; i++){
        //    u.s2x[i] = u.s2x[i-1]+T.s_S.s2x[i];
        //}
        cout<<"B2A end"<<endl;
        F_shuffle(role, bitlen, Idx, T.s_K, k_res, T.len);
        cout<<"K shuffle"<<endl;
        F_shuffle(role, bitlen, Idx, u, s_res, u.len);
        cout<<"u shuffle"<<endl;

        //Sh2 O(m);

        

        R.s_S.s2x[0] = s_res.s2x[0];
        R.s_K.s2x[0] = k_res.s2x[0];
        for(int i=1; i<m; i++){
            R.s_S.s2x[i] = s_res.s2x[i] - s_res.s2x[i-1];
            R.s_K.s2x[i] = k_res.s2x[i];
        }

        party = new ABYParty((e_role)(role-1), address, 7766, seclvl, 32, nthreads, mt_alg);
        sharings = party->GetSharings();
        circ = sharings[S_BOOL]->GetCircuitBuildRoutine();
        ac = sharings[S_ARITH]->GetCircuitBuildRoutine();
        yc = sharings[S_YAO]->GetCircuitBuildRoutine();

        Data* r1 = new Data[R.len];
        for(int i=0; i<R.len; i++){
            randombytes_buf(&r1[i], 0);
        }
       
        
        val = ac->PutSharedSIMDINGate(R.len, R.s_S.s2x, bitlen);
        rl = circ->PutSharedSIMDINGate(R.len, r1, bitlen);
        val = circ->PutA2BGate(val, yc);
        val = circ->PutXORGate(val, rl);
        //val = circ->PutA2BGate(val, yc);
        share* s_out1 = circ->PutOUTGate(val, ALL);
        party->ExecCircuit();
        std::cout << party->GetTiming(P_SETUP) << "\t" << party->GetTiming(P_ONLINE) << "\t" << party->GetTiming(P_TOTAL) << std::endl;

        Data* sval1 = new Data[R.len];
        s_out1->get_clear_value_vec(&sval1, &bitlen, (uint32_t*)&R.len);

        boolSh2 tmp(R.len);
        tmp.share(role, 2*role%3, 1, sval1);
        for(int i=0; i<R.len; i++){
            R.s_S.s2x[i] = tmp.s2x[i]^r1[i];
        }
        delete party;
        cout<<"A2B"<<endl;
    
    }else{
        //P0
        
        int flag = 0;
        
        for(int i=1; i<Ind.size(); i++){
            if(Ind[i] == 0)
                Idx[flag++] = i-1;
            else{
                int index = rand()%(Ind.size()-m)+m;
                while(Idx[index] != 0)
                    index = rand()%(Ind.size()-m)+m;

                Idx[index] = i-1;
                
                }
        }
        Idx[flag] = Ind.size()-1;          //!!!
        F_shuffle(role, bitlen, Idx, T.s_K, k_res, T.len);
        F_shuffle(role, bitlen, Idx, u, s_res, u.len);

    }
}







void aggre_m(int role, int operation, vector<uint32_t>& Ind, boolSh2& Rs, boolSh2& O, int m, crypto* crypt){
//operation = 1: max
    
    
    
    O.init(m);
    int n = Rs.len;
    int rows = n/2, cols = log2(n);
    boolSh2 sInd(cols);
    //Data* a1, *b1, *a2, *b2;
    Data* a, *b1, *a2, *b2;
    a = new Data[rows+1];
    b1 = new Data[rows+1];
    a2 = new Data[rows+1];
    b2 = new Data[rows+1];

    //vector<vector<share*>> Ma(rows+1, vector<share*>(cols+1)), Mb(rows+1, vector<share*>(cols+1));
    Data** Ma, **Mb;
    Ma = new Data*[rows+1];
    Mb = new Data*[rows+1];
    for(int i=0; i<rows+1; i++){
        Ma[i] = new Data[cols+1];
        Mb[i] = new Data[cols+1];
    }
    /*share*** Ma, *** Mb;
    Ma = new share**[rows+1];
    Mb = new share**[rows+1];
    for(int i=0; i<rows+1; i++){
        Ma[i] = new share*[cols+1];
        Mb[i] = new share*[cols+1];
    }*/


    boolSh2 L(Rs.len);boolSh2 R(L.len);
    vector<uint32_t> Idx(Ind.size());
    uint32_t bitlen = 64;
    string address = "127.0.0.1";
	e_mt_gen_alg mt_alg = MT_OT;
    uint32_t nthreads = 1;
	seclvl seclvl = crypt->get_seclvl();
	ABYParty* party;
    boolSh2 s_Ind(Ind.size());
	
//P1,P2
    if(role != 0){
        //cout<<"P"<<role<<endl;
        
        if(role == 1){
            //party = new ABYParty((e_role)0, address, 7766, seclvl, 32, nthreads, mt_alg);
            //RecvFromParty(&sInd.len, 1, (role+2)%3, XSHARETAG);           //没必要
            for(int i=1; i<rows+1; i++){
                boolRecvFromParty(sInd.s2x, sInd.len, 0, XSHARETAG);
                for(int j=1; j<cols+1; j++){
                    Ma[i][j] = sInd.s2x[j-1];
                    //cout<<"P1 "<<i<<"-"<<j<<" : "<<Ma[i][j]<<endl;
                }
            }
            boolRecvFromParty(s_Ind.s2x, s_Ind.len, 0, XSHARETAG);
       
            

            cout<<"P1 recv end"<<endl;
        }else{
            //party = new ABYParty((e_role)1, address, 7766, seclvl, 32, nthreads, mt_alg);
            for(int i=1; i<rows+1; i++){
                sInd.share(role, 0, 0, Ind.data());
                for(int j=1; j<cols+1; j++){
                    Ma[i][j] = sInd.s2x[j-1];
                    //cout<<"P2 "<<i<<"-"<<j<<" : "<<Ma[i][j]<<" "<<sInd.s2x[i-1]<<endl;
                }
                    
                
            }
            s_Ind.share(role, 0, 0, Ind.data());
            
            //cout<<"share end"<<endl;
            cout<<"P2 share end"<<endl;
           
        }
        
        //cout<<"P"<<role<<" 11111"<<endl;
       
        //vector<Sharing*>& sharings = party->GetSharings();

        //Circuit* bc = sharings[S_BOOL]->GetCircuitBuildRoutine();


        for(int i=1; i<=rows; i++){

            b1[i] = Rs.s2x[2*i-2];

            a2[i] = s_Ind.s2x[2*i-1];
            b2[i] = Rs.s2x[2*i-1];
            a[i] = Ma[i][1];

        }
        //cout<<"第一列"<<endl;
        mpc_agg(role, operation, bitlen, a, a2, b1, b2, Mb, rows, 1, crypt);
        //cout<<"P"<<role<<"第一列结束"<<rows<<" "<<cols<<endl;
        //cout<<"111"<<endl;
        //cout<<"p12==========================================="<<endl;
        //验这里dddd
        int flag = 0;
        for(int j=2; j<=cols; j++){
            //cout<<"j: "<<j<<endl;
            for(int i=1; i<=rows; i++){
                //cout<<", i: "<<i<<endl;
                
                int k = (i-1)/pow(2, j-1);
                int q = (i-1)%(2*(j-1));

                int t = pow(2, j-2);
                
                a[i] = Ma[i][j];
                //a1[i] = Ma[(2*k+1)*t][j-1];
                b1[i] = Mb[(2*k+1)*t][j-1];
                //auto p1 = make_pair(a1, b1);
                int indexb1 = (2*k+1)*pow(2, j-2);
                //cout<<"i: "<<i<<",j: "<<j<<", b1 index: "<<indexb1<<endl;
                if((2*k+1)*pow(2, j-2)>rows){
                    cout<<"wrong"<<endl;
                }
                
                if(q == 0){
                    int ind = (2*k+1)*pow(2, j-1);
                    //a2 = bc->PutSharedINGate(sInd.s2x[ind], bitlen);
                    //b2 = bc->PutSharedINGate(Rs.s2x[ind], bitlen);

                    a2[i] = s_Ind.s2x[ind];
                    b2[i] = Rs.s2x[ind];
                   
                   //cout<<"q==0 i: "<<i<<",j: "<<j<<", b2 index: "<<ind<<endl;
                   if(ind>n){
                    cout<<"wrong"<<endl;
                   }

                }else{
                    int inda = (k+1)*pow(2, j-1) - pow(2, j-2) + q+1 - pow(2, ceil(log2(q+1))-1);
                    int indb = ceil(log2(q+1));
                    
                    
                    
                    a2[i] = Ma[inda][indb];
                    b2[i] = Mb[inda][indb];
                    //cout<<"i: "<<i<<",j: "<<j<<", b2 index: "<<inda<<endl;
                    if(inda>rows || indb>cols)
                        cout<<"wrong"<<endl;

                }

            }
            mpc_agg(role, operation, bitlen, a, a2, b1, b2, Mb, rows, j, crypt);        
        }

        cout<<"P"<<role<<"========================================="<<endl;
        
        
        Data* sst = new Data[L.len];
        //sst[0] = bc->PutSharedINGate(Rs.s2x[0], bitlen);
        L.s2x[0] = Rs.s2x[0];
        //circ->PutPrintValueGate(sst[0], "sst");
        for(int i=2; i<=L.len; i++){
            L.s2x[i-1] = Mb[static_cast<int>(i-pow(2, floor(log2(i-1))))][static_cast<int>(ceil(log2(i)))];
            //circ->PutPrintValueGate(sst[i-1], "sst");
        }
        

        F_shuffle(role, bitlen, Idx, L, R, L.len);

        for(int i=0; i<m; i++){
            O.s2x[i] = R.s2x[i];
            //cout<<"O: "<<O.s2x[i]<<endl;
        }

        cout<<"P"<<role<<" O.len: "<<O.len<<endl;
      
        
        //return O;
    }else{
        //P0    
        //cout<<"P"<<role<<endl;
        for(int i=1; i<rows+1; i++){
            Ma[i][1] = Ind[2*i-1]&Ind[2*i-2];
        }
        for(int j=2; j<cols+1; j++){
            for(int i=1; i<rows+1; i++){
                int k = (i-1)/pow(2, j-1);
                int q = (i-1)%(2*(j-1));

                int t = pow(2, j-2);
                Ma[i][j] = Ma[(2*k+1)*t][j-1];
                

                if(q == 0){
                    int ind = (2*k+1)*pow(2, j-1);
                    Ma[i][j] &= Ind[ind];
                }else{
                    int inda = (k+1)*pow(2, j-1) - pow(2, j-2) + q+1 - pow(2, ceil(log2(q+1))-1);
                    int indb = ceil(log2(q+1));
                    Ma[i][j] &= Ma[inda][indb];
                }
                           
            }
        }
        //sInd.init(cols);
        Data* tmp = new Data[cols];
        for(int i=1; i<rows+1; i++){
            
            for(int j=1; j<cols+1; j++){
                tmp[j-1] = Ma[i][j];
            }
            sInd.share(role, 2, 0, tmp);

            boolSendToParty(sInd.s2x, sInd.len, 1, XSHARETAG);
        }
        s_Ind.share(role, 2, 0, Ind.data());
        boolSendToParty(s_Ind.s2x, s_Ind.len, 1, XSHARETAG);

        
        

        //cout<<"send end"<<endl;

        cout<<"P0 send end"<<endl;
        





        int flag = 0;

        //cout<<"Ind.size(): "<<Ind.size()<<endl;
        
        int reverse_flag = Ind.size()-1;
        for(int i=1; i<Ind.size(); i++){
            if(Ind[i] == 0)
                Idx[flag++] = i-1;
            else{
                /*int index = rand()%(Ind.size()-m)+m;
                while(Idx[index] != 0)
                    index = rand()%(Ind.size()-m)+m;
                */
                Idx[reverse_flag--] = i-1;
                
                }
        }
        Idx[flag] = Ind.size()-1;

        /*cout<<"Idx: =============="<<endl;
        for(int i=0; i<Idx.size(); i++){
            cout<<Idx[i]<<endl;
        }*/

        //cout<<"222"<<endl;
        cout<<"Idx size: "<<Idx.size()<<" "<<flag<<endl;
        for(int i=0; i<Idx.size(); i++){
            //cout<<"Idx: "<<Idx[i]<<endl;
        }
        
        F_shuffle(role, bitlen, Idx, L, R, L.len);
        //cout<<"shuffle"<<endl;
        //cout<<"333"<<endl;
    cout<<"P0 end"<<endl;
        //return O;

    }

}



#endif
