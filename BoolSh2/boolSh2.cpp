#include "boolSh2.h"
//int buffersize = 20370;
    //class Sh2;
int buffernum = 20000;
    boolSh2::boolSh2(int length){
        s2x = new Data[length];
        len = length;
    };
    void boolSh2::init(int length){
        s2x = new Data[length];
        len = length;
    };
// Sets random seed
void init_sharing() {
  if (sodium_init() == -1) {
    printf("Cannot initialize sodium library.\n");
    exit(-1);
  }
}
template<typename T>
    void boolSh2::share(int role, int dest_role, int data_owner, T* data){
        
        
        if(role == data_owner){
            //vector<uint32_t> r(len);
            //cout<<len<<endl;
            Data* r = new Data[len];
            //cout<<"111"<<endl;
            for(int i=0; i<len; i++){
                randombytes_buf(&r[i], sizeof(Data));
                s2x[i] = data[i]^r[i];
                //cout<<"P1: "<<i<<" "<<r[i]<<" "<<s2x[i]<<endl;
            }
            Data l = len;
           // cout<<"P1 l: "<<l<<endl;
            boolSendToParty(&l, 1, dest_role, 0);
            
            //Send(comm, r, len, 0);
            boolSendToParty(r, len, dest_role, 1);
    

        }else{
            uint64_t l;
            boolRecvFromParty(&l, 1, dest_role, 0);
            //len = lenth;
            //cout<<"l: "<<l<<endl;
            len = l;
            init(len);

            //uint32_t* data = new uint32_t[len];

            boolRecvFromParty(s2x, len, dest_role, 1);  
           

        }
    }
    void boolSh2::restruct(int role, Data* r_data, int dest_role){

       
        //r_data = new uint32_t[len];
        
        if(role == 1){
            //Send(comm, s2x, len, 0);
            boolSendToParty(s2x, len, dest_role, 3);
            //Recv(comm, r_data, len, 0);
            boolRecvFromParty(r_data, len, dest_role, 4);

            /*for(int i=0; i<len; i++){
                cout<<r_data[i]<<endl;
            }*/

        }else{
            //Recv(comm, r_data, len, 0);
            boolRecvFromParty(r_data, len, dest_role, 3);
            //Send(comm, s2x, len, 0);
            boolSendToParty(s2x, len, dest_role, 4);
        }


        for(int i=0; i<len; i++){
            r_data[i]^=s2x[i];
        }
    }