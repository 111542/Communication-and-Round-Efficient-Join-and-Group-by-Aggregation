#include "Sh2.h"

//int buffersize = 20370;
    //class Sh2;
//int buffernum = 20000;
    Sh2::Sh2(int length){
        s2x = new uint64_t[length];
        len = length;
    };
    void Sh2::init(int length){
        s2x = new uint64_t[length];
        len = length;
    };

    void Sh2::share(int role, int dest_role, int data_owner, Data* data){
        
        
        if(role == data_owner){
            //vector<uint64_t> r(len);
            //cout<<len<<endl;
            uint64_t* r = new uint64_t[len];
            //cout<<"111"<<endl;
            for(int i=0; i<len; i++){
                r[i] = rand()%(data[i]+1-1000)+1000;              //为了后面的重随机化时不出现负数
                //cout<<r[i]<<endl;
            }

            boolSendToParty(&len, 1, dest_role, 0);
            
            //Send(comm, r, len, 0);
            boolSendToParty(r, len, dest_role, 1);
    
            for(int i=0; i<len; i++){
                s2x[i] = data[i]-r[i];
                //cout<<s2x[i]<<endl;
            }

        }else{
            boolRecvFromParty(&len, 1, dest_role, 0);
            init(len);

            //uint64_t* data = new uint64_t[len];

            boolRecvFromParty(s2x, len, dest_role, 1);  

        }
    }
    void Sh2::restruct(int role, uint64_t* r_data, int dest_role){

       
        //r_data = new uint64_t[len];
        
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
            r_data[i]+=s2x[i];
        }
    }