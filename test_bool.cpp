#include<iostream>
#include <cstdio>
#include <cstdlib>
#include <unistd.h>
#include <fcntl.h>

#include<vector>
#include <algorithm>
//#include "BoolSh2/boolSh2.cpp"
//#include "LowMC/lowmc.h"
#include "Table/Table.h"
//#include "exp.h"
#include <sys/time.h>
#include <sodium.h>
#include <thread>
#include <chrono>
#include <mutex>


#define INTERFACE_NAME "lo"
#define INTERVAL_SECONDS 1




using namespace std;
mutex mtx;
uint64_t total_bytes = 0;




class Table{
	public:
	vector<uint32_t> K, S;
	int len;

	void print(){
		for(int i=0; i<len; i++){
			cout<<K[i]<<" "<<S[i]<<endl;
		}
	}
};



//uint32_t blocksize = 64;
//uint32_t sboxes = 13, keysize = 128, rounds = 12;

int main(int argc, char** argv){

	uint32_t bitlen = 64;
	int role, data_owner=1;
	uint32_t num, dis_num;
	Table hp1, hp2, cs;


	num = atol(argv[1]);
	dis_num = atol(argv[2]);


	//cin>>role;
	MPI_Init(&argc, &argv);
    init_sharing();

	//role = *argv[1] - '0';
	//cout<<role<<endl;
	//uint16_t port = 8080+role;

	int secparam = 128;
	crypto* crypt = new crypto(secparam, (uint8_t*)const_seed);
	LowMCParams param = {sboxes, keysize, blocksize, keysize, rounds};

	int ssize;
	MPI_Comm_size(MPI_COMM_WORLD, &ssize);

	//cout<<"ssize: "<<ssize<<endl;

	
	
	MPI_Comm_rank(MPI_COMM_WORLD, &role);
	//cout<<role<<endl;

	boolSh2 s_xk(num), s_xs(num), s_yk(dis_num), s_ys(dis_num);
	//cout<<"s_xk.len: "<<s_xk.len<<endl;
	Data* r_data = new Data[num];
//cout<<"ppppppppppppppppppppppppppppppppppppppppppp"<<endl;

	/*if(role == 0){
		cout<<"P"<<role<<":"<<endl;
	}*/

	



	int dest_role;

	
	//cout<<"len: "<<sX.len<<"-------------"<<endl;
	Data* K1, *S1, *K2, *S2;
    K1 = new Data[num];
    S1 = new Data[num];
    K2 = new Data[dis_num];
    S2 = new Data[dis_num];

	struct timeval t1, t2;



	if(role == 1){


        for(int i=0; i<num; i++){
            Data* r1, r2, r3, r4;
            randombytes_buf(&K1[i], sizeof(Data));
            randombytes_buf(&K2[i], sizeof(Data));
            randombytes_buf(&S1[i], sizeof(Data));
            randombytes_buf(&S2[i], sizeof(Data));

        }
		//init(num, dis_num, id_low, dis_low, lev_num, hp1, cs);
		
		

		//hp1.print();
		//cs.print();
		dest_role = 2;
gettimeofday(&t1, 0);


		s_xk.share(role, dest_role, 1, K1);   //
			//cout<<"xk"<<endl;
			//cout<<"------------share end----------------"<<endl;
		s_xs.share(role, dest_role, 1, S1);
			//cout<<"xs"<<endl;
		s_yk.share(role, dest_role, 1, K2);
			cout<<"yk"<<endl;
		s_ys.share(role, dest_role, 1, S2);

		//s_xk.restruct(role, r_data, dest_role);
gettimeofday(&t2, 0);
cout<<"ttttttttttttttttttttttttttttt: "<<t2.tv_sec-t1.tv_sec+(t2.tv_usec-t1.tv_usec)*1e-6<<endl;

		m_vRandomBits.Create(2 * param.blocksize * param.blocksize * param.nrounds + param.nrounds * param.blocksize, crypt);
            //m_vRandomBits.Print(0, m_vRandomBits.GetSize());
		MPI_Send(m_vRandomBits.GetArr(), m_vRandomBits.GetSize(), MPI_UNSIGNED_CHAR, dest_role, XCHANGE_MSG_TAG, MPI_COMM_WORLD);
		
            //send(comm_1, m_vRandomBits.GetArr(), m_vRandomBits.GetSize(), 0);

		        int exp_key_bitlen = param.blocksize * (param.nrounds+1);
        //Use a dummy key for benchmark reasons
	    key.Create(exp_key_bitlen, crypt);

		
		 cout<<"P1 init"<<endl;

			
	}else if(role == 2){

		dest_role = 1;
//这里可能出问题，要测一下
		 
		//cout<<"P"<<role<<":"<<endl;
		cout<<"sxk.len: "<<s_xk.len<<endl;
		s_xk.share(role, dest_role, 1, K1);  
			//cout<<"xk"<<endl;
		s_xs.share(role, dest_role, 1, S1);
			//cout<<"xs"<<endl;
			
			
		s_yk.share(role, dest_role, 1, K2);
			//cout<<"yk"<<endl;
		s_ys.share(role, dest_role, 1, S2);
		//cout<<"ys"<<endl;
		//s_xk.restruct(role, r_data, dest_role);


		
		BYTE* output =  new BYTE[(2 * param.blocksize * param.blocksize * param.nrounds + param.nrounds * param.blocksize)/8];
            //if(recv(comm_0, output, 2 * param.blocksize * param.blocksize * param.nrounds + param.nrounds * param.blocksize, 0) == -1)
                //cout<<"wrong"<<endl;;
		MPI_Recv(output, 2 * param.blocksize * param.blocksize * param.nrounds + param.nrounds * param.blocksize/8, MPI_UNSIGNED_CHAR, dest_role, XCHANGE_MSG_TAG, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
        m_vRandomBits.AttachBuf(output, (uint64_t)(2 * param.blocksize * param.blocksize * param.nrounds + param.nrounds * param.blocksize)/8);
        //m_vRandomBits.Print(0, m_vRandomBits.GetSize());

		        int exp_key_bitlen = param.blocksize * (param.nrounds+1);
        //Use a dummy key for benchmark reasons
	    key.Create(exp_key_bitlen, crypt);	

		
		 cout<<"P2 init"<<endl;
	
	}

	STable sX(num), sY(dis_num), res;
	sX.s_K = s_xk;	
	sX.s_S = s_xs;
	sY.s_K = s_yk;
	sY.s_S = s_ys;
	
	struct timeval begin_join, end_join, begin_groupby, end_groupby, begin_sum, end_sum, begin_max, end_max, begin_avg, end_avg;
    long seconds_join, seconds_groupby, seconds_sum, seconds_max, seconds_avg, micro_join, micro_groupby, micro_sum, micro_max, micro_avg;
    double elapsed_join, elapsed_groupby, elapsed_sum, elapsed_max, elapsed_avg;

    // start timer
	
    gettimeofday(&begin_join, 0);


//	join(role, bitlen, 1, crypt, param, sX, sY, res);
	

cout<<"P"<<role<<" K: "<<res.len<<endl;
	// stop timer
    gettimeofday(&end_join, 0);
    seconds_join = end_join.tv_sec - begin_join.tv_sec;
    micro_join = end_join.tv_usec - begin_join.tv_usec;
    elapsed_join = seconds_join + micro_join * 1e-6;


	string aname = "join";
    //printf("%s \tpassed in:\t\t%f\n", aname.c_str(), elapsed_join);
	cout<<"P"<<role<<" join  setup_time: "<<setup_time<<", total_time: "<<elapsed_join<<endl;
    std::cout <<"P"<<role<< " Total traffic on interface-send_bytes: " <<getsendbytes()<<", recv_bytes: "<<getrecvbytes()  << std::endl;
	cout<<"P"<<role<<"ABY setup send comm: "<<send_setup_comm<<" , recv comm: "<<recv_setup_comm<<endl;
	cout<<"P"<<role<<"ABY online send comm: "<<send_online_comm<<" , recv comm: "<<recv_online_comm<<endl;

	
	vector<uint32_t> Ind(sX.len);
	STable R(sX.len);

	
		Data* rK, *rS;
		rK = new Data[sX.len];
		rS = new Data[sX.len];
	
		
	
	
	/*if(role == 1){

		res.s_S.s2x[0] +=1;
		
		uint32_t tmp = res.s_K.s2x[0];
		res.s_K.s2x[0] = res.s_K.s2x[res.len-2];
		res.s_K.s2x[res.len-2] = tmp;

		//res.s_K.restruct(role, rK, 2);
		//res.s_S.restruct(role, rS, 2);

		//for(int i=0; i<res.len; i++){
		//	cout<<rK[i]<<" "<<rS[i]<<endl;
		//}

		
	}else if(role == 2){
		
		uint32_t tmp = res.s_K.s2x[0];
		res.s_K.s2x[0] = res.s_K.s2x[res.len-2];
		res.s_K.s2x[res.len-2] = tmp;

		//res.s_K.restruct(role, rK, 1);
		//res.s_S.restruct(role, rS, 1);
	}*/

	gettimeofday(&begin_groupby, 0);
	uint32_t m;
	m = group_by(role, bitlen, Ind, &param, crypt, sX, R);
	gettimeofday(&end_groupby, 0);
	seconds_groupby = end_groupby.tv_sec - begin_groupby.tv_sec;
    micro_groupby = end_groupby.tv_usec - begin_groupby.tv_usec;
    elapsed_groupby = seconds_groupby + micro_groupby * 1e-6;


	string name = "group_by";
    printf("%s \tpassed in:\t\t%f\n", name.c_str(), elapsed_groupby);

	/*cout<<"Ind.size(): "<<Ind.size()<<endl;
	if(role == 0)
		for(int i=0; i<Ind.size(); i++){
			cout<<Ind[i]<<endl;
		}
	*/
	//pwd(role, bitlen, res, R, crypt, &param);

	int m_pad=m;
	
//cout<<"P"<<role<<" m: "<<m<<endl;
	if(!isPowerOfTwo(m)){
		m_pad = pow(2, ceil(log2(m)));
	}
	//cout<<"m_pad: "<<m_pad<<endl;

	STable R_pad(m_pad);
	if(role != 0){
		for(int i=0; i<m; i++){
		R_pad.s_K.s2x[i] = R.s_K.s2x[i];
		R_pad.s_S.s2x[i] = R.s_S.s2x[i];

		}
		for(int i=m; i<m_pad; i++){
			R_pad.s_K.s2x[i] = 10000+rand()%1000;
			R_pad.s_S.s2x[i] = rand()%1000+1000;
		}
	}

	

	STable Rsum(m), Rmax(m_pad);
	//Sh2 smax;
	//Sh2 smax(m);
	gettimeofday(&begin_sum, 0);
	aggre_sum(role, bitlen, Ind, R, Rsum, m, crypt);
	
	
	
	gettimeofday(&end_sum, 0);
	seconds_sum = end_sum.tv_sec - begin_sum.tv_sec;
        micro_sum = end_sum.tv_usec - begin_sum.tv_usec;
    elapsed_sum = seconds_sum + micro_sum * 1e-6;

std::cout <<"P"<<role<< " Total traffic on interface-send_bytes: " <<getsendbytes()<<", recv_bytes: "<<getrecvbytes()  << std::endl;
        cout<<"P"<<role<<"ABY setup send comm: "<<send_setup_comm<<" , recv comm: "<<recv_setup_comm<<endl;
        cout<<"P"<<role<<"ABY online send comm: "<<send_online_comm<<" , recv comm: "<<recv_online_comm<<endl;




	name = "sum";
     printf("%s \tpassed in:\t\t%f\n", name.c_str(), elapsed_sum);

	//if(role == 0)
	    //aggre_m(role, 1, Ind, R.s_S, m, crypt);
	//else
	
	vector<uint32_t> Ind_mad(Ind.size()+m_pad-m);
	for(int i=0; i<Ind.size(); i++){
		Ind_mad[i] = Ind[i];
	}
	for(int i=Ind.size(); i<Ind_mad.size(); i++){
		Ind_mad[i] = 0;
	}
	
	gettimeofday(&begin_max, 0);
//	aggre_m(role, 1, Ind_mad, R_pad.s_S, Rmax.s_S, m_pad, crypt);
	gettimeofday(&end_max, 0);
	//cout<<"=============================================================="<<endl;
	seconds_max = end_max.tv_sec - begin_max.tv_sec;
    micro_max = end_max.tv_usec - begin_max.tv_usec;
    elapsed_max = seconds_max + micro_max * 1e-6;


	name = "max";
    printf("%s \tpassed in:\t\t%f\n", name.c_str(), elapsed_max);


    std::cout <<"P"<<role<< " Total traffic on interface-send_bytes: " <<getsendbytes()<<", recv_bytes: "<<getrecvbytes()  << std::endl;
        cout<<"P"<<role<<"ABY setup send comm: "<<send_setup_comm<<" , recv comm: "<<recv_setup_comm<<endl;
        cout<<"P"<<role<<"ABY online send comm: "<<send_online_comm<<" , recv comm: "<<recv_online_comm<<endl;


//cout<<"P"<<role<<"Rmax.s_S.len: "<<Rmax.s_S.len<<endl;

	//cout<<"================="<<endl;
	//uint32_t* max_S = new uint32_t[m];
	/*if(role == 1){




		R.s_K.restruct(role, rK, 2);
		R.s_S.restruct(role, rS, 2);
		//Rmax.s_S.restruct(role, max_S, 2);

		//cout<<"P1 smax.len: "<<smax.len<<endl;

		for(int i=0; i<R.len; i++){
			cout<<rK[i]<<" "<<rS[i]<<endl;
		}

		cout<<"====================="<<endl;
		for(int i=0; i<Rmax.len; i++){
			//cout<<max_S[i]<<endl;
		}
	}else if(role == 2){
		R.s_K.restruct(role, rK, 1);
		R.s_S.restruct(role, rS, 1);
		//Rmax.s_S.restruct(role, max_S, 1);
		//cout<<"P2 smax.len: "<<smax.len<<endl;
	}*/

	MPI_Finalize();
}




//因为join里比较循环里加了break，重复的不会出现，但为了测试方便
