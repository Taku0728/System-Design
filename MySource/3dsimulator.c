#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <math.h>
#include "MT.h" //メルセンヌ・ツイスター//

#ifdef __unix__
#include <unistd.h>
#elif defined _WIN32
#include <windows.h>
#define sleep(x) Sleep(1000 * x) //sleep関数の定義//
#endif

#define Tf (44) //分裂1回あたりの遊走距離(20ミクロン×Tf)//
#define Tf_m (44) //中皮細胞
#define P_sleep (0.05) //休眠確率
#define Dir_val (3.0) //中皮細胞の増殖・遊走の方向バイアス
// #define P_up 0.5 //前進確率
// #define P_down 0.1 //後進確率(%)//
// #define P_side 0.2 //横進確率(%)//
#define P_spring (50) //1hあたりの湧き出し確率(%)//

#define P_f_division (0.5) //線維芽細胞の分裂確率
#define P_m_division (1.0) //中皮細胞の分裂確率
#define P_f_migration (0.3) //線維芽細胞の遊走確率
#define P_m_migration (0.7) //中皮細胞の遊走確率
#define Survival_cond (5) //生存可能な周囲細胞数

#define R_sqrt3 (0.577350) //１分のルート3の高速化
#define R_sqrt2 (0.707107) //１分のルート2の高速化
#define P_room (19.10408) //8/sqrt(3) + 12/sqrt(2) + 6の高速化

#define N (101) //傷の大きさ
#define H (25) //組織の距離
#define TMPFILE "tempfile.tmp" //一時ファイル//
#define GNUPLOT "gnuplot" //gnuplotの場所//
#define INIT_INTERVAL (2) //初期待ち時間(s)//
#define INTERVAL (0.8) //待ち時間(s)//

int world[N][N][H] = {0}; //セルの状態//
int nextworld[N][N][H] = {0}; //次のセルの状態//
int number[N][N][H] = {0}; //細胞の状態//

void fputworld(); //gnuplot出力//
void initworld();
void nextt(int t_count); //状態更新//
void calcnext(int i, int j, int k);	//ルール適用//
void spring();
void completion(int i, int j, int k);	//空きの補完	

int f_action(int i, int j, int k); 				//行動（分裂・遊走・静止）の決定
int isViable(int i, int j, int k, int type);	//隣合う遷移細胞が存在するか
double d_value(int i0, int j0, int k0, int i1, int j1, int j3);			//向かう方向の期待値
// void top_count(world,t);
// void sq_count(world,t); 

int main(int argc,char *argv[]){
    int t = 0; //経過時間(h)//
	int i,j,k;
	int m_cell = 0; //中皮細胞の個数//
	int f_cell = 0; //線維芽細胞の個数//
	int t_count; //ステップ//
	int scale;
	int MAXT;
	int fission_total = 0;
	FILE *pipe;
	FILE *file;
	FILE *fp;

	init_genrand((unsigned)time(NULL)); //乱数の初期化//
	
	scale = Tf/22; //1hあたりの遊走距離(20ミクロン×scale)//
	MAXT = scale*300; //7日間のシミュレーション//

    //初期条件//
	printf("t = 0 h\n");
	initworld();
	fputworld();

    if((pipe = popen(GNUPLOT " -persist","w")) == NULL){
		fprintf(stderr," Cannot open the pipe! \n");
		exit(1);
	}
	
	//gnuplotの設定//
	fprintf(pipe,"unset key\n");
	fprintf(pipe,"set xrange [-10:%d]\n",N+10);
	fprintf(pipe,"set yrange [-10:%d]\n",N+10);
	fprintf(pipe,"set zrange [0:%d]\n",H);
	fprintf(pipe, "set size rectangular\n");
	fprintf(pipe, "unset xtics\n");
	fprintf(pipe, "unset ytics\n");
	fprintf(pipe, "unset ztics\n");
	
	//グラフに出力//
	fprintf(pipe, "set title 't = %d h'\n",t);
	fprintf(pipe, "set title font 'Arial,15'\n");
	fprintf(pipe,"splot \"" TMPFILE "\" index 0 w p ps 0.1 pt 65 lt 5, \"" TMPFILE "\" index 1 w p ps 0.1 pt 65 lt 2\n");
	fprintf(pipe,"name='move%d'\n load 'savegif.gp'\n",t);
	fflush(pipe);
	sleep(INIT_INTERVAL);
	
	//細胞数をカウント//
	for(i=0;i<=N-1;++i){
		for(j=0;j<=N-1;++j){
			for(k=0;k<=H-1;++k){
				if(world[i][j][k] == 1){
					++m_cell;
				}
				else if (world[i][j][k] == 2){
					++f_cell;
				}
			}
		}
	}
	
	//細胞数を出力//
	file = fopen("count_m_cell.txt","w");
	fprintf(file,"%d\n",m_cell);
	fclose(file);
	file = fopen("count_f_cell.txt","w");
	fprintf(file,"%d\n",f_cell);
	fclose(file);
	

	for(t_count=1;t_count<=MAXT;++t_count){
		if(t_count % scale == 0){
			t += 1;  //グラフ時間表示更新//
			// top_count(world,t); //治癒速度測定①//
			// sq_count(world,t); //治癒速度測定②//
			if((genrand_int32()%100 + 1) <= P_spring){
				spring(world,number); //湧き出し//
			}
			spring();
		}
		//状態更新//
		nextt(t_count);

		//グラフに出力//
		printf("t = %d h\n",t);
		fputworld();
		fprintf(pipe, "set title 't = %d h'\n",t);
		fprintf(pipe, "set title font 'Arial,15'\n");
		fprintf(pipe,"splot \"" TMPFILE "\" index 0 w p ps 0.1 pt 65 lt 5, \"" TMPFILE "\" index 1 w p ps 0.1 pt 65 lt 2\n");
		fprintf(pipe,"name='move%d'\n load 'savegif.gp'\n",t);
		fflush(pipe);
		sleep(INTERVAL);
	}
	return 0;
}

void fputworld(){
	int i,j,k;
	FILE *fp;

	if((fp = fopen(TMPFILE,"w")) == NULL){ 
		fprintf(stderr," Cannot open the file! \n");
		exit(1);
	}
	
	//線維芽細胞の位置を出力//
	for(i=0;i<=N-1;++i){
		for(j=0;j<=N-1;++j){
            for (k=0;k<=H-1;++k){
			    if(world[i][j][k] == 1){ 
	  	            fprintf(fp,"%d %d %d\n",i,j,k);
                }
			}
		}
	}
	
	fprintf(fp,"\n\n");
	
	
	//中皮細胞の位置を出力//
	for(i=0;i<=N-1;++i){
		for(j=0;j<=N-1;++j){
            for (k=0;k<=H-1;++k){
			    if(world[i][j][k] == 2){ 
	  	            fprintf(fp,"%d %d %d\n",i,j,k);
                }
			}
		}
	}
	
	fclose(fp);
}

void initworld(){
  int i;
  int j;
  int k;

	//初期条件//
    for(i=0;i<=N-1;++i){
        for (j=0;j<=N-1;++j){
            world[i][j][0] = 1;
		    number[i][j][0] = genrand_int32()%Tf + 1;
    	    world[i][j][H-1] = 1;
    	    number[i][j][H-1] = genrand_int32()%Tf + 1;
        }
    }
	for(i=0;i<=N-1;++i){
		world[i][0][1] = 2;
		number[i][0][1] = genrand_int32()%Tf_m + 1;
		world[0][i][1] = 2;
		number[0][i][1] = genrand_int32()%Tf_m + 1;
		world[i][N-1][1] = 2;
		number[i][N-1][1] = genrand_int32()%Tf_m + 1;
		world[N-1][i][1] = 2;
		number[N-1][i][1] = genrand_int32()%Tf_m + 1;
		world[i][0][H-2] = 2;
		number[i][0][H-2] = genrand_int32()%Tf_m + 1;
		world[0][i][H-2] = 2;
		number[0][i][H-2] = genrand_int32()%Tf_m + 1;
		world[i][N-1][H-2] = 2;
		number[i][N-1][H-2] = genrand_int32()%Tf_m + 1;
		world[N-1][i][H-2] = 2;
		number[N-1][i][H-2] = genrand_int32()%Tf_m + 1;
	}
}

void spring(){
	int site_x,site_y;
	long i,j,k;
	int ztop[N][N];
	int zbot[N][N];
	int bort;
	
	//天井と床を特定
	for(i=0;i<=N-1;++i){
		for(j=0;j<=N-1;++j){
			k = 0;
			zbot[i][j] = 0;
			ztop[i][j] = N-1;
			while(k<=H-1){
				if(world[i][j][k] == 0){
					zbot[i][j] = k;
					break;
				}
				++k;
			}
			k = N-1;
			while(k>=0){
				if(world[i][j][k] == 0){
					ztop[i][j] = k;
					break;
				}
				--k;
			}
		}
	}

	//湧き出しの位置を決定
	for(i=1;i<=(10000/ (N*N));++i){
		site_x = genrand_int32()%(N-2) + 1;
		site_y = genrand_int32()%(N-2) + 1;
		bort = genrand_int32()%2;
		if (bort == 0){
			if(world[site_x][site_y] == 0){ //空きセルに湧き出し//
				world[site_x][site_y][zbot[site_x][site_y]] = 2;
				number[site_x][site_y][zbot[site_x][site_y]] = genrand_int32()%Tf_m + 1;
			}
		}
		else{
			if(world[site_x][site_y] == 0){ //空きセルに湧き出し//
				world[site_x][site_y][ztop[site_x][site_y]] = 2;
				number[site_x][site_y][ztop[site_x][site_y]] = genrand_int32()%Tf_m + 1;
			}
		}
	}
}

void nextt(int t_count){
	int i,j,k;
	
	// for(i=0;i<=N-1;++i){
	// 	for(j=0;j<=N-1;++j){
	// 		for(k=0;k<=H-1;++k){
	// 			nextworld[i][j][k] = world[i][j][k]; //前状態の引継//
	// 		}
	// 	}
	// }
	//状態更新を開始する正方形の頂点を決定//
	for (i=0;i<=N-1;++i) {
		for (j=0;j<=N-1;++j) {
			for (k=0;k<=H-1;++k) {
				if (number[i][j][k] != 0) {
					switch (t_count%4){
						case 1:
							calcnext(i, j, k);
							break;
						case 2:
							calcnext(N - 1 - i, j, k);
							break;
						case 3:
							calcnext(i, N - 1 - j, k);
							break;
						case 4:
							calcnext(N - 1 - i, N - 1 - j, k);
							break;
						case 5:
							calcnext(i, j, H - 1 - k);
							break;
						case 6:
							calcnext(N - 1 - i, j, H - 1 - k);
							break;
						case 7:
							calcnext(i, N - 1 - j, H - 1 - k);
							break;
						default:
							calcnext(N - 1 - i, N - 1 - j, H - 1 - k);
							break;
					}
				}
			}
		}
	}	
	
	//空きの補完
	for (i=0;i<=N-1;++i) {
		for (j=0;j<=N-1;++j) {
			for (k=0;k<=H-1;++k) {
				if (number[i][j][k] == 0) {
					switch (t_count%4){
						case 1:
							completion(i, j, k);
							break;
						case 2:
							completion(N - 1 - i, j, k);
							break;
						case 3:
							completion(i, N - 1 - j, k);
							break;
						case 4:
							completion(N - 1 - i, N - 1 - j, k);
							break;
						case 5:
							completion(i, j, H - 1 - k);
							break;
						case 6:
							completion(N - 1 - i, j, H - 1 - k);
							break;
						case 7:
							completion(i, N - 1 - j, H - 1 - k);
							break;
						default:
							completion(N - 1 - i, N - 1 - j, H - 1 - k);
							break;
					}
				}
			}
		}
	}	



	//境界条件・中皮細胞//
	for(i=0;i<=N-1;i++){
		if(world[i][0][1] != 2){
			world[i][0][1] = 2;
			number[i][0][1] = genrand_int32()%Tf_m + 1; //供給細胞の状態//
		}
		if(world[i][N-1][1] != 2){
			world[i][N-1][1] = 2;
			number[i][N-1][1] = genrand_int32()%Tf_m + 1; //供給細胞の状態//
		}
		if(world[i][0][H-2] != 2){
			world[i][0][H-2] = 2;
			number[i][0][H-2] = genrand_int32()%Tf_m + 1; //供給細胞の状態//
		}
		if(world[i][N-1][H-2] != 2){
			world[i][N-1][H-2] = 2;
			number[i][N-1][H-2] = genrand_int32()%Tf_m + 1; //供給細胞の状態//
		}
		if(world[0][i][1] != 2){
			world[0][i][1] = 2;
			number[0][i][1] = genrand_int32()%Tf_m + 1; //供給細胞の状態//
		}
		if(world[N-1][i][1] != 2){
			world[N-1][i][1] = 2;
			number[N-1][i][1] = genrand_int32()%Tf_m + 1; //供給細胞の状態//
		}
		if(world[0][i][H-2] != 2){
			world[0][i][H-2] = 2;
			number[0][i][N-2] = genrand_int32()%Tf_m + 1; //供給細胞の状態//
		}
		if(world[N-1][i][H-2] != 2){
			world[N-1][i][H-2] = 2;
			number[N-1][i][H-2] = genrand_int32()%Tf_m + 1; //供給細胞の状態//
		}
	}

	//境界条件・線維芽細胞//
	for(i=0;i<=N-1;++i){
		for(j=0;j<=N-1;++j){
			if(world[i][j][0] != 1){
				world[i][j][0] = 1;
				number[i][j][0] = genrand_int32()%Tf + 1; //供給細胞の状態//
			}
			if(world[i][j][H-1] != 1){
				world[i][j][H-1] = 1;
				number[i][j][H-1] = genrand_int32()%Tf + 1; //供給細胞の状態//
			}
		}
	}
	
	//状態更新//
	// for(i=0;i<=N-1;++i){
	// 	for(j=0;j<=N-1;++j){
	// 		for (k=0;k<=H-1;++k){
	// 			world[i][j][k] = nextworld[i][j][k];
	// 		}
	// 	}
	// }
	
}

void calcnext(int i,int j,int k){
	int action;
	int type = world[i][j][k];
	int i0, j0, k0, i1, j1, k1;
	int si, sj, sk;
	int di, dj, dk;
	double p = 0;				//確率用変数
	double p_dir = 1;			//方向による確率用変数
	double p_tem = 0;			//確率用変数
	double p_all = 0;			//確率係数の合計
	int count = 0;				//カウント用変数
	double val = 0.0;			//隣接カウント用変数
	int Room[3][3][3] = {0};	//周辺の状態を記録
	int Dir_P[3][3][3] = {0};	//方向バイアスの記録
	int cond = 1;				//条件用変数

	//周辺状態を確認
	i0 = max(0, i - 1);
	j0 = max(0, j - 1);
	k0 = max(0, k - 1);
	i1 = min(N - 1, i + 1);
	j1 = min(N - 1, j + 1);
	k1 = min(H - 1, k + 1);
	for (si = i0; si <= i1; ++si) {
		for (sj = j0; sj <= j1; ++sj) {
			for (sk = k0; sk <= k1; ++sk) {
				if (world[si][sj][sk] != 0) {
					count = 0;
					if(si == i) ++count;
					if(sj == j) ++count;
					if(sk == k) ++count;
					switch (count) {
						case 0:							//頂点
							val +=  R_sqrt3;			//確率係数 = 1/sqrt(3)
							break;
						case 1:							//辺心
							val += R_sqrt2;				//確率係数 = 1/sqrt(2)
							break;
						default:						//面心
							val += 1.0;					//確率係数 = 1.0
							break;
					}
				}
				else{
					if (isViable(si, sj, sk, type)) {	//生存可能か
						Room[si - i0][sj - j0][sk - k0] = 1;
						p_dir = 1;
						if (type == 2){
							p_dir = d_value(i, j, k, si, sj, sk);
							Dir_P[si - i0][sj - j0][sk - k0] = p_dir;
						}
						count = 0;
						if(si == i) ++count;
						if(sj == j) ++count;
						if(sk == k) ++count;
						switch (count) {
							case 0:							//頂点
								p_all += p_dir * R_sqrt3;	//確率係数 = 1/sqrt(3)
								break;
							case 1:							//辺心
								p_all += p_dir * R_sqrt2;	//確率係数 = 1/sqrt(2)
								break;
							default:						//面心
								p_all += p_dir;				//確率係数 = 1.0
								break;
						}
					}
				}
			}
		}
	}

	//生存条件より小さければ消滅
	if (val < Survival_cond) {
		world[i][j][k] = 0;
		number[i][j][k] = 0;
		return;
	}
	if (number[i][j][k] == -1){
		return;
	}
	if (p_all == 0) {
		//番号の更新
		++number[i][j][k];
		if (number[i][j][k] > Tf) {
			if(genrand_real1() <=  P_sleep){
    			number[i][j][k] = -1;
    		}
    		else{
    			number[i][j][k] = 1;
    		}
		}
		return;
	}
	//空いているセルがない<<

	//空いているセルがある>>
	//行動を決定
	action = f_action(i, j, k);
	if (action == 0) {
		++number[i][j][k];
		return;
	}
	p_tem = 0;
	cond = 0;
	count = 0;
	//分裂
	if (action == 1) {
		p = genrand_real1() * p_all;			//正規化;
		for (si = i0; si <= i1; ++si) {
			for (sj = j0; sj <= j1; ++sj) {
				for (sk = k0; sk <= k1; ++sk) {
					if (world[si][sj][sk] != 0 || cond) {
						continue;
					}
					//生存制約
					if (Room[si - i0][sj - j0][sk - k0] != 1) {
						continue;
					}
					p_dir = 1;
					//中皮細胞の方向バイアス
					if (type == 2){
						p_dir = Dir_P[si - i0][sj - j0][sk - k0];
					}
					count = 0;
					if (si == i) ++count;
					if (sj == j) ++count;
					if (sk == k) ++count;
					switch (count) {
						case 0:							//頂点
							p_tem += p_dir * R_sqrt3;	//確率係数 = 1/sqrt(3)
							break;
						case 1:							//辺心
							p_tem += p_dir * R_sqrt2;	//確率係数 = 1/sqrt(2)
							break;
						case 2:							//面心
							p_tem += p_dir;				//確率係数 = 1.0
							break;
						default:						//中心
							break;						//確率係数 = 0
					}
					if (p_tem >= p) {					//分裂方向決定
						di = si;
						dj = sj;
						dk = sk;
						cond = 1;
						break;
					}
				}
			}
		}
		//番号の更新//
		if (cond){
			world[di][dj][dk] = type;
    		number[di][dj][dk] = 1;
			if (genrand_real1() <=  P_sleep) {
    			number[i][j][k] = -1;
    		}
    		else{
    			number[i][j][k] = 1;
    		}
		}
	}

	//遊走//
    else if (action == 2) {
		// p = genrand_real1() * P_room;
		p = genrand_real1() * p_all;
		for (si = i0; si <= i1; ++si) {
			for (sj = j0; sj <= j1; ++sj) {
				for (sk = k0; sk <= k1; ++sk) {
					if (world[si][sj][sk] != 0 || cond) {
						continue;				
					}
					//中皮細胞の増殖制約
					if (Room[si - i0][sj - j0][sk - k0] != 1) {
						continue;
					}
					p_dir = 1;
					//中皮細胞の方向バイアス
					if (type == 2){
						p_dir = Dir_P[si - i0][sj - j0][sk - k0];
					}
					count = 0;
					if (si == i) ++count;
					if (sj == j) ++count;
					if (sk == k) ++count;
					switch (count) {
						case 0:							//頂点
							p_tem += p_dir * R_sqrt3;	//確率係数 = 1/sqrt(3)
							break;
						case 1:							//辺心
							p_tem += p_dir * R_sqrt2;	//確率係数 = 1/sqrt(2)
							break;
						case 2:							//面心
							p_tem += p_dir;				//確率係数 = 1.0
							break;
						default:						//中心
							break;						//確率係数 = 0
					}
					if (p_tem >= p) {
						di = si;
						dj = sj;
						dk = sk;
						cond = 1;
						break;
					}
				}
			}
		}
		if (cond){
			world[di][dj][dk] = type;
			number[di][dj][dk] = number[i][j][k] + 1;
			world[i][j][k] = 0;
			number[i][j][k] = 0;
		}
    }
}

int f_action(int i, int j, int k){
	int type = world[i][j][k];
	double p_div;
	int Tf_tem;
	if (type == 1) {
		p_div = P_f_division;				//線維芽細胞
		Tf_tem = Tf;
	}	
	else {
		p_div = P_m_division;				//中皮細胞
		Tf_tem = Tf_m;
	}
	//分裂周期(M)かどうか
	if(number[i][j][k] == Tf_tem) {
		if (genrand_real1() <= p_div){
			return 1;						//分裂
		}
		else {
			number[i][j][k] = 1;
		}					
	}
	if (type == 1){
		if (genrand_real1() <= P_f_migration) {
			return 2;							//遊走
		}
	}
	else {
		if (genrand_real1() <= P_m_migration) {
			return 2;							//遊走
		}
	}
	return 0;								//静止
}

int isViable(int i, int j, int k, int type){
	double val = 0.0;
	double val2 = 0.0;
	int count = 0;
	int i0 = max(0, i - 1);
	int j0 = max(0, j - 1);
	int k0 = max(0, k - 1);

	int i1 = min(N - 1, i + 1);
	int j1 = min(N - 1, j + 1);
	int k1 = min(H - 1, k + 1);
	int si, sj, sk;
	
	//線維芽細胞
	if (type == 1){
		//中皮細胞の上では増殖できない
		if (k < (H - 1)/2) {
			if (world[i][j][k0] == 2 || world[i][j][max(0, k0 - 1)] == 2){
				return 0;
			}
		}
		else {
			if (world[i][j][k1] == 2 || world[i][j][min(H - 1, k1 + 1)] == 2){
				return 0;
			}
		}
		
		for (si = i0; si <= i1; ++si) {
			for (sj = j0; sj <= j1; ++sj) {
				for (sk = k0; sk <= k1; ++sk) {
					if (world[si][sj][sk] == 0) {
						continue;
					}
					count = 0;
					if(si == i) ++count;
					if(sj == j) ++count;
					if(sk == k) ++count;
					switch (count) {
						case 0:						//頂点
							val += R_sqrt3;			//確率係数 = 1/sqrt(3)
							break;
						case 1:						//辺心
							val += R_sqrt2;			//確率係数 = 1/sqrt(2)
							break;
						default:					//面心
							val += 1.0;				//確率係数 = 1.0
							break;
					}
				}
			}
		}
		if (val >= Survival_cond) {
			return 1;
		}
		else {
			return 0;
		}
	}

	//中皮細胞
	else if (type == 2){
		//線維芽細胞の下では増殖できない
		if (k < (H - 1)/2) {
			if (world[i][j][k1] == 1 || world[i][j][min(H - 1, k1 + 1)] == 1){
				return 0;
			}
		}
		else {
			if (world[i][j][k0] == 1 || world[i][j][max(0, k0 - 1)] == 1){
				return 0;
			}
		}
		for (si = i0; si <= i1; ++si) {
			for (sj = j0; sj <= j1; ++sj) {
				for (sk = k0; sk <= k1; ++sk) {
					if (world[si][sj][sk] == 0) {
						continue;
					}
					count = 0;
					if(si == i) ++count;
					if(sj == j) ++count;
					if(sk == k) ++count;
					switch (count) {
						case 0:						//頂点
							val += R_sqrt3;			//確率係数 = 1/sqrt(3)
							break;
						case 1:						//辺心
							val += R_sqrt2;			//確率係数 = 1/sqrt(2)
							break;
						default:					//面心
							val += 1.0;				//確率係数 = 1.0
							break;
					}
					//線維芽細胞が隣接する
					if (world[si][sj][sk] == 1) {
						switch (count) {
							case 0:						//頂点
								val2 += R_sqrt3;		//確率係数 = 1/sqrt(3)
								break;
							case 1:						//辺心
								val2 += R_sqrt2;		//確率係数 = 1/sqrt(2)
								break;
							default:					//面心
								val2 += 1.0;			//確率係数 = 1.0
								break;
						}
					}
				}
			}
		}
		if (val >= Survival_cond && val2 >= 2.0) {
			return 1;
		}
		else {
			return 0;
		}
	}
	
}

double d_value(int i0, int j0, int k0, int i1, int j1, int k1){
	double dist0;		//出発地点と中心の距離
	double dist1;		//目的地点と中心の距離
	int count;
	double p_tem = 0;
	dist0 = sqrt(pow(i0 - (N - 1)/2.0, 2.0) + pow(j0 - (N - 1)/2.0, 2.0));
	dist1 = sqrt(pow(i1 - (N - 1)/2.0, 2.0) + pow(j1 - (N - 1)/2.0, 2.0) + pow(k1 - k0, 2.0));
	if (i1 == i0) ++count;
	if (j1 == j0) ++count;
	if (k1 == k0) ++count;
	switch (count) {
		case 0:							//頂点
			p_tem = R_sqrt3;			//確率係数 = 1/sqrt(3)
			break;
		case 1:							//辺心
			p_tem = R_sqrt2;			//確率係数 = 1/sqrt(2)
			break;
		case 2:							//面心
			p_tem = 1.0;				//確率係数 = 1.0
			break;
		default:						//中心
			break;						//確率係数 = 0
	}
	return pow(Dir_val, (dist0 - dist1)* p_tem);
}

void completion(int i, int j, int k) {
	int i0 = max(0, i - 1);
	int j0 = max(0, j - 1);
	int k0 = max(0, k - 1);
	int i1 = min(N - 1, i + 1);
	int j1 = min(N - 1, j + 1);
	int k1 = min(H - 1, k + 1);
	int si, sj, sk;
	int count0 = 0;
	int count1 = 0;
	int val = 0;

	for (si = i0; si <= i1; ++si) {
		for (sj = j0; sj <= j1; ++sj) {
			for (sk = k0; sk <= k1; ++sk) {
				if (world[i][j][k] == 2) {
					count0 = 0;
					if(si == i) ++count0;
					if(sj == j) ++count0;
					if(sk == k) ++count0;
					switch (count0) {
						case 0:						//頂点
							val += R_sqrt3;			//確率係数 = 1/sqrt(3)
							break;
						case 1:						//辺心
							val += R_sqrt2;			//確率係数 = 1/sqrt(2)
							break;
						default:					//面心
							val += 1.0;				//確率係数 = 1.0
							break;
					}
					break;							//同じ高さでは１度のみ
				}
				else if (world[i][j][k] == 1 && !count1){
					count0 = 0;
					if (si == i) ++count0;
					if (sj == j) ++count0;
					if (sk == k) ++count0;
					if (count0 = 2) {
						count1 = 1;
					}
				}
			}
		}
	}
	if (val > 4.5 && count1 > 0){
		world[i][j][k] = 2;
		number[i][j][k] = genrand_int32()%Tf_m + 1;
	}
}