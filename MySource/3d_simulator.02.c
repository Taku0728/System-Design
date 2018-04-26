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
#define Survival_cond2 (3) //中皮細胞の生存に必要な周囲の線維芽細胞数
#define Completion_coef (4.5) //中皮細胞の補填に必要な周囲の中皮細胞数
#define Disp_cond (3.0) //線と判断する分散

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

//行動（分裂・遊走・静止）の決定
int f_action(int i, int j, int k); 				
//隣合う遷移細胞が存在するか
int isViable(int i, int j, int k, int type);	
//向かう方向の期待値
double getDvalue(int i0, int j0, int k0, int i1, int j1, int j3);			
//[i0][j0][k0] -> [i1][j1][k1] の分裂・遊走確率係数
double getPcoef(int i0, int j0, int k0, int i1, int j1, int k1);
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
				if (number[i][j][k] != 2) {
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
	int i0, j0, k0, i1, j1, k1, i2, j2, k2;

	//確率用基準変数
	double Pbdy = 0;
	//確率係数用一時変数
	double Ptem = 0;
	//増殖確率係数の合計
	double SumPdiv = 0;
	//遊走確率係数の合計
	double SumPmig = 0;
	//増殖確率係数の記録
	double Pdiv[3][3][3] = {0};
	//遊走確率係数の記録
	double Pmig[3][3][3] = {0};

	//生存条件を満たさなければ消滅
	if (!isViable(i, j, k, type)){
		world[i][j][k] = 0;
		number[i][j][k] = 0;
		return;
	}
	if (number[i][j][k] == -1){
		return;
	}

	//周辺状態を確認
	i0 = max(0, i - 1);
	j0 = max(0, j - 1);
	k0 = max(0, k - 1);
	i1 = min(N - 1, i + 1);
	j1 = min(N - 1, j + 1);
	k1 = min(H - 1, k + 1);
	for (i2 = i0; i2 <= i1; ++i2) {
		for (j2 = j0; j2 <= j1; ++j2) {
			for (k2 = k0; k2 <= k1; ++k2) {
				//確率係数とその累積
				Ptem = getPcoef(i, j, k, i2, j2, k2);
				//中皮細胞の場合、方向バイアスをかける
				if (type == 2) {
					Ptem *= getDvalue(i0, j0, k0, i2, j2, k2); 
				}
				//遊走確率係数の累積
				Pmig[i2 - i0][j2 - j0][k2 - k0] = Ptem;
				SumPmig += Ptem;
				//空セルなら増殖確率係数の累積
				if (world[i2][j2][k2] == 0 && isViable(i2, j2, k2, type)) {
					Pdiv[i2 - i0][j2 - j0][k2 - k0] = Ptem;
					SumPdiv += Ptem;
				}
			}
		}
	}
	
	action = f_action(i, j, k);
	if (action == 0) {
		++number[i][j][k];
		return;
	}

	Ptem = 0;
	//分裂の場合
	if (action == 1) {
		Pbdy = genrand_real1() * SumPdiv;
		if (SumPdiv == 0) {
			if (genrand_real1() <=  P_sleep) {
    			number[i][j][k] = -1;
    		}
    		else{
    			number[i][j][k] = 1;
    		}
			return;
		}
		for (i2 = i0; i2 <= i1; ++i2) {
			for (j2 = j0; j2 <= j1; ++j2) {
				for (k2 = k0; k2 <= k1; ++k2) {
					Ptem += Pdiv[i2 - i0][j2 - j0][k2 - k0];
					if (Ptem < Pbdy) {
						continue;
					}
					world[i2][j2][k2] = type;
					if (genrand_real1() <=  P_sleep) {
    					number[i][j][k] = -1;
    				}
    				else{
    					number[i][j][k] = 1;
    				}
					if (genrand_real1() <=  P_sleep) {
    							number[i2][j2][k2] = -1;
    				}
    				else{
    					number[i2][j2][k2] = 1;
    				}
					return;
				}
			}
		}
	}
	//遊走の場合
	else if (action == 2) {
		Pbdy = genrand_real1() * SumPmig;
		if (SumPmig == 0) {
			++number[i][j][k];
			return;
		}
		for (i2 = i0; i2 <= i1; ++i2) {
			for (j2 = j0; j2 <= j1; ++j2) {
				for (k2 = k0; k2 <= k1; ++k2) {
					Ptem += Pmig[i2 - i0][j2 - j0][k2 - k0];
					if (Ptem >= Pbdy) {
						if (world[i2][j2][k2] != 0) {
							++number[i][j][k];
							return;
						}
						world[i2][j2][k2] = type;
						number[i2][j2][k2] = number[i][j][k] + 1;
						world[i][j][k] = 0;
						number[i][j][k] = 0;
						return;
					}

				}
			}
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

int isViable (int i, int j, int k, int type) {
	double tem = 0;
	double val = 0;
	double val2 = 0;
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
						//頂点の確率係数 = 1/sqrt(3)
						case 0:						
							val += R_sqrt3;			
							break;
						//辺心nの確率係数 = 1/sqrt(2)
						case 1:						
							val += R_sqrt2;			
							break;
						//面心の確率係数 = 1.0
						case 2:						
							val += 1.0;				
							break;
						//中心確率係数 = 0.0
						default:					
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
					tem = 0;
					count = 0;
					if(si == i) ++count;
					if(sj == j) ++count;
					if(sk == k) ++count;
					switch (count) {
						//頂点の確率係数 = 1/sqrt(3)
						case 0:						
							tem = R_sqrt3;			
							break;
						//辺心nの確率係数 = 1/sqrt(2)
						case 1:						
							tem = R_sqrt2;			
							break;
						//面心の確率係数 = 1.0
						case 2:						
							tem = 1.0;				
							break;
						//中心確率係数 = 0.0
						default:					
							break;				
					}
					val += tem;
					//線維芽細胞が隣接する			
					if (world[si][sj][sk] == 1) {
						val2 += tem;
					}
				}
			}
		}
		if (val >= Survival_cond && val2 >= Survival_cond2) {
			return 1;
		}
		else {
			return 0;
		}
	}
	
}

double getDvalue(int i0, int j0, int k0, int i1, int j1, int k1){
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
	int i2, j2, k2;
	int count0 = 0;
	int count1 = 0;
	int val = 0;
	//分散
	double disp[3] = {0};

	for (i2 = i0; i2 <= i1; ++i2) {
		for (j2 = j0; j2 <= j1; ++j2) {
			for (k2 = k0; k2 <= k1; ++k2) {
				if (world[i][j][k] == 2) {
					count0 = 0;
					if(i2 == i) ++count0;
					if(j2 == j) ++count0;
					if(k2 == k) ++count0;
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
					disp[0] += pow(i - i2, 2.0) / 9;
					disp[1] += pow(j - j2, 2.0) / 9;
					disp[2] += pow(k - k2, 2.0) / 9;
				}
				else if (world[i][j][k] == 1 && !count1){
					count0 = 0;
					if (i2 == i) ++count0;
					if (j2 == j) ++count0;
					if (k2 == k) ++count0;
					if (count0 = 2) {
						count1 = 1;
					}
				}
			}
		}
	}
	count0 = 0;
	for (i2 = 0; i2 <= 2; ++i2) {
		if (disp[i] < Disp_cond && val > Completion_coef && count1 != 0) {
			world[i][j][k] = 2;
			number[i][j][k] = genrand_int32()%Tf_m + 1;
			return;
		}
	}
}

double getPcoef (int i0, int j0, int k0, int i1, int j1, int k1) {
	int count = 0;
	if(i1 == i0) ++count;
	if(j1 == j0) ++count;
	if(k1 == k0) ++count;
	switch (count) {
		case 0:							//頂点
			return R_sqrt3;				//確率係数 = 1/sqrt(3)
		case 1:							//辺心
			return R_sqrt2;				//確率係数 = 1/sqrt(2)
		default:						//面心
			return 1.0;				//確率係数 = 1.0
	}
}