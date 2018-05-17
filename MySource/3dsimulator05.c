#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <math.h>
#include "MT.h" //メルセンヌ・ツイスター//


#ifdef _WIN32
#include <windows.h>
#define SLEEP(x) Sleep(1000 * x) //sleep関数の定義//
#else
#include <unistd.h>
#define SLEEP(x) usleep(x * 1000000)

int max(int a, int b) {
	return a < b ? b : a;
}

int min(int a, int b) {
	return a > b ? b : a;
}
#endif

#define Tf (44) //分裂1回あたりの遊走距離(20ミクロン×Tf)//
#define Tf_m (28) //中皮細胞
#define P_sleep (0.01) //休眠確率

#define P_f_division (0.1) //線維芽細胞の分裂確率
#define P_m_division (1.0) //中皮細胞の分裂確率
#define P_f_migration (0.2) //線維芽細胞の遊走確率
#define P_m_migration (0.8) //中皮細胞の遊走確率
#define Survival_cond (4.0) //生存可能な周囲細胞数
#define Survival_cond2 (1.0) //中皮細胞の生存に必要な周囲の線維芽細胞数
#define Dir_val (3.0) //中皮細胞の増殖・遊走の方向バイアス
#define Completion_min (4) //中皮細胞の補填に必要な周囲の中皮細胞数の最小
#define Completion_max (8) //中皮細胞の補填に必要な周囲の中皮細胞数の最大

#define N (51) //傷の大きさ
#define H (21) //組織の距離
#define GNUPLOT "gnuplot" //gnuplotの場所//
#define INIT_INTERVAL (5.0) //初期待ち時間(s)//
#define INTERVAL (2.0) //待ち時間(s)//

//ACCELARATION OF 1/SQRT(3)
const double R_sqrt3 = 0.57735026919;
//ACCELARATION OF 1/SQRT(2)
const double R_sqrt2 = 0.70710678118;

//FREQUENCY OF POLUMER DIFFUSION PER t_count
const int FREQ = 1;
//DIFFUSION COEFFICIENT OF FIBRIN
const double COE_DIF = 0.15;
//FIBRINOLYSIS REACTION RARE
const double RATE_FIB = 0.15;

//PROBABILITY OF SPRING [t_count*vacancy]
const double P_spring = 0.0001;
//TEMP FILE FOR OUTPUT
const char *TMPFILE = "tempfile.tmp";
const char *TMPFILE2 = "tempfile2.tmp";
//DISPLAY SIZE OF CELLS
const float FONTSIZE = 1.1; //細胞の表示サイズ

int world[N][N][H] = {0}; //セルの状態//
int prevworld[N][N][H] = {0}; //次のセルの状態//
int number[N][N][H] = {0}; //細胞の状態//

//CONCENTRATION OF FIBRIN
double CFibrin[N][N][H] = {0};
double NextCFibrin[N][N][H] = {0};
//CONCENTRATION OF PLASMIN
double CPlasmin[N][N][H] = {0};
double NextCPlasmin[N][N][H] = {0};

//gnuplot出力//
void showworld(FILE *pipe, FILE *pipe2, int t, int state);
void initworld();
void nextt(int t_count); //状態更新//
void calcnext(int i, int j, int k);        //ルール適用//
void spring();
void completion(int i, int j, int k);        //空きの補完
//配列のシャッフル
void randsortarray(int *a,int len);

//DIFFUSION AND FIBRINOLYSIS OF FIBRIN
void actFibrin(int i, int j, int k);

//行動（分裂・遊走・静止）の決定
int f_action(int i, int j, int k);                                 
//生存可能かの判断
int isViable(int i, int j, int k, int type);
//存続可能かの判断
int isSurvival(int i, int j, int k, int type);
//向かう方向の期待値
double getDvalue(int i0, int j0, int k0, int i1, int j1, int j3);          
//[i0][j0][k0] -> [i1][j1][k1] の分裂・遊走確率係数
double getPcoef(int i0, int j0, int k0, int i1, int j1, int k1);


int main(int argc,char *argv[]){
  int t = 0; //経過時間(h)//
	int i,j,k;
	int e;
	int m_cell = 0; //中皮細胞の個数//
	int f_cell = 0; //線維芽細胞の個数//
	int t_count; //ステップ//
	int scale;
	int MAXT;
	int fission_total = 0;
	FILE *pipe;
	FILE *pipe2;
	FILE *file;
	FILE *fp;


	init_genrand((unsigned)time(NULL)); //乱数の初期化//
	
	scale = Tf/22; //1hあたりの遊走距離(20ミクロン×scale)//
	MAXT = scale*1000; //7日間のシミュレーション//

	if((pipe = popen(GNUPLOT " -persist","w")) == NULL){
		fprintf(stderr," Cannot open the pipe! \n");
		exit(1);
	}
	//gnuplotの設定//
	fprintf(pipe,"unset key\n");
	fprintf(pipe,"set xrange [-10:%d]\n",N+10);
	fprintf(pipe,"set yrange [-10:%d]\n",N+10);
	fprintf(pipe,"set zrange [0:%d]\n",H);
	fprintf(pipe, "unset xtics\n");
	fprintf(pipe, "unset ytics\n");
	fprintf(pipe, "unset ztics\n");

	if((pipe2 = popen(GNUPLOT " -persist","w")) == NULL){
		fprintf(stderr," Cannot open the pipe! \n");
		exit(1);
	}
	//gnuplotの設定//
	fprintf(pipe2,"unset key\n");
	fprintf(pipe2,"set xrange [-10:%d]\n",N+10);
	fprintf(pipe2,"set yrange [-10:%d]\n",N+10);
	fprintf(pipe2,"set zrange [0:%d]\n",H);
	fprintf(pipe2, "unset xtics\n");
	fprintf(pipe2, "unset ytics\n");
	fprintf(pipe2, "unset ztics\n");
	
  	//初期条件//
	printf("t = 0 h\n");
	initworld();

	//グラフに出力//
	showworld(pipe, pipe2, t, 1);
	SLEEP(INIT_INTERVAL);
	
	//細胞数をカウント//
	// for(i=0;i<=N-1;++i){
	// 	for(j=0;j<=N-1;++j){
	// 		for(k=0;k<=H-1;++k){
	// 			if(world[i][j][k] == 1){
	// 				++m_cell;
	// 			}
	// 			else if (world[i][j][k] == 2){
	// 				++f_cell;
	// 			}
	// 		}
	// 	}
	// }
	
	//細胞数を出力//
	// file = fopen("count_m_cell.txt","w");
	// fprintf(file,"%d\n",m_cell);
	// fclose(file);
	// file = fopen("count_f_cell.txt","w");
	// fprintf(file,"%d\n",f_cell);
	// fclose(file);
	


	for(t_count=1;t_count<=MAXT;++t_count){
		if(t_count % scale == 0){
			t += 1;  //グラフ時間表示更新//
			// top_count(world,t); //治癒速度測定①//
			// sq_count(world,t); //治癒速度測定②//
		}
		for(i=0;i<=N-1;++i){
		    for(j=0;j<=N-1;++j){
		        for (k=0;k<=H-1;++k){
					prevworld[i][j][k] = world[i][j][k];
		        }
		    }
		}
		spring();

		//状態更新//
		nextt(t_count);
		//終了条件//
		e = N*N/100;
		for(i=0;i<=N-1;++i){
		    for(j=0;j<=N-1;++j){
		        for (k=0;k<=H-1;++k){
					if (e > 0 && prevworld[i][j][k] == 0 && world[i][j][k] != prevworld[i][j][k]){
						e -= 1;
					}
		        }
		    }
		}
		if (e > 0){
			printf("SIMULATION OVER\n");
			showworld(pipe, pipe2, t, 0);
			break;
		}
		
		//グラフに出力//
		printf("t = %d h\n",t);
		showworld(pipe, pipe2, t, 1);
		SLEEP(INTERVAL);
	}
	return 0;
}

void showworld(FILE *pipe, FILE *pipe2, int t, int state){
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
	
	fprintf(fp,"\n\n");
	
	//コラーゲンの位置を出力//
	for(i=0;i<=N-1;++i){
		for(j=0;j<=N-1;++j){
			for (k=0;k<=H-1;++k){
			  	if(world[i][j][k] == 3){ 
					fprintf(fp,"%d %d %d\n",i,j,k);
				}
			}
		}
	}

    fprintf(fp,"\n\n");
	fclose(fp);

	if((fp = fopen(TMPFILE2,"w")) == NULL){ 
		fprintf(stderr," Cannot open the file! \n");
		exit(1);
	}

    //フィブリンの位置を3段階で出力
    //フィブリンの位置を出力 1//
	for(i=0;i<=N-1;++i){
		for(j=0;j<=N-1;++j){
	  		for (k=0;k<=H-1;++k){
			  	if(CFibrin[i][j][k] > 0.75){ 
					fprintf(fp,"%d %d %d\n",i,j,k);
				}
			}
		}
	}

    fprintf(fp,"\n\n");
    
    //フィブリンの位置を出力 2//
	for(i=0;i<=N-1;++i){
		for(j=0;j<=N-1;++j){
	  		for (k=0;k<=H-1;++k){
			  	if(CFibrin[i][j][k] < 0.75 && CFibrin[i][j][k] > 0.5 ){ 
					fprintf(fp,"%d %d %d\n",i,j,k);
				}
			}
		}
	}

    fprintf(fp,"\n\n");
    
    //フィブリンの位置を出力 3//
	for(i=0;i<=N-1;++i){
		for(j=0;j<=N-1;++j){
	 		for (k=0;k<=H-1;++k){
			  	if(CFibrin[i][j][k] < 0.5 && CFibrin[i][j][k] > 0.25 ){ 
					fprintf(fp,"%d %d %d\n",i,j,k);
				}
			}
		}
	}

	fclose(fp);
	

	if (state){
		fprintf(pipe, "set title 't = %d h'\n",t);
	}
	else {
		fprintf(pipe, "set title 't = %d h OVER'\n",t);
	}
	//線維芽細胞
	fprintf(pipe, "splot \"%s\" index 0 w p ps %f pt 4 lt 5, ", TMPFILE, FONTSIZE);
	//中皮細胞
	fprintf(pipe, "\"%s\" index 1 w p ps %f pt 4 lt 2, ", TMPFILE, FONTSIZE);
	//コラーゲン
	fprintf(pipe, "\"%s\" index 2 w p ps %f pt 4 lt 3\n", TMPFILE, FONTSIZE);
	fflush(pipe);

	if (state){
		fprintf(pipe2, "set title 't = %d h'\n",t);
	}
	else {
		fprintf(pipe2, "set title 't = %d h OVER'\n",t);
	}
	//フィブリン1
	fprintf(pipe2, "splot \"%s\" index 0 w p ps %f pt 4 lt 1, ", TMPFILE2, 0.3);
	//フィブリン2
	fprintf(pipe2, "\"%s\" index 1 w p ps %f pt 4 lt 1, ", TMPFILE2, 0.2);
	//フィブリン3
	fprintf(pipe2, "\"%s\" index 2 w p ps %f pt 4 lt 1\n", TMPFILE2, 0.1);
	fflush(pipe2);
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
		for (j=0;j<=N-1;++j){
			if (genrand_real1() < 0.01){
	  			world[i][j][1] = 3;
			}
			if (genrand_real1() < 0.01){
	  			world[i][j][H - 2] = 3;
			}
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

	//フィブリンで充満
    for(i=0;i<=N-1;++i){
	    for(j=0;j<=N-1;++j){
		    for(k=0;k<=H-1;++k){
                if(world[i][j][k] == 0){
				    CFibrin[i][j][k] = 1;
                }
            }
    	}
	}
}

void spring(){
	long i,j,k;

	for (i=1; i<N; ++i){
		for (j=1; j<N; ++j){
			for (k=1; k<H; ++k){
				if (world[i][j][k] == 0 && genrand_real1() < P_spring && isViable(i,j,k,2)){
					if (genrand_real1() < 0.5) {
						world[i][j][k] = 2;
						number[i][j][k] = genrand_int32()%Tf_m + 1;
					}
					else {
						world[i][j][k] = 2;
						number[i][j][k] = genrand_int32()%Tf_m + 1;
					}				
				}
			}
		}
	}
}

void nextt(int t_count){
	int i,j,k;
	
	//細胞の更新順をランダムに
	int n;
	int a[N*N*H];
	for (n = 0; n <= N*N*H - 1; ++n) {
		a[n] = n;
	}
	randsortarray(a,N*N*H);
	//ルール適用
	for (n = 0; n <= N*N*H - 1; ++n) {
		i = a[n]/(N*H);
		j = (a[n] - i*N*H)/H;
		k = a[n] - i*N*H - j*H;
        if (world[i][j][k] == 1 || world[i][j][k] == 2) {
		    calcnext(i, j, k);
        }
	}

	for (n=0; n < FREQ; ++n) {
		//プラスミンがフィブリンのあるところに生成
    	for(i=0;i<=N-1;++i){
			for(j=0;j<=N-1;++j){
				for(k=0;k<=H-1;++k){
    	            CPlasmin[i][j][k] = CFibrin[i][j][k];
    	        }
    	    }
    	}
    	//フィブリンの拡散と線溶
    	for(i=0;i<=N-1;++i){
			for(j=0;j<=N-1;++j){
				for(k=0;k<=H-1;++k){
					if (world[i][j][k] == 0){
    	           		actFibrin(i,j,k);
					}
    	        }
    	    }
    	}
    	for(i=0;i<=N-1;++i){
			for(j=0;j<=N-1;++j){
				for(k=0;k<=H-1;++k){
    	            CFibrin[i][j][k] = NextCFibrin[i][j][k];
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
			number[0][i][H-2] = genrand_int32()%Tf_m + 1; //供給細胞の状態//
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

	//壁でフィブリンが消失//
    for(i=0;i<=N-1;++i){
	    for(j=0;j<=N-1;++j){
		    for(k=0;k<=H-1;++k){
			    if(i == 0) {
				    CFibrin[i][j][k] = 0;
                }
                if(i == N-1) {
				    CFibrin[i][j][k] = 0;
                }
                if(j == 0) {
				    CFibrin[i][j][k] = 0;
                }
                if(j == N-1) {
				    CFibrin[i][j][k] = 0;
                }
		    }
	    }
	}
}

void calcnext(int i,int j,int k){
	int action;
	int type = world[i][j][k];
	int i0, j0, k0, i1, j1, k1, i2, j2, k2;

	//確率用基準変数
	double Pbdy = 0;
	//確率係数用一時変数
	double Ptem = 0;
	//確率係数の合計
	double SumPval = 0;
	//増殖確率係数の合計
	double SumPdiv = 0;
	//遊走確率係数の合計
	double SumPmig = 0;
	//増殖確率係数の記録
	double Pdiv[3][3][3] = {0};
	//遊走確率係数の記録
	double Pmig[3][3][3] = {0};

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
				if (world[i2][j2][k2] != 0) {
					SumPval += Ptem;
				}
				//中皮細胞の場合、方向バイアスをかける
				if (type == 2) {
					Ptem *= getDvalue(i, j, k, i2, j2, k2);
				}
				if (world[i2][j2][k2] == 0 || (world[i2][j2][k2] == 3 && type == 1)) {
					//遊走確率係数の累積
					Pmig[i2 - i0][j2 - j0][k2 - k0] = Ptem;
					SumPmig += Ptem;
					if (isViable(i2, j2, k2, type)) {
						//増殖確率係数の累積
						Pdiv[i2 - i0][j2 - j0][k2 - k0] = Ptem;
						SumPdiv += Ptem;
					}
				}
			}
		}
	}

	//生存条件を満たさなければ消滅
	// if (!isSurvival(i, j, k, type)){
	if (SumPval < Survival_cond) {
		world[i][j][k] = 0;
		number[i][j][k] = 0;
		return;
	}
	if (number[i][j][k] == -1){
		return;
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
						if (Pdiv[i2 - i0][j2 - j0][k2 - k0] > 0) {
							world[i2][j2][k2] = type;
							number[i2][j2][k2] = number[i][j][k] + 1;
							if(type == 1){
								world[i][j][k] = 3;
							}
							else if (type == 2){
								world[i][j][k] = 0;
							}
							number[i][j][k] = 0;
						}
						else {
							++number[i][j][k];
						}
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
		p_div = P_f_division;                                //線維芽細胞
		Tf_tem = Tf;
	}        
	else {
		p_div = P_m_division;                                //中皮細胞
		Tf_tem = Tf_m;
	}
	//分裂周期(M)かどうか
	if(number[i][j][k] == Tf_tem) {
		if (genrand_real1() <= p_div){
			return 1;                                                //分裂
		}
		else {
			number[i][j][k] = 1;
		}
	}
	if (type == 1){
		if (genrand_real1() <= P_f_migration) {
			return 2;                                                        //遊走
		}
	}
	else {
		if (genrand_real1() <= P_m_migration) {
			return 2;                                                        //遊走
		}
	}
	return 0;                                                                //静止
}

int isViable(int i, int j, int k, int type) {
	double tem = 0;
	double val = 0;
	double val2 = 0;
	int i0 = max(0, i - 1);
	int j0 = max(0, j - 1);
	int k0 = max(0, k - 1);
	int i1 = min(N - 1, i + 1);
	int j1 = min(N - 1, j + 1);
	int k1 = min(H - 1, k + 1);
	int i2, j2, k2;
	
	//線維芽細胞
	if (type == 1){
		if (CFibrin[i][j][k] < 0.1) {
			return 0;
		}
		for (i2 = i0; i2 <= i1; ++i2) {
			for (j2 = j0; j2 <= j1; ++j2) {
				for (k2 = k0; k2 <= k1; ++k2) {
					if (world[i2][j2][k2] == 0) {
						continue;
					}
					if (world[i2][j2][k2] == 2) {
						val -= 0.5 * getPcoef(i, j, k, i2, j2, k2);
					}
					else {
						val += getPcoef(i, j, k, i2, j2, k2);
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
	else if (type == 2) {
		for (i2 = i0; i2 <= i1; ++i2) {
			for (j2 = j0; j2 <= j1; ++j2) {
				for (k2 = k0; k2 <= k1; ++k2) {
					if (world[i2][j2][k2] == 0) {
						continue;
					}
					tem = getPcoef(i, j, k, i2, j2, k2);
					val += tem;
					//線維芽細胞が隣接する
					if (world[i2][j2][k2] == 1 || world[i2][j2][k2] == 3) {
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
	else{
		return 0;
	}
}

int isSurvival(int i, int j, int k, int type) {
	double tem = 0;
	double val = 0;
	double val2 = 0;
	int i0 = max(0, i - 1);
	int j0 = max(0, j - 1);
	int k0 = max(0, k - 1);
	int i1 = min(N - 1, i + 1);
	int j1 = min(N - 1, j + 1);
	int k1 = min(H - 1, k + 1);
	int i2, j2, k2;
	
	for (i2 = i0; i2 <= i1; ++i2) {
		for (j2 = j0; j2 <= j1; ++j2) {
			for (k2 = k0; k2 <= k1; ++k2) {
				if (world[i2][j2][k2] == 0) {
					continue;
				}
				tem = getPcoef(i, j, k, i2, j2, k2);
				val += tem;
				if (world[i2][j2][k2] == 1) {
					val2 += tem;
				}
			}
		}
	}
	if (val >= Survival_cond) {
		if (type == 1) {
			return 1;
		}
		else if (val2 >= Survival_cond2) {
			return 1;
		}
	}
	return 0;
}

double getDvalue(int i0, int j0, int k0, int i1, int j1, int k1){
	double dist0;                //出発地点と中心の距離
	double dist1;                //目的地点と中心の距離
	int count;
	double p_tem = 0;
	dist0 = sqrt(pow(i0 - (N - 1)/2.0, 2.0) + pow(j0 - (N - 1)/2.0, 2.0));
	dist1 = sqrt(pow(i1 - (N - 1)/2.0, 2.0) + pow(j1 - (N - 1)/2.0, 2.0) + pow(k1 - k0, 2.0));
	p_tem = getPcoef(i0, j0, k0, i1, j1, k1);
	return pow(Dir_val, (dist0 - dist1) * p_tem);
}

void completion(int i, int j, int k) {
	int i0 = max(0, i - 1);
	int j0 = max(0, j - 1);
	int k0 = max(0, k - 1);
	int i1 = min(N - 1, i + 1);
	int j1 = min(N - 1, j + 1);
	int k1 = min(H - 1, k + 1);
	int i2, j2, k2;
	int val = 0;

	for (i2 = i0; i2 <= i1; ++i2) {
		for (j2 = j0; j2 <= j1; ++j2) {
			for (k2 = k0; k2 <= k1; ++k2) {
				if (world[i][j][k] == 2) {
					val += getPcoef(i, j, k, i2, j2, k2);
				}
			}
		}
	}
	
	if (Completion_min < val && val < Completion_max && isViable(i, j, k, 2)) {
		world[i][j][k] = 2;
		number[i][j][k] = genrand_int32()%Tf_m + 1;
		return;
	}
}

double getPcoef(int i0, int j0, int k0, int i1, int j1, int k1) {
	int count = 0;
	if(i1 == i0) {
		++count;
	}
	if(j1 == j0) {
		++count;
	}
	if(k1 == k0) {
		++count;
	}
	switch (count) {
		case 0:                                            //頂点
			return R_sqrt3;                                //確率係数 = 1/sqrt(3)
		case 1:                                            //辺心
			return R_sqrt2;                                //確率係数 = 1/sqrt(2)
		case 2:                                            //面心
			return 1.0;                                    //確率係数 = 1.0
		default:                                           //中心
			return 0;                                      //確率係数 = 0
	}
}

void randsortarray(int *a,int len) {
	int i, j, tem;
	for (i = 0; i <= len - 1; ++i) {
		j = genrand_int32()%len;
		tem = a[i];
		a[i] = a[j];
		a[j] = tem;
	}
}

void actFibrin(int i, int j, int k) {
	int i0 = max(0, i - 1);
	int j0 = max(0, j - 1);
	int k0 = max(0, k - 1);
	int i1 = min(N - 1, i + 1);
	int j1 = min(N - 1, j + 1);
	int k1 = min(H - 1, k + 1);
    NextCFibrin[i][j][k] = CFibrin[i][j][k] + COE_DIF * (- 6 * CFibrin[i][j][k] 
    + CFibrin[i0][j][k] + CFibrin[i1][j][k] + CFibrin[i][j0][k] + CFibrin[i][j1][k] 
    + CFibrin[i][j][k0] + CFibrin[i][j][k1] )
    - RATE_FIB * CFibrin[i][j][k] * CPlasmin[i][j][k];
    
    NextCFibrin[i][j][k] = min(1, NextCFibrin[i][j][k]);
    NextCFibrin[i][j][k] = max(0, NextCFibrin[i][j][k]);

}
