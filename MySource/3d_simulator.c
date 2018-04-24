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

#define Tf 44 //分裂1回あたりの遊走距離(20ミクロン×Tf)//
#define P_sleep 5 //休眠確率(%)//
#define P_up 50 //前進確率(%)//
#define P_down 10 //後進確率(%)//
#define P_side 20 //横進確率(%)//
#define P_spring 50 //1hあたりの湧き出し確率(%)//

#define P_f_division 0.2 //線維芽細胞の分裂確率
#define P_m_division 0.6 //中皮細胞の分裂確率
#define P_migration 0.4 //細胞の遊走確率

#define R_sqrt3 0.57735 //１分のルート3の高速化
#define R_sqrt2 0.70711 //１分のルート2の高速化

#define N 101 //傷の大きさ//
#define H 50 //組織の距離
#define TMPFILE "tempfile.tmp" //一時ファイル//
#define GNUPLOT "gnuplot" //gnuplotの場所//
#define INIT_INTERVAL 20 //初期待ち時間(s)//
#define INTERVAL 1 //待ち時間(s)//

int world[N][N][H] = {0}; //セルの状態//
int nextworld[N][N][H] = {0}; //次のセルの状態//
int number[N][N][H] = {0}; //細胞の状態//

void fputworld(); //gnuplot出力//
void initworld();
void nextt(int t_count); //状態更新//
void calcnext(int i,int j,int k); //ルール適用//
void spring();

int f_action(int i, int j, int k); 	//行動（分裂・遊走・静止）の決定
int on_fcell(int i, int j, int k);		//隣合う遷移細胞が存在するか
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
	MAXT = scale*1000; //7日間のシミュレーション//

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
	
	//分裂由来細胞を記録//
	fp = fopen("fission.txt","w");
	fprintf(fp,"%d\n",0);
	fclose(fp);

	// for(t_count=1;t_count<=MAXT;++t_count){
	// 	if(t_count % scale == 0){
	// 		t += 1;  //グラフ時間表示更新//
	// 		// top_count(world,t); //治癒速度測定①//
	// 		// sq_count(world,t); //治癒速度測定②//
	// 		if((genrand_int32()%100 + 1) <= P_spring){
	// 			spring(world,number); //湧き出し//
	// 		}
	// 	}
	// 	//状態更新//
	// 	nextt(t_count);

	// 	//グラフに出力//
	// 	printf("t = 0 h\n");
	// 	fputworld();
	// 	fprintf(pipe, "set title 't = %d h'\n",t);
	// 	fprintf(pipe, "set title font 'Arial,15'\n");
	// 	fprintf(pipe,"splot \"" TMPFILE "\" index 0 w p ps 0.1 pt 65 lt 5,\"" TMPFILE "\" index 1 w p ps 0.1 pt 65\n");
	// 	fprintf(pipe,"name='move%d'\n load 'savegif.gp'\n",t);
	// 	fflush(pipe);
	// 	sleep(INTERVAL);
	// }

	
	return 0;
}

void fputworld(){
	int i,j,k;
	FILE *fp;

	if((fp = fopen(TMPFILE,"w")) == NULL){ 
		fprintf(stderr," Cannot open the file! \n");
		exit(1);
	}
	
	//細胞の位置を出力//
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
	
	//細胞の位置を出力//
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
	
	//細胞の位置を出力//
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
    	    world[i][j][H-1] = 1;
        }
    }
	for(i=0;i<=N-1;++i){
		world[i][0][1] = 2;
		world[0][i][1] = 2;
		world[i][N-1][1] = 2;
		world[N-1][i][1] = 2;
		world[i][0][H-2] = 2;
		world[0][i][H-2] = 2;
		world[i][N-1][H-2] = 2;
		world[N-1][i][H-2] = 2;
	}

    for(i=0;i<=N-1;++i){
        for (j=0;j<=N-1;++j){
		    number[i][j][0] = genrand_int32()%Tf + 1;
    	    number[i][j][H-1] = genrand_int32()%Tf + 1;
        }
    }
	for(i=0;i<=N-1;++i){
		number[i][0][1] = genrand_int32()%Tf + 1;
		number[0][i][1] = genrand_int32()%Tf + 1;
		number[i][N-1][1] = genrand_int32()%Tf + 1;
		number[N-1][i][1] = genrand_int32()%Tf + 1;
		number[i][0][H-2] = genrand_int32()%Tf + 1;
		number[0][i][H-2] = genrand_int32()%Tf + 1;
		number[i][N-1][H-2] = genrand_int32()%Tf + 1;
		number[N-1][i][H-2] = genrand_int32()%Tf + 1;
	}
}

void spring(){
	int site_x,site_y;
	long i,j,k;
	int ztop[N][N];
	int zbot[N][N];
	int bort;
	
	//湧き出しの位置を決定//

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

	// for(i=0;i<=10000000;++i){
		site_x = genrand_int32()%(N-2) + 1;
		site_y = genrand_int32()%(N-2) + 1;
		bort = genrand_int32()%2;
		
		if (bort == 0){
			if(world[site_x][site_y] == 0){ //空きセルに湧き出し//
				world[site_x][site_y][zbot[site_x][site_y]] = 2;
				number[site_x][site_y][zbot[site_x][site_y]] = genrand_int32()%Tf + 1;
			}
		}
		else{
			if(world[site_x][site_y] == 0){ //空きセルに湧き出し//
				world[site_x][site_y][ztop[site_x][site_y]] = 2;
				number[site_x][site_y][ztop[site_x][site_y]] = genrand_int32()%Tf + 1;
			}
		}
	// }
}

void nextt(int t_count){
	int i,j,k;
	
	for(i=0;i<=N-1;++i){
		for(j=0;j<=N-1;++j){
			for(k=0;k<=H-1;++k){
				nextworld[i][j][k] = world[i][j][k]; //前状態の引継//
			}
		}
	}
	//状態更新を開始する正方形の頂点を決定//
	if(t_count%4 == 1){
		for(i=0;i<=N-1;++i){
			for(j=0;j<=N-1;++j){
				for (k=0;k<=H-1;++k){
					if(number[i][j][k] > 0){
						calcnext(i,j,k); //ルール適用//
					}
				}
			}
		}
	}
	else if(t_count%4 == 2){
		for(i=N-1;i>=0;--i){
			for(j=0;j<=N-1;++j){
				for (k=0;k<=H-1;++k){
					if(number[i][j][k] > 0){
						calcnext(i,j,k); //ルール適用//
					}
				}
			}
		}
	}
	else if(t_count%4 == 3){
		for(i=0;i<=N-1;++i){
			for(j=N-1;j>=0;--j){
				for (k=0;k<=H-1;++k){
					if(number[i][j][k] > 0){
						calcnext(i,j,k); //ルール適用//
					}
				}
			}
		}
	}
	else{
		for(i=N-1;i>=0;--i){
			for(j=N-1;j>=0;--j){
				for (k=0;k<=H-1;++k){
					if(number[i][j][k] > 0){
						calcnext(i,j,k); //ルール適用//
					}
				}
			}
		}
	}
	
	
	//境界条件・中皮細胞//
	for(i=0;i<=N-1;i++){
		if(nextworld[i][0][1] < 1){
			nextworld[i][0][1] = 2;
			number[i][0][1] = genrand_int32()%Tf + 1; //供給細胞の状態//
		}
		if(nextworld[i][N-1][1] < 1){
			nextworld[i][N-1][1] = 2;
			number[i][N-1][1] = genrand_int32()%Tf + 1; //供給細胞の状態//
		}
		if(nextworld[i][0][H-2] < 1){
			nextworld[i][0][H-2] = 2;
			number[i][0][H-2] = genrand_int32()%Tf + 1; //供給細胞の状態//
		}
		if(nextworld[i][N-1][H-2] < 1){
			nextworld[i][N-1][H-2] = 2;
			number[i][N-1][H-2] = genrand_int32()%Tf + 1; //供給細胞の状態//
		}
		if(nextworld[0][i][1] < 1){
			nextworld[0][i][1] = 2;
			number[0][i][1] = genrand_int32()%Tf + 1; //供給細胞の状態//
		}
		if(nextworld[N-1][i][1] < 1){
			nextworld[N-1][i][1] = 2;
			number[N-1][i][1] = genrand_int32()%Tf + 1; //供給細胞の状態//
		}
		if(nextworld[0][i][H-2] < 1){
			nextworld[0][i][H-2] = 2;
			number[0][i][N-2] = genrand_int32()%Tf + 1; //供給細胞の状態//
		}
		if(nextworld[N-1][i][H-2] < 1){
			nextworld[N-1][i][H-2] = 2;
			number[N-1][i][H-2] = genrand_int32()%Tf + 1; //供給細胞の状態//
		}
	}

	//境界条件・線維芽細胞//
	for(i=0;i<N-1;++i){
		for(j=0;j<=N-1;++j){
			if(nextworld[i][j][0] < 1){
				nextworld[i][j][0] = 1;
				number[i][j][0] = genrand_int32()%Tf + 1; //供給細胞の状態//
			}
			if(nextworld[i][j][H-1] < 1){
				nextworld[i][j][H-1] = 1;
				number[i][j][H-1] = genrand_int32()%Tf + 1; //供給細胞の状態//
			}
		}
	}
	
	//状態更新//
	for(i=0;i<=N-1;++i){
		for(j=0;j<=N-1;++j){
			for (k=0;k<=H-1;++k){
				world[i][j][k] = nextworld[i][j][k];
			}
		}
	}
	
}

void calcnext(int i,int j,int k){
	int action;
	int type = world[i][j][k];
	int i0, j0, k0, i1, j1, k1;
	int si, sj, sk;
	int di, dj, dk;
	double p = 0;			//確率用変数
	double p_tem = 0;		//確率用変数
	double p_all = 0;		//確率係数の合計
	int count = 0;			//カウント用変数
	int count1 = 0;
	int Room[3][3][3] = {};	//周辺の状態を記録
	int cond = 1;		//条件用変数
	
	//行動を決定
	action = f_action(i, j, k);
	if (action == 0){
		if(genrand_int32()%100 + 1 <=  P_sleep){
    		number[i][j][k] = -1;
    	}
    	else{
    		number[i][j][k] = 1;
    	}
		return;
	}

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
					++count1;
					continue;							//スキップ
				}
				if (type == 2) {						//中皮細胞の場合
					if(on_fcell(si,sj,sk)){	
						Room[si - i0][sj - j0][sk - k0] == 1;
					}
					else{
						cond = 0;
					}
				}
				if (cond) {
					if(si == i) ++count;
					if(sj == j) ++count;
					if(sk == k) ++count;
					switch (count) {
						case 0:							//頂点
							p_all += R_sqrt3;			//確率係数 = 1/sqrt(3)
							break;
						case 1:							//辺心
							p_all += R_sqrt2;			//確率係数 = 1/sqrt(2)
							break;
						default:						//面心
							p_all += 1.0;				//確率係数 = 1.0
							break;
					}
					count = 0;
				}
				cond = 1;
			}
		}
	}
	//近傍に細胞が存在しなければ消滅
	if (count1 == 0){
		world[i][j][k] = 0;
		number[i][j][k] = 0;
	}
	if (p_all == 0){
		if(genrand_int32()%100 + 1 <=  P_sleep){
    		number[i][j][k] = -1;
    	}
    	else{
    		number[i][j][k] = 1;
    	}
		return;
	}
	//空いているセルがない<<

	//空いているセルがある>>
	p_tem = 0;
	//分裂
	if (action == 1) {
		cond = 0;
		p = genrand_real1() * p_all;			//正規化;
		for (si = i0; si <= i1; ++si) {
			for (sj = j0; sj <= j1; ++sj) {
				for (sk = k0; sk <= k1; ++sk) {
					if (world[si][sj][sk] != 0 || cond) {
						continue;				//スキップ
					}
					//中皮細胞の増殖制約
					if (type == 2 && Room[si - i0][sj - j0][sk - k0] != 1) {
						continue;
					}
					if (si == i) ++count;
					if (sj == j) ++count;
					if (sk == k) ++count;
					switch (count) {
						case 0:					//頂点
							p_tem += R_sqrt3;	//確率係数 = 1/sqrt(3)
							break;
						case 1:					//辺心
							p_tem += R_sqrt2;	//確率係数 = 1/sqrt(2)
							break;
						case 2:					//面心
							p_tem += 1.0;		//確率係数 = 1.0
							break;
						default:				//中心
							break;				//確率係数 = 0
					}
					if (p_tem > p) {			//分裂方向決定
						di = si;
						dj = sj;
						dk = sk;
						cond = 1;
					}
					count = 0;
				}
			}
		}
	}
	//遊走//
    else if (action == 2) {
		cond = 0;
		p = genrand_real1();
		for (si = i0; si <= i1; ++si) {
			for (sj = j0; sj <= j1; ++sj) {
				for (sk = k0; sk <= k1; ++sk) {
					if (world[si][sj][sk] != 0 || cond) {
						continue;				//スキップ
					}
					//中皮細胞の増殖制約
					if (type == 2 && Room[si - i0][sj - j0][sk - k0] != 1) {
						continue;
					}
					p_tem += 1 / 26;
					if (p_tem > p) {
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
			world[i][j][k] = 0;
			number[i][j][k] = 0;
		}
    }

	//番号の更新//
	if (cond) {
		world[di][dj][dk] = type;
    	if (genrand_int32()%100 + 1 <=  P_sleep) {
    		number[di][dj][dk] = -1;
    	}
    	else{
    		number[di][dj][dk] = 1;
    	}
	}

}

int f_action(int i, int j, int k){
	int type = world[i][j][k];
	double p = genrand_real1();
	double p_div;
	if(type == 1) {
		p_div = P_f_division;				//線維芽細胞
	}	
	else {
		p_div = P_m_division;				//中皮細胞
	}				
	//分裂周期(M)かどうか
	if(number[i][j][k] == Tf && p < p_div) {
		return 1;							//分裂
	}
	p = genrand_real1();
	if(p < P_migration / (1 - p_div)){
		return 2;							//遊走
	}
	return 0;								//静止
}

int on_fcell(int i, int j, int k){
	if(world[max(0, i - 1)][j][k] == 1 ||
		world[min(N - 1, i + 1)][j][k] == 1 ||
		world[i][max(0, j - 1)][k] == 1 ||
		world[i][min(N - 1, j + 1)][k] == 1 ||
		world[i][j][max(0, k - 1)] == 1 ||
		world[i][j][min(H - 1, k + 1)] == 1)
	{
		return 1;
	}
	else {
		return 0;
	}
}