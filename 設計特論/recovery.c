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

#define N 101 //傷の大きさ//
#define BUFSIZE 256 //バッファサイズ//
#define INIT_INTERVAL 2 //初期待ち時間(s)//
#define INTERVAL 0.5 //待ち時間(s)//
#define TMPFILE "tempfile.tmp" //一時ファイル//
#define GNUPLOT "gnuplot" //gnuplotの場所//

void fputworld(int world[][N],int fission[][N]); //gnuplot出力//
void nextt(int world[][N],int number[][N],int t_count,int fission[][N]); //状態更新//
void calcnext(int world[][N],int nextworld[][N],int number[][N],int mig_count[][N],int fission[][N],int i,int j); //ルール適用//
void initworld(int world[][N],int number[][N]); //初期状態//
void spring(int world[][N],int number[][N]); //湧き出し//
void top_count(int world[][N], int t); //先頭細胞速度//
void sq_count(int world[][N],int t); //治癒速度//


int main(int argc,char *argv[])
{
	int t = 0; //経過時間(h)//
	int i,j,k;
	int world[N][N] = {0}; //セルの状態//
	int number[N][N] = {0}; //細胞の状態//
	int fission[N][N] = {0}; //分裂由来細胞//
	int cell = 0; //細胞の個数//
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
	initworld(world,number);
	top_count(world,t);
	sq_count(world,t);
	fputworld(world,fission);
	
	if((pipe = popen(GNUPLOT " -persist","w")) == NULL){
		fprintf(stderr," Cannot open the pipe! \n");
		exit(1);
	}
	
	//gnuplotの設定//
	fprintf(pipe,"unset key\n");
	fprintf(pipe,"set xrange [0:%d]\n",N);
	fprintf(pipe,"set yrange [0:%d]\n",N);
	fprintf(pipe, "set size square\n");
	fprintf(pipe, "unset xtics\n");
	fprintf(pipe, "unset ytics\n");
	
	//グラフに出力//
	fprintf(pipe, "set title 't = %d h'\n",t);
	fprintf(pipe, "set title font 'Arial,15'\n");
	fprintf(pipe,"plot \"" TMPFILE "\" index 0 w p ps 1 pt 4 lt 3\n");
	fprintf(pipe,"name='move%d'\n load 'savegif.gp'\n",t);
	fflush(pipe);
	sleep(INIT_INTERVAL);
	
	//細胞数をカウント//
		for(i=0;i<=N-1;++i){
	    		for(j=0;j<=N-1;++j){
	    			if(world[i][j] == 1){
	    				cell += 1;
	    			}
	    		}
	    	}
	
	//細胞数を出力//
	file = fopen("count_cell.txt","w");
	fprintf(file,"%d\n",cell);
	fclose(file);
	
	//分裂由来細胞を記録//
	fp = fopen("fission.txt","w");
	fprintf(fp,"%d\n",0);
	fclose(fp);
	
	//ステップ更新//
	for(t_count=1;t_count<=MAXT;++t_count){
		cell = 0;
		fission_total = 0;
		
		if(t_count % scale == 0){
			t += 1;  //グラフ時間表示更新//
			top_count(world,t); //治癒速度測定①//
			sq_count(world,t); //治癒速度測定②//
			if((genrand_int32()%100 + 1) <= P_spring){
				spring(world,number); //湧き出し//
			}
		}
		
		//状態更新//
		nextt(world,number,t_count,fission);
		
		//グラフに出力//
		printf("t = %d h\n",t);
		fputworld(world,fission);
		fprintf(pipe, "set title 't = %d h'\n",t);
		fprintf(pipe, "plot \"" TMPFILE "\" index 0 w p ps 1 pt 4 lt 3, \"" TMPFILE "\" index 1 w p ps 1 pt 4 lt 7\n");
		fprintf(pipe,"name='move%d'\n load 'savegif.gp'\n",t_count);
		fflush(pipe);
		sleep(INTERVAL);

		
		//細胞数をカウント//
		for(i=0;i<=N-1;++i){
	    		for(j=0;j<=N-1;++j){
	    			if(world[i][j] == 1){
	    				cell += 1;
	    			}
	    			if(fission[i][j] == 1){
	    				fission_total += 1;
	    			}
	    		}
	    	}
	    
		//細胞数を出力//
		file = fopen("count_cell.txt","a");
		fprintf(file,"%d\n",cell);
		fclose(file);
		
		//分裂由来細胞を記録//
		fp = fopen("fission.txt","a");
	    fprintf(fp,"%d\n",fission_total);
	    fclose(fp);
		
	    if(cell >= N * N) break; //完全治癒//
		
	}
	
	return 0;
}


void nextt(int world[][N],int number[][N],int t_count,int fission[][N])
{
	int i,j;
	int nextworld[N][N] = {0}; //次のセルの状態//
	int mig_count[N][N] = {0}; //同一セルの遊走状況確認//
	
	
	for(i=0;i<=N-1;++i){
		for(j=0;j<=N-1;++j){
				nextworld[i][j] = world[i][j]; //前状態の引継//
			}
		}
	//状態更新を開始する正方形の頂点を決定//
	if(t_count%4 == 1){
		for(i=0;i<=N-1;++i){
			for(j=0;j<=N-1;++j){
				if(number[i][j] > 0){
					calcnext(world,nextworld,number,mig_count,fission,i,j); //ルール適用//
				}
			}
		}
	}
	
	else if(t_count%4 == 2){
		for(i=N-1;i>=0;--i){
			for(j=0;j<=N-1;++j){
				if(number[i][j] > 0){
					calcnext(world,nextworld,number,mig_count,fission,i,j); //ルール適用//
				}
			}
		}
	}
	
	else if(t_count%4 == 3){
		for(i=0;i<=N-1;++i){
			for(j=N-1;j>=0;--j){
				if(number[i][j] > 0){
					calcnext(world,nextworld,number,mig_count,fission,i,j); //ルール適用//
				}
			}
		}
	}
	else{
		for(i=N-1;i>=0;--i){
			for(j=N-1;j>=0;--j){
				if(number[i][j] > 0){
					calcnext(world,nextworld,number,mig_count,fission,i,j); //ルール適用//
				}
			}
		}
	}
	
	
	//境界条件//
	for(i=0;i<=N-1;i++){
		if(nextworld[i][0] < 1){
			nextworld[i][0] = 1;
			number[i][0] = genrand_int32()%Tf + 1; //供給細胞の状態//
		}
		if(nextworld[i][N-1] < 1){
			nextworld[i][N-1] = 1;
			number[i][N-1] = genrand_int32()%Tf + 1; //供給細胞の状態//
		}
	}
    for(i=1;i<=N-2;i++){
    	if(nextworld[0][i] < 1){
    		nextworld[0][i] = 1;
    		number[0][i] = genrand_int32()%Tf + 1; //供給細胞の状態//
    	}
    	if(nextworld[N-1][i] < 1){
    		nextworld[N-1][i] = 1;
    		number[N-1][i] = genrand_int32()%Tf + 1; //供給細胞の状態//
    	}
    }
	
	//状態更新//
	for(i=0;i<=N-1;++i){
		for(j=0;j<=N-1;++j){
				world[i][j] = nextworld[i][j];
			}
		}
	
}


void calcnext(int world[][N],int nextworld[][N],int number[][N],int mig_count[][N],int fission[][N],int i,int j)
{
	int rule;
	int state[4] = {0};
	int k;
	int x, y;
	int a,b;
	int mig_check = 0;
	int rand[4] = {100};
	int rand_sum;
	int state_sum = 0;
	int position[4] = {0};
	
	
	//分裂・遊走を決定//
	if(number[i][j] == Tf){
		rule = 1;
	}
	else{
		rule = 2;
	}
	
	        //考慮する周囲のセル//
	/*
	            |  (i, j+1)  |
    --------------------------------------
      (i-1, j)  |  (i, j  )  |  (i+1, j)
    --------------------------------------
                |  (i, j-1)  |
                                        */
	
	//周辺状態を確認//
	if((i==0 && j==0)||(i==0 && j==N-1)||(i==N-1 && j==N-1)||(i==N-1 && j==0)){
		for(k=0;k<=3;k++){
			state[k] = 1;
		}
	}
	else if(i==0){
		state[0] = 1;
        state[1] = 1;
        state[2] = world[i+1][j];
        state[3] = 1;
	}
	else if(i==N-1){
		state[0] = world[i-1][j];
        state[1] = 1;
        state[2] = 1;
        state[3] = 1;
	}
	else if(j==0){
		state[0] = 1;
        state[1] = world[i][j+1];
        state[2] = 1;
        state[3] = 1;
	}
    else if(j==N-1){
		state[0] = 1;
        state[1] = 1;
        state[2] = 1;
        state[3] = world[i][j-1];
	}
	else{
		state[0] = world[i-1][j];
        state[1] = world[i][j+1];
        state[2] = world[i+1][j];
        state[3] = world[i][j-1];
	}
	
	for(k=0;k<=3;++k){
		state_sum += state[k];
	}
	
	//四方を囲まれたら次のステップで埋まる//
	if(state_sum==4 && i!=0 && j!=0 && i!=N-1 && j!=N-1 && world[i][j]==0){
		nextworld[i][j] = 1;
		number[i][j] = genrand_int32()%Tf + 1;
		fission[i][j] = 0;
	}
	
	//座標変換//
	a = i - (N-1)/2; 
	b = j - (N-1)/2; 
	
	//中心方向を決定//
	if(b>=a && b>=-a){
		position[3] = 1;
		position[1] = -1;
	}
	else if(b<=a && b>=-a){
		position[0] = 1;
		position[2] = -1;
	}
	else if(b<=a && b<=-a){
		position[1] = 1;
		position[3] = -1;
	}
	else{
		position[2] = 1;
		position[0] = -1;
	}

    //分裂//
    if (rule == 1){
    	if(world[i][j] == 1){
      		do{
      			rand_sum = 0;
      			
      			for(k=0;k<=3;++k){
      				if(state[k] == 0){
      					if(position[k] > 0){
      						rand[k] = genrand_int32()%100 + 1; //空きセルに対して乱数を付与//
      						if(rand[k] <= P_up){
      							rand_sum += 1; //分裂回数をカウント//
      						}
      					}
      					else if(position[k] < 0){
      						rand[k] = genrand_int32()%100 + 1; //空きセルに対して乱数を付与//
      						if(rand[k] <= P_down){
      							rand_sum += 1; //分裂回数をカウント//
      						}
      					}
      					else{
      						rand[k] = genrand_int32()%100 + 1; //空きセルに対して乱数を付与//
      						if(rand[k] <= P_side){
      							rand_sum += 1; //分裂回数をカウント//
      						}
      					}
      				}
      			}
      		}while(rand_sum > 1); //分裂回数制限//
      		
      	    for(k=0;k<=3;++k){ //分裂方向を決定//
      	    	if (k == 0){
      	    		x = i-1;
                    y = j;
      	    	}
                else if (k == 1){
                    x = i;
                    y = j+1;
                }
                else if (k == 2){
                    x = i+1;
                    y = j;
                }
                else{
                    x = i;
                    y = j-1;
                }
      	    	if(position[k] > 0){
      	    		if (state[k] == 0 && rand[k] <= P_up){ //空きセルにP_upの確率で分裂//
      	    			nextworld[i][j] = 1;
      		            nextworld[x][y] = 1;
      	    		    number[x][y] = 1;
      	    			fission[x][y] = 1;
      	    		}
      	    	}
      	    	else if(position[k] < 0){
      	    		if (state[k] == 0 && rand[k] <= P_down){ //空きセルにP_downの確率で分裂//
      	    			nextworld[i][j] = 1;
      		            nextworld[x][y] = 1;
      	    		    number[x][y] = 1;
      	    			fission[x][y] = 1;
      	    		}
      	    	}
      	    	else{
      	    		if (state[k] == 0 && rand[k] <= P_side){ //空きセルにP_sideの確率で分裂//
      	    			nextworld[i][j] = 1;
      		            nextworld[x][y] = 1;
      	    		    number[x][y] = 1;
      	    			fission[x][y] = 1;
      	    		}
      	    	}
      	    	}
      	    }
    	
    	//番号の更新//
    	if(genrand_int32()%100 + 1 <=  P_sleep){
    		number[i][j] = -1;
    	}
    	else{
    		number[i][j] = 1;
    	}
    }
    
	//遊走//
    else if (rule == 2) {
    	if(world[i][j] == 1){
      		do{
      			rand_sum = 0;
      			
      			for(k=0;k<=3;++k){
      				if(state[k] == 0){
      					if(position[k] > 0){
      						rand[k] = genrand_int32()%100 + 1; //空きセルに対して乱数を付与//
      						if(rand[k] <= P_up){
      							rand_sum += 1; //遊走回数をカウント//
      						}
      					}
      					else if(position[k] < 0){
      						rand[k] = genrand_int32()%100 + 1; //空きセルに対して乱数を付与//
      						if(rand[k] <= P_down){
      							rand_sum += 1; //遊走回数をカウント//
      						}
      					}
      					else{
      						rand[k] = genrand_int32()%100 + 1; //空きセルに対して乱数を付与//
      						if(rand[k] <= P_side){
      							rand_sum += 1; //遊走回数をカウント//
      						}
      					}
      				}
      			}
      		}while(rand_sum > 1); //遊走回数制限//
      		
      	    for(k=0;k<=3;++k){ //遊走方向を決定//
      	    	if (k == 0){
      	    		x = i-1;
                    y = j;
      	    	}
                else if (k == 1){
                    x = i;
                    y = j+1;
                }
                else if (k == 2){
                    x = i+1;
                    y = j;
                }
                else{
                    x = i;
                    y = j-1;
                }
      	    	if(position[k] > 0){
      	    		if (state[k] == 0 && rand[k] <= P_up){ //空きセルにP_upの確率で遊走//
      	    			if(mig_count[x][y] == 0){ //同一セルへの遊走を禁止//
      	    				nextworld[i][j] = 0;
      				        nextworld[x][y] = 1;
      				        mig_count[x][y] = 1;
      	    			    number[x][y] = number[i][j] + 1;
      	    			    mig_check = 1;
      	    				if(fission[i][j] == 1){ //分裂情報の引継ぎ//
      	    					fission[i][j] = 0;
      	    					fission[x][y] = 1;
      	    				}
      	    			}
      	    		}
      	    	}
      	    	else if(position[k] < 0){
      	    		if (state[k] == 0 && rand[k] <= P_down){ //空きセルにP_downの確率で遊走//
      	    			if(mig_count[x][y] == 0){ //同一セルへの遊走を禁止//
      	    				nextworld[i][j] = 0;
      				        nextworld[x][y] = 1;
      				        mig_count[x][y] = 1;
      	    			    number[x][y] = number[i][j] + 1;
      	    			    mig_check = 1;
      	    				if(fission[i][j] == 1){ //分裂情報の引継ぎ//
      	    					fission[i][j] = 0;
      	    					fission[x][y] = 1;
      	    				}
      	    		    }
      	    		}
      	    	}
      	    	else{
      	    		if (state[k] == 0 && rand[k] <= P_side){ //空きセルにP_sideの確率で遊走//
      	    			if(mig_count[x][y] == 0){ //同一セルへの遊走を禁止//
      	    				nextworld[i][j] = 0;
      				        nextworld[x][y] = 1;
      				        mig_count[x][y] = 1;
      	    			    number[x][y] = number[i][j] + 1;
      	    			    mig_check = 1;
      	    				if(fission[i][j] == 1){ //分裂情報の引継ぎ//
      	    					fission[i][j] = 0;
      	    					fission[x][y] = 1;
      	    				}
      	    		    }
      	    		}
      	    	}
      	    }
    	if(mig_check < 1){
    		number[i][j] += 1;
    	}
    }
    }
}


void fputworld(int world[][N],int fission[][N])
{
	int i,j;
	FILE *fp;

	if((fp = fopen(TMPFILE,"w")) == NULL){ 
		fprintf(stderr," Cannot open the file! \n");
		exit(1);
	}
	
	//細胞の位置を出力//
	for(i=0;i<=N-1;++i){
		for(j=0;j<=N-1;++j){
			if(world[i][j] == 1 && fission[i][j] == 0){ 
	  	fprintf(fp,"%d %d\n",i,j);
			}
		}
	}
	
	fprintf(fp,"\n\n");
	
	//細胞の位置を出力//
	for(i=0;i<=N-1;++i){
		for(j=0;j<=N-1;++j){
			if(world[i][j] == 1 && fission[i][j] == 1){ 
	  	fprintf(fp,"%d %d\n",i,j);
			}
		}
	}
	
	fclose(fp); 
	
}


void initworld(int world[][N],int number[][N])
{
  int i;
  int j;

	//初期条件//
    for(i=0;i<=N-1;i++){
        world[i][0] = 1;
    	world[i][N-1] = 1;
    }
    for(j=1;j<=N-2;j++){
        world[0][j] = 1;
        world[N-1][j] = 1;
    }
	
	for(i=0;i<=N-1;i++){
		number[i][0] = genrand_int32()%Tf + 1;
    	number[i][N-1] = genrand_int32()%Tf + 1;
    }
    for(j=1;j<=N-2;j++){
        number[0][j] = genrand_int32()%Tf + 1;
        number[N-1][j] = genrand_int32()%Tf + 1;
    }
}


void spring(int world[][N],int number[][N])
{
	int site_x,site_y;
	long i;
	
	//湧き出しの位置を決定//
	for(i=0;i<=10000000;i++){	
		site_x = genrand_int32()%(N-2) + 1;
		site_y = genrand_int32()%(N-2) + 1;
		
		if(world[site_x][site_y] == 0){ //空きセルに湧き出し//
			world[site_x][site_y] = 1;
			number[site_x][site_y] = genrand_int32()%Tf + 1;
			break;
		}
	}
}


void top_count(int world[][N], int t)
{
	int x,y;
	int i,j,s,k,l;
	double D[5] = {0};
	int e[5];
	int e_count = 0;
	FILE *fp;
	
	//e[]の初期化//
	for(i=1;i<=4;i++){
		e[i] = -1;
	}

	//中心//
	s = 0;
	x = (N-1)/2;
	y = (N-1)/2;
	if(world[x][y]==1){
		for(i=1;i<=4;i++){
			e[i] = 0;
		}
	}
	
	//中心以外//
	for(s=1;s<=(N-1)/2;s++){
		for(i=-s;i<=s;i++){
			for(j=-s;j<=s;j++){
				x = i + (N-1)/2;
				y = j + (N-1)/2;
				if(world[x][y] == 1){
					if(j>=i && j>=-i && e[1] == -1){
						e[1] = s;
					}
					else if(j<=i && j>=-i && e[2] == -1){
						e[2] = s;
					}
					else if(j<=i && j<=-i && e[3] == -1){
						e[3] = s;
					}
					else if(j>=i && j<=-i && e[4] == -1){
						e[4] = s;
					}
				}
			}
		}
	}
	
	//治癒状況を記録//
	if(t==0){
		fp = fopen("recovery_top.txt","w");
	}
	else{
		fp = fopen("recovery_top.txt","a");
	}
	
	if(t<10){
		fprintf(fp, " %d   ",t);
	}
	else if(t<100){
		fprintf(fp, " %d  ",t);
	}
	else{
		fprintf(fp, " %d ",t);
	}
	
	for(i=1;i<=4;++i){ //中心からの距離を記録//
		if(e[i]<10){
			fprintf(fp, "  %d ",e[i]);
		}
		else{
			fprintf(fp, " %d ",e[i]);
		}
		e_count += e[i];
	}
	
	if(e_count/4 < 10){ //中心からの平均距離を記録//
		fprintf(fp,"  %d \n",e_count/4);
	}
	else{
		fprintf(fp," %d \n",e_count/4);
	}
	
	fclose(fp);
	
}

void sq_count(int world[][N],int t)
{
	int e[4];
	int cell;
	int i,j;
	int e_count = 0;
	FILE *fp;
	
	//行または列の細胞数をカウント//
	for(i=1;i<=(N-1)/2;++i){
		cell = 0;
		for(j=1;j<=N-2;++j){
			if(world[i][j] == 1){
				cell += 1;
			}
		}
		if(cell < N-2){
			e[0] = (N-1)/2 - (i-1); //治癒の進行状況//
			break;
		}
		if(i == (N-1)/2 && cell == N-2){
			e[0] = 0; //中心まで治癒//
		}
	}
	
	for(j=1;j<=(N-1)/2;++j){
		cell = 0;
		for(i=1;i<=N-2;++i){
			if(world[i][j] == 1){
				cell += 1;
			}
		}
		if(cell < N-2){
			e[3] = (N-1)/2 - (j-1); //治癒の進行状況//
			break;
		}
		if(j == (N-1)/2 && cell == N-2){
			e[3] = 0; //中心まで治癒//
		}
	}

	for(i=N-2;i>=(N-1)/2;--i){
		cell = 0;
		for(j=1;j<=N-2;++j){
			if(world[i][j] == 1){
				cell += 1;
			}
		}
		if(cell < N-2){
			e[2] = (i+1) - (N-1)/2; //治癒の進行状況//
			break;
		}
		if(i == (N-1)/2 && cell == N-2){
			e[2] = 0; //中心まで治癒//
		}
	}
	
	for(j=N-2;j>=(N-1)/2;--j){
		cell = 0;
		for(i=1;i<=N-2;++i){
			if(world[i][j] == 1){
				cell += 1;
			}
		}
		if(cell < N-2){
			e[1] = (j+1) - (N-1)/2; //治癒の進行状況//
			break;
		}
		if(j == (N-1)/2 && cell == N-2){
			e[1] = 0; //中心まで治癒//
		}
	}
	
	
	//治癒状況を記録//
	if(t==0){
		fp = fopen("recovery_sq.txt","w");
	}
	else{
		fp = fopen("recovery_sq.txt","a");
	}
	
	if(t<10){
		fprintf(fp, " %d   ",t);
	}
	else if(t<100){
		fprintf(fp, " %d  ",t);
	}
	else{
		fprintf(fp, " %d ",t);
	}
	
	for(i=0;i<=3;++i){ //中心からの距離を記録//
		if(e[i]<10){
			fprintf(fp, "  %d ",e[i]);
		}
		else{
			fprintf(fp, " %d ",e[i]);
		}
		e_count += e[i];
	}
	
	if(e_count/4 < 10){ //中心からの平均距離を記録//
		fprintf(fp,"  %d \n",e_count/4);
	}
	else{
		fprintf(fp," %d \n",e_count/4);
	}
	
	fclose(fp);
	
}

