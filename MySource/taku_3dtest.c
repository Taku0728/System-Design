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

#define N 101 //傷の大きさ//
#define TMPFILE "tempfile.tmp" //一時ファイル//
#define GNUPLOT "gnuplot" //gnuplotの場所//
#define INIT_INTERVAL 2 //初期待ち時間(s)//
#define INTERVAL 0.5 //待ち時間(s)//

void fputworld(int world[][N][N],int fission[][N][N]); //gnuplot出力//
void initworld(int world[][N][N],int number[][N][N]);

int main(int argc,char *argv[]){
    int t = 0; //経過時間(h)//
	int i,j,k;
	int world[N][N][N] = {0}; //セルの状態//
	int number[N][N][N] = {0}; //細胞の状態//
	int fission[N][N][N] = {0}; //分裂由来細胞//
	FILE *pipe;

	init_genrand((unsigned)time(NULL)); //乱数の初期化//

    //初期条件//
	printf("t = 0 h\n");
	initworld(world,number);	
	fputworld(world,fission);

    if((pipe = popen(GNUPLOT " -persist","w")) == NULL){
		fprintf(stderr," Cannot open the pipe! \n");
		exit(1);
	}
	
	//gnuplotの設定//
	fprintf(pipe,"unset key\n");
	fprintf(pipe,"set xrange [0:%d]\n",N);
	fprintf(pipe,"set yrange [0:%d]\n",N);
	fprintf(pipe,"set zrange [0:%d]\n",N);
	// fprintf(pipe, "set size square\n");
	fprintf(pipe, "set ticslevel 0\n");
	// fprintf(pipe, "unset xtics\n");
	// fprintf(pipe, "unset ytics\n");
	// fprintf(pipe, "unset ztics\n");
	
	//グラフに出力//
	fprintf(pipe, "set title 't = %d h'\n",t);
	fprintf(pipe, "set title font 'Arial,15'\n");
	fprintf(pipe,"splot \"" TMPFILE "\" index 0 w p ps 1 pt 4 lt 3\n");
	fprintf(pipe,"name='move%d'\n load 'savegif.gp'\n",t);
	fflush(pipe);
	// sleep(INIT_INTERVAL);
}

void fputworld(int world[][N][N],int fission[][N][N])
{
	int i,j,k;
	FILE *fp;

	if((fp = fopen(TMPFILE,"w")) == NULL){ 
		fprintf(stderr," Cannot open the file! \n");
		exit(1);
	}
	
	//細胞の位置を出力//
	for(i=0;i<=N-1;++i){
		for(j=0;j<=N-1;++j){
            for (k=0;k<=N-1;++k){
			    if(world[i][j][k] == 1 && fission[i][j][k] == 0){ 
	  	            fprintf(fp,"%d %d %d\n",i,j,k);
                }
			}
		}
	}
	
	fprintf(fp,"\n\n");
	
	//細胞の位置を出力//
	for(i=0;i<=N-1;++i){
		for(j=0;j<=N-1;++j){
            for (k=0;k<=N-1;++k){
			    if(world[i][j][k] == 1 && fission[i][j][k] == 1){ 
	  	            fprintf(fp,"%d %d %d\n",i,j,k);
                }
			}
		}
	}
	
	fclose(fp); 
	
}

void initworld(int world[][N][N],int number[][N][N])
{
  int i;
  int j;
  int k;

	//初期条件//
    // for(j=0;j<=N-1;++j){
    //     for (k=0;k<=N-1;++k){
    //         world[0][j][k] = 1;
    //         world[N-1][j][k] = 1;
    //     }
    // }
    // for(i=0;i<=N-1;++i){
    //     for (k=0;k<=N-1;++k){
    //         world[i][0][k] = 1;
    //         world[i][N-1][k] = 1;
    //     }
    // }
    for(i=0;i<=N-1;++i){
        for (j=0;j<=N-1;++j){
            world[i][j][0] = 1;
    	    world[i][j][N-1] = 1;
        }
    }

    for(j=0;j<=N-1;++j){
        for (k=0;k<=N-1;++k){
		    number[0][j][k] = genrand_int32()%Tf + 1;
    	    number[N-1][j][k] = genrand_int32()%Tf + 1;
        }
    }
    for(i=0;i<=N-1;++i){
        for (k=0;k<=N-1;++k){
		    number[i][0][k] = genrand_int32()%Tf + 1;
    	    number[i][N-1][k] = genrand_int32()%Tf + 1;
        }
    }
    for(i=0;i<=N-1;++i){
        for (j=0;j<=N-1;++j){
		    number[i][j][0] = genrand_int32()%Tf + 1;
    	    number[i][j][N-1] = genrand_int32()%Tf + 1;
        }
    }

}