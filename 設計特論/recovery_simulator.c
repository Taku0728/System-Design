#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <math.h>
#include "MT.h" //�����Z���k�E�c�C�X�^�[//

#ifdef __unix__
#include <unistd.h>
#elif defined _WIN32
#include <windows.h>
#define sleep(x) Sleep(1000 * x) //sleep�֐��̒�`//
#endif

#define Tf 44 //����1�񂠂���̗V������(20�~�N�����~Tf)//
#define P_sleep 5 //�x���m��(%)//
#define P_up 50 //�O�i�m��(%)//
#define P_down 10 //��i�m��(%)//
#define P_side 20 //���i�m��(%)//
#define P_spring 50 //1h������̗N���o���m��(%)//

#define N 101 //���̑傫��//
#define BUFSIZE 256 //�o�b�t�@�T�C�Y//
#define INIT_INTERVAL 2 //�����҂�����(s)//
#define INTERVAL 0.5 //�҂�����(s)//
#define TMPFILE "tempfile.tmp" //�ꎞ�t�@�C��//
#define GNUPLOT "gnuplot" //gnuplot�̏ꏊ//

void fputworld(int world[][N],int fission[][N]); //gnuplot�o��//
void nextt(int world[][N],int number[][N],int t_count,int fission[][N]); //��ԍX�V//
void calcnext(int world[][N],int nextworld[][N],int number[][N],int mig_count[][N],int fission[][N],int i,int j); //���[���K�p//
void initworld(int world[][N],int number[][N]); //�������//
void spring(int world[][N],int number[][N]); //�N���o��//
void top_count(int world[][N], int t); //�擪�זE���x//
void sq_count(int world[][N],int t); //�������x//


int main(int argc,char *argv[])
{
	int t = 0; //�o�ߎ���(h)//
	int i,j,k;
	int world[N][N] = {0}; //�Z���̏��//
	int number[N][N] = {0}; //�זE�̏��//
	int fission[N][N] = {0}; //����R���זE//
	int cell = 0; //�זE�̌�//
	int t_count; //�X�e�b�v//
	int scale;
	int MAXT;
	int fission_total = 0;
	FILE *pipe;
	FILE *file;
	FILE *fp;
	
	init_genrand((unsigned)time(NULL)); //�����̏�����//
	
	scale = Tf/22; //1h������̗V������(20�~�N�����~scale)//
	MAXT = scale*168; //7���Ԃ̃V�~�����[�V����//
	
	//������//
	printf("t = 0 h\n");
	initworld(world,number);
	top_count(world,t);
	sq_count(world,t);
	fputworld(world,fission);
	
	if((pipe = popen(GNUPLOT " -persist","w")) == NULL){
		fprintf(stderr," Cannot open the pipe! \n");
		exit(1);
	}
	
	//gnuplot�̐ݒ�//
	fprintf(pipe,"unset key\n");
	fprintf(pipe,"set xrange [0:%d]\n",N);
	fprintf(pipe,"set yrange [0:%d]\n",N);
	fprintf(pipe, "set size square\n");
	fprintf(pipe, "unset xtics\n");
	fprintf(pipe, "unset ytics\n");
	fprintf(pipe, "set colorsequence classic\n");
	fprintf(pipe, "unset border\n");
	
	//�O���t�ɏo��//
	fprintf(pipe, "set title 't = %d h'\n",t);
	fprintf(pipe, "set title font 'Arial,15'\n");
	fprintf(pipe,"plot \"" TMPFILE "\" index 0 w p ps 1 pt 4 lt 3\n");
	fflush(pipe);
	sleep(INIT_INTERVAL);
	
	//�זE����J�E���g//
		for(i=0;i<=N-1;++i){
	    		for(j=0;j<=N-1;++j){
	    			if(world[i][j] == 1){
	    				cell += 1;
	    			}
	    		}
	    	}
	
	//�זE����o��//
	file = fopen("count_cell.txt","w");
	fprintf(file,"%d\n",cell);
	fclose(file);
	
	//����R���זE��L�^//
	fp = fopen("fission.txt","w");
	fprintf(fp,"%d\n",0);
	fclose(fp);
	
	//�X�e�b�v�X�V//
	for(t_count=1;t_count<=MAXT;++t_count){
		cell = 0;
		fission_total = 0;
		
		if(t_count % scale == 0){
			t += 1;  //�O���t���ԕ\���X�V//
			top_count(world,t); //�������x����@//
			sq_count(world,t); //�������x����A//
			if((genrand_int32()%100 + 1) <= P_spring){
				spring(world,number); //�N���o��//
			}
		}
		
		//��ԍX�V//
		nextt(world,number,t_count,fission);
		
		//�O���t�ɏo��//
		printf("t = %d h\n",t);
		fputworld(world,fission);
		fprintf(pipe, "set title 't = %d h'\n",t);
		fprintf(pipe, "plot \"" TMPFILE "\" index 0 w p ps 1 pt 4 lt 3, \"" TMPFILE "\" index 1 w p ps 1 pt 4 lt 1\n");
		fflush(pipe);
		sleep(INTERVAL);

		
		//�זE����J�E���g//
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
	    
		//�זE����o��//
		file = fopen("count_cell.txt","a");
		fprintf(file,"%d\n",cell);
		fclose(file);
		
		//����R���זE��L�^//
		fp = fopen("fission.txt","a");
	    fprintf(fp,"%d\n",fission_total);
	    fclose(fp);
		
	    if(cell >= N * N) break; //���S����//
		
	}
	
	return 0;
}


void nextt(int world[][N],int number[][N],int t_count,int fission[][N])
{
	int i,j;
	int nextworld[N][N] = {0}; //���̃Z���̏��//
	int mig_count[N][N] = {0}; //����Z���̗V���󋵊m�F//
	
	
	for(i=0;i<=N-1;++i){
		for(j=0;j<=N-1;++j){
				nextworld[i][j] = world[i][j]; //�O��Ԃ̈��p//
			}
		}
	//��ԍX�V��J�n���鐳���`�̒��_���//
	if(t_count%4 == 1){
		for(i=0;i<=N-1;++i){
			for(j=0;j<=N-1;++j){
				if(number[i][j] > 0){
					calcnext(world,nextworld,number,mig_count,fission,i,j); //���[���K�p//
				}
			}
		}
	}
	
	else if(t_count%4 == 2){
		for(i=N-1;i>=0;--i){
			for(j=0;j<=N-1;++j){
				if(number[i][j] > 0){
					calcnext(world,nextworld,number,mig_count,fission,i,j); //���[���K�p//
				}
			}
		}
	}
	
	else if(t_count%4 == 3){
		for(i=0;i<=N-1;++i){
			for(j=N-1;j>=0;--j){
				if(number[i][j] > 0){
					calcnext(world,nextworld,number,mig_count,fission,i,j); //���[���K�p//
				}
			}
		}
	}
	else{
		for(i=N-1;i>=0;--i){
			for(j=N-1;j>=0;--j){
				if(number[i][j] > 0){
					calcnext(world,nextworld,number,mig_count,fission,i,j); //���[���K�p//
				}
			}
		}
	}
	
	
	//���E��//
	for(i=0;i<=N-1;i++){
		if(nextworld[i][0] < 1){
			nextworld[i][0] = 1;
			number[i][0] = genrand_int32()%Tf + 1; //�����זE�̏��//
		}
		if(nextworld[i][N-1] < 1){
			nextworld[i][N-1] = 1;
			number[i][N-1] = genrand_int32()%Tf + 1; //�����זE�̏��//
		}
	}
    for(i=1;i<=N-2;i++){
    	if(nextworld[0][i] < 1){
    		nextworld[0][i] = 1;
    		number[0][i] = genrand_int32()%Tf + 1; //�����זE�̏��//
    	}
    	if(nextworld[N-1][i] < 1){
    		nextworld[N-1][i] = 1;
    		number[N-1][i] = genrand_int32()%Tf + 1; //�����זE�̏��//
    	}
    }
	
	//��ԍX�V//
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
	
	
	//����E�V�����//
	if(number[i][j] == Tf){
		rule = 1;
	}
	else{
		rule = 2;
	}
	
	        //�l��������͂̃Z��//
	/*
	            |  (i, j+1)  |
    --------------------------------------
      (i-1, j)  |  (i, j  )  |  (i+1, j)
    --------------------------------------
                |  (i, j-1)  |
                                        */
	
	//���ӏ�Ԃ�m�F//
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
	
	//�l����͂܂ꂽ�玟�̃X�e�b�v�Ŗ��܂�//
	if(state_sum==4 && i!=0 && j!=0 && i!=N-1 && j!=N-1 && world[i][j]==0){
		nextworld[i][j] = 1;
		number[i][j] = genrand_int32()%Tf + 1;
		fission[i][j] = 0;
	}
	
	//���W�ϊ�//
	a = i - (N-1)/2; 
	b = j - (N-1)/2; 
	
	//���S�������//
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

    //����//
    if (rule == 1){
    	if(world[i][j] == 1){
      		do{
      			rand_sum = 0;
      			
      			for(k=0;k<=3;++k){
      				if(state[k] == 0){
      					if(position[k] > 0){
      						rand[k] = genrand_int32()%100 + 1; //�󂫃Z���ɑ΂��ė�����t�^//
      						if(rand[k] <= P_up){
      							rand_sum += 1; //����񐔂�J�E���g//
      						}
      					}
      					else if(position[k] < 0){
      						rand[k] = genrand_int32()%100 + 1; //�󂫃Z���ɑ΂��ė�����t�^//
      						if(rand[k] <= P_down){
      							rand_sum += 1; //����񐔂�J�E���g//
      						}
      					}
      					else{
      						rand[k] = genrand_int32()%100 + 1; //�󂫃Z���ɑ΂��ė�����t�^//
      						if(rand[k] <= P_side){
      							rand_sum += 1; //����񐔂�J�E���g//
      						}
      					}
      				}
      			}
      		}while(rand_sum > 1); //����񐔐���//
      		
      	    for(k=0;k<=3;++k){ //����������//
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
      	    		if (state[k] == 0 && rand[k] <= P_up){ //�󂫃Z����P_up�̊m���ŕ���//
      	    			nextworld[i][j] = 1;
      		            nextworld[x][y] = 1;
      	    		    number[x][y] = 1;
      	    			fission[x][y] = 1;
      	    		}
      	    	}
      	    	else if(position[k] < 0){
      	    		if (state[k] == 0 && rand[k] <= P_down){ //�󂫃Z����P_down�̊m���ŕ���//
      	    			nextworld[i][j] = 1;
      		            nextworld[x][y] = 1;
      	    		    number[x][y] = 1;
      	    			fission[x][y] = 1;
      	    		}
      	    	}
      	    	else{
      	    		if (state[k] == 0 && rand[k] <= P_side){ //�󂫃Z����P_side�̊m���ŕ���//
      	    			nextworld[i][j] = 1;
      		            nextworld[x][y] = 1;
      	    		    number[x][y] = 1;
      	    			fission[x][y] = 1;
      	    		}
      	    	}
      	    	}
      	    }
    	
    	//�ԍ��̍X�V//
    	if(genrand_int32()%100 + 1 <=  P_sleep){
    		number[i][j] = -1;
    	}
    	else{
    		number[i][j] = 1;
    	}
    }
    
	//�V��//
    else if (rule == 2) {
    	if(world[i][j] == 1){
      		do{
      			rand_sum = 0;
      			
      			for(k=0;k<=3;++k){
      				if(state[k] == 0){
      					if(position[k] > 0){
      						rand[k] = genrand_int32()%100 + 1; //�󂫃Z���ɑ΂��ė�����t�^//
      						if(rand[k] <= P_up){
      							rand_sum += 1; //�V���񐔂�J�E���g//
      						}
      					}
      					else if(position[k] < 0){
      						rand[k] = genrand_int32()%100 + 1; //�󂫃Z���ɑ΂��ė�����t�^//
      						if(rand[k] <= P_down){
      							rand_sum += 1; //�V���񐔂�J�E���g//
      						}
      					}
      					else{
      						rand[k] = genrand_int32()%100 + 1; //�󂫃Z���ɑ΂��ė�����t�^//
      						if(rand[k] <= P_side){
      							rand_sum += 1; //�V���񐔂�J�E���g//
      						}
      					}
      				}
      			}
      		}while(rand_sum > 1); //�V���񐔐���//
      		
      	    for(k=0;k<=3;++k){ //�V���������//
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
      	    		if (state[k] == 0 && rand[k] <= P_up){ //�󂫃Z����P_up�̊m���ŗV��//
      	    			if(mig_count[x][y] == 0){ //����Z���ւ̗V����֎~//
      	    				nextworld[i][j] = 0;
      				        nextworld[x][y] = 1;
      				        mig_count[x][y] = 1;
      	    			    number[x][y] = number[i][j] + 1;
      	    			    mig_check = 1;
      	    				if(fission[i][j] == 1){ //������̈��p��//
      	    					fission[i][j] = 0;
      	    					fission[x][y] = 1;
      	    				}
      	    			}
      	    		}
      	    	}
      	    	else if(position[k] < 0){
      	    		if (state[k] == 0 && rand[k] <= P_down){ //�󂫃Z����P_down�̊m���ŗV��//
      	    			if(mig_count[x][y] == 0){ //����Z���ւ̗V����֎~//
      	    				nextworld[i][j] = 0;
      				        nextworld[x][y] = 1;
      				        mig_count[x][y] = 1;
      	    			    number[x][y] = number[i][j] + 1;
      	    			    mig_check = 1;
      	    				if(fission[i][j] == 1){ //������̈��p��//
      	    					fission[i][j] = 0;
      	    					fission[x][y] = 1;
      	    				}
      	    		    }
      	    		}
      	    	}
      	    	else{
      	    		if (state[k] == 0 && rand[k] <= P_side){ //�󂫃Z����P_side�̊m���ŗV��//
      	    			if(mig_count[x][y] == 0){ //����Z���ւ̗V����֎~//
      	    				nextworld[i][j] = 0;
      				        nextworld[x][y] = 1;
      				        mig_count[x][y] = 1;
      	    			    number[x][y] = number[i][j] + 1;
      	    			    mig_check = 1;
      	    				if(fission[i][j] == 1){ //������̈��p��//
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
	
	//�זE�̈ʒu��o��//
	for(i=0;i<=N-1;++i){
		for(j=0;j<=N-1;++j){
			if(world[i][j] == 1 && fission[i][j] == 0){ 
	  	fprintf(fp,"%d %d\n",i,j);
			}
		}
	}
	
	fprintf(fp,"\n\n");
	
	//�זE�̈ʒu��o��//
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

	//������//
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
	
	//�N���o���̈ʒu���//
	for(i=0;i<=10000000;i++){
		site_x = genrand_int32()%99 + 1;
		site_y = genrand_int32()%99 + 1;
		
		if(world[site_x][site_y] == 0){ //�󂫃Z���ɗN���o��//
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
	
	//e[]�̏�����//
	for(i=1;i<=4;i++){
		e[i] = -1;
	}

	//���S//
	s = 0;
	x = (N-1)/2;
	y = (N-1)/2;
	if(world[x][y]==1){
		for(i=1;i<=4;i++){
			e[i] = 0;
		}
	}
	
	//���S�ȊO//
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
	
	//�����󋵂�L�^//
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
	
	for(i=1;i<=4;++i){ //���S����̋�����L�^//
		if(e[i]<10){
			fprintf(fp, "  %d ",e[i]);
		}
		else{
			fprintf(fp, " %d ",e[i]);
		}
		e_count += e[i];
	}
	
	if(e_count/4 < 10){ //���S����̕��ϋ�����L�^//
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
	
	//�s�܂��͗�̍זE����J�E���g//
	for(i=1;i<=(N-1)/2;++i){
		cell = 0;
		for(j=1;j<=N-2;++j){
			if(world[i][j] == 1){
				cell += 1;
			}
		}
		if(cell < N-2){
			e[0] = (N-1)/2 - (i-1); //�����̐i�s��//
			break;
		}
		if(i == (N-1)/2 && cell == N-2){
			e[0] = 0; //���S�܂Ŏ���//
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
			e[3] = (N-1)/2 - (j-1); //�����̐i�s��//
			break;
		}
		if(j == (N-1)/2 && cell == N-2){
			e[3] = 0; //���S�܂Ŏ���//
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
			e[2] = (i+1) - (N-1)/2; //�����̐i�s��//
			break;
		}
		if(i == (N-1)/2 && cell == N-2){
			e[2] = 0; //���S�܂Ŏ���//
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
			e[1] = (j+1) - (N-1)/2; //�����̐i�s��//
			break;
		}
		if(j == (N-1)/2 && cell == N-2){
			e[1] = 0; //���S�܂Ŏ���//
		}
	}
	
	
	//�����󋵂�L�^//
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
	
	for(i=0;i<=3;++i){ //���S����̋�����L�^//
		if(e[i]<10){
			fprintf(fp, "  %d ",e[i]);
		}
		else{
			fprintf(fp, " %d ",e[i]);
		}
		e_count += e[i];
	}
	
	if(e_count/4 < 10){ //���S����̕��ϋ�����L�^//
		fprintf(fp,"  %d \n",e_count/4);
	}
	else{
		fprintf(fp," %d \n",e_count/4);
	}
	
	fclose(fp);
	
}

