#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <math.h>
#include "MT.h"

#ifdef __unix__
#include <unistd.h>
#elif defined _WIN32
#include <windows.h>
#define sleep(x) Sleep(1000 * x)
#endif

#define TYPES 2                             // #simulation types
#define M 100                               // the width of the world
#define N 120                               // the height of the world
#define WOUND M                             // the size of the wound
#define SHEET 10                            // the thickness of the sheet
#define P_REST 0.05                         // the probability of rest
#define P_DEG 0.01                          // the probability of degradation of a sheet
#define P_FRONT 0.50                        // the probability of action on the cell in the front of cell (i,j)
#define P_SIDE 0.20                         // the probability of action on the cell in the side of cell (i,j)
#define P_REAR 0.10                         // the probability of action on the cell in the rear of cell (i,j)
#define N_STEPS_PER_H 3                     // #steps / hour
#define RECOVERY_RATE (1*N_STEPS_PER_H)     // the period of cell recovery [hour]
#define CYCLE_PERIOD 22*N_STEPS_PER_H       // the period of cells (22 hours)
#define MAX_T RECOVERY_RATE*(WOUND/2)       // max #steps
#define INIT_INTERVAL 5                     // initial interval for plotting
#define INTERVAL 0.25                       // interval for plotting
#define CENTER_SOLID 0                      // 1: use sheet with solid center
#define PIC_FLAG 0                          // 1: save gif files (interval -> 0.5)
#define TMPFILE "tempfile.tmp"              // the name of a temporary file
#define GNUPLOT "gnuplot"

void InitializeWorld(int world[][N]);
void InitializeCellTime(int cell_time[][N]);
void InitializeCluster(int cluster[]);
void CellTimeProgress(int cell_time[][N]);
void NextT(int world[][N], int cell_time[][N], int cluster[], int t);
void SheetDegradation(int world[][N], int next_world[][N], int t);
void WoundRecovery(int next_world[][N], int next_cell_time[][N], int t);
void CalcNext(int world[][N], int next_world[][N], int cell_time[][N], int next_cell_time[][N], int cluster[], int i, int j);
void StateCheck(int world[][N], int state[], int i, int j);
int NumberOfStates(void);
void StateCount(int state[], int n_state[]);
void Fission(int world[][N], int next_world[][N], int next_cell_time[][N], int state[], int n_state[], int cluster[], int i, int j);
void Migration(int world[][N], int next_world[][N], int cell_time[][N], int next_cell_time[][N], int state[], int n_state[], int cluster[], int i, int j);
int DirectionX(int i, int k);
int DirectionY(int j, int k);
void ProbabilityGenerator(double prob[], int i, int j);
void PutWorld(int world[][N]);
void Draw(int t);
void GifDraw(int t);
void Copy(int from[][N], int to[][N]);
void PrintState(int s, int world[][N], FILE *fp);
void PlotConfig(void);
void MergeCluster(int cluster[], int i1, int j1, int i2, int j2);
int ClusterNumber(int cluster[], int i, int j);
int Dim1Index(int i, int j);
int CellTimeGenerator(void);

/*
  Global variables
*/
int type_number;
int state_size; // #states
FILE *pipe; // communication pipe to gnuplot

/*
  main function: Simulation of Lifegame, a 2-dimensional cellular automaton
*/
int main(int argc, char *argv[]){
  int t;              // time
  int adhesion_flag = -1;
  int world[M][N] = {0};  // the array of cells
  int cell_time[M][N] = {0};  // the state of cells
  int cluster[M*N] = {0};  // cluster number of cells

  if (argc < 2){     // if an initial state is not assigned, then print error message
    fprintf(stderr, "Usage: cell_simulator (Simulation type number)\ninitial state type:\n1: Adhesion (no sheet)\n2: Adhesion (with sheet)\n");
    exit(1);
  }

  type_number = atoi(argv[1]);
  state_size = NumberOfStates();

  if (type_number > TYPES || type_number <= 0) {
    fprintf(stderr, "Error!! Mistaken usage: (type number) in {1, 2}.\n");
    exit(1);
  }

  init_genrand((unsigned)time(NULL));

  InitializeWorld(world);
  InitializeCellTime(cell_time);
  InitializeCluster(cluster);

  // open the pipe
  if ((pipe = popen(GNUPLOT " -persist", "w")) == NULL){
    fprintf(stderr, "Error!! Can't open the pipe.\n");
    exit(1);
  }

  // initial setting of gnuplot
  PlotConfig();

  for (t = 0; t <= MAX_T + 1; t++){
    PutWorld(world);

    fprintf(pipe, "set title 't = %f h'\n", (double)t/(double)N_STEPS_PER_H);
    if (PIC_FLAG == 1){
      GifDraw(t);
    }
    else {
      Draw(t);
    }

    if (t == 0){
      sleep(INIT_INTERVAL);
    }
    else {
      sleep(INTERVAL);
    }

    // check whether adhesion occurred or not
    if (ClusterNumber(cluster, 0, 0) == ClusterNumber(cluster, 0, N-1)){
      adhesion_flag = 1;
      break;  // stop simulation
    }

    CellTimeProgress(cell_time);
    NextT(world, cell_time, cluster, t);
  }

  fprintf(pipe, "exit\n");
  fflush(pipe);
  pclose(pipe);

  if (adhesion_flag == 1){
    printf("Adhesion time: %f [h]\n", (double)t / (double)N_STEPS_PER_H);
  }
  else {
    printf("Adhesion did not occur.\n");
  }

  return 0;
}

/* ************************************************** Functions ************************************************** */

/*
  function InitializeWorld: load the initial state in the array
*/
void InitializeWorld(int world[][N]){
  int i, j;

  if (type_number == 1 || type_number == 2) { // Initial state is above + below
    // fibrin
    for (i = 0; i < M; i++){
      for (j = 0; j < N; j++){
        world[i][j] = 3;
      }
    }

    // fibroblast
    for (i = 0; i < M; i++){
      world[i][0] = 2;
      world[i][N-1] = 2;
    }

    // mesothelium
    for (i = 0; i < M; i++){
      if (i < (M/2 - WOUND/2) || i >= (M/2 + WOUND/2)){
        world[i][1] = 1;
        world[i][N-2] = 1;
      }
    }

    if (type_number == 2){
      // sheet
      for (i = 0; i < M; i++){
        for (j = 0; j < SHEET; j++){
          world[i][2+j] = 5;
        }
      }
    }
  }

  return;
}

/*
  function InitializeCellTime: load the initial cell time in the array
*/
void InitializeCellTime(int cell_time[][N]){
  int i, j;

  if (type_number == 1 || type_number == 2){
    // boundary (above + below)
    for(i = 0; i < M; i++){
      cell_time[i][0] = CellTimeGenerator();
      cell_time[i][N-1] = CellTimeGenerator();
    }
  }

  return;
}

/*
  function InnitializeCluster: load the initial cluster number in the array
*/
void InitializeCluster(int cluster[]){
  int i;

  if (type_number == 1 || type_number == 2){
    for (i = 0; i <= Dim1Index(M-1, N-1); i++){
      if (i > Dim1Index(0, 0) && i <= Dim1Index(M-1, 0)){
        cluster[i] = Dim1Index(0, 0);
      }
      else if (i > Dim1Index(0, N-1) && i <= Dim1Index(M-1, N-1)){
        cluster[i] = Dim1Index(0, N-1);
      }
      else {
        cluster[i] = -1;
      }
    }
  }

  return;
}

/*
  function CellTimeProgress: update all the cell time
*/
void CellTimeProgress(int cell_time[][N]){
  int i, j;

  // time progress (no progress on edges) ***************************************
  for (i = 1; i < M-1; i++){
    for (j = 1; j < N-1; j++){
      if (cell_time[i][j] != -1){
        if (cell_time[i][j] == CYCLE_PERIOD - 1){
          if (genrand_real2() < P_REST){
            cell_time[i][j] = -1; // cell stop
          }
          else{
            cell_time[i][j] = 0;
          }
        }
        else{
          cell_time[i][j]++;
        }
      }
    }
  }

  return;
}

/*
  function NextT: update all the cell states
*/
void NextT(int world[][N], int cell_time[][N], int cluster[], int t){
  int next_world[M][N] = {0}; // an array for containing the next state of cells
  int next_cell_time[M][N] = {0}; // an array for containing the next time of cells
  int i, j;

  // copy the array world / cell_time to the array next_world / next_cell_time
  Copy(world, next_world);
  Copy(cell_time, next_cell_time);

  // apply the rule
  for (i = 0; i < M; i++){
    for (j = 0; j < N; j++){
      CalcNext(world, next_world, cell_time, next_cell_time, cluster, i, j);
    }
  }

  // the rules are not applied to the cells on the boundary
  for(i = 0; i < M; i++){
    next_world[i][0] = world[i][0];
    next_world[i][N-1] = world[i][N-1];
    next_cell_time[i][0] = CellTimeGenerator();
    next_cell_time[i][N-1] = CellTimeGenerator();
  }

  if (type_number == 2){
    SheetDegradation(world, next_world, t);
  }

  if (type_number == 1 || type_number == 2){
    WoundRecovery(next_world, next_cell_time, t);
  }

  // update the array world
  Copy(next_world, world);
  Copy(next_cell_time, cell_time);

  return;
}

/*
  function SheetDestruction: degradation of the sheet
*/
void SheetDegradation(int world[][N], int next_world[][N], int t){
  int wound_size;
  int i, j;
  double C;

  wound_size = M - (2 * t / RECOVERY_RATE);

  // sheet degradation (1 row degradation)
  for (i = 0; i < M; i++){
    if (world[i][2] == 5){
      if (genrand_real2() < P_DEG){
        if (CENTER_SOLID == 1){
          if (i < (M/2 - wound_size/2) || i >= (M/2 + wound_size/2)){ // sheet over the mesothelium
            for (j = 0; j < SHEET; j++){
              next_world[i][2+j] = 0; // sheet degradation occurred
            }
          }
          else { // sheet over the wound
            C = ((double)wound_size / 2)*((double)wound_size / 2);  // normalization constant
            if (genrand_real2() < 0.5*((double)i - (double)M/2)*((double)i - (double)M/2) / C){
              for (j = 0; j < SHEET; j++){
                next_world[i][2+j] = 3; // sheet degradation occurred
              }
            }
            else {
              for (j = 0; j < SHEET; j++){
                next_world[i][2+j] = 5; // sheet degradation didn't occur
              }
            }
          }
        }
        else {
          for (j = 0; j < SHEET; j++){
            if (i < (M/2 - wound_size/2) || i >= (M/2 + wound_size/2)){
              next_world[i][2+j] = 0; // sheet degradation occurred
            }
            else {
              next_world[i][2+j] = 3; // sheet degradation occurred
            }
          }
        }
      }
    }
  }

  return;
}

/*
  function WoundRecovery: Recover the wound and fibrinolysis occurs
*/
void WoundRecovery(int next_world[][N], int next_cell_time[][N], int t){
  int wound_size;
  int i, j;

  if (t % RECOVERY_RATE == 0){
    wound_size = M - (2 * t / RECOVERY_RATE);
    for (i = 0; i < M; i++){
      if (i < (M/2 - wound_size/2) || i >= (M/2 + wound_size/2)){
        next_world[i][1] = 1;
        next_world[i][N-2] = 1;

        // fibrinolysis
        for (j = 2; j < N-2; j++){
          if (next_world[i][j] == 3){
            next_world[i][j] = 0;
          }
        }
      }
    }
  }

  return;
}

/*
  function CalcNext: update a cell state
*/
void CalcNext(int world[][N], int next_world[][N], int cell_time[][N], int next_cell_time[][N], int cluster[], int i, int j){
  int rule;     // the type of rule applied in this step
  int state[4]; // for checking the state of cells around the cell (i,j)
  int n_state[state_size];  // the value of n_state[l] represents how many cells are in the state l

  StateCheck(world, state, i, j);
  StateCount(state, n_state);

  if (cell_time[i][j] == 0){
    rule = 1; // fission
  }
  else if (cell_time[i][j] != -1) {
    rule = 2; // migration
  }

  if (type_number == 1 || type_number == 2){
    if (rule == 1){
      Fission(world, next_world, next_cell_time, state, n_state, cluster, i, j);
    }
    else if (rule == 2) {
      Migration(world, next_world, cell_time, next_cell_time, state, n_state, cluster, i, j);
    }
  }

  return;
}

/*
  function StateCheck(world, state, i, j): store the state of cells around the cell (i,j)
*/
void StateCheck(int world[][N], int state[], int i, int j){
  /*
    the position of cells:
                      |  (i, j+1)  |
          -------------------------------------
            (i-1, j)  |  (i, j  )  |  (i+1, j)
          -------------------------------------
                      |  (i, j-1)  |
  */

  if (i == 0){
    state[0] = 0;   // left
    state[2] = world[i+1][j];   // right
  }
  else if (i == M-1){
    state[0] = world[i-1][j];   // left
    state[2] = 0;   // right
  }
  else {
    state[0] = world[i-1][j];   // left
    state[2] = world[i+1][j];   // right
  }

  if (j == 0){
    state[1] = world[i][j+1];   // up
    state[3] = 0;   // down
  }
  else if (j == N-1){
    state[1] = 0;   // up
    state[3] = world[i][j-1];   // down
  }
  else{
    state[1] = world[i][j+1];   // up
    state[3] = world[i][j-1];   // down
  }

  return;
}

/*
  function NumberOfStates: returns #states
*/
int NumberOfStates(void){
  if (type_number == 1){
    return 5;           // #state = 5: State = {0,1,2,3,4}
  }

  else if (type_number == 2) {
    return 6;           // #state = 6: State = {0,1,2,3,4,5}
  }
}


/*
  function StateCount: stores #cells with state l around the cell (i,j) in n_state[l]
*/
void StateCount(int state[], int n_state[]){
  int k, l;

  // initialize the array n_state[]
  for (k = 0; k < state_size; k++){
    n_state[k] = 0;
  }

  for (k = 0; k < 4; k++){  // check the state of cells around the cell (i,j)
    for (l = 0; l < state_size; l++){
      if (state[k] == l) {  // if the state is found to be l, increase the value in n_state[l] by 1.
        n_state[l]++;
        break;
      }
    }
  }

  return;
}

void Fission(int world[][N], int next_world[][N], int next_cell_time[][N], int state[], int n_state[], int cluster[], int i, int j){
  int k;
  int x, y;
  int n_area;
  double random;
  double cumulative_p = 0;
  double prob[4] = {0};  // probability of action on cells around the cell (i,j)

  if (world[i][j] == 2){
    n_area = n_state[3] + n_state[4];
    if (n_area <= 0) {  // fission does not occur
      return;
    }

    ProbabilityGenerator(prob, i, j);

    random = genrand_real2();

    for (k = 0; k < 4; k++){
      // the cell (x, y) is the candidate for the destination of migration
      x = DirectionX(i, k);
      y = DirectionY(j, k);

      if (x == -1 || x == M || y == -1 || y == N){  // migration (i,j) -> (x,y) does not occur due to the boundary
        continue;
      }
      if (state[k] != 3 && state[k] != 4){
        continue;
      }
      cumulative_p += prob[k];
      if (random < cumulative_p){
        next_world[x][y] = 2;
        next_cell_time[x][y] = 1;
        MergeCluster(cluster, i, j, x, y);
        break;  // a trial is not needed any more once migration occurred
      }
    }
  }

  return;
}

void Migration(int world[][N], int next_world[][N], int cell_time[][N], int next_cell_time[][N], int state[], int n_state[], int cluster[], int i, int j){
  int k;
  int x, y;
  int n_area;
  double random;
  double cumulative_p = 0;
  double prob[4] = {0};  // probability of action on cells around the cell (i,j)

  if (world[i][j] == 2){
    n_area = n_state[3] + n_state[4];
    if (n_area <= 0) {  // migration does not occur
      return;
    }

    ProbabilityGenerator(prob, i, j);

    random = genrand_real2();

    for (k = 0; k < 4; k++){
      // the cell (x, y) is the candidate for the destination of migration
      x = DirectionX(i, k);
      y = DirectionY(j, k);

      if (x == -1 || x == M || y == -1 || y == N){  // migration (i,j) -> (x,y) does not occur due to the boundary
        continue;
      }
      if (state[k] != 3 && state[k] != 4){
        continue;
      }
      cumulative_p += prob[k];
      if (random < cumulative_p){
        // swap the state between the cell (i,j) and the cell (x,y)
        next_world[i][j] = 4;
        next_world[x][y] = 2;
        next_cell_time[x][y] = cell_time[i][j];
        MergeCluster(cluster, i, j, x, y);
        break;  // a trial is not needed any more once migration occurred
      }
    }
  }

  return;
}

/*
  function DirectionX(i, k): returns the x-coordinate of state[k]
*/
int DirectionX(int i, int k){
  /*
    the position of cells:
                      |  (i, j+1)  |
          -------------------------------------
            (i-1, j)  |  (i, j  )  |  (i+1, j)
          -------------------------------------
                      |  (i, j-1)  |

    state[0] = (i-1, j)
    state[1] = (i, j+1)
    state[2] = (i+1, j)
    state[3] = (i, j-1)
  */

  if (k == 0){
    return i - 1;
  }
  else if (k == 1){
    return i;
  }
  else if (k == 2){
    return i + 1;
  }
  else if (k == 3){
    return i;
  }
}

/*
  function DirectionY(j, k): returns the y-coordinate of state[k]
*/
int DirectionY(int j, int k){
  /*
    the position of cells:
                      |  (i, j+1)  |
          -------------------------------------
            (i-1, j)  |  (i, j  )  |  (i+1, j)
          -------------------------------------
                      |  (i, j-1)  |

    state[0] = (i-1, j)
    state[1] = (i, j+1)
    state[2] = (i+1, j)
    state[3] = (i, j-1)
  */

  if (k == 0){
    return j;
  }
  else if (k == 1){
    return j + 1;
  }
  else if (k == 2){
    return j;
  }
  else if (k == 3){
    return j - 1;
  }
}

/*
  function DirectionX(i, k): returns the x-coordinate of state[k]
*/
void ProbabilityGenerator(double prob[], int i, int j){
  /*
    the position of cells:
                      |  (i, j+1)  |
          -------------------------------------
            (i-1, j)  |  (i, j  )  |  (i+1, j)
          -------------------------------------
                      |  (i, j-1)  |

    prob[0] = probability of action on the cell (i-1, j)
    prob[1] = probability of action on the cell (i, j+1)
    prob[2] = probability of action on the cell (i+1, j)
    prob[3] = probability of action on the cell (i, j-1)
  */

  if (j < N / 2){
    prob[0] = P_SIDE;
    prob[1] = P_FRONT;
    prob[2] = P_SIDE;
    prob[3] = P_REAR;
  }
  else {
    prob[0] = P_SIDE;
    prob[1] = P_REAR;
    prob[2] = P_SIDE;
    prob[3] = P_FRONT;
  }

  return;
}

/*
  function PutWorld: print the current state of each cell
*/
void PutWorld(int world[][N]){
  int i, j;
  FILE *fp;  // file pointer of a temporary file

  // open a temporary file
  if ((fp = fopen(TMPFILE, "w")) == NULL){
    fprintf(stderr, "File can't open.\n");
    exit(1);
  }

  if (type_number == 1 || type_number == 2) {
    // index 0: State 1 (mesothelium)
    PrintState(1, world, fp);
    fprintf(fp, "-100 -100\n"); // dummy plot (for fixing color of cells)
    fprintf(fp, "\n\n");

    // index 1: State 2 (fibroblast)
    PrintState(2, world, fp);
    fprintf(fp, "\n\n");

    // index 2: State 3 (fibrin)
    PrintState(3, world, fp);
    fprintf(fp, "-150 -150\n"); // dummy plot (for fixing color of cells)
    fprintf(fp, "\n\n");

    // index 3: State 4 (collagen)
    PrintState(4, world, fp);
    fprintf(fp, "-200 -200\n"); // dummy plot (for fixing color of cells)
    fprintf(fp, "\n\n");

    if (type_number == 2){
      // index 4: State 5 (sheet)
      PrintState(5, world, fp);
      fprintf(fp, "-250 -250\n"); // dummy plot (for fixing color of cells)
      fprintf(fp, "\n\n");
    }
  }

  fclose(fp);
}

/*
  function Draw: draw with gnuplot
*/
void Draw(int t){
  if (type_number == 1){
    // State 1 (mesothelium) in index 0: blue
    // State 2 (fibroblast) in index 1: red
    // State 3 (fibrin) in index 2: yellow
    // State 4 (collagen) in index 3: pink
    fprintf(pipe, "plot \"" TMPFILE "\" index 0 w p ps 1 pt 5 lt 3, \"" TMPFILE "\" index 1 w p ps 1 pt 5 lt 1, \"" TMPFILE "\" index 2 w p ps 1 pt 4 lt 6, \"" TMPFILE "\" index 3 w p ps 1 pt 4 lt 4\n");
  }

  else if (type_number == 2){
    // State 1 (mesothelium) in index 0: blue
    // State 2 (fibroblast) in index 1: red
    // State 3 (fibrin) in index 2: yellow
    // State 4 (collagen) in index 3: pink
    // State 5 (sheet) in index 4: green
    fprintf(pipe, "plot \"" TMPFILE "\" index 0 w p ps 1 pt 5 lt 3, \"" TMPFILE "\" index 1 w p ps 1 pt 5 lt 1, \"" TMPFILE "\" index 2 w p ps 1 pt 4 lt 6, \"" TMPFILE "\" index 3 w p ps 1 pt 4 lt 4, \"" TMPFILE "\" index 4 w p ps 1 pt 4 lt 2\n");
  }

  fflush(pipe);
}

/*
  function GifDraw: draw with gnuplot (for exporting gif files)
*/
void GifDraw(int t){
  if (type_number == 1){
    // State 1 (mesothelium) in index 0: blue
    // State 2 (fibroblast) in index 1: red
    // State 3 (fibrin) in index 2: yellow
    // State 4 (collagen) in index 3: pink
    fprintf(pipe, "plot \"" TMPFILE "\" index 0 w p ps 1 pt 5 lt 3, \"" TMPFILE "\" index 1 w p ps 1 pt 5 lt 1, \"" TMPFILE "\" index 2 w p ps 1 pt 4 lt 18, \"" TMPFILE "\" index 3 w p ps 1 pt 4 lt 15\n");
  }

  else if (type_number == 2){
    // State 1 (mesothelium) in index 0: blue
    // State 2 (fibroblast) in index 1: red
    // State 3 (fibrin) in index 2: yellow
    // State 4 (collagen) in index 3: pink
    // State 5 (sheet) in index 4: green
    fprintf(pipe, "plot \"" TMPFILE "\" index 0 w p ps 1 pt 5 lt 3, \"" TMPFILE "\" index 1 w p ps 1 pt 5 lt 1, \"" TMPFILE "\" index 2 w p ps 1 pt 4 lt 18, \"" TMPFILE "\" index 3 w p ps 1 pt 4 lt 15, \"" TMPFILE "\" index 4 w p ps 1 pt 4 lt 2\n");
  }

  fprintf(pipe,"name='move%03d'\n load 'savegif.gp'\n", t);

  fflush(pipe);
}

void Copy(int from[][N], int to[][N]){
  int i, j;

  for (i = 0; i < M; i++){
    for (j = 0; j < N; j++){
      to[i][j] = from[i][j];
    }
  }

  return;
}

void PrintState(int s, int world[][N], FILE *fp){
  int i, j;

  for (i = 0; i < M; i++){
    for (j = 0; j < N; j++){
      if (world[i][j] == s){
        fprintf(fp, "%d %d\n", i, j);
      }
    }
  }

  return;
}

/*
  function PlotConfig: for configuration of gnuplot
*/
void PlotConfig(void){

  fprintf(pipe, "set colorsequence classic\n");
  fprintf(pipe, "set size ratio %f\n", (double)N / (double)M);  // canvas shape: height = (N/M)*width
  fprintf(pipe, "set xrange [0:%d]\n", M);
  fprintf(pipe, "set yrange [0:%d]\n", N);
  fprintf(pipe, "unset xtics\n");       // no x-axis
  fprintf(pipe, "unset ytics\n");       // no y-axis
  fprintf(pipe, "unset border\n");      // no frame
  fprintf(pipe, "unset key\n");         // no legend

  return;
}

void MergeCluster(int cluster[], int i1, int j1, int i2, int j2){
  int id1, id2;

  id1 = ClusterNumber(cluster, i1, j1);
  id2 = ClusterNumber(cluster, i2, j2);

  if (id1 != id2){
    cluster[id2] = id1;
  }

  return;
}

int ClusterNumber(int cluster[], int i, int j){
  int id;
  int k, l;

  id = Dim1Index(i, j);

  if (cluster[id] == -1){
    return id;
  }
  else {
    k = cluster[id] % M;
    l = cluster[id] / M;

    return ClusterNumber(cluster, k, l);
  }
}

/*
  function Dim1Index: returns the 1-dim. index of array[i][j]
*/
int Dim1Index(int i, int j){
  /*
    1-dim. indices of elements in a 2-dim. array are assigned as below:

      |  (N-1)M  | (N-1)M+1 |    ...   |   NM-1   |
      ---------------------------------------------
      |    :     |    :     |    ...   |    :     |
      ---------------------------------------------
      |    M     |   M+1    |    ...   |   2M-1   |
      ---------------------------------------------
      |    0     |    1     |    ...   |    M-1   |

    j (y-axis)
    ^
    |--> i  (x-axis)
  */

  return i + j*M;
}

/*
  function CellTimeGenerator: generate an integer in [0, CYCLE_PERIOD-1] uniformly at random
*/
int CellTimeGenerator(void){
  return (int)floor(CYCLE_PERIOD*genrand_real2());
}
