#include <math.h>
#include <complex.h>
#include <time.h>
#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <sys/stat.h>
#include <unistd.h>
#include <mpfr.h>

#include "MTGenerator.h"

/*State of optional functions:  (remember to always update text below after changing this file)
 -probDensDistFunMode - INACTIVE
 -FREQUENT UPDATES OF NEIGHBOUR LIST - INACTIVE
 -Rectangular volume moves (not romboidal) - INACTIVE
 -consfStepByStep - INACTIVE
 -tworzenie IDENTYCZNEGO rozkładu binarnego jak przy gaussie - INACTIVE
 -oddziaływanie homogeniczne=sigma a heterogeniczne=sigma+deltaDiameter - INACTIVE
 -oddziaływanie homogeniczne=sigma+deltaDiameter a heterogeniczne=sigma-deltaDiameter - INACTIVE
 -exchangeParticlesMoves - INACTIVE; type: los#1:spośród N; los#2:spośród particle[los#1].neighCounter
 -ADDED-NO NECESSARY
*/

int N,gaps,activeN,loadedConfiguration,loadType=0,loadedSetStartGenerator,loadedSetGenerator,iterationsNumber,
    growing,skipFirstIteration,saveConfigurations,
    useSpecificDirectory,useFileToIterate,fileIterateIterationsNumber=0,actIteration=0,multiplyArgument,
    onlyMath[2]={0,0},initMode,neighUpdatingFrequency;
long cyclesOfEquilibration,cyclesOfMeasurement,timeEq=0,timeMe=0,timeMath=0,intervalSampling,intervalOutput,intervalResults,savedConfigurationsInt,generatorStartPoint=0;
double maxDeltaR,desiredAcceptanceRatioR,desiredAcceptanceRatioV,
       startMinPacFrac,startMaxPacFrac,minArg,maxArg,loadedArg,
       intervalMin[10],intervalMax[10],intervalDelta[10],
       startArg[2],deltaR,deltaV=0.1, //*sigma(=1)
       sigma=1,deltaDiameter,randomStartStep[2],VcpPerParticle,
       neighRadius,neighRadius2,multiplyFactor,pressureRealOfNotFluid,
       iterationTable[1000][3],dr[3];
int alpha,beta,parM,parN,parXCells,parYCells;  //parametry sterujące
char buffer[200]="",bufferN[20],bufferGaps[20],bufferG[5],bufferD[20],bufferInitMode[5],bufferFolderIndex[5],bufferAlpha[20],bufferBeta[20],
     resultsFileName[200]="ResultsSummary.txt",
     excelResultsFileName[200]="ExcelResultsSummary.txt",
     configurationsFileName[200]="Configurations",
     configurationsListFileName[200]="ConfigurationsList.txt",
     loadConfigurationsFileName[200]="Configurations",
     probDensDistFunFileName[200]="ProbDensDistFun",              //probDensDistFunMode
     probDensDistFunResultsFileName[200]="ProbDensDistFunRes",
     loadedJOBID[50]="j-none";
/////////////////  PARTICLE functions{
typedef struct particle {
    double r[2], normR[2];  //x,y
    int neighbours[20], neighCounter;
    int type;           //type: 0-black (large), 1-white (small),  tutaj typ w sumie nie byłby potrzebny, ale jak już jest zrobiony, to zostawiam - nie będzie trzeba dodatkowych warunków na rysowanie czarno/białych dysków w mathematice
    double diameter;    //układy binarne: diameter=type?sigma-deltaDiameter:sigma+deltaDiameter; w przypadku initMode=7(random) i beta=1 średnice zadane są rozkładem Gaussa - wówczas pole 'type' przyjmuje 1(white) dla diameter<sigma i 0(black) dla diameter>sigma
} particle;


void getParticlesDistanceSquared (particle *p1, particle *p2, double boxMatrix[2][2]) {
    double normalizedRX=p1->normR[0]-p2->normR[0],
           normalizedRY=p1->normR[1]-p2->normR[1],
           rx=p1->r[0]-p2->r[0],
           ry=p1->r[1]-p2->r[1];
    rx-=round(normalizedRX)*boxMatrix[0][0]+round(normalizedRY)*boxMatrix[0][1];
    ry-=round(normalizedRX)*boxMatrix[1][0]+round(normalizedRY)*boxMatrix[1][1];
    dr[0]=rx; dr[1]=ry; dr[2]=rx*rx+ry*ry;
}

void adjustNeighRadius (double volume) {
    neighRadius=1.6*sqrt(volume)/sqrt(N);
    neighRadius2=neighRadius*neighRadius;
}

void updateNeighbourList (particle *particles, double boxMatrix[2][2], double volume) {
    adjustNeighRadius(volume);
    for (int i=0;i<activeN;i++) particles[i].neighCounter=0;
    for (int i=0;i<activeN-1;i++) for (int j=i+1;j<activeN;j++) {
        getParticlesDistanceSquared(&particles[i],&particles[j],boxMatrix);
        if (dr[2]<neighRadius2) {
            particles[i].neighbours[particles[i].neighCounter++]=j;
            particles[j].neighbours[particles[j].neighCounter++]=i;
        }
    }
}

void checkSinglePeriodicBoundaryConditions (particle *p, double boxMatrix[2][2]) {
    for (int j=0;j<2;j++) {
        for (int i=0;i<2;i++) p->r[i]-=floor(p->normR[j])*boxMatrix[i][j];
        p->normR[j]=fmod(p->normR[j],1); if (p->normR[j]<0) p->normR[j]++;
    }
}

void checkPeriodicBoundaryConditions (particle *particles, double boxMatrix[2][2]) {
    for (int i=0;i<activeN;i++)
        checkSinglePeriodicBoundaryConditions(&particles[i],boxMatrix);
}

int createRandomGaps (particle *particles, double boxMatrix[2][2], double volume) {
    printf("Creating %d random gaps... ",gaps); fflush(stdout);
    adjustNeighRadius(volume);
    int gapsIndexes[gaps], attempt=0, innerAttempt=0, allReady=0;
    do {
        attempt++;
        for (int i=0;i<gaps;i++) {
            gapsIndexes[i] = (int)(MTRandom0to1(randomStartStep)*N);
            for (int j=0;j<i;j++) {
                getParticlesDistanceSquared(&particles[gapsIndexes[i]],&particles[gapsIndexes[j]],boxMatrix);
                if (dr[2]<neighRadius2) {
                    i--;
                    innerAttempt++;
                    break;
                }
                if (j+1==i) innerAttempt=0;
            }
            if (innerAttempt>100000) break;
            if (innerAttempt==0 && i+1==gaps) allReady=1;
        }
    } while (!allReady && attempt<10000000);

    if (attempt>=10000000) {
        printf("ERROR: Couldn't create %d random gaps in %d steps.\n",gaps,attempt);
        return 0;
    } else {
        bool change; do {
            change=false;
            for (int i=gaps-1;i>0;i--) {
                if (gapsIndexes[i-1]>gapsIndexes[i]) {
                    int buffer=gapsIndexes[i];
                    gapsIndexes[i]=gapsIndexes[i-1]; gapsIndexes[i-1]=buffer;
                    change=true;
                }
            }
        } while (change);
        int actualGapIndex=0;
        for (int i=0;i<activeN;i++) {
            if (actualGapIndex<gaps) while (i+actualGapIndex==gapsIndexes[actualGapIndex]) {
                actualGapIndex++;
                if (actualGapIndex==gaps) break;
            }
            for (int j=0;j<2;j++) {
                particles[i].r[j]=particles[i+actualGapIndex].r[j];
                particles[i].normR[j]=particles[i+actualGapIndex].normR[j];
            }
        }
        printf("done\n"); fflush(stdout);
        return 1;
    }
}

double getInitDiameterBinarySystem (int type) {
    double diameter;
    switch (alpha) {
        case 0: diameter=type?sigma-deltaDiameter:sigma+deltaDiameter; break;
        case 1: diameter=type?sigma-deltaDiameter:sigma; break;
    }
    return diameter;
}

double getInitDiameterByDistributionFunction (double (*function)(double), double randomMin, double randomMax) {  //funkcja(double) MUSI być przeskalowana tak, by jej (gęstości prawdopobieństwa) MAX był równy 1
    double diameter;
    do diameter=randomMin+(randomMax-randomMin)*MTRandom0to1(randomStartStep); while (MTRandom0to1(randomStartStep)>function(diameter));
    return diameter;
}

double gaussianFunction (double x) {
    double my=sigma, delta=deltaDiameter;  //my-wartosc srednia rozkladu, delta-odchylenie standardowe
    return /* 1/delta/sqrt(2*M_PI)* */exp(-pow(x-my,2)/(2*delta*delta));  //wykomentowany człon to amplituda funkcji w punkcie MAX - w taki sposób funkcja jest skalowana do maxValue=1
}

double cosinusFunction (double x) {
    return cos((x-sigma)/deltaDiameter);
}

void invertTable (int table[50][50], int size[2]) {
    for (int i=0;i<size[0];i++) for (int j=0;j<size[1];j++)
        if (table[i][j]==0) table[i][j]=1; else table[i][j]=0;
}

int initPositions (particle *particles, double boxMatrix[2][2], double detBoxMatrix, double matrixOfParticlesSize[2], int n[2], double matrixCellXY[6][6][2], double pacFrac, double volume) {
    if (generatorStartPoint==0) {
        printf("Setting start position of p-random number generator to actual CPU time (for INIT PROCEDURES)...\n");
        InitRandomMT();
    } else {
        printf("Setting start position of p-random number generator to %ld (for INIT PROCEDURES)...\n",generatorStartPoint);
        InitMT((unsigned int)generatorStartPoint);
    }

    double mod=sqrt(N/matrixOfParticlesSize[0]/matrixOfParticlesSize[1]), interval[2][2], actualPosition[2];
    for (int i=0;i<2;i++) for (int j=0;j<2;j++) interval[i][j]=boxMatrix[i][j]/matrixOfParticlesSize[j]/mod*n[j];
    int rowCounter=0, columnCounter=0;

    //szablony struktur HD
    int pTC[50][50], typeN[2], horizontalShift;                 //typeN w obu wymiarach [0-x,1-y] musi byc <=50 (wymiar pTC, ew. mozna powiekszyc)
    for (int i=0;i<50;i++) for (int j=0;j<50;j++) pTC[i][j]=0;  //pTC(particleTypeCell) stanowi elementarna komorke typow nadawanych czastkom, horizontalShift - przesuniecie w lewo (wyrazone w liczbie czastek) komorki elementarnej typow
    switch (initMode) {
        case 0: {typeN[0]=1; typeN[1]=1; horizontalShift=0;} break;    //HD
        case 1: case 2: {typeN[0]=3; typeN[1]=2; horizontalShift=0;    //S1
            pTC[2][0]=1;
            pTC[0][1]=1;} break;
        case 3: case 4: {typeN[0]=7; typeN[1]=2; horizontalShift=2;    //S3
            pTC[1][0]=1; pTC[2][0]=1; pTC[6][0]=1;
            pTC[1][1]=1; pTC[3][1]=1; pTC[4][1]=1;} break;
        case 5: case 6: {typeN[0]=6; typeN[1]=2; horizontalShift=3;    //S6
            pTC[4][0]=1; pTC[5][0]=1;
            pTC[0][1]=1; pTC[1][1]=1; pTC[2][1]=1; pTC[4][1]=1;} break;
        case 8: case 9: {typeN[0]=6; typeN[1]=2; horizontalShift=3;    //S7
            pTC[1][0]=1; pTC[2][0]=1; pTC[4][0]=1; pTC[5][0]=1;
            pTC[0][1]=1; pTC[1][1]=1; pTC[2][1]=1;} break;
        case 10: case 11: {typeN[0]=9; typeN[1]=6; horizontalShift=0;    //S19
            pTC[1][0]=1; pTC[2][0]=1; pTC[3][0]=1; pTC[5][0]=1; pTC[6][0]=1; pTC[7][0]=1; pTC[8][0]=1;
            pTC[0][1]=1; pTC[1][1]=1; pTC[2][1]=1; pTC[3][1]=1; pTC[5][1]=1; pTC[6][1]=1; pTC[7][1]=1;
            pTC[0][2]=1; pTC[1][2]=1; pTC[2][2]=1; pTC[3][2]=1; pTC[4][2]=1;
            pTC[0][3]=1; pTC[1][3]=1; pTC[2][3]=1; pTC[3][3]=1; pTC[5][3]=1; pTC[6][3]=1; pTC[7][3]=1;
            pTC[1][4]=1; pTC[2][4]=1; pTC[3][4]=1; pTC[5][4]=1; pTC[6][4]=1; pTC[7][4]=1; pTC[8][4]=1;
            pTC[4][5]=1; pTC[5][5]=1; pTC[6][5]=1; pTC[7][5]=1; pTC[8][5]=1;} break;
        case 12: case 13: {typeN[0]=5; typeN[1]=10; horizontalShift=0;    //D7-1
            pTC[1][0]=1; pTC[2][0]=1;
            pTC[0][1]=1; pTC[1][1]=1; pTC[2][1]=1;
            pTC[1][2]=1; pTC[2][2]=1;
            pTC[3][5]=1; pTC[4][5]=1;
            pTC[0][6]=1; pTC[3][6]=1; pTC[4][6]=1;
            pTC[3][7]=1; pTC[4][7]=1;} break;
        case 14: case 15: {typeN[0]=21; typeN[1]=2; horizontalShift=12;    //D7-2
            pTC[1][0]=1; pTC[2][0]=1; pTC[13][0]=1; pTC[14][0]=1; pTC[17][0]=1; pTC[18][0]=1; pTC[19][0]=1;
            pTC[0][1]=1; pTC[1][1]=1; pTC[2][1]=1; pTC[5][1]=1; pTC[6][1]=1; pTC[17][1]=1; pTC[18][1]=1;} break;
        case 16: case 17: {typeN[0]=19; typeN[1]=2; horizontalShift=15;    //D7-3
            pTC[1][0]=1; pTC[2][0]=1; pTC[8][0]=1; pTC[9][0]=1; pTC[10][0]=1; pTC[16][0]=1; pTC[17][0]=1;
            pTC[0][1]=1; pTC[1][1]=1; pTC[2][1]=1; pTC[8][1]=1; pTC[9][1]=1; pTC[12][1]=1; pTC[13][1]=1;} break;
        case 18: case 19: {typeN[0]=5; typeN[1]=10; horizontalShift=0;    //D10  (D7-1[x3])
            pTC[1][0]=1; pTC[2][0]=1;
            pTC[0][1]=1; pTC[1][1]=1; pTC[2][1]=1;
            pTC[1][2]=1; pTC[2][2]=1; pTC[4][2]=1;
            pTC[3][3]=1; pTC[4][3]=1;
            pTC[3][5]=1; pTC[4][5]=1;
            pTC[0][6]=1; pTC[3][6]=1; pTC[4][6]=1;
            pTC[1][7]=1; pTC[3][7]=1; pTC[4][7]=1;
            pTC[1][8]=1; pTC[2][8]=1;} break;
        case 20: case 21: {typeN[0]=5; typeN[1]=10; horizontalShift=0;    //D13  (D7-1[x6])
            pTC[1][0]=1; pTC[2][0]=1; pTC[4][0]=1;
            pTC[0][1]=1; pTC[1][1]=1; pTC[2][1]=1;
            pTC[1][2]=1; pTC[2][2]=1; pTC[4][2]=1;
            pTC[3][3]=1; pTC[4][3]=1;
            pTC[1][4]=1; pTC[2][4]=1;
            pTC[1][5]=1; pTC[3][5]=1; pTC[4][5]=1;
            pTC[0][6]=1; pTC[3][6]=1; pTC[4][6]=1;
            pTC[1][7]=1; pTC[3][7]=1; pTC[4][7]=1;
            pTC[1][8]=1; pTC[2][8]=1;
            pTC[3][9]=1; pTC[4][9]=1;} break;
        case 22: case 23: {typeN[0]=21; typeN[1]=2; horizontalShift=12;    //D8  (D7-2[x3])
            pTC[1][0]=1; pTC[2][0]=1; pTC[11][0]=1; pTC[13][0]=1; pTC[14][0]=1; pTC[17][0]=1; pTC[18][0]=1; pTC[19][0]=1;
            pTC[0][1]=1; pTC[1][1]=1; pTC[2][1]=1; pTC[5][1]=1; pTC[6][1]=1; pTC[15][1]=1; pTC[17][1]=1; pTC[18][1]=1;} break;
        case 24: case 25: {typeN[0]=21; typeN[1]=2; horizontalShift=12;    //D9  (D7-2[x6])
            pTC[1][0]=1; pTC[2][0]=1; pTC[4][0]=1; pTC[11][0]=1; pTC[13][0]=1; pTC[14][0]=1; pTC[17][0]=1; pTC[18][0]=1; pTC[19][0]=1;
            pTC[0][1]=1; pTC[1][1]=1; pTC[2][1]=1; pTC[5][1]=1; pTC[6][1]=1; pTC[8][1]=1; pTC[15][1]=1; pTC[17][1]=1; pTC[18][1]=1;} break;
        case 26: case 27: {typeN[0]=7; typeN[1]=14; horizontalShift=0;    //D19-1
            pTC[1][0]=1; pTC[2][0]=1; pTC[3][0]=1;
            pTC[0][1]=1; pTC[1][1]=1; pTC[2][1]=1; pTC[3][1]=1;
            pTC[0][2]=1; pTC[1][2]=1; pTC[2][2]=1; pTC[3][2]=1; pTC[4][2]=1;
            pTC[0][3]=1; pTC[1][3]=1; pTC[2][3]=1; pTC[3][3]=1;
            pTC[1][4]=1; pTC[2][4]=1; pTC[3][4]=1;
            pTC[4][7]=1; pTC[5][7]=1; pTC[6][7]=1;
            pTC[0][8]=1; pTC[4][8]=1; pTC[5][8]=1; pTC[6][8]=1;
            pTC[0][9]=1; pTC[3][9]=1; pTC[4][9]=1; pTC[5][9]=1; pTC[6][9]=1;
            pTC[0][10]=1; pTC[4][10]=1; pTC[5][10]=1; pTC[6][10]=1;
            pTC[4][11]=1; pTC[5][11]=1; pTC[6][11]=1;} break;
        case 28: case 29: {typeN[0]=43; typeN[1]=2; horizontalShift=30;    //D19-2
            pTC[1][0]=1; pTC[2][0]=1; pTC[3][0]=1; pTC[18][0]=1; pTC[19][0]=1; pTC[20][0]=1; pTC[24][0]=1; pTC[25][0]=1; pTC[26][0]=1; pTC[27][0]=1; pTC[30][0]=1; pTC[31][0]=1; pTC[32][0]=1; pTC[33][0]=1; pTC[34][0]=1; pTC[37][0]=1; pTC[38][0]=1; pTC[39][0]=1; pTC[40][0]=1;
            pTC[0][1]=1; pTC[1][1]=1; pTC[2][1]=1; pTC[3][1]=1; pTC[7][1]=1; pTC[8][1]=1; pTC[9][1]=1; pTC[24][1]=1; pTC[25][1]=1; pTC[26][1]=1; pTC[30][1]=1; pTC[31][1]=1; pTC[32][1]=1; pTC[33][1]=1; pTC[36][1]=1; pTC[37][1]=1; pTC[38][1]=1; pTC[39][1]=1; pTC[40][1]=1;} break;
        case 30: case 31: {typeN[0]=39; typeN[1]=2; horizontalShift=33;    //D19-3
            pTC[1][0]=1; pTC[2][0]=1; pTC[3][0]=1; pTC[11][0]=1; pTC[12][0]=1; pTC[13][0]=1; pTC[14][0]=1; pTC[17][0]=1; pTC[18][0]=1; pTC[19][0]=1; pTC[20][0]=1; pTC[28][0]=1; pTC[29][0]=1; pTC[30][0]=1; pTC[33][0]=1; pTC[34][0]=1; pTC[35][0]=1; pTC[36][0]=1; pTC[37][0]=1;
            pTC[0][1]=1; pTC[1][1]=1; pTC[2][1]=1; pTC[3][1]=1; pTC[11][1]=1; pTC[12][1]=1; pTC[13][1]=1; pTC[16][1]=1; pTC[17][1]=1; pTC[18][1]=1; pTC[19][1]=1; pTC[20][1]=1; pTC[23][1]=1; pTC[24][1]=1; pTC[25][1]=1; pTC[33][1]=1; pTC[34][1]=1; pTC[35][1]=1; pTC[36][1]=1;} break;
        case 32: case 33: {typeN[0]=37; typeN[1]=2; horizontalShift=21;    //D19-4
            pTC[1][0]=1; pTC[2][0]=1; pTC[3][0]=1; pTC[6][0]=1; pTC[7][0]=1; pTC[8][0]=1; pTC[11][0]=1; pTC[12][0]=1; pTC[13][0]=1; pTC[14][0]=1; pTC[21][0]=1; pTC[22][0]=1; pTC[23][0]=1; pTC[24][0]=1; pTC[25][0]=1; pTC[32][0]=1; pTC[33][0]=1; pTC[34][0]=1; pTC[35][0]=1;
            pTC[0][1]=1; pTC[1][1]=1; pTC[2][1]=1; pTC[3][1]=1; pTC[10][1]=1; pTC[11][1]=1; pTC[12][1]=1; pTC[13][1]=1; pTC[14][1]=1; pTC[21][1]=1; pTC[22][1]=1; pTC[23][1]=1; pTC[24][1]=1; pTC[27][1]=1; pTC[28][1]=1; pTC[29][1]=1; pTC[32][1]=1; pTC[33][1]=1; pTC[34][1]=1;} break;
        case 34: case 35: {typeN[0]=7; typeN[1]=14; horizontalShift=0;    //D22-1  (D19-1[x3])
            pTC[1][0]=1; pTC[2][0]=1; pTC[3][0]=1;
            pTC[0][1]=1; pTC[1][1]=1; pTC[2][1]=1; pTC[3][1]=1;
            pTC[0][2]=1; pTC[1][2]=1; pTC[2][2]=1; pTC[3][2]=1; pTC[4][2]=1;
            pTC[0][3]=1; pTC[1][3]=1; pTC[2][3]=1; pTC[3][3]=1;
            pTC[1][4]=1; pTC[2][4]=1; pTC[3][4]=1; pTC[5][4]=1; pTC[6][4]=1;
            pTC[5][5]=1;
            pTC[4][7]=1; pTC[5][7]=1; pTC[6][7]=1;
            pTC[0][8]=1; pTC[4][8]=1; pTC[5][8]=1; pTC[6][8]=1;
            pTC[0][9]=1; pTC[3][9]=1; pTC[4][9]=1; pTC[5][9]=1; pTC[6][9]=1;
            pTC[0][10]=1; pTC[4][10]=1; pTC[5][10]=1; pTC[6][10]=1;
            pTC[1][11]=1; pTC[2][11]=1; pTC[4][11]=1; pTC[5][11]=1; pTC[6][11]=1;
            pTC[2][12]=1;} break;
        case 36: case 37: {typeN[0]=7; typeN[1]=14; horizontalShift=0;    //D25-1  (D19-1[x6])
            pTC[1][0]=1; pTC[2][0]=1; pTC[3][0]=1; pTC[5][0]=1; pTC[6][0]=1;
            pTC[0][1]=1; pTC[1][1]=1; pTC[2][1]=1; pTC[3][1]=1;
            pTC[0][2]=1; pTC[1][2]=1; pTC[2][2]=1; pTC[3][2]=1; pTC[4][2]=1;
            pTC[0][3]=1; pTC[1][3]=1; pTC[2][3]=1; pTC[3][3]=1;
            pTC[1][4]=1; pTC[2][4]=1; pTC[3][4]=1; pTC[5][4]=1; pTC[6][4]=1;
            pTC[5][5]=1;
            pTC[2][6]=1;
            pTC[1][7]=1; pTC[2][7]=1; pTC[4][7]=1; pTC[5][7]=1; pTC[6][7]=1;
            pTC[0][8]=1; pTC[4][8]=1; pTC[5][8]=1; pTC[6][8]=1;
            pTC[0][9]=1; pTC[3][9]=1; pTC[4][9]=1; pTC[5][9]=1; pTC[6][9]=1;
            pTC[0][10]=1; pTC[4][10]=1; pTC[5][10]=1; pTC[6][10]=1;
            pTC[1][11]=1; pTC[2][11]=1; pTC[4][11]=1; pTC[5][11]=1; pTC[6][11]=1;
            pTC[2][12]=1;
            pTC[5][13]=1;} break;
        case 38: case 39: {typeN[0]=7; typeN[1]=14; horizontalShift=0;    //D25-2  (D19-1[x3])
            pTC[1][0]=1; pTC[2][0]=1; pTC[3][0]=1;
            pTC[0][1]=1; pTC[1][1]=1; pTC[2][1]=1; pTC[3][1]=1;
            pTC[0][2]=1; pTC[1][2]=1; pTC[2][2]=1; pTC[3][2]=1; pTC[4][2]=1;
            pTC[0][3]=1; pTC[1][3]=1; pTC[2][3]=1; pTC[3][3]=1; pTC[5][3]=1;
            pTC[1][4]=1; pTC[2][4]=1; pTC[3][4]=1; pTC[5][4]=1; pTC[6][4]=1;
            pTC[4][5]=1; pTC[5][5]=1; pTC[6][5]=1;
            pTC[4][7]=1; pTC[5][7]=1; pTC[6][7]=1;
            pTC[0][8]=1; pTC[4][8]=1; pTC[5][8]=1; pTC[6][8]=1;
            pTC[0][9]=1; pTC[3][9]=1; pTC[4][9]=1; pTC[5][9]=1; pTC[6][9]=1;
            pTC[0][10]=1; pTC[2][10]=1; pTC[4][10]=1; pTC[5][10]=1; pTC[6][10]=1;
            pTC[1][11]=1; pTC[2][11]=1; pTC[4][11]=1; pTC[5][11]=1; pTC[6][11]=1;
            pTC[1][12]=1; pTC[2][12]=1; pTC[3][12]=1;} break;
        case 40: case 41: {typeN[0]=7; typeN[1]=14; horizontalShift=0;    //D31  (D19-1[x6])
            pTC[1][0]=1; pTC[2][0]=1; pTC[3][0]=1; pTC[5][0]=1; pTC[6][0]=1;
            pTC[0][1]=1; pTC[1][1]=1; pTC[2][1]=1; pTC[3][1]=1; pTC[5][1]=1;
            pTC[0][2]=1; pTC[1][2]=1; pTC[2][2]=1; pTC[3][2]=1; pTC[4][2]=1;
            pTC[0][3]=1; pTC[1][3]=1; pTC[2][3]=1; pTC[3][3]=1; pTC[5][3]=1;
            pTC[1][4]=1; pTC[2][4]=1; pTC[3][4]=1; pTC[5][4]=1; pTC[6][4]=1;
            pTC[4][5]=1; pTC[5][5]=1; pTC[6][5]=1;
            pTC[1][6]=1; pTC[2][6]=1; pTC[3][6]=1;
            pTC[1][7]=1; pTC[2][7]=1; pTC[4][7]=1; pTC[5][7]=1; pTC[6][7]=1;
            pTC[0][8]=1; pTC[2][8]=1; pTC[4][8]=1; pTC[5][8]=1; pTC[6][8]=1;
            pTC[0][9]=1; pTC[3][9]=1; pTC[4][9]=1; pTC[5][9]=1; pTC[6][9]=1;
            pTC[0][10]=1; pTC[2][10]=1; pTC[4][10]=1; pTC[5][10]=1; pTC[6][10]=1;
            pTC[1][11]=1; pTC[2][11]=1; pTC[4][11]=1; pTC[5][11]=1; pTC[6][11]=1;
            pTC[1][12]=1; pTC[2][12]=1; pTC[3][12]=1;
            pTC[4][13]=1; pTC[5][13]=1; pTC[6][13]=1;} break;
        case 42: case 43: {typeN[0]=43; typeN[1]=2; horizontalShift=30;    //D22-2  (D19-2[x3])
            pTC[1][0]=1; pTC[2][0]=1; pTC[3][0]=1; pTC[15][0]=1; pTC[16][0]=1; pTC[18][0]=1; pTC[19][0]=1; pTC[20][0]=1; pTC[22][0]=1; pTC[24][0]=1; pTC[25][0]=1; pTC[26][0]=1; pTC[27][0]=1; pTC[30][0]=1; pTC[31][0]=1; pTC[32][0]=1; pTC[33][0]=1; pTC[34][0]=1; pTC[37][0]=1; pTC[38][0]=1; pTC[39][0]=1; pTC[40][0]=1;
            pTC[0][1]=1; pTC[1][1]=1; pTC[2][1]=1; pTC[3][1]=1; pTC[7][1]=1; pTC[8][1]=1; pTC[9][1]=1; pTC[21][1]=1; pTC[22][1]=1; pTC[24][1]=1; pTC[25][1]=1; pTC[26][1]=1; pTC[28][1]=1; pTC[30][1]=1; pTC[31][1]=1; pTC[32][1]=1; pTC[33][1]=1; pTC[36][1]=1; pTC[37][1]=1; pTC[38][1]=1; pTC[39][1]=1; pTC[40][1]=1;} break;
        case 44: case 45: {typeN[0]=43; typeN[1]=2; horizontalShift=30;    //D25-3  (D19-2[x6])
            pTC[1][0]=1; pTC[2][0]=1; pTC[3][0]=1; pTC[5][0]=1; pTC[6][0]=1; pTC[15][0]=1; pTC[16][0]=1; pTC[18][0]=1; pTC[19][0]=1; pTC[20][0]=1; pTC[22][0]=1; pTC[24][0]=1; pTC[25][0]=1; pTC[26][0]=1; pTC[27][0]=1; pTC[30][0]=1; pTC[31][0]=1; pTC[32][0]=1; pTC[33][0]=1; pTC[34][0]=1; pTC[37][0]=1; pTC[38][0]=1; pTC[39][0]=1; pTC[40][0]=1; pTC[42][0]=1;
            pTC[0][1]=1; pTC[1][1]=1; pTC[2][1]=1; pTC[3][1]=1; pTC[5][1]=1; pTC[7][1]=1; pTC[8][1]=1; pTC[9][1]=1; pTC[11][1]=1; pTC[12][1]=1; pTC[21][1]=1; pTC[22][1]=1; pTC[24][1]=1; pTC[25][1]=1; pTC[26][1]=1; pTC[28][1]=1; pTC[30][1]=1; pTC[31][1]=1; pTC[32][1]=1; pTC[33][1]=1; pTC[36][1]=1; pTC[37][1]=1; pTC[38][1]=1; pTC[39][1]=1; pTC[40][1]=1;} break;
        case 46: case 47: {typeN[0]=39; typeN[1]=2; horizontalShift=33;    //D20  (D19-3[x3])
            pTC[1][0]=1; pTC[2][0]=1; pTC[3][0]=1; pTC[9][0]=1; pTC[11][0]=1; pTC[12][0]=1; pTC[13][0]=1; pTC[14][0]=1; pTC[17][0]=1; pTC[18][0]=1; pTC[19][0]=1; pTC[20][0]=1; pTC[28][0]=1; pTC[29][0]=1; pTC[30][0]=1; pTC[33][0]=1; pTC[34][0]=1; pTC[35][0]=1; pTC[36][0]=1; pTC[37][0]=1;
            pTC[0][1]=1; pTC[1][1]=1; pTC[2][1]=1; pTC[3][1]=1; pTC[11][1]=1; pTC[12][1]=1; pTC[13][1]=1; pTC[16][1]=1; pTC[17][1]=1; pTC[18][1]=1; pTC[19][1]=1; pTC[20][1]=1; pTC[23][1]=1; pTC[24][1]=1; pTC[25][1]=1; pTC[31][1]=1; pTC[33][1]=1; pTC[34][1]=1; pTC[35][1]=1; pTC[36][1]=1;} break;
        case 48: case 49: {typeN[0]=39; typeN[1]=2; horizontalShift=33;    //D21  (D19-3[x6])
            pTC[1][0]=1; pTC[2][0]=1; pTC[3][0]=1; pTC[9][0]=1; pTC[11][0]=1; pTC[12][0]=1; pTC[13][0]=1; pTC[14][0]=1; pTC[17][0]=1; pTC[18][0]=1; pTC[19][0]=1; pTC[20][0]=1; pTC[22][0]=1; pTC[28][0]=1; pTC[29][0]=1; pTC[30][0]=1; pTC[33][0]=1; pTC[34][0]=1; pTC[35][0]=1; pTC[36][0]=1; pTC[37][0]=1;
            pTC[0][1]=1; pTC[1][1]=1; pTC[2][1]=1; pTC[3][1]=1; pTC[5][1]=1; pTC[11][1]=1; pTC[12][1]=1; pTC[13][1]=1; pTC[16][1]=1; pTC[17][1]=1; pTC[18][1]=1; pTC[19][1]=1; pTC[20][1]=1; pTC[23][1]=1; pTC[24][1]=1; pTC[25][1]=1; pTC[31][1]=1; pTC[33][1]=1; pTC[34][1]=1; pTC[35][1]=1; pTC[36][1]=1;} break;
        case 50: case 51: {typeN[0]=12; typeN[1]=4; horizontalShift=6;    //S37
            invertTable(pTC,typeN);
            pTC[1][0]=0; pTC[6][0]=0;
            pTC[0][1]=0; pTC[6][1]=0;
            pTC[0][2]=0; pTC[7][2]=0;
            pTC[7][3]=0; pTC[8][3]=0; pTC[9][3]=0; pTC[10][3]=0; pTC[11][3]=0;} break;
        case 52: case 53: {typeN[0]=15; typeN[1]=10; horizontalShift=0;    //S61
            invertTable(pTC,typeN);
            pTC[1][0]=0; pTC[7][0]=0;
            pTC[0][1]=0; pTC[7][1]=0;
            pTC[0][2]=0; pTC[8][2]=0;
            pTC[8][3]=0; pTC[14][3]=0;
            pTC[9][4]=0; pTC[10][4]=0; pTC[11][4]=0; pTC[12][4]=0; pTC[13][4]=0; pTC[14][4]=0;
            pTC[8][5]=0; pTC[14][5]=0;
            pTC[0][6]=0; pTC[8][6]=0;
            pTC[0][7]=0; pTC[7][7]=0;
            pTC[1][8]=0; pTC[7][8]=0;
            pTC[1][9]=0; pTC[2][9]=0; pTC[3][9]=0; pTC[4][9]=0; pTC[5][9]=0; pTC[6][9]=0;} break;
        case 54: case 55: {typeN[0]=3; typeN[1]=6; horizontalShift=0;    //D1-1
            pTC[0][0]=1;
            pTC[1][3]=1;} break;
        case 56: case 57: {typeN[0]=7; typeN[1]=2; horizontalShift=2;    //D1-2
            pTC[0][0]=1;
            pTC[2][1]=1;} break;
        case 58: case 59: {typeN[0]=3; typeN[1]=6; horizontalShift=0;    //D2  (D1-1[x3])
            pTC[0][0]=1;
            pTC[1][1]=1;
            pTC[1][3]=1;
            pTC[0][4]=1;} break;
        case 60: case 61: {typeN[0]=parN+parM; typeN[1]=2*(parN-parM)-4; horizontalShift=0; invertTable(pTC,typeN);    //re-entrant
            //szkielet (jedno-dyskowy) struktury
            for (int col=0;col<parN;col++) pTC[col][0]=0;
            for (int row=1;parN-2-row>parM;row++) {
                int skeletonLeftSide=row/2;
                pTC[skeletonLeftSide][row]=0; pTC[skeletonLeftSide+(parN-1-row)][row]=0;
            }
            for (int col=0;col<=(parN-2-parM)/2;col++) pTC[col][parN-2-parM]=0;
            for (int col=(parN-2-parM)/2+parM+1;col<typeN[0];col++) pTC[col][parN-2-parM]=0;
            for (int row=parN-1-parM;row<typeN[1];row++) {
                int skeletonLeftSide=row/2-(row-(parN-2-parM));
                pTC[skeletonLeftSide][row]=0; pTC[skeletonLeftSide+(parM+1+row-(parN-2-parM))][row]=0;
            }
            } break;
        case 62: case 63: {typeN[0]=15; typeN[1]=10; horizontalShift=0;    //S61-zigzag[+1], w nawiasie jest 'poszerzenie' poziomej 'belki' (z jednej strony, wyrażone w atomach) w stosunku do S61
            invertTable(pTC,typeN);
            pTC[0][0]=0; pTC[8][0]=0;
            pTC[0][1]=0; pTC[7][1]=0;
            pTC[0][2]=0; pTC[8][2]=0;
            pTC[0][3]=0; pTC[7][3]=0;
            pTC[0][4]=0; for (int i=8;i<=14;i++) pTC[i][4]=0;
            pTC[0][5]=0; pTC[7][5]=0;
            pTC[0][6]=0; pTC[8][6]=0;
            pTC[0][7]=0; pTC[7][7]=0;
            pTC[0][8]=0; pTC[8][8]=0;
            for (int i=0;i<=7;i++) pTC[i][9]=0; } break;
        case 64: case 65: {typeN[0]=15; typeN[1]=10; horizontalShift=0;    //S61-zigzag[+2]
            invertTable(pTC,typeN);
            pTC[0][0]=0; pTC[8][0]=0;
            pTC[0][1]=0; pTC[7][1]=0;
            pTC[0][2]=0; pTC[8][2]=0;
            pTC[0][3]=0; pTC[7][3]=0;
            pTC[0][4]=0; pTC[1][4]=0; for (int i=7;i<=14;i++) pTC[i][4]=0;
            pTC[0][5]=0; pTC[7][5]=0;
            pTC[0][6]=0; pTC[8][6]=0;
            pTC[0][7]=0; pTC[7][7]=0;
            pTC[0][8]=0; pTC[8][8]=0;
            for (int i=0;i<=8;i++) pTC[i][9]=0; pTC[14][9]=0;} break;
        case 66: case 67: {typeN[0]=18; typeN[1]=12; horizontalShift=0;    //S91
            invertTable(pTC,typeN);
            pTC[3][0]=0; pTC[10][0]=0;
            pTC[2][1]=0; pTC[10][1]=0;
            pTC[2][2]=0; pTC[11][2]=0;
            pTC[1][3]=0; pTC[11][3]=0;
            pTC[1][4]=0; pTC[12][4]=0;
            pTC[0][5]=0; for (int i=12;i<=17;i++) pTC[i][5]=0;
            pTC[1][6]=0; pTC[12][6]=0;
            pTC[1][7]=0; pTC[11][7]=0;
            pTC[2][8]=0; pTC[11][8]=0;
            pTC[2][9]=0; pTC[10][9]=0;
            pTC[3][10]=0; pTC[10][10]=0;
            for (int i=3;i<=9;i++) pTC[i][11]=0;} break;
        case 68: case 69: {typeN[0]=18; typeN[1]=12; horizontalShift=0;    //S91-zigzag[+1]
            invertTable(pTC,typeN);
            pTC[2][0]=0; pTC[11][0]=0;
            pTC[1][1]=0; pTC[11][1]=0;
            pTC[2][2]=0; pTC[11][2]=0;
            pTC[2][3]=0; pTC[10][3]=0;
            pTC[2][4]=0; pTC[11][4]=0;
            for (int i=0;i<=1;i++) pTC[i][5]=0; for (int i=11;i<=17;i++) pTC[i][5]=0;
            pTC[2][6]=0; pTC[11][6]=0;
            pTC[2][7]=0; pTC[10][7]=0;
            pTC[2][8]=0; pTC[11][8]=0;
            pTC[1][9]=0; pTC[11][9]=0;
            pTC[2][10]=0; pTC[11][10]=0;
            for (int i=2;i<=10;i++) pTC[i][11]=0;} break;
        case 70: case 71: {typeN[0]=18; typeN[1]=12; horizontalShift=0;    //S91-zigzag[+2]
            invertTable(pTC,typeN);
            pTC[2][0]=0; pTC[11][0]=0;
            pTC[2][1]=0; pTC[10][1]=0;
            pTC[2][2]=0; pTC[11][2]=0;
            pTC[1][3]=0; pTC[11][3]=0;
            pTC[2][4]=0; pTC[11][4]=0;
            for (int i=0;i<=2;i++) pTC[i][5]=0; for (int i=10;i<=17;i++) pTC[i][5]=0;
            pTC[2][6]=0; pTC[11][6]=0;
            pTC[1][7]=0; pTC[11][7]=0;
            pTC[2][8]=0; pTC[11][8]=0;
            pTC[2][9]=0; pTC[10][9]=0;
            pTC[2][10]=0; pTC[11][10]=0;
            for (int i=1;i<=11;i++) pTC[i][11]=0;} break;
        case 72: case 73: {typeN[0]=15; typeN[1]=10; horizontalShift=0;    //S61-zigzag[+3]
            invertTable(pTC,typeN);
            pTC[14][0]=0; pTC[9][0]=0;
            pTC[14][1]=0; pTC[8][1]=0;
            pTC[1][2]=0; pTC[7][2]=0;
            pTC[1][3]=0; pTC[6][3]=0;
            for (int i=0;i<=2;i++) pTC[i][4]=0; for (int i=6;i<=14;i++) pTC[i][4]=0;
            pTC[1][5]=0; pTC[6][5]=0;
            pTC[1][6]=0; pTC[7][6]=0;
            pTC[14][7]=0; pTC[8][7]=0;
            pTC[14][8]=0; pTC[9][8]=0;
            for (int i=0;i<=9;i++) pTC[i][9]=0; pTC[13][9]=0; pTC[14][9]=0;} break;
    }
    if ((initMode>0 && initMode<=6 && initMode%2==0) || (initMode>7 && initMode%2==1)) invertTable(pTC,typeN);

    for (int i=0;i<N;i++) {
        int cellNumber[2][2]={{columnCounter/n[0],rowCounter/n[1]},{columnCounter%n[0],rowCounter%n[1]}}; //cellNumber[0/1][X/Y]: 0-numer komorki, 1-kolumna/rzad W komorce
        actualPosition[0]=cellNumber[0][0]*interval[0][0]+cellNumber[0][1]*interval[0][1]+matrixCellXY[cellNumber[1][0]][cellNumber[1][1]][0]*sqrt(pacFrac);
        actualPosition[1]=cellNumber[0][1]*interval[1][1]+cellNumber[0][0]*interval[1][0]+matrixCellXY[cellNumber[1][0]][cellNumber[1][1]][1]*sqrt(pacFrac);

        for (int j=0;j<2;j++) particles[i].r[j]=actualPosition[j];
        particles[i].normR[0]=(boxMatrix[1][1]*particles[i].r[0]-boxMatrix[0][1]*particles[i].r[1])/detBoxMatrix;
        particles[i].normR[1]=-(boxMatrix[1][0]*particles[i].r[0]-boxMatrix[0][0]*particles[i].r[1])/detBoxMatrix;

        if (initMode==7) {
            if (beta==0) { //układ binarny dysków BLACK/WHITE w stosunku 50% (za co odpowiada późniejsza część algorytmu)
                //particles[i].type=getInitDiameterByDistributionFunction(gaussianFunction,sigma-6*deltaDiameter,sigma+6*deltaDiameter)<sigma?1:0;  //tworzenie IDENTYCZNEGO rozkładu binarnego jak przy gaussie (mniejsze/większe od sigma) jeżeli start generatora będzie ten sam - do tego dochodzi jeszcze poniżej zakomentowanie kodu sprawdzającego udział 50/50
                if (MTRandom0to1(randomStartStep)<0.5) particles[i].type=0; else particles[i].type=1;

                /*//rzędy/kolumny czarno-białe
                if (rowCounter%2==0) particles[i].type=0; else particles[i].type=1;
                //if (columnCounter==6) if (rowCounter%2==0) particles[i].type=0; else particles[i].type=1;
                */

                particles[i].diameter=getInitDiameterBinarySystem(particles[i].type);
            } else if (beta==1) { //rozkład Gaussa średnic dysków o wartości średniej \my=sigma i odchyleniu standardowym \delta=deltaDiameter
                particles[i].diameter=getInitDiameterByDistributionFunction(gaussianFunction,sigma-6*deltaDiameter,sigma+6*deltaDiameter);
                particles[i].type=particles[i].diameter<sigma?1:0;
            } else if (beta==2) { //rozkład Cossinusowy z MAX w sigma i okresie 2*Pi*deltaDiameter
                particles[i].diameter=getInitDiameterByDistributionFunction(cosinusFunction,sigma-M_PI/2*deltaDiameter,sigma+M_PI/2*deltaDiameter);
                particles[i].type=particles[i].diameter<sigma?1:0;
            }
        } else {
            particles[i].type=pTC[(columnCounter+rowCounter/typeN[1]*horizontalShift)%typeN[0]][rowCounter%typeN[1]];
            if (!(initMode==60 || initMode==61)) particles[i].diameter=getInitDiameterBinarySystem(particles[i].type);
        }

        columnCounter++;
        if (columnCounter*1.000001>=matrixOfParticlesSize[0]*mod) {
            rowCounter++;
            columnCounter=0;
        }
    }

    if (initMode==7 && beta==0) {  //doprowadzenie układu binarnego BLACK/WHITE do udziału po 50% każdego typu
        int blackCounter=0;
        for (int i=0;i<N;i++) if (particles[i].type==0) blackCounter++;
        while (blackCounter>N/2) {
            int randIndex = (int)(MTRandom0to1(randomStartStep)*N);
            if (particles[randIndex].type==0) {particles[randIndex].type=1; particles[randIndex].diameter=getInitDiameterBinarySystem(particles[randIndex].type); blackCounter--;}
        }
        while (blackCounter<N/2) {
            int randIndex = (int)(MTRandom0to1(randomStartStep)*N);
            if (particles[randIndex].type==1) {particles[randIndex].type=0; particles[randIndex].diameter=getInitDiameterBinarySystem(particles[randIndex].type); blackCounter++;}
        }
    }

    if (initMode==60 || initMode==61) {  //pogrubienie szkieletu struktury re-entrant
        adjustNeighRadius(volume);
        int typeBuffer[N], skeletonType=initMode%2;
        for (int i=0;i<N;i++) typeBuffer[i]=particles[i].type;
        for (int betaIterator=0;betaIterator<beta;betaIterator++) {
            for (int i=0;i<N-1;i++) for (int j=i+1;j<N;j++) {
                getParticlesDistanceSquared(&particles[i],&particles[j],boxMatrix);
                if (dr[2]<neighRadius2 && particles[i].type!=particles[j].type) {
                    if (particles[i].type!=skeletonType) typeBuffer[i]=skeletonType;
                    else typeBuffer[j]=skeletonType;
                }
            }
            for (int i=0;i<N;i++) particles[i].type=typeBuffer[i];
        }
        for (int i=0;i<N;i++) particles[i].diameter=getInitDiameterBinarySystem(particles[i].type);
    }

    if (gaps>0) return createRandomGaps(particles,boxMatrix,volume);
    else return 1;
}

int getEnergyAll (particle *particles, double boxMatrix[2][2]) {
    int energy=0;
    for (int i=0;i<activeN-1;i++) for (int j=i+1;j<activeN;j++) {
        getParticlesDistanceSquared(&particles[i],&particles[j],boxMatrix);
        double drRoot=sqrt(dr[2]);
        if (drRoot<0.5*(particles[i].diameter+particles[j].diameter)) {
        //if (drRoot<(particles[i].type==particles[j].type?sigma:sigma+deltaDiameter)) {  //oddziaływanie homogeniczne=sigma a heterogeniczne=sigma+deltaDiameter 1/3
        //if (drRoot<(particles[i].type==particles[j].type?sigma+deltaDiameter:sigma-deltaDiameter)) {  //oddziaływanie homogeniczne=sigma+deltaDiameter a heterogeniczne=sigma-deltaDiameter 1/3
            energy=1;
            j=activeN; i=activeN; break;
        }
    }
    return energy;
}

int getEnergy (particle *particles, particle *dispPart, int index, double boxMatrix[2][2]) {
    int energy=0;
    for (int i=0;i<particles[index].neighCounter;i++) {
        getParticlesDistanceSquared(&particles[particles[index].neighbours[i]],dispPart,boxMatrix);
        double drRoot=sqrt(dr[2]);
        if (drRoot<0.5*(particles[index].diameter+particles[particles[index].neighbours[i]].diameter)) {
        //if (drRoot<(particles[index].type==particles[particles[index].neighbours[i]].type?sigma:sigma+deltaDiameter)) {  //oddziaływanie homogeniczne=sigma a heterogeniczne=sigma+deltaDiameter 2/3
        //if (drRoot<(particles[index].type==particles[particles[index].neighbours[i]].type?sigma+deltaDiameter:sigma-deltaDiameter)) {  //oddziaływanie homogeniczne=sigma+deltaDiameter a heterogeniczne=sigma-deltaDiameter 2/3
            energy=1;
            i=particles[index].neighCounter; break;
        }
    }
    return energy;
}

int attemptToDisplaceAParticle (particle *particles, int index, double boxMatrix[2][2], double detBoxMatrix) {
    int result=1;
    particle displacedParticle;
    for (int i=0;i<2;i++) displacedParticle.r[i]=particles[index].r[i]+(MTRandom0to1(randomStartStep)-0.5)*deltaR;
    displacedParticle.normR[0]=(boxMatrix[1][1]*displacedParticle.r[0]-boxMatrix[0][1]*displacedParticle.r[1])/detBoxMatrix;
    displacedParticle.normR[1]=-(boxMatrix[1][0]*displacedParticle.r[0]-boxMatrix[0][0]*displacedParticle.r[1])/detBoxMatrix;
    double newEnPotVicinity=getEnergy(particles,&displacedParticle,index,boxMatrix);

    if (newEnPotVicinity==0) {
        for (int i=0;i<2;i++) {
            particles[index].r[i]=displacedParticle.r[i];
            particles[index].normR[i]=displacedParticle.normR[i];
        }
        checkSinglePeriodicBoundaryConditions(&particles[index],boxMatrix);
    } else result=0;
    return result;
}

int attemptToExchangeParticles(particle *particles, int index1, int index2, double boxMatrix[2][2]) {
    int result=1;
    int bufferType=particles[index1].type; double bufferDiameter=particles[index1].diameter;
    particles[index1].type=particles[index2].type; particles[index1].diameter=particles[index2].diameter;
    particles[index2].type=bufferType; particles[index2].diameter=bufferDiameter;

    double newEnPotVicinity=getEnergy(particles,&particles[index1],index1,boxMatrix);
    if (newEnPotVicinity==0) newEnPotVicinity=getEnergy(particles,&particles[index2],index2,boxMatrix);
    if (newEnPotVicinity!=0) {
        result=0;
        particles[index2].type=particles[index1].type; particles[index2].diameter=particles[index1].diameter;
        particles[index1].type=bufferType; particles[index1].diameter=bufferDiameter;
    }
    return result;
}

void cloneParticlesForSpecificBoxMatrix (particle* clonedParticles, particle *particles, double boxMatrix[2][2]) {
    for (int i=0;i<activeN;i++) {
        for (int j=0;j<2;j++) clonedParticles[i].normR[j]=particles[i].normR[j];
        for (int j=0;j<2;j++) clonedParticles[i].r[j]=boxMatrix[j][0]*particles[i].normR[0]+boxMatrix[j][1]*particles[i].normR[1];
    }
}

int attemptToChangeVolume (particle *particles, double pressure, double boxMatrix[2][2], double *detBoxMatrix, double *volume) {
    int result=1;
    double newBoxMatrix[2][2];
    if ((*volume)/VcpPerParticle/N<pressureRealOfNotFluid) {
    //if (pressure>pressureRealOfNotFluid) {  //dozwolone zmiany ksztaltu pudla (faza stala)
        newBoxMatrix[0][0]=boxMatrix[0][0]+(MTRandom0to1(randomStartStep)-0.5)*deltaV;
        newBoxMatrix[1][1]=boxMatrix[1][1]+(MTRandom0to1(randomStartStep)-0.5)*deltaV;
        newBoxMatrix[0][1]=boxMatrix[0][1]+(MTRandom0to1(randomStartStep)-0.5)*deltaV;  //Rectangular volume moves (not romboidal) /**/ around '+(MTRandom0to1(randomStartStep)-0.5)*deltaV'
    } else {    //NIEdozwolone zmiany ksztaltu pudla (faza plynu)
        newBoxMatrix[0][0]=boxMatrix[0][0]+(MTRandom0to1(randomStartStep)-0.5)*deltaV;
        double modifier=newBoxMatrix[0][0]/boxMatrix[0][0];
        newBoxMatrix[1][1]=boxMatrix[1][1]*modifier;
        newBoxMatrix[0][1]=boxMatrix[0][1]*modifier;
    }
    newBoxMatrix[1][0]=newBoxMatrix[0][1];
    double newDetBoxMatrix=newBoxMatrix[0][0]*newBoxMatrix[1][1]-newBoxMatrix[1][0]*newBoxMatrix[0][1],
           newVolume=fabs(newDetBoxMatrix);

    particle particlesInNewBox[activeN];
    cloneParticlesForSpecificBoxMatrix(particlesInNewBox,particles,newBoxMatrix);
    for (int i=0;i<activeN;i++) for (int j=0;j<particles[i].neighCounter;j++) {
        if (i<particles[i].neighbours[j]) {
            getParticlesDistanceSquared(&particlesInNewBox[i],&particlesInNewBox[particles[i].neighbours[j]],newBoxMatrix);
            double drRoot=sqrt(dr[2]);
            if (drRoot<0.5*(particles[i].diameter+particles[particles[i].neighbours[j]].diameter)) {
            //if (drRoot<(particles[i].type==particles[particles[i].neighbours[j]].type?sigma:sigma+deltaDiameter)) {  //oddziaływanie homogeniczne=sigma a heterogeniczne=sigma+deltaDiameter 3/3
            //if (drRoot<(particles[i].type==particles[particles[i].neighbours[j]].type?sigma+deltaDiameter:sigma-deltaDiameter)) {  //oddziaływanie homogeniczne=sigma+deltaDiameter a heterogeniczne=sigma-deltaDiameter 3/3
                result=0;
                i=activeN; j=particles[i].neighCounter; break;
            }
        }
    }
    if (result) {
        double arg=-(pressure*(newVolume-(*volume))-(((double)N)*log(newVolume/(*volume))+log((newBoxMatrix[0][0]+newBoxMatrix[1][1])/(boxMatrix[0][0]+boxMatrix[1][1]))));
        if (MTRandom0to1(randomStartStep)>exp(arg)) result=0;
        else {
            *volume=newVolume;
            for (int i=0;i<2;i++) for (int j=0;j<2;j++) boxMatrix[i][j]=newBoxMatrix[i][j];
            for (int i=0;i<activeN;i++) for (int j=0;j<2;j++) particles[i].r[j]=particlesInNewBox[i].r[j];
            *detBoxMatrix=newDetBoxMatrix;
        }
    }
    return result;
}
/////////////////  } PARTICLE functions

int createIterationTable () {
    char startArguments[50]; FILE *fileStartArguments = fopen("startArguments.txt","rt");
    if (fileStartArguments==NULL) {printf("Missing file: startArguments.txt\n"); return 1;}
    while (fgets(startArguments,50,fileStartArguments)!=NULL) {
        sscanf(startArguments,"%c",startArguments); char *pEnd;
        iterationTable[fileIterateIterationsNumber][0]=strtod(startArguments,&pEnd);
        iterationTable[fileIterateIterationsNumber][1]=strtod(pEnd,&pEnd);
        iterationTable[fileIterateIterationsNumber++][2]=strtod(pEnd,NULL);
    }
    fclose(fileStartArguments);
    return 0;
}

void addAppendix (char *fileName, char *JOBID, bool jobIdOn) {
    strcpy(buffer,"2D_N-"); strncat(buffer,bufferN,20);
    strncat(buffer,"_gaps-",10); strncat(buffer,bufferGaps,20);
    strncat(buffer,"_G-",5); strncat(buffer,bufferG,5);
    strncat(buffer,"_badanie-",10); strncat(buffer,bufferFolderIndex,5);
    strncat(buffer,"_D-",5); strncat(buffer,bufferD,20);
    strncat(buffer,"_inM-",6); strncat(buffer,bufferInitMode,5);
    strncat(buffer,"_a-",5); strncat(buffer,bufferAlpha,20);
    strncat(buffer,"_b-",5); strncat(buffer,bufferBeta,20);

    mkdir(buffer,S_IRWXU);
    strncat(buffer,"/",2);
    if (jobIdOn) {
        strncat(buffer,JOBID,50);
        strncat(buffer,"_",2);
    }
    strncat(buffer,fileName,200);
    strcpy(fileName,buffer);
}

void getNextArgument (double prevArg[2], bool countIterations) {
    if (countIterations) if (--iterationsNumber==0) growing=-1;
    if (useFileToIterate) {
        if (++actIteration<fileIterateIterationsNumber) {
            for (int i=0;i<2;i++) prevArg[i]=growing?iterationTable[actIteration][i+1]:iterationTable[fileIterateIterationsNumber-1-actIteration][i+1];
            if (growing) startMinPacFrac=iterationTable[actIteration][0];
            else startMaxPacFrac=iterationTable[fileIterateIterationsNumber-1-actIteration][0];
        } else growing=-1;
    } else if (growing==1) {
        if (multiplyArgument) prevArg[0]*=multiplyFactor;
        else for (int i=0;i<10;i++)
            if (round(prevArg[0]*10000)>=round(intervalMin[i]*10000) && round(prevArg[0]*10000)<round(intervalMax[i]*10000)) {
                double newArg=round(prevArg[0]*10000)+round(intervalDelta[i]*10000);
                prevArg[0]=newArg/10000.0;
                break;
            }
        if (round(prevArg[0]*10000)>round(maxArg*10000)) growing=-1;
    } else if (growing==0) {
        if (multiplyArgument) prevArg[0]/=multiplyFactor;
        else for (int i=0;i<10;i++)
            if (round(prevArg[0]*10000)>round(intervalMin[i]*10000) && round(prevArg[0]*10000)<=round(intervalMax[i]*10000)) {
                double newArg=round(prevArg[0]*10000)-round(intervalDelta[i]*10000);
                prevArg[0]=newArg/10000.0;
                break;
            }
        if (round(prevArg[0]*10000)<round(minArg*10000)) growing=-1;
    }
}

bool isLineCorrect(char linia[4096]) {
    sscanf(linia,"%c",linia);
    int actIndex=0, dataIndex=0; while (dataIndex<3) {
        char data[50]="";
        int licznik=0, dotCounter=0;
        while (linia[actIndex]!='\t' && (dotCounter=linia[actIndex]=='.'?dotCounter+1:dotCounter)<=1 && licznik<50) data[licznik++]=linia[actIndex++];
        if (dotCounter>1 || licznik>=50) {dataIndex=10; continue;} //test of single dot in a number
        actIndex++; dataIndex++;
        if (dataIndex<10 && ((data[0]!='-' && data[1]!='.') || (data[0]=='-' && data[2]!='.'))) dataIndex=10; //test of dot position after first digit
    } if (dataIndex<10 && linia[actIndex]!='\n') dataIndex=10; //test of max 3 numbers in a row

    return (dataIndex<10);
}

void updateTableAndGetActualMean (double table[100], double & mean, int const & changeIndex, double const & changeValue) {
    mean-=table[changeIndex]*0.01; table[changeIndex]=changeValue; mean+=changeValue*0.01;
}

int main(int argumentsNumber, char **arguments) {
/////////////////////////////////////////////// DANE WEJSCIOWE
    int testValue; do {
        char config[800];
        FILE *fileConfig = fopen("config.txt","rt");
        if (fileConfig==NULL) {
            printf("Missing file: config.txt\n");
            return 0;
        }
        int dataIndex=0,intervalLicznik=0;
        while(fgets(config,800,fileConfig)!=NULL) {
            sscanf(config,"%c",config);
            int actIndex=0,licznik=0;
            char data[20];
            while (config[actIndex]!='=') actIndex++;
            actIndex++;
            while (config[actIndex]!=';') data[licznik++]=config[actIndex++];
            switch (dataIndex) {
                case 0:testValue=strtol(data,NULL,10);break;
                case 1:N=strtol(data,NULL,10);break;
                case 2:gaps=strtol(data,NULL,10);break;
                case 3:initMode=strtol(data,NULL,10);break;
                case 4:deltaDiameter=strtod(data,NULL);break;
                case 5:alpha=strtol(data,NULL,10);break;
                case 6:beta=strtol(data,NULL,10);break;
                case 7:parM=strtol(data,NULL,10);break;
                case 8:parN=strtol(data,NULL,10);break;
                case 9:parXCells=strtol(data,NULL,10);break;
                case 10:parYCells=strtol(data,NULL,10);break;
                case 11:pressureRealOfNotFluid=strtod(data,NULL);break;
                case 12:growing=strtol(data,NULL,10);break;
                case 13:loadedConfiguration=strtol(data,NULL,10);break;
                case 14:loadedArg=strtod(data,NULL);break;
                case 15:{strcpy(loadedJOBID,"j-"); strncat(loadedJOBID,data,licznik);}break;
                case 16:loadedSetStartGenerator=strtol(data,NULL,10);break;
                case 17:loadedSetGenerator=strtol(data,NULL,10);break;
                case 18:iterationsNumber=strtol(data,NULL,10);break;
                case 19:intervalSampling=strtol(data,NULL,10);break;
                case 20:intervalOutput=strtol(data,NULL,10);break;
                case 21:saveConfigurations=strtol(data,NULL,10);break;
                case 22:savedConfigurationsInt=strtol(data,NULL,10);break;
                case 23:neighUpdatingFrequency=strtol(data,NULL,10);break;
                case 24:skipFirstIteration=strtol(data,NULL,10);break;
                case 25:useSpecificDirectory=strtol(data,NULL,10);break;
                case 26:cyclesOfEquilibration=strtol(data,NULL,10);break;
                case 27:cyclesOfMeasurement=strtol(data,NULL,10);break;
                case 28:intervalResults=strtol(data,NULL,10);break;
                case 29:maxDeltaR=strtod(data,NULL);break;
                case 30:desiredAcceptanceRatioR=strtod(data,NULL);break;
                case 31:desiredAcceptanceRatioV=strtod(data,NULL);break;
                case 32:useFileToIterate=strtol(data,NULL,10);break;
                case 33:startMinPacFrac=strtod(data,NULL);break;
                case 34:startMaxPacFrac=strtod(data,NULL);break;
                case 35:minArg=strtod(data,NULL);break;
                case 36:maxArg=strtod(data,NULL);break;
                case 37:multiplyArgument=strtol(data,NULL,10);break;
                case 38:multiplyFactor=strtod(data,NULL);break;
                default:
                    switch ((dataIndex-39)%3) {
                        case 0: intervalMin[intervalLicznik++/3]=strtod(data,NULL);break;
                        case 1: intervalMax[intervalLicznik++/3]=strtod(data,NULL);break;
                        case 2: intervalDelta[intervalLicznik++/3]=strtod(data,NULL);break;
                    }
                    break;
            }
            dataIndex++;
            for (int i=0;i<20;i++) data[i]=' ';
        }
        fclose(fileConfig);
    } while (testValue!=12345);

    //zlecanie parametrow z poziomu wiersza polecen:
    char JOBID[50]="j-"; int pointNumber=0;
    if (argumentsNumber==1) {
        strncat(JOBID,"none",50);
        if (useFileToIterate) if(createIterationTable()) return 0;
    } else {
        int correctNumberOfArguments=1;
        switch (strtol(arguments[1],NULL,10)) {
            case 0: //ustaw JOBID
                if (argumentsNumber==3) {
                    strncat(JOBID,arguments[2],50); strcpy(loadedJOBID,JOBID);
                    if (useFileToIterate) if(createIterationTable()) return 0;
                } else correctNumberOfArguments=0; break;
            case 1: //ustaw JOBID, singleRun dla parametrow zadanych bezposrednio
                if (argumentsNumber==14) {
                    useFileToIterate=0;
                    strncat(JOBID,arguments[2],50); strcpy(loadedJOBID,JOBID);
                    growing=strtol(arguments[9],NULL,10);
                    if (growing) {
                        startMinPacFrac=strtod(arguments[3],NULL); minArg=strtod(arguments[4],NULL);
                    } else {
                        startMaxPacFrac=strtod(arguments[3],NULL); maxArg=strtod(arguments[4],NULL);
                    }
                    N=strtol(arguments[5],NULL,10);
                    gaps=strtol(arguments[6],NULL,10);
                    deltaDiameter=strtod(arguments[7],NULL);
                    initMode=strtol(arguments[8],NULL,10);
                    iterationsNumber=strtol(arguments[10],NULL,10);
                    useSpecificDirectory=strtol(arguments[11],NULL,10);
                    alpha=strtol(arguments[12],NULL,10);
                    beta=strtol(arguments[13],NULL,10);
                    //generatorStartPoint=7;  //homo/hetero - 3: 44/80, 5: 41/86, 7:42/84 - generowanie konfiguracji binarnych losowych dla N=56
                } else correctNumberOfArguments=0; break;
            case 2: //ustaw JOBID, run z najistotniejszymi parametrami z 'config.txt' nadpisanymi z poziomu wywolania
                if (argumentsNumber==15) {
                    if (useFileToIterate) if(createIterationTable()) return 0;
                    strncat(JOBID,arguments[2],50);
                    N=strtol(arguments[3],NULL,10);
                    gaps=strtol(arguments[4],NULL,10);
                    deltaDiameter=strtod(arguments[5],NULL);
                    initMode=strtol(arguments[6],NULL,10);
                    growing=strtol(arguments[7],NULL,10);
                    iterationsNumber=strtol(arguments[8],NULL,10);
                    useSpecificDirectory=strtol(arguments[9],NULL,10);
                    skipFirstIteration=strtol(arguments[10],NULL,10);
                    pointNumber=strtol(arguments[11],NULL,10);
                    generatorStartPoint=strtol(arguments[12],NULL,10);
                    alpha=strtol(arguments[13],NULL,10);
                    beta=strtol(arguments[14],NULL,10);
                } else correctNumberOfArguments=0; break;
            case 3: //ustaw JOBID, tryb loadowany #1 od zadanego argumentu w odpowiednim folderze i trybie
                if (argumentsNumber==15) {
                    if (useFileToIterate) if(createIterationTable()) return 0;
                    strncat(JOBID,arguments[2],50);
                    strcpy(loadedJOBID,"j-"); strncat(loadedJOBID,arguments[3],50);
                    N=strtol(arguments[4],NULL,10);
                    gaps=strtol(arguments[5],NULL,10);
                    deltaDiameter=strtod(arguments[6],NULL);
                    initMode=strtol(arguments[7],NULL,10);
                    growing=strtol(arguments[8],NULL,10);
                    loadedConfiguration=1;
                    loadedArg=strtod(arguments[9],NULL);
                    iterationsNumber=strtol(arguments[10],NULL,10);
                    useSpecificDirectory=strtol(arguments[11],NULL,10);
                    skipFirstIteration=strtol(arguments[12],NULL,10);
                    alpha=strtol(arguments[13],NULL,10);
                    beta=strtol(arguments[14],NULL,10);
                } else correctNumberOfArguments=0; break;
            case 4: //ustaw JOBID, tryb loadowany #2 od zadanego numeru punktu (0->startArg) w odpowiednim folderze i trybie
                if (argumentsNumber==15) {
                    if (useFileToIterate) if(createIterationTable()) return 0;
                    strncat(JOBID,arguments[2],50);
                    strcpy(loadedJOBID,"j-"); strncat(loadedJOBID,arguments[3],50);
                    N=strtol(arguments[4],NULL,10);
                    gaps=strtol(arguments[5],NULL,10);
                    deltaDiameter=strtod(arguments[6],NULL);
                    initMode=strtol(arguments[7],NULL,10);
                    growing=strtol(arguments[8],NULL,10);
                    loadedConfiguration=1; loadType=1;
                    pointNumber=strtol(arguments[9],NULL,10);
                    iterationsNumber=strtol(arguments[10],NULL,10);
                    useSpecificDirectory=strtol(arguments[11],NULL,10);
                    skipFirstIteration=strtol(arguments[12],NULL,10);
                    alpha=strtol(arguments[13],NULL,10);
                    beta=strtol(arguments[14],NULL,10);
                } else correctNumberOfArguments=0; break;
            case 5: //ustaw JOBID, zrob tryb ONLYMATH, gdzie argument wskazuje ile poczatkowych linii Results ma byc pominietych
                if (argumentsNumber==14) {
                    if (useFileToIterate) if(createIterationTable()) return 0;
                    strncat(JOBID,arguments[2],50); strcpy(loadedJOBID,JOBID);
                    onlyMath[0]=1;
                    onlyMath[1]=strtol(arguments[3],NULL,10);
                    N=strtol(arguments[4],NULL,10);
                    gaps=strtol(arguments[5],NULL,10);
                    deltaDiameter=strtod(arguments[6],NULL);
                    initMode=strtol(arguments[7],NULL,10);
                    growing=strtol(arguments[8],NULL,10);
                    loadedConfiguration=0;
                    pointNumber=strtol(arguments[9],NULL,10);
                    iterationsNumber=strtol(arguments[10],NULL,10);
                    useSpecificDirectory=strtol(arguments[11],NULL,10);
                    skipFirstIteration=0;
                    alpha=strtol(arguments[12],NULL,10);
                    beta=strtol(arguments[13],NULL,10);
                } else correctNumberOfArguments=0; break;
            default: {
                printf("Wrong type of run! (0-6)\n");
                return 0;
            } break;
        }
        if (!correctNumberOfArguments) {
            printf("Wrong number of arguments for this type of run!\n");
            printf("If type of run is '0', next arguments: $JOBID\n");
            printf("If type of run is '1', next arguments: $JOBID, startMinPacFrac, minArg, N, gaps, deltaDiameter, initMode, growing, iterationsNumber, useSpecificDirectory, alpha, beta\n");
            printf("If type of run is '2', next arguments: $JOBID, N, gaps, deltaDiameter, initMode, growing, iterationsNumber, useSpecificDirectory, skipFirstIteration, pointNumber, generatorStartPoint, alpha, beta\n");
            printf("If type of run is '3', next arguments: $JOBID, JOBID of configuration to load, N, gaps, deltaDiameter, initMode, growing, loadedArg, iterationsNumber, useSpecificDirectory, skipFirstIteration, alpha, beta\n");
            printf("If type of run is '4', next arguments: $JOBID, JOBID of configuration to load, N, gaps, deltaDiameter, initMode, growing, pointNumber, iterationsNumber, useSpecificDirectory, skipFirstIteration, alpha, beta\n");
            printf("If type of run is '5', next arguments: $JOBID, lines to skip from Results, N, gaps, deltaDiameter, initMode, growing, pointNumber, iterationsNumber, useSpecificDirectory, alpha, beta\n");
            return 0;
        }
    }

    //ostatnie konfiguracje
    if (useFileToIterate) {
        if (growing) {
            startMinPacFrac=iterationTable[0][0];
            for (int i=0;i<2;i++) startArg[i]=iterationTable[0][i+1];
        } else {
            startMaxPacFrac=iterationTable[fileIterateIterationsNumber-1][0];
            for (int i=0;i<2;i++) startArg[i]=iterationTable[fileIterateIterationsNumber-1][i+1];
        }
    } else {
        startArg[0]=growing?minArg:maxArg;
        startArg[1]=0;
    }
    pressureRealOfNotFluid/=(sigma*sigma);
    deltaR=maxDeltaR*sigma; deltaV*=sigma;
    for (int i=0;i<pointNumber;i++) getNextArgument(startArg,false);
    if (loadedConfiguration && loadType) loadedArg=startArg[0];

    if (initMode==60 || initMode==61) {
        if (parM<1) {printf("ERROR: Wrong parameter parM: %d (should be greater or equal to 1).\n",parM); return 0;}
        else if (parN<parM+3) {printf("ERROR: Wrong parameter parN: %d (should be greater or equal to parM+3).\n",parN); return 0;}
        else if (parXCells<1 || parYCells<1) {printf("ERROR: Wrong parameter(s) parX/YCells (both should be greater or equal to 1).\n"); return 0;}
        else if (parM-beta*2<1) {printf("ERROR: Wrong parameter beta: %d (parM-beta*2 should be greater or equal to 1).\n",beta); return 0;}
        else if (parN-beta*3<parM-beta*2+3) {printf("ERROR: Wrong parameter beta: %d (parN-beta*3 should be greater or equal to parM-beta*2+3).\n",beta); return 0;}
        N=(parN+parM)*(2*(parN-parM)-4)*parXCells*parYCells;  //przy strukturze re-entrant N jest sterowane z poziomu par-N/M/XCells/YCells
    }
    activeN=N-gaps;

    if (sqrt(N/168.0)!=floor(sqrt(N/168.0)) && sqrt(N/196.0)!=floor(sqrt(N/196.0)) &&
        sqrt(N/48.0)!=floor(sqrt(N/48.0)) && sqrt(N/56.0)!=floor(sqrt(N/56.0)) &&
        sqrt(N/108.0)!=floor(sqrt(N/108.0)) && sqrt(N/300.0)!=floor(sqrt(N/300.0)) &&
        sqrt(N/588.0)!=floor(sqrt(N/588.0)) && sqrt(N/1444.0)!=floor(sqrt(N/1444.0)) &&
        sqrt(N/7396.0)!=floor(sqrt(N/7396.0)) && sqrt(N/2028.0)!=floor(sqrt(N/2028.0)) &&
        sqrt(N/5476.0)!=floor(sqrt(N/5476.0)) && sqrt(N/36.0)!=floor(sqrt(N/36.0)) &&
        fabs(sqrt(N)-floor(sqrt(N)))>0.000001 && N%780!=0 && !(initMode==60 || initMode==61) &&
        sqrt(N/432.0)!=floor(sqrt(N/432.0))) {
        printf("ERROR: Not supported N: %d.\n",N);
        return 0;
    } else if (((initMode==1 || initMode==2) && sqrt(N/168.0)!=floor(sqrt(N/168.0))) ||
               ((initMode==3 || initMode==4 || initMode==56 || initMode==57) && sqrt(N/196.0)!=floor(sqrt(N/196.0))) ||
               ((initMode==5 || initMode==6 || initMode==8 || initMode==9) && sqrt(N/48.0)!=floor(sqrt(N/48.0))) ||
               ((initMode==10 || initMode==11) && sqrt(N/108.0)!=floor(sqrt(N/108.0))) ||
               ((initMode==12 || initMode==13 || initMode==18 || initMode==19 || initMode==20 || initMode==21 || initMode==52 || initMode==53 || initMode==62 || initMode==63 || initMode==64 || initMode==65 || initMode==72 || initMode==73) && sqrt(N/300.0)!=floor(sqrt(N/300.0))) ||
               ((initMode==14 || initMode==15 || initMode==22 || initMode==23 || initMode==24 || initMode==25 || initMode==26 || initMode==27 || initMode==34 || initMode==35 || initMode==36 || initMode==37 || initMode==38 || initMode==39 || initMode==40 || initMode==41) && sqrt(N/588.0)!=floor(sqrt(N/588.0))) ||
               ((initMode==16 || initMode==17) && sqrt(N/1444.0)!=floor(sqrt(N/1444.0))) ||
               ((initMode==28 || initMode==29 || initMode==42 || initMode==43 || initMode==44 || initMode==45) && sqrt(N/7396.0)!=floor(sqrt(N/7396.0))) ||
               ((initMode==30 || initMode==31 || initMode==46 || initMode==47 || initMode==48 || initMode==49) && sqrt(N/2028.0)!=floor(sqrt(N/2028.0))) ||
               ((initMode==32 || initMode==33) && sqrt(N/5476.0)!=floor(sqrt(N/5476.0))) ||
               ((initMode==50 || initMode==51) && sqrt(N/192.0)!=floor(sqrt(N/192.0))) ||
               ((initMode==54 || initMode==55 || initMode==58 || initMode==59) && sqrt(N/36.0)!=floor(sqrt(N/36.0))) ||
               ((initMode==66 || initMode==67 || initMode==68 || initMode==69 || initMode==70 || initMode==71) && sqrt(N/432.0)!=floor(sqrt(N/432.0)))) {  //a initMode=0(pure HD) i 7(random) niech będą możliwe na dowolnych spośród zdefiniowanych wyżej
        printf("ERROR: Not supported initMode (%d) for N: %d.\n",initMode,N);
        return 0;
    }

    //nazwy folderow na podstawie parametrow programu
    sprintf(bufferG,"%d",growing); sprintf(bufferN,"%d",N); sprintf(bufferGaps,"%d",gaps); sprintf(bufferD,"%.3E",deltaDiameter);
    sprintf(bufferInitMode,"%d",initMode); sprintf(bufferAlpha,"%d",alpha); sprintf(bufferBeta,"%d",beta);

    int folderIndex=useSpecificDirectory, checkNext;
    char bufferCheckFolderExisting[200];
    FILE *checkFolderExisting;
    if (!folderIndex) do {
        sprintf(bufferFolderIndex,"%d",++folderIndex);
        strcpy(bufferCheckFolderExisting,resultsFileName);
        addAppendix(bufferCheckFolderExisting,JOBID,false);
        checkFolderExisting = fopen(bufferCheckFolderExisting,"rt");
        if (checkFolderExisting!=NULL) {
            fclose(checkFolderExisting);
            checkNext=1;
        } else checkNext=0;
    } while (checkNext);
    sprintf(bufferFolderIndex,"%d",folderIndex);
    addAppendix(resultsFileName,JOBID,false);
    addAppendix(excelResultsFileName,JOBID,false);
    addAppendix(configurationsFileName,JOBID,true);
    addAppendix(loadConfigurationsFileName,loadedJOBID,true); strncat(loadConfigurationsFileName,"_arg-",6); sprintf(buffer,"%.4E",loadedArg); strncat(loadConfigurationsFileName,buffer,100); strncat(loadConfigurationsFileName,".txt",5);
    addAppendix(configurationsListFileName,JOBID,false);
    addAppendix(probDensDistFunFileName,JOBID,true); addAppendix(probDensDistFunResultsFileName,JOBID,false);  //probDensDistFunMode

    particle particles[N];
    double args[10];

    FILE *fileResults, *fileExcelResults, *fileConfigurations, *fileSavedConfigurations, *fileConfigurationsList, *fileAllResults, *fileProbDensDistFun, *fileProbDensDistFunResults;    //probDensDistFunMode
    fileResults = fopen(resultsFileName,"rt"); if (fileResults==NULL) {
        fileResults = fopen(resultsFileName,"a");
        fprintf(fileResults,"Cycles\tPressure*\tVolume\tdVolume\tBoxMatrix[0][0]\tdBoxMatrix[0][0]\tBoxMatrix[1][1]\tdBoxMatrix[1][1]\tBoxMatrix[1][0]([0][1])\tdBoxMatrix[1][0]([0][1])\tRho\tdRho\tV/V_cp\tdV/V_cp\tS1111\tdS1111\tS1122\tdS1122\tS1212\tdS1212\tS2222\tdS2222\tS1112\tdS1112\tS1222\tdS1222\tavNu\tdAvNu\tavNu2\tdAvNu2\tavB*\tdAvB*\tavMy*\tdAvMy*\tavE*\tdAvE*\tavNuDir1\tdAvNuDir1\tavNuDir2\tdAvNuDir2\tavBDir1*\tdAvBDir1*\tavBDir2*\tdAvBDir2*\tavMyDir1*\tdAvMyDir1*\tavMyDir2*\tdAvMyDir2*\n");
        fclose(fileResults);
    }
    fileExcelResults = fopen(excelResultsFileName,"rt"); if (fileExcelResults==NULL) {
        fileExcelResults = fopen(excelResultsFileName,"a");
        fprintf(fileExcelResults,"Pressure*\tV/V_cp\tavNu\tavNu2\tavB*\tavMy*\tavE*\n");
        fclose(fileExcelResults);
    }




/////////////////////////////////////////////// WARUNKI POCZATKOWE

    double arg[2]={startArg[0],startArg[1]}, oldBoxMatrix[2][2];
    while (growing>=0) {
        double pressureReduced=arg[0], pressure=pressureReduced/sigma/sigma,            //pressureReduced=\tau*\sigma^2/kT, kT=1[hard] - Dwa punkty widzenia: 1) zmniejszenie sigma ZWIĘKSZA JEDNOSTKE p^*, zatem ten sam STAN FIZYCZNY jak przy sigma=1 bedzie przy mniejszym p^*. pReal to tak naprawde pReducedAdjusted. Inny, równoważny punkt widzenia, to 2) pReal redukuje objetosc, ktora NIE jest wyrazana w jednostkach sigma. Objetosc jest obliczana z boxMatrix, ktory jest inicjowany z czynnikiem *sigma, zatem MA jednostkę, a NIE jest zredukowany. Przeciez gdyby sigma=2, to boxMatrix bylby 2x wiekszy, a 'w jednostkach sigma' (zredukowany) powinien pozostac identyczny
               boxMatrix[2][2],detBoxMatrix,matrixOfParticlesSize[2],unitCellAtCP[2],   //obydwa sprowadzają się do tego, że przy liczeniu prawdopodobieństwa ma być jednostka zredukowana (bezwymiarowa): 1) zakłada, że volume jest zredukowane, więc dostosowuje pReduced do stanu fizycznego; 2) zakłada, że pReduced już jest OK (w końcu jest zredukowane), tylko po prostu objętość NIE jest zredukowana, i trzeba ją zredukować dzieląc przez sigma^2
               matrixCellXY[6][6][2];
        int n[2]; //n[X/Y], matrixCell[n[XMax]][n[YMax]][x/y], zatem: n[X/Y](max)=6
        unitCellAtCP[0]=sigma; unitCellAtCP[1]=sigma*sqrt(3);
        n[0]=1; n[1]=2;
        matrixCellXY[0][0][0]=0; matrixCellXY[0][0][1]=0;
        matrixCellXY[0][1][0]=unitCellAtCP[0]/2.0; matrixCellXY[0][1][1]=unitCellAtCP[1]/2.0;
        VcpPerParticle=unitCellAtCP[0]*unitCellAtCP[1]/(double)n[0]/(double)n[1];
        if (initMode==60 || initMode==61) {matrixOfParticlesSize[0]=(parN+parM)*parXCells; matrixOfParticlesSize[1]=(2*(parN-parM)-4)*parYCells;}     //struktury: re-entrant
        else if (sqrt(N/168.0)==floor(sqrt(N/168.0))) {matrixOfParticlesSize[0]=12; matrixOfParticlesSize[1]=14;}      //struktury: S1
        else if (sqrt(N/196.0)==floor(sqrt(N/196.0))) {matrixOfParticlesSize[0]=14; matrixOfParticlesSize[1]=14;} //struktury: S3, D1-2
        else if (sqrt(N/48.0)==floor(sqrt(N/48.0))) {matrixOfParticlesSize[0]=6; matrixOfParticlesSize[1]=8;}     //struktury: S6, S7, S37(Nx4)
        else if (sqrt(N/56.0)==floor(sqrt(N/56.0))) {matrixOfParticlesSize[0]=7; matrixOfParticlesSize[1]=8;}     //struktury: random i pure HD (one zadziałają również na wszystkich innych)
        else if (sqrt(N/108.0)==floor(sqrt(N/108.0))) {matrixOfParticlesSize[0]=9; matrixOfParticlesSize[1]=12;}     //struktury: S19
        else if (sqrt(N/300.0)==floor(sqrt(N/300.0))) {matrixOfParticlesSize[0]=15; matrixOfParticlesSize[1]=20;}     //struktury: D7-1 (D10, D13), S61, S61-zigzag[dowolny]
        else if (sqrt(N/588.0)==floor(sqrt(N/588.0))) {matrixOfParticlesSize[0]=21; matrixOfParticlesSize[1]=28;}     //struktury: D7-2 (D8, D9), D19-1 (D22-1, D25-1, D25-2, D31)
        else if (sqrt(N/1444.0)==floor(sqrt(N/1444.0))) {matrixOfParticlesSize[0]=38; matrixOfParticlesSize[1]=38;}     //struktury: D7-3
        else if (sqrt(N/7396.0)==floor(sqrt(N/7396.0))) {matrixOfParticlesSize[0]=86; matrixOfParticlesSize[1]=86;}     //struktury: D19-2 (D22-2, D25-3)
        else if (sqrt(N/2028.0)==floor(sqrt(N/2028.0))) {matrixOfParticlesSize[0]=39; matrixOfParticlesSize[1]=52;}     //struktury: D19-3 (D20, D21)
        else if (sqrt(N/5476.0)==floor(sqrt(N/5476.0))) {matrixOfParticlesSize[0]=74; matrixOfParticlesSize[1]=74;}     //struktury: D19-4
        else if (sqrt(N/36.0)==floor(sqrt(N/36.0))) {matrixOfParticlesSize[0]=6; matrixOfParticlesSize[1]=6;}     //struktury: D1-1 (D2)
        else if (sqrt(N/432.0)==floor(sqrt(N/432.0))) {matrixOfParticlesSize[0]=18; matrixOfParticlesSize[1]=24;}     //struktury: S91, S91-zigzag[dowolny]
        else if (floor(sqrt(N))==sqrt(N)) {matrixOfParticlesSize[0]=matrixOfParticlesSize[1]=sqrt(N);}     //struktury: random i pure HD (one zadziałają również na wszystkich innych)
        else if (N%780==0) {matrixOfParticlesSize[0]=26; matrixOfParticlesSize[1]=30;}     //struktury: random i pure HD (one zadziałają również na wszystkich innych)
        double NLinearMod = sqrt(N/matrixOfParticlesSize[0]/matrixOfParticlesSize[1]);
        if (startArg[0]==arg[0]) {
            for (int i=0;i<2;i++) boxMatrix[i][i]=matrixOfParticlesSize[i]*unitCellAtCP[i]/(double)n[i]*NLinearMod*(growing?sqrt(startMinPacFrac):sqrt(startMaxPacFrac));
            boxMatrix[1][0]=0.0; boxMatrix[0][1]=0.0;
        } else {
            for (int i=0;i<2;i++) for (int j=0;j<2;j++) boxMatrix[i][j]=oldBoxMatrix[i][j];
        }
        detBoxMatrix=boxMatrix[0][0]*boxMatrix[1][1]-boxMatrix[1][0]*boxMatrix[0][1];
        double volume=fabs(detBoxMatrix), rho=N/volume, pacFrac=1.0/VcpPerParticle/rho;

        if (!onlyMath[0]) {
            if (arg[0]==startArg[0] && !loadedConfiguration) {
                printf("INIT POS.- N: %d, gaps: %d, growing: %d, StartPressRed: %.7E (StartDen: %.7E, startPacFrac: %.7E, alpha: %d, beta: %d), deltaDiameter: %.3E\n",N,gaps,growing,startArg[0],rho,pacFrac,alpha,beta,deltaDiameter);
                if (!initPositions(particles,boxMatrix,detBoxMatrix,matrixOfParticlesSize,n,matrixCellXY,pacFrac,volume)) return 0;
                updateNeighbourList(particles,boxMatrix,volume);
                printf("Checking overlaps in inited configuration... "); fflush(stdout);
                if (getEnergyAll(particles,boxMatrix)==1) {
                    printf("Configuration's init parameters generate overlap(s) [energy=1].\n");
                    return 0;
                } else  {printf("done\n"); fflush(stdout);}
            } else if (loadedConfiguration) {
                char configurations[4096];
                FILE *fileCTL = fopen(loadConfigurationsFileName,"rt");
                if (fileCTL==NULL) {
                    printf("Missing file (configuration): %s\n",loadConfigurationsFileName);
                    return 0;
                }
                for (int i=0;i<3;i++) fgets(configurations,4096,fileCTL);
                int character,dataType=0,pIndex=-1; char data[50]=""; int actIndex=0;
                while (dataType<11) {
                    character=fgetc(fileCTL); //character is in int, but it can work as char
                    if (dataType<10) { //stage #1 configuration parameters
                        if (character==' ') {
                            data[actIndex++]=' '; //end of data without clearing the entire array
                            args[dataType++]=strtod(data,NULL);
                            if (dataType==10) {
                                boxMatrix[0][0]=args[5]; boxMatrix[1][1]=args[6]; boxMatrix[1][0]=args[7]; boxMatrix[0][1]=args[7];
                                deltaR=args[8]; deltaV=args[9];
                                arg[0]=args[3]; pressureReduced=arg[0]; pressure=pressureReduced/sigma/sigma;

                                detBoxMatrix=boxMatrix[0][0]*boxMatrix[1][1]-boxMatrix[1][0]*boxMatrix[0][1];
                                volume=fabs(detBoxMatrix); rho=N/volume; pacFrac=1.0/VcpPerParticle/rho;
                            }
                            actIndex=0;
                            continue;
                        }
                        data[actIndex++]=character;
                    } else { //stage #2 configuration parameters (coordinates of particles)
                        if (character=='b' || character=='w') {pIndex++; particles[pIndex].type=(character=='b')?0:1; continue;}
                        if (pIndex>=0) {
                            for (int i=0;i<3;i++) {
                                character=fgetc(fileCTL);
                                while (character!=',' && character!=']') {
                                    data[actIndex++]=character;
                                    character=fgetc(fileCTL);
                                } data[actIndex++]=' '; actIndex=0;
                                if (i<2) particles[pIndex].r[i]=strtod(data,NULL);
                                else particles[pIndex].diameter=strtod(data,NULL);
                            }
                            particles[pIndex].normR[0]=(boxMatrix[1][1]*particles[pIndex].r[0]-boxMatrix[0][1]*particles[pIndex].r[1])/detBoxMatrix;
                            particles[pIndex].normR[1]=-(boxMatrix[1][0]*particles[pIndex].r[0]-boxMatrix[0][0]*particles[pIndex].r[1])/detBoxMatrix;
                            while (character!=']') character=fgetc(fileCTL); fgetc(fileCTL); //next to read: 'b' or 'w'
                            if (pIndex>=activeN-1) dataType++;
                        }
                    }
                }
                fclose(fileCTL);
                printf("LOADING POS.- N: %d, gaps: %d, growing: %d, StartPressRed: %.7E (startDen: %.7E, startPacFrac: %.7E, alpha: %d, beta: %d), RandStart: %.1f, RandStep: %.1f, Cycles: %ld, DeltaR: %.4E, DeltaV: %.4E\n",N,gaps,growing,args[3],args[2],pacFrac,alpha,beta,args[0],args[1],(long)args[4],args[8],args[9]);
                //for (int i=0;i<2;i++) for (int j=0;j<2;j++) printf("boxMatrix[%d][%d]=%.17E\n",i,j,boxMatrix[i][j]);
                //for (int i=0;i<activeN;i++) printf("%d: %.17E,  %.17E,  %d\n",i,particles[i].r[0],particles[i].r[1],particles[i].type);return 0;
                updateNeighbourList(particles,boxMatrix,volume);
                printf("Checking overlaps in loaded file... "); fflush(stdout);
                if (getEnergyAll(particles,boxMatrix)==1) {
                    printf("Configuration from loaded file contains overlap(s) [energy=1].\n");
                    return 0;
                } else  {printf("done\n"); fflush(stdout);}
            }
        }

        if (skipFirstIteration) {
            printf("Skipping first iteration...\n");
            skipFirstIteration=0;
        } else {
            if (!onlyMath[0]) {
                if (loadedConfiguration) {
                    if (loadedSetStartGenerator) {
                        printf("Setting start position of p-random number generator to position from file...\n");
                        InitMT((unsigned int)args[0]);
                        randomStartStep[0]=args[0];
                    } else {
                        printf("Setting start position of p-random number generator to position from file - DISABLED\n");
                        if (generatorStartPoint==0) {
                            generatorStartPoint=time(0);
                            printf("Setting start position of p-random number generator to actual CPU time...\n");
                        } else printf("Setting start position of p-random number generator to %ld...\n",generatorStartPoint);
                        InitMT((unsigned int)generatorStartPoint);
                        randomStartStep[0]=generatorStartPoint;
                    }
                    randomStartStep[1]=0;
                    if (loadedSetGenerator) {
                        printf("Setting p-random number generator to last position from file...\n");
                        for (double i=0;i<args[1];i++) MTGenerate(randomStartStep);
                    } else printf("Setting p-random number generator to last position from file - DISABLED\n");
                } else {//generatorStartPoint=0;//ADDED-NO NECESSARY
                    if (generatorStartPoint==0) {
                        generatorStartPoint=time(0);
                        printf("Setting start position of p-random number generator to actual CPU time...\n");
                    } else printf("Setting start position of p-random number generator to %ld...\n",generatorStartPoint);
                    InitMT((unsigned int)generatorStartPoint);
                    randomStartStep[0]=generatorStartPoint;
                    randomStartStep[1]=0;
                }
                printf("Start of equilibration at reduced pressure: %.7E (startDen: %.7E, startPacFrac: %.7E)... (%ld cycles)\n",pressureReduced,rho,pacFrac,cyclesOfEquilibration);
            } else printf("Start of mathOnly mode for: N: %d, gaps: %d, growing: %d, pressRed: %.7E, deltaDiameter: %.3E\n",N,gaps,growing,pressureReduced,deltaDiameter);
            fflush(stdout);




/////////////////////////////////////////////// RDZEN MC

            long volumeMoveChance=(int)ceil(activeN/sqrt(activeN)),
                exchangeMoveChance=0/*(int)ceil(activeN/4.0)*/,   //exchangeParticlesMoves (exchangeMoveChance=0 for INACTIVE, exchangeMoveChance=activeN/4 for ACTIVE)
                fullCycle=activeN+volumeMoveChance+exchangeMoveChance,
                fullCycleWithoutVMoves=activeN+exchangeMoveChance,
                cycle=0,   //UWAGA cycle LONG, nie moze byc za duzo cykli
                timeStart,timeEquilibration=0,timeEnd,
                attemptedNumberR=0, displacedNumberR=0,
                attemptedNumberV=0, displacedNumberV=0,
                cyclesOfMeasurementBuffer=arg[1]==0?cyclesOfMeasurement:0;
            double deltaRTable[100], deltaRMean=deltaR, deltaVTable[100], deltaVMean=deltaV;
            for (int i=0;i<100;i++) {deltaRTable[i]=deltaRMean; deltaVTable[i]=deltaVMean;}
            int simulationStage=cyclesOfEquilibration>0?0:cyclesOfMeasurementBuffer>0?1:2;  //0-equilibration, 1-measurement, 2-end
            int volumeMove=0, cycleCounter=0, indexScanned=(matrixOfParticlesSize[0]*round(matrixOfParticlesSize[1]*NLinearMod/2.0)-round(matrixOfParticlesSize[0]/2.0))*NLinearMod;

            char allResultsFileName[200],bufferConfigurations[200],bufferSavedConfigurations[200],bufferProbDensDistFun[200],bufferProbDensDistFunResults[200],bufferPressure[100];    //probDensDistFunMode
            strcpy(allResultsFileName,configurationsFileName); strcpy(bufferConfigurations,configurationsFileName);
            strcpy(bufferProbDensDistFun,probDensDistFunFileName); strcpy(bufferProbDensDistFunResults,probDensDistFunResultsFileName); //probDensDistFunMode
            sprintf(bufferPressure,"%.4E",pressureReduced);
            strncat(allResultsFileName,"_arg-",6); strncat(allResultsFileName,bufferPressure,100); strncat(allResultsFileName,"_Results.txt",13);
            strncat(bufferConfigurations,"_arg-",6); strncat(bufferConfigurations,bufferPressure,100); strcpy(bufferSavedConfigurations,bufferConfigurations); strncat(bufferConfigurations,".txt",5); strncat(bufferSavedConfigurations,"_transient.txt",15);
            strncat(bufferProbDensDistFun,"_arg-",6); strncat(bufferProbDensDistFun,bufferPressure,100); strncat(bufferProbDensDistFun,".txt",5); //probDensDistFunMode
            strncat(bufferProbDensDistFunResults,"_arg-",6); strncat(bufferProbDensDistFunResults,bufferPressure,100);    //probDensDistFunMode
            strncat(bufferProbDensDistFunResults,"_N-",4); sprintf(bufferPressure,"%d",N); strncat(bufferProbDensDistFunResults,bufferPressure,100);
            strncat(bufferProbDensDistFunResults,"_inM-",6); sprintf(bufferPressure,"%d",initMode); strncat(bufferProbDensDistFunResults,bufferPressure,100);
            strncat(bufferProbDensDistFunResults,".txt",5);

            fileAllResults = fopen(allResultsFileName,"a");
            if (saveConfigurations) fileSavedConfigurations = fopen(bufferSavedConfigurations,"a");
            //fileProbDensDistFun = fopen(bufferProbDensDistFun,"a");  //probDensDistFunMode - comment if not desired 1/7
            //char confsStepByStepFileName[200]="confsStepByStep.txt"; addAppendix(confsStepByStepFileName,JOBID,false); FILE *fileConfsStepByStep = fopen(confsStepByStepFileName,"a");  //consfStepByStep - comment if not desired 1/3
            if (onlyMath[0]) simulationStage=2;

            timeStart=time(0);
            while (simulationStage<2) {
                int randIndex;
                if (volumeMove) {
                    randIndex = (int)(MTRandom0to1(randomStartStep)*fullCycleWithoutVMoves);
                    volumeMove=0;
                } else randIndex = (int)(MTRandom0to1(randomStartStep)*fullCycle);
                if (randIndex<activeN) {
                    attemptedNumberR++;
                    if (attemptToDisplaceAParticle(particles,randIndex,boxMatrix,detBoxMatrix))
                        displacedNumberR++;
                } else if (randIndex<fullCycleWithoutVMoves) {
                    //type: los#1:spośród N; los#2:spośród N, ale tak długo, aż typ będzie inny niż w #1
                    /*int randIndex1=(int)(MTRandom0to1(randomStartStep)*activeN),
                        randIndex2=(int)(MTRandom0to1(randomStartStep)*activeN);
                    while (particles[randIndex1].type==particles[randIndex2].type) randIndex2=(int)(MTRandom0to1(randomStartStep)*activeN);*/
                    //

                    //type: los#1:spośród N; los#2:spośród particle[los#1].neighCounter
                    int randIndex1=(int)(MTRandom0to1(randomStartStep)*activeN),
                        randIndex2=particles[randIndex1].neighbours[(int)(MTRandom0to1(randomStartStep)*particles[randIndex1].neighCounter)];
                    if (particles[randIndex1].neighCounter>0 && particles[randIndex1].type!=particles[randIndex2].type)
                    //

                    attemptToExchangeParticles(particles,randIndex1,randIndex2,boxMatrix);
                } else {
                    volumeMove=1;
                    attemptedNumberV++;
                    if (attemptToChangeVolume(particles,pressure,boxMatrix,&detBoxMatrix,&volume))
                        displacedNumberV++;
                }

                /*bool saveConf=false;    //consfStepByStep - comment if not desired 2/3
                if (iStep<fullCycle*10 || ((long)iStep)%fullCycle==0) saveConf=true;
                if (saveConf) {
                    fprintf(fileConfsStepByStep,"multimers[x_,y_,kI_]:={{Opacity[If[x==0 && y==0,0.4,0]],Red,Polygon[{{0,0},{%.12E,%.12E},{%.12E,%.12E},{%.12E,%.12E}}]},Opacity[If[x==0 && y==0,1,0.3]]",
                            boxMatrix[0][0],boxMatrix[1][0],boxMatrix[0][0]+boxMatrix[0][1],boxMatrix[1][0]+boxMatrix[1][1],boxMatrix[0][1],boxMatrix[1][1]);
                    for (int i=0;i<activeN;i++) {
                        fprintf(fileConfsStepByStep,",%c[%.12E+x,%.12E+y,%.17E]",particles[i].type==0?'b':'w',particles[i].r[0],particles[i].r[1],particles[i].diameter);
                    }
                    fprintf(fileConfsStepByStep,"};\nAppendTo[VvsSTEP,{%ld,%.12E}];configurationsList=Append[configurationsList,g[%ld,multimers[%.12E,%.12E,kolorIndex=1],multimers[%.12E,%.12E,kolorIndex=1],multimers[%.12E,%.12E,kolorIndex=1],multimers[0,0,kolorIndex=1]]];\n",
                            (long)iStep,fabs(boxMatrix[0][0]*boxMatrix[1][1]-boxMatrix[1][0]*boxMatrix[0][1]),(long)iStep,boxMatrix[0][0],boxMatrix[1][0],boxMatrix[0][0]+boxMatrix[0][1],boxMatrix[1][0]+boxMatrix[1][1],boxMatrix[0][1],boxMatrix[1][1]);
                }*/

                cycleCounter++;
                if (cycleCounter>=fullCycle) {
                    cycleCounter=0;
                    cycle++;

                    // !!!FREQUENT UPDATES OF NEIGHBOUR LIST DURING EQUILIBRATION!!! - na potrzeby silnego ściskania z relatywnie mało-gęstych konfiguracji początkowych
                    //if (cycle%2==0) updateNeighbourList(particles,boxMatrix,volume);

                    if (cycle%intervalSampling==0) {
                        if (simulationStage==0 && cycle>cyclesOfEquilibration) {
                            simulationStage=1;
                            printf("Equilibration finished after: %ld cycles (%ldsec).\n",cyclesOfEquilibration,time(0)-timeStart);
                            fflush(stdout);
                        }
                        double acceptanceRatioR = displacedNumberR/(double)attemptedNumberR,
                               acceptanceRatioV = displacedNumberV/(double)attemptedNumberV;
                        if (cycle%neighUpdatingFrequency==0) updateNeighbourList(particles,boxMatrix,volume);

                        //ADDED-NO NECESSARY/////////
                        /*if (cycle%20000==0) {
                            printf("Checking overlaps in transient configuration... "); fflush(stdout);
                            if (getEnergyAll(particles,boxMatrix)==1) {
                                printf("Transient configuration contains overlap(s) [energy=1]. ");
                                char allResultsErrorFileName[200],configurationErrorFileName[200];
                                strcpy(allResultsErrorFileName,allResultsFileName); strncat(allResultsErrorFileName,".err",5);
                                strcpy(configurationErrorFileName,bufferConfigurations); strncat(configurationErrorFileName,".err",5);
                                if (rename(allResultsFileName,allResultsErrorFileName)==0) printf("Results file successfully renamed (.err). ");
                                else printf("Error renaming results file (.err). ");
                                if (rename(bufferConfigurations,configurationErrorFileName)==0) printf("Configuration file successfully renamed (.err).\n");
                                else printf("Error renaming configuration file (.err).\n");
                                return 0;
                            } else  {printf("done\n"); fflush(stdout);}
                        }*/

                        if (simulationStage==1) {
                            if (timeEquilibration==0) {
                                timeEquilibration=time(0);

                                printf("Checking overlaps in equilibrated configuration... "); fflush(stdout);
                                if (getEnergyAll(particles,boxMatrix)==1) {
                                    printf("Equilibrated configuration contains overlap(s) [energy=1]. ");
                                    char allResultsErrorFileName[200];
                                    strcpy(allResultsErrorFileName,allResultsFileName); strncat(allResultsErrorFileName,".err",5);
                                    if (rename(allResultsFileName,allResultsErrorFileName)==0) printf("Results file successfully renamed (.err).\n");
                                    else printf("Error renaming results file (.err).\n");
                                    return 0;
                                } else  {printf("done\n"); fflush(stdout);}
                            }

                            rho=N/volume; pacFrac=1.0/VcpPerParticle/rho;
                            if (cycle%intervalResults==0) {
                                fprintf(fileAllResults,"%.17E\t%.17E\t%.17E\t\n",boxMatrix[0][0],boxMatrix[1][1],boxMatrix[0][1]);

                                //Probability Density Distribution Function
                                //probDensDistFunMode - comment if not desired 2/7
                                /*for (int i=0;i<activeN;i++) for (int j=0;j<particles[i].neighCounter;j++) {
                                    if (i<particles[i].neighbours[j]) {
                                        getParticlesDistanceSquared(&particles[i],&particles[particles[i].neighbours[j]],boxMatrix);
                                        double drRoot=sqrt(dr[2]);

                                        fprintf(fileProbDensDistFun,"{%d,%d,%.17E},",i,particles[i].neighbours[j],drRoot-0.5*(particles[i].diameter+particles[particles[i].neighbours[j]].diameter));
                                    }
                                }*/
                            }

                            if (saveConfigurations && cycle%savedConfigurationsInt==0) {
                                fprintf(fileSavedConfigurations,"%ld\t%.12E\t%.12E\t%.12E\t{",(cycle+(long)args[4]),boxMatrix[0][0],boxMatrix[1][1],boxMatrix[0][1]);
                                int activeNMinus1=activeN-1;
                                for (int i=0;i<activeNMinus1;i++)
                                    fprintf(fileSavedConfigurations,"%c[%.17E,%.17E,%.17E],",particles[i].type==0?'b':'w',particles[i].r[0],particles[i].r[1],particles[i].diameter);
                                fprintf(fileSavedConfigurations,"%c[%.17E,%.17E,%.17E]}",particles[activeNMinus1].type==0?'b':'w',particles[activeNMinus1].r[0],particles[activeNMinus1].r[1],particles[activeNMinus1].diameter);
                                fprintf(fileSavedConfigurations,"\n");
                            }

                            if (cycle%intervalOutput==0) {
                                printf("Cycle: %ld, ",(cycle+(long)args[4]));
                                printf("simulation time: full-%ldsec, measurement-%ldsec\n",time(0)-timeStart,time(0)-timeEquilibration);
                                printf("   AccRatR: %.4E, dR: %.4E, AccRatV: %.4E, dV: %.4E\n",acceptanceRatioR,deltaR,acceptanceRatioV,deltaV);
                                printf("   Dens: %.4E, V/V_cp: %.4E, PressRed: %.4E\n",rho,pacFrac,pressureReduced);
                                printf("   box00: %.8E, box11: %.8E, box01(10): %.8E\n",boxMatrix[0][0],boxMatrix[1][1],boxMatrix[0][1]);
                                fflush(stdout);

                                fileConfigurations = fopen(bufferConfigurations,"w");
                                fprintf(fileConfigurations,"Rho: %.12E\tV/V_cp: %.12E\tPressureRed: %.12E\tRandStart: %.1f\tRandSteps: %.1f\tCycles: %ld\tEquilTime: %ldsec\tMeasuTime: %ldsec\n\n",rho,pacFrac,
                                        pressureReduced,randomStartStep[0],randomStartStep[1],(cycle+(long)args[4]),(timeEquilibration-timeStart),(timeEnd-timeEquilibration));
                                fprintf(fileConfigurations,"=====================\tConfigurations (data to load)\t=====================\n");
                                fprintf(fileConfigurations,"%.1f %.1f %.17E %.17E %ld %.17E %.17E %.17E %.17E %.17E {",randomStartStep[0],randomStartStep[1],rho,pressureReduced,(cycle+(long)args[4]),boxMatrix[0][0],boxMatrix[1][1],boxMatrix[0][1],deltaR,deltaV);
                                for (int i=0;i<activeN;i++)
                                    fprintf(fileConfigurations,"%c[%.17E,%.17E,%.17E],",particles[i].type==0?'b':'w',particles[i].r[0],particles[i].r[1],particles[i].diameter);
                                fprintf(fileConfigurations,"{Opacity[0.2],Red,Polygon[{{0,0},{%.12E,%.12E},{%.12E,%.12E},{%.12E,%.12E}}]},{Opacity[0.2],Green,Disk[{%.12E,%.12E},%.12E]}}",boxMatrix[0][0],boxMatrix[1][0],boxMatrix[0][0]+boxMatrix[0][1],boxMatrix[1][0]+boxMatrix[1][1],boxMatrix[0][1],boxMatrix[1][1],particles[indexScanned].r[0],particles[indexScanned].r[1],neighRadius);
                                fprintf(fileConfigurations,"\n==========================================\n\n\nboxMatrix[0][0]=%.12E, boxMatrix[1][1]=%.12E, boxMatrix[0][1]=boxMatrix[1][0]=%.12E",boxMatrix[0][0],boxMatrix[1][1],boxMatrix[0][1]);
                                fclose(fileConfigurations);
                            }
                        } //else {//dostosowywanie delt w trakcie pomiaru
                            if (acceptanceRatioR>desiredAcceptanceRatioR) deltaR*=1.05; else deltaR*=0.95;
                            if (deltaR>maxDeltaR*sigma) deltaR=maxDeltaR*sigma;
                            if (acceptanceRatioV>desiredAcceptanceRatioV) deltaV*=1.05; else deltaV*=0.95;

                            int sampleNumberMod100=(cycle/intervalSampling)%100;
                            updateTableAndGetActualMean(deltaRTable,deltaRMean,sampleNumberMod100,deltaR); deltaR=deltaRMean;
                            updateTableAndGetActualMean(deltaVTable,deltaVMean,sampleNumberMod100,deltaV); deltaV=deltaVMean;
                        //}
                        attemptedNumberR=0; displacedNumberR=0;
                        attemptedNumberV=0; displacedNumberV=0;
                    }
                    if (simulationStage==1 && cycle-cyclesOfEquilibration>=cyclesOfMeasurementBuffer) simulationStage=2;
                }
            }
            fclose(fileAllResults);
            if (saveConfigurations) fclose(fileSavedConfigurations);
            //fclose(fileProbDensDistFun); //probDensDistFunMode - comment if not desired 3/7
            //fclose(fileConfsStepByStep);    //consfStepByStep - comment if not desired 3/3
            if (timeEquilibration==0) timeEquilibration=time(0);
            printf("Checking overlaps in final configuration... "); fflush(stdout);
            if (!onlyMath[0] && getEnergyAll(particles,boxMatrix)==1) {
                printf("Final configuration contains overlap(s) [energy=1]. ");
                char allResultsErrorFileName[200],configurationErrorFileName[200];
                strcpy(allResultsErrorFileName,allResultsFileName); strncat(allResultsErrorFileName,".err",5);
                strcpy(configurationErrorFileName,bufferConfigurations); strncat(configurationErrorFileName,".err",5);
                if (rename(allResultsFileName,allResultsErrorFileName)==0) printf("Results file successfully renamed (.err). ");
                else printf("Error renaming results file (.err). ");
                if (rename(bufferConfigurations,configurationErrorFileName)==0) printf("Configuration file successfully renamed (.err).\n");
                else printf("Error renaming configuration file (.err).\n");
                return 0;
            } else {printf("done\n"); fflush(stdout);}
            timeEnd=time(0);




/////////////////////////////////////////////// OBLICZENIE WYNIKOW

            mpfr_prec_t prec=1000;
            printf("Start of calculation of results...\n");

            //obliczenie liczby linii danych (potrzebne do podziału na zespoły i obliczenia średnich błędów)
            printf("Calculation of data lines... "); fflush(stdout);
            fileAllResults=fopen(allResultsFileName,"rt");
            char linia[4096]; double dataLicznik=0; int faultyLines=0, onlyMathLinesBuffer=onlyMath[1];
            if (onlyMath[0]) for (int i=0;i<onlyMathLinesBuffer;i++) {
                if (fgets(linia,300,fileAllResults)!=NULL && !isLineCorrect(linia)) onlyMathLinesBuffer++;
            }
            while(fgets(linia,300,fileAllResults)!=NULL) {
                if (isLineCorrect(linia)) dataLicznik++; else faultyLines++;
            }
            fclose(fileAllResults);
            printf("done (Found %ld data lines [%d faulty lines occurred].",(long)dataLicznik,faultyLines);
            if ((long)dataLicznik%10>0) printf(" Last %ld lines won't be considered, due to calculations of averages in 10 sets.)\n",(long)dataLicznik%10); else printf(")\n");
            dataLicznik-=(long)dataLicznik%10;

            //obliczenie srednich wartosci mierzonych wielkosci
            printf("Calculation of averages... "); fflush(stdout);

            mpfr_t avVolumeSet[10], avBoxMatrixSet[10][3], avRhoSet[10], avPacFracSet[10]; //wyniki dzielone są na 10 zespołów (obliczane są nieskorelowane wzajemnie średnie "lokalne", do obliczenia błędu średniej "globalnej")
            for (int i=0;i<10;i++) {
                mpfr_init2(avVolumeSet[i],prec); mpfr_set_si(avVolumeSet[i],0,MPFR_RNDN);
                for (int j=0;j<3;j++) {
                    mpfr_init2(avBoxMatrixSet[i][j],prec);
                    mpfr_set_si(avBoxMatrixSet[i][j],0,MPFR_RNDN);
                }
                mpfr_init2(avRhoSet[i],prec); mpfr_set_si(avRhoSet[i],0,MPFR_RNDN);
                mpfr_init2(avPacFracSet[i],prec); mpfr_set_si(avPacFracSet[i],0,MPFR_RNDN);
            }
            fileAllResults=fopen(allResultsFileName,"rt");
            if (onlyMath[0]) for (int i=0;i<onlyMathLinesBuffer;i++) fgets(linia,300,fileAllResults);
            double lineCounter=0;
            while(fgets(linia,300,fileAllResults)!=NULL && lineCounter<dataLicznik) {
                sscanf(linia,"%c",linia);
                int actIndex=0;
                int dataIndex=0; double dataD[3]; while (dataIndex<3) {
                    char data[50]="";
                    int licznik=0, dotCounter=0;
                    while (linia[actIndex]!='\t' && (dotCounter=linia[actIndex]=='.'?dotCounter+1:dotCounter)<=1 && licznik<50) data[licznik++]=linia[actIndex++];
                    if (dotCounter>1 || licznik>=50) {dataIndex=10; continue;}
                    actIndex++;
                    if (dataIndex<10 && ((data[0]!='-' && data[1]!='.') || (data[0]=='-' && data[2]!='.'))) dataIndex=10;
                    else dataD[dataIndex++]=strtod(data,NULL);
                } if (dataIndex<10 && linia[actIndex]!='\n') dataIndex=10;
                if (dataIndex<10) {
                    int setIndex=(int)(lineCounter/dataLicznik*10.0);
                    for (int i=0;i<3;i++) mpfr_add_d(avBoxMatrixSet[setIndex][i],avBoxMatrixSet[setIndex][i],dataD[i],MPFR_RNDN);
                    lineCounter++;
                }
            }
            fclose(fileAllResults);
            mpfr_t avVolume, avBoxMatrix[3], avRho, avPacFrac;
            mpfr_init2(avVolume,prec); mpfr_set_si(avVolume,0,MPFR_RNDN);
            for (int j=0;j<3;j++) {mpfr_init2(avBoxMatrix[j],prec); mpfr_set_si(avBoxMatrix[j],0,MPFR_RNDN);}
            mpfr_init2(avRho,prec); mpfr_set_si(avRho,0,MPFR_RNDN);
            mpfr_init2(avPacFrac,prec); mpfr_set_si(avPacFrac,0,MPFR_RNDN);
            mpfr_t buff1, buff2; mpfr_init2(buff1,prec); mpfr_init2(buff2,prec);
            for (int i=0;i<10;i++) {
                for (int j=0;j<3;j++) {
                    mpfr_div_d(avBoxMatrixSet[i][j],avBoxMatrixSet[i][j],dataLicznik*0.1,MPFR_RNDN);
                    mpfr_add(avBoxMatrix[j],avBoxMatrix[j],avBoxMatrixSet[i][j],MPFR_RNDN);
                }
                mpfr_mul(buff1,avBoxMatrixSet[i][0],avBoxMatrixSet[i][1],MPFR_RNDN);
                mpfr_mul(buff2,avBoxMatrixSet[i][2],avBoxMatrixSet[i][2],MPFR_RNDN);
                mpfr_sub(buff1,buff1,buff2,MPFR_RNDN);
                mpfr_abs(avVolumeSet[i],buff1,MPFR_RNDN); mpfr_add(avVolume,avVolume,avVolumeSet[i],MPFR_RNDN);
                mpfr_si_div(avRhoSet[i],N,avVolumeSet[i],MPFR_RNDN); mpfr_add(avRho,avRho,avRhoSet[i],MPFR_RNDN);
                mpfr_d_div(avPacFracSet[i],1.0/VcpPerParticle,avRhoSet[i],MPFR_RNDN); mpfr_add(avPacFrac,avPacFrac,avPacFracSet[i],MPFR_RNDN);
            }
            mpfr_clear(buff1); mpfr_clear(buff2);
            mpfr_div_si(avVolume,avVolume,10,MPFR_RNDN); mpfr_div_si(avRho,avRho,10,MPFR_RNDN); mpfr_div_si(avPacFrac,avPacFrac,10,MPFR_RNDN);
            for (int i=0;i<3;i++) mpfr_div_si(avBoxMatrix[i],avBoxMatrix[i],10,MPFR_RNDN);

            //obliczenie bledow mierzonych wielkosci
            mpfr_t dAvVolume, dAvBoxMatrix[3], dAvRho, dAvPacFrac;
            mpfr_init2(dAvVolume,prec); mpfr_set_si(dAvVolume,0,MPFR_RNDN);
            for (int j=0;j<3;j++) {mpfr_init2(dAvBoxMatrix[j],prec); mpfr_set_si(dAvBoxMatrix[j],0,MPFR_RNDN);}
            mpfr_init2(dAvRho,prec); mpfr_set_si(dAvRho,0,MPFR_RNDN);
            mpfr_init2(dAvPacFrac,prec); mpfr_set_si(dAvPacFrac,0,MPFR_RNDN);
            mpfr_init2(buff1,prec);
            for (int i=0;i<10;i++) {mpfr_sub(buff1,avVolume,avVolumeSet[i],MPFR_RNDN); mpfr_mul(buff1,buff1,buff1,MPFR_RNDN); mpfr_add(dAvVolume,dAvVolume,buff1,MPFR_RNDN);} mpfr_div_d(dAvVolume,dAvVolume,90.0,MPFR_RNDN); mpfr_sqrt(dAvVolume,dAvVolume,MPFR_RNDN); //10*9 (n(n-1))
            for (int j=0;j<3;j++) {for (int i=0;i<10;i++) {mpfr_sub(buff1,avBoxMatrix[j],avBoxMatrixSet[i][j],MPFR_RNDN); mpfr_mul(buff1,buff1,buff1,MPFR_RNDN); mpfr_add(dAvBoxMatrix[j],dAvBoxMatrix[j],buff1,MPFR_RNDN);} mpfr_div_d(dAvBoxMatrix[j],dAvBoxMatrix[j],90.0,MPFR_RNDN); mpfr_sqrt(dAvBoxMatrix[j],dAvBoxMatrix[j],MPFR_RNDN);}
            for (int i=0;i<10;i++) {mpfr_sub(buff1,avRho,avRhoSet[i],MPFR_RNDN); mpfr_mul(buff1,buff1,buff1,MPFR_RNDN); mpfr_add(dAvRho,dAvRho,buff1,MPFR_RNDN);} mpfr_div_d(dAvRho,dAvRho,90.0,MPFR_RNDN); mpfr_sqrt(dAvRho,dAvRho,MPFR_RNDN);
            for (int i=0;i<10;i++) {mpfr_sub(buff1,avPacFrac,avPacFracSet[i],MPFR_RNDN); mpfr_mul(buff1,buff1,buff1,MPFR_RNDN); mpfr_add(dAvPacFrac,dAvPacFrac,buff1,MPFR_RNDN);} mpfr_div_d(dAvPacFrac,dAvPacFrac,90.0,MPFR_RNDN); mpfr_sqrt(dAvPacFrac,dAvPacFrac,MPFR_RNDN);
            mpfr_clear(buff1);
            printf("done\n");

            //obliczenie srednich iloczynow elementow tensora odkształceń
            printf("Calculation of average values of products of strain tensor's elements... "); fflush(stdout);
            mpfr_t e1111Set[10], e1122Set[10], e1212Set[10], e2222Set[10], e1112Set[10], e1222Set[10], eVSQSet[10], eXXPeYYSQSet[10], eXXMeYYSQSet[10],
                   HxyHyx,HxxHyy,HxxHxy,Hxx2,Hyy2,mod0,mod1;
            for (int i=0;i<10;i++) {
                mpfr_init2(e1111Set[i],prec); mpfr_set_si(e1111Set[i],0,MPFR_RNDN); mpfr_init2(e1122Set[i],prec); mpfr_set_si(e1122Set[i],0,MPFR_RNDN);
                mpfr_init2(e1212Set[i],prec); mpfr_set_si(e1212Set[i],0,MPFR_RNDN); mpfr_init2(e2222Set[i],prec); mpfr_set_si(e2222Set[i],0,MPFR_RNDN);
                mpfr_init2(e1112Set[i],prec); mpfr_set_si(e1112Set[i],0,MPFR_RNDN); mpfr_init2(e1222Set[i],prec); mpfr_set_si(e1222Set[i],0,MPFR_RNDN);
                mpfr_init2(eVSQSet[i],prec); mpfr_set_si(eVSQSet[i],0,MPFR_RNDN); mpfr_init2(eXXPeYYSQSet[i],prec); mpfr_set_si(eXXPeYYSQSet[i],0,MPFR_RNDN);
                mpfr_init2(eXXMeYYSQSet[i],prec); mpfr_set_si(eXXMeYYSQSet[i],0,MPFR_RNDN);}
            mpfr_init2(HxyHyx,prec); mpfr_init2(HxxHyy,prec); mpfr_init2(HxxHxy,prec); mpfr_init2(Hxx2,prec); mpfr_init2(Hyy2,prec); mpfr_init2(mod0,prec); mpfr_init2(mod1,prec);
            fileAllResults=fopen(allResultsFileName,"rt");
            if (onlyMath[0]) for (int i=0;i<onlyMathLinesBuffer;i++) fgets(linia,300,fileAllResults);
            lineCounter=0; int setIndex, oldSetIndex=-1;
            mpfr_t hxyhyx, hxxPhyy, hxx2, hyy2, e11, e22, e12; mpfr_init2(hxyhyx,prec); mpfr_init2(hxxPhyy,prec); mpfr_init2(hxx2,prec); mpfr_init2(hyy2,prec); mpfr_init2(e11,prec); mpfr_init2(e22,prec); mpfr_init2(e12,prec);
            mpfr_init2(buff1,prec); mpfr_init2(buff2,prec); mpfr_t buff3; mpfr_init2(buff3,prec);
            while(fgets(linia,300,fileAllResults)!=NULL && lineCounter<dataLicznik) {
                setIndex=(int)(lineCounter/dataLicznik*10.0);
                if (setIndex!=oldSetIndex) {
                    mpfr_mul(HxyHyx,avBoxMatrixSet[setIndex][2],avBoxMatrixSet[setIndex][2],MPFR_RNDN);
                    mpfr_mul(HxxHyy,avBoxMatrixSet[setIndex][0],avBoxMatrixSet[setIndex][1],MPFR_RNDN);
                    mpfr_mul(HxxHxy,avBoxMatrixSet[setIndex][0],avBoxMatrixSet[setIndex][2],MPFR_RNDN);
                    mpfr_mul(Hxx2,avBoxMatrixSet[setIndex][0],avBoxMatrixSet[setIndex][0],MPFR_RNDN);
                    mpfr_mul(Hyy2,avBoxMatrixSet[setIndex][1],avBoxMatrixSet[setIndex][1],MPFR_RNDN);
                    mpfr_sub(mod0,HxyHyx,HxxHyy,MPFR_RNDN);
                    mpfr_mul(mod1,mod0,mod0,MPFR_RNDN); mpfr_d_div(mod1,0.5,mod1,MPFR_RNDN);
                    oldSetIndex=setIndex;
                }
                sscanf(linia,"%c",linia);
                int actIndex=0;
                double h[3],h11,h22,h12;
                int dataIndex=0; while (dataIndex<3) {
                    char data[50]="";
                    int licznik=0, dotCounter=0;
                    while (linia[actIndex]!='\t' && (dotCounter=linia[actIndex]=='.'?dotCounter+1:dotCounter)<=1 && licznik<50) data[licznik++]=linia[actIndex++];
                    if (dotCounter>1 || licznik>=50) {dataIndex=10; continue;}
                    actIndex++;
                    if (dataIndex<10) {
                        if ((data[0]!='-' && data[1]!='.') || (data[0]=='-' && data[2]!='.')) dataIndex=10;
                        else h[dataIndex++]=strtod(data,NULL);
                    }
                } if (dataIndex<10 && linia[actIndex]!='\n') dataIndex=10;
                if (dataIndex<10) {
                    h11=h[0]; h22=h[1]; h12=h[2];

                    mpfr_set_d(hxyhyx,h12,MPFR_RNDN); mpfr_mul(hxyhyx,hxyhyx,hxyhyx,MPFR_RNDN);
                    mpfr_set_d(hxxPhyy,h11,MPFR_RNDN); mpfr_add_d(hxxPhyy,hxxPhyy,h22,MPFR_RNDN);
                    mpfr_set_d(hxx2,h11,MPFR_RNDN); mpfr_mul(hxx2,hxx2,hxx2,MPFR_RNDN);
                    mpfr_set_d(hyy2,h22,MPFR_RNDN); mpfr_mul(hyy2,hyy2,hyy2,MPFR_RNDN);

                    //e11
                    mpfr_mul_d(buff1,avBoxMatrixSet[setIndex][2],h11,MPFR_RNDN); mpfr_mul_d(buff1,buff1,h12,MPFR_RNDN);
                    mpfr_mul(buff2,avBoxMatrixSet[setIndex][0],HxyHyx,MPFR_RNDN); mpfr_sub(buff1,buff1,buff2,MPFR_RNDN);
                    mpfr_mul_d(buff2,avBoxMatrixSet[setIndex][2],h12,MPFR_RNDN); mpfr_mul_d(buff2,buff2,h22,MPFR_RNDN);
                    mpfr_add(buff1,buff1,buff2,MPFR_RNDN); mpfr_mul(buff1,buff1,avBoxMatrixSet[setIndex][1],MPFR_RNDN); mpfr_mul_si(buff1,buff1,2,MPFR_RNDN);
                    mpfr_sub(buff2,hxyhyx,Hxx2,MPFR_RNDN); mpfr_add(buff2,buff2,hxx2,MPFR_RNDN); mpfr_mul(buff2,buff2,Hyy2,MPFR_RNDN);
                    mpfr_sub(buff1,buff2,buff1,MPFR_RNDN); mpfr_sub(buff2,hyy2,HxyHyx,MPFR_RNDN); mpfr_add(buff2,buff2,hxyhyx,MPFR_RNDN);
                    mpfr_mul(buff2,buff2,HxyHyx,MPFR_RNDN); mpfr_add(buff1,buff1,buff2,MPFR_RNDN); mpfr_mul(e11,mod1,buff1,MPFR_RNDN);
                    //e22
                    mpfr_mul_d(buff1,HxxHxy,h12,MPFR_RNDN); mpfr_mul(buff1,buff1,hxxPhyy,MPFR_RNDN); mpfr_mul(buff2,HxxHyy,HxyHyx,MPFR_RNDN);
                    mpfr_sub(buff1,buff1,buff2,MPFR_RNDN); mpfr_mul_si(buff1,buff1,2,MPFR_RNDN); mpfr_add(buff2,hxyhyx,hyy2,MPFR_RNDN);
                    mpfr_sub(buff2,buff2,Hyy2,MPFR_RNDN); mpfr_mul(buff2,buff2,Hxx2,MPFR_RNDN); mpfr_sub(buff1,buff2,buff1,MPFR_RNDN);
                    mpfr_add(buff2,hxx2,hxyhyx,MPFR_RNDN); mpfr_sub(buff2,buff2,HxyHyx,MPFR_RNDN); mpfr_mul(buff2,buff2,HxyHyx,MPFR_RNDN);
                    mpfr_add(buff1,buff1,buff2,MPFR_RNDN); mpfr_mul(e22,mod1,buff1,MPFR_RNDN);
                    //e12
                    mpfr_mul_d(buff1,avBoxMatrixSet[setIndex][2],h12,MPFR_RNDN); mpfr_mul(buff1,buff1,hxxPhyy,MPFR_RNDN); mpfr_add(buff2,hxx2,hxyhyx,MPFR_RNDN);
                    mpfr_mul(buff2,buff2,avBoxMatrixSet[setIndex][1],MPFR_RNDN); mpfr_sub(buff1,buff1,buff2,MPFR_RNDN); mpfr_mul(buff1,buff1,avBoxMatrixSet[setIndex][2],MPFR_RNDN);
                    mpfr_add(buff2,hxyhyx,hyy2,MPFR_RNDN); mpfr_mul(buff2,buff2,HxxHxy,MPFR_RNDN); mpfr_mul_d(buff3,HxxHyy,h12,MPFR_RNDN);
                    mpfr_mul(buff3,buff3,hxxPhyy,MPFR_RNDN); mpfr_add(buff1,buff1,buff3,MPFR_RNDN); mpfr_sub(buff1,buff1,buff2,MPFR_RNDN);
                    mpfr_mul(e12,mod1,buff1,MPFR_RNDN);

                    mpfr_mul(buff1,e11,e11,MPFR_RNDN); mpfr_add(e1111Set[setIndex],e1111Set[setIndex],buff1,MPFR_RNDN);
                    mpfr_mul(buff1,e11,e22,MPFR_RNDN); mpfr_add(e1122Set[setIndex],e1122Set[setIndex],buff1,MPFR_RNDN);
                    mpfr_mul(buff1,e12,e12,MPFR_RNDN); mpfr_add(e1212Set[setIndex],e1212Set[setIndex],buff1,MPFR_RNDN);
                    mpfr_mul(buff1,e22,e22,MPFR_RNDN); mpfr_add(e2222Set[setIndex],e2222Set[setIndex],buff1,MPFR_RNDN);
                    mpfr_mul(buff1,e11,e12,MPFR_RNDN); mpfr_add(e1112Set[setIndex],e1112Set[setIndex],buff1,MPFR_RNDN);
                    mpfr_mul(buff1,e12,e22,MPFR_RNDN); mpfr_add(e1222Set[setIndex],e1222Set[setIndex],buff1,MPFR_RNDN);

                    //obliczenia bezposrednie B i My (a nie poprzez podatnosci)
                    mpfr_set_d(buff1,h11,MPFR_RNDN); mpfr_mul_d(buff1,buff1,h22,MPFR_RNDN);
                    mpfr_set_d(buff2,h12,MPFR_RNDN); mpfr_mul_d(buff2,buff2,h12,MPFR_RNDN);
                    mpfr_sub(buff1,buff1,buff2,MPFR_RNDN); mpfr_abs(buff1,buff1,MPFR_RNDN); mpfr_sub(buff1,buff1,avVolumeSet[setIndex],MPFR_RNDN);
                    mpfr_pow_si(buff1,buff1,2,MPFR_RNDN); mpfr_add(eVSQSet[setIndex],eVSQSet[setIndex],buff1,MPFR_RNDN);
                    mpfr_add(buff1,e11,e22,MPFR_RNDN); mpfr_pow_si(buff1,buff1,2,MPFR_RNDN); mpfr_add(eXXPeYYSQSet[setIndex],eXXPeYYSQSet[setIndex],buff1,MPFR_RNDN);
                    mpfr_sub(buff1,e11,e22,MPFR_RNDN); mpfr_pow_si(buff1,buff1,2,MPFR_RNDN); mpfr_add(eXXMeYYSQSet[setIndex],eXXMeYYSQSet[setIndex],buff1,MPFR_RNDN);

                    lineCounter++;
                }
            }
            mpfr_clear(buff1); mpfr_clear(buff2); mpfr_clear(buff3);
            mpfr_clear(hxyhyx); mpfr_clear(hxxPhyy); mpfr_clear(hxx2); mpfr_clear(hyy2); mpfr_clear(e11); mpfr_clear(e22); mpfr_clear(e12);
            fclose(fileAllResults);
            mpfr_t e1111, e1122, e1212, e2222, e1112, e1222, eVSQ, eXXPeYYSQ, eXXMeYYSQ;
            mpfr_init2(e1111,prec); mpfr_set_si(e1111,0,MPFR_RNDN); mpfr_init2(e1122,prec); mpfr_set_si(e1122,0,MPFR_RNDN); mpfr_init2(e1212,prec); mpfr_set_si(e1212,0,MPFR_RNDN);
            mpfr_init2(e2222,prec); mpfr_set_si(e2222,0,MPFR_RNDN); mpfr_init2(e1112,prec); mpfr_set_si(e1112,0,MPFR_RNDN); mpfr_init2(e1222,prec); mpfr_set_si(e1222,0,MPFR_RNDN);
            mpfr_init2(eVSQ,prec); mpfr_set_si(eVSQ,0,MPFR_RNDN); mpfr_init2(eXXPeYYSQ,prec); mpfr_set_si(eXXPeYYSQ,0,MPFR_RNDN); mpfr_init2(eXXMeYYSQ,prec); mpfr_set_si(eXXMeYYSQ,0,MPFR_RNDN);
            for (int i=0;i<10;i++) {
                mpfr_div_d(e1111Set[i],e1111Set[i],dataLicznik*0.1,MPFR_RNDN); mpfr_add(e1111,e1111,e1111Set[i],MPFR_RNDN);
                mpfr_div_d(e1122Set[i],e1122Set[i],dataLicznik*0.1,MPFR_RNDN); mpfr_add(e1122,e1122,e1122Set[i],MPFR_RNDN);
                mpfr_div_d(e1212Set[i],e1212Set[i],dataLicznik*0.1,MPFR_RNDN); mpfr_add(e1212,e1212,e1212Set[i],MPFR_RNDN);
                mpfr_div_d(e2222Set[i],e2222Set[i],dataLicznik*0.1,MPFR_RNDN); mpfr_add(e2222,e2222,e2222Set[i],MPFR_RNDN);
                mpfr_div_d(e1112Set[i],e1112Set[i],dataLicznik*0.1,MPFR_RNDN); mpfr_add(e1112,e1112,e1112Set[i],MPFR_RNDN);
                mpfr_div_d(e1222Set[i],e1222Set[i],dataLicznik*0.1,MPFR_RNDN); mpfr_add(e1222,e1222,e1222Set[i],MPFR_RNDN);
                mpfr_div_d(eVSQSet[i],eVSQSet[i],dataLicznik*0.1,MPFR_RNDN); mpfr_add(eVSQ,eVSQ,eVSQSet[i],MPFR_RNDN);
                mpfr_div_d(eXXPeYYSQSet[i],eXXPeYYSQSet[i],dataLicznik*0.1,MPFR_RNDN); mpfr_add(eXXPeYYSQ,eXXPeYYSQ,eXXPeYYSQSet[i],MPFR_RNDN);
                mpfr_div_d(eXXMeYYSQSet[i],eXXMeYYSQSet[i],dataLicznik*0.1,MPFR_RNDN); mpfr_add(eXXMeYYSQ,eXXMeYYSQ,eXXMeYYSQSet[i],MPFR_RNDN);
            }
            mpfr_div_si(e1111,e1111,10,MPFR_RNDN); mpfr_div_si(e1122,e1122,10,MPFR_RNDN); mpfr_div_si(e1212,e1212,10,MPFR_RNDN);
            mpfr_div_si(e2222,e2222,10,MPFR_RNDN); mpfr_div_si(e1112,e1112,10,MPFR_RNDN); mpfr_div_si(e1222,e1222,10,MPFR_RNDN);
            mpfr_div_si(eVSQ,eVSQ,10,MPFR_RNDN); mpfr_div_si(eXXPeYYSQ,eXXPeYYSQ,10,MPFR_RNDN); mpfr_div_si(eXXMeYYSQ,eXXMeYYSQ,10,MPFR_RNDN);
            //obliczenie bledow iloczynow elementow tensora odkształceń
            mpfr_t dE1111, dE1122, dE1212, dE2222, dE1112, dE1222, dEVSQ, dEXXPeYYSQ, dEXXMeYYSQ;
            mpfr_init2(dE1111,prec); mpfr_set_si(dE1111,0,MPFR_RNDN); mpfr_init2(dE1122,prec); mpfr_set_si(dE1122,0,MPFR_RNDN); mpfr_init2(dE1212,prec); mpfr_set_si(dE1212,0,MPFR_RNDN);
            mpfr_init2(dE2222,prec); mpfr_set_si(dE2222,0,MPFR_RNDN); mpfr_init2(dE1112,prec); mpfr_set_si(dE1112,0,MPFR_RNDN); mpfr_init2(dE1222,prec); mpfr_set_si(dE1222,0,MPFR_RNDN);
            mpfr_init2(dEVSQ,prec); mpfr_set_si(dEVSQ,0,MPFR_RNDN); mpfr_init2(dEXXPeYYSQ,prec); mpfr_set_si(dEXXPeYYSQ,0,MPFR_RNDN); mpfr_init2(dEXXMeYYSQ,prec); mpfr_set_si(dEXXMeYYSQ,0,MPFR_RNDN);
            mpfr_init2(buff1,prec);
            for (int i=0;i<10;i++) {mpfr_sub(buff1,e1111,e1111Set[i],MPFR_RNDN); mpfr_mul(buff1,buff1,buff1,MPFR_RNDN); mpfr_add(dE1111,dE1111,buff1,MPFR_RNDN);} mpfr_div_d(dE1111,dE1111,90.0,MPFR_RNDN); mpfr_sqrt(dE1111,dE1111,MPFR_RNDN); //10*9 (n(n-1))
            for (int i=0;i<10;i++) {mpfr_sub(buff1,e1122,e1122Set[i],MPFR_RNDN); mpfr_mul(buff1,buff1,buff1,MPFR_RNDN); mpfr_add(dE1122,dE1122,buff1,MPFR_RNDN);} mpfr_div_d(dE1122,dE1122,90.0,MPFR_RNDN); mpfr_sqrt(dE1122,dE1122,MPFR_RNDN);
            for (int i=0;i<10;i++) {mpfr_sub(buff1,e1212,e1212Set[i],MPFR_RNDN); mpfr_mul(buff1,buff1,buff1,MPFR_RNDN); mpfr_add(dE1212,dE1212,buff1,MPFR_RNDN);} mpfr_div_d(dE1212,dE1212,90.0,MPFR_RNDN); mpfr_sqrt(dE1212,dE1212,MPFR_RNDN);
            for (int i=0;i<10;i++) {mpfr_sub(buff1,e2222,e2222Set[i],MPFR_RNDN); mpfr_mul(buff1,buff1,buff1,MPFR_RNDN); mpfr_add(dE2222,dE2222,buff1,MPFR_RNDN);} mpfr_div_d(dE2222,dE2222,90.0,MPFR_RNDN); mpfr_sqrt(dE2222,dE2222,MPFR_RNDN);
            for (int i=0;i<10;i++) {mpfr_sub(buff1,e1112,e1112Set[i],MPFR_RNDN); mpfr_mul(buff1,buff1,buff1,MPFR_RNDN); mpfr_add(dE1112,dE1112,buff1,MPFR_RNDN);} mpfr_div_d(dE1112,dE1112,90.0,MPFR_RNDN); mpfr_sqrt(dE1112,dE1112,MPFR_RNDN);
            for (int i=0;i<10;i++) {mpfr_sub(buff1,e1222,e1222Set[i],MPFR_RNDN); mpfr_mul(buff1,buff1,buff1,MPFR_RNDN); mpfr_add(dE1222,dE1222,buff1,MPFR_RNDN);} mpfr_div_d(dE1222,dE1222,90.0,MPFR_RNDN); mpfr_sqrt(dE1222,dE1222,MPFR_RNDN);
            for (int i=0;i<10;i++) {mpfr_sub(buff1,eVSQ,eVSQSet[i],MPFR_RNDN); mpfr_mul(buff1,buff1,buff1,MPFR_RNDN); mpfr_add(dEVSQ,dEVSQ,buff1,MPFR_RNDN);} mpfr_div_d(dEVSQ,dEVSQ,90.0,MPFR_RNDN); mpfr_sqrt(dEVSQ,dEVSQ,MPFR_RNDN);
            for (int i=0;i<10;i++) {mpfr_sub(buff1,eXXPeYYSQ,eXXPeYYSQSet[i],MPFR_RNDN); mpfr_mul(buff1,buff1,buff1,MPFR_RNDN); mpfr_add(dEXXPeYYSQ,dEXXPeYYSQ,buff1,MPFR_RNDN);} mpfr_div_d(dEXXPeYYSQ,dEXXPeYYSQ,90.0,MPFR_RNDN); mpfr_sqrt(dEXXPeYYSQ,dEXXPeYYSQ,MPFR_RNDN);
            for (int i=0;i<10;i++) {mpfr_sub(buff1,eXXMeYYSQ,eXXMeYYSQSet[i],MPFR_RNDN); mpfr_mul(buff1,buff1,buff1,MPFR_RNDN); mpfr_add(dEXXMeYYSQ,dEXXMeYYSQ,buff1,MPFR_RNDN);} mpfr_div_d(dEXXMeYYSQ,dEXXMeYYSQ,90.0,MPFR_RNDN); mpfr_sqrt(dEXXMeYYSQ,dEXXMeYYSQ,MPFR_RNDN);
            mpfr_clear(buff1);
            printf("done\n");

            //obliczenie podatnosci, wspolczynnika Poissona i modulow sprezystosci
            printf("Calculation of compliances, Poisson's ratio and elastic moduli... "); fflush(stdout);
            //Eijkl - bezwymiarowe [strain], volume - [sigma^2], kT=1[hard]
            mpfr_t s1111, dS1111, s1122, dS1122, s1212, dS1212, s2222, dS2222, s1112, dS1112, s1222, dS1222;
            mpfr_init2(buff1,prec); mpfr_init2(buff2,prec);
            mpfr_init2(s1111,prec); mpfr_init2(dS1111,prec); mpfr_mul(s1111,e1111,avVolume,MPFR_RNDN); mpfr_mul(buff1,e1111,dAvVolume,MPFR_RNDN); mpfr_abs(buff1,buff1,MPFR_RNDN); mpfr_mul(buff2,dE1111,avVolume,MPFR_RNDN); mpfr_abs(buff2,buff2,MPFR_RNDN); mpfr_add(dS1111,buff1,buff2,MPFR_RNDN); //Sijkl=Vp*<EijEkl>/(kT)
            mpfr_init2(s1122,prec); mpfr_init2(dS1122,prec); mpfr_mul(s1122,e1122,avVolume,MPFR_RNDN); mpfr_mul(buff1,e1122,dAvVolume,MPFR_RNDN); mpfr_abs(buff1,buff1,MPFR_RNDN); mpfr_mul(buff2,dE1122,avVolume,MPFR_RNDN); mpfr_abs(buff2,buff2,MPFR_RNDN); mpfr_add(dS1122,buff1,buff2,MPFR_RNDN);
            mpfr_init2(s1212,prec); mpfr_init2(dS1212,prec); mpfr_mul(s1212,e1212,avVolume,MPFR_RNDN); mpfr_mul(buff1,e1212,dAvVolume,MPFR_RNDN); mpfr_abs(buff1,buff1,MPFR_RNDN); mpfr_mul(buff2,dE1212,avVolume,MPFR_RNDN); mpfr_abs(buff2,buff2,MPFR_RNDN); mpfr_add(dS1212,buff1,buff2,MPFR_RNDN);
            mpfr_init2(s2222,prec); mpfr_init2(dS2222,prec); mpfr_mul(s2222,e2222,avVolume,MPFR_RNDN); mpfr_mul(buff1,e2222,dAvVolume,MPFR_RNDN); mpfr_abs(buff1,buff1,MPFR_RNDN); mpfr_mul(buff2,dE2222,avVolume,MPFR_RNDN); mpfr_abs(buff2,buff2,MPFR_RNDN); mpfr_add(dS2222,buff1,buff2,MPFR_RNDN);
            mpfr_init2(s1112,prec); mpfr_init2(dS1112,prec); mpfr_mul(s1112,e1112,avVolume,MPFR_RNDN); mpfr_mul(buff1,e1112,dAvVolume,MPFR_RNDN); mpfr_abs(buff1,buff1,MPFR_RNDN); mpfr_mul(buff2,dE1112,avVolume,MPFR_RNDN); mpfr_abs(buff2,buff2,MPFR_RNDN); mpfr_add(dS1112,buff1,buff2,MPFR_RNDN);
            mpfr_init2(s1222,prec); mpfr_init2(dS1222,prec); mpfr_mul(s1222,e1222,avVolume,MPFR_RNDN); mpfr_mul(buff1,e1222,dAvVolume,MPFR_RNDN); mpfr_abs(buff1,buff1,MPFR_RNDN); mpfr_mul(buff2,dE1222,avVolume,MPFR_RNDN); mpfr_abs(buff2,buff2,MPFR_RNDN); mpfr_add(dS1222,buff1,buff2,MPFR_RNDN);

            //TODO: double->mpfr_t
            double S11=(mpfr_get_d(s1111,MPFR_RNDN)+mpfr_get_d(s2222,MPFR_RNDN))*0.5, dS11=fabs((mpfr_get_d(dS1111,MPFR_RNDN)+mpfr_get_d(dS2222,MPFR_RNDN))*0.5),
                   S66=4.0*mpfr_get_d(s1212,MPFR_RNDN),
                   avNu=-mpfr_get_d(s1122,MPFR_RNDN)/S11, dAvNu=fabs(mpfr_get_d(dS1122,MPFR_RNDN)/S11)+fabs(dS11*mpfr_get_d(s1122,MPFR_RNDN)/S11/S11), //nu obliczane z Sxxyy
                   avNu2=S66/S11*0.5-1, dAvNu2=fabs(0.5/S11*4.0*mpfr_get_d(dS1212,MPFR_RNDN))+fabs(0.5*S66/S11/S11*dS11), //nu obliczane z Sxyxy (inny rodzaj scinania, przy izotropowych ukladach powinno byc tyle samo co nu1)

                   //\lambdaReduced=\lambda*\sigma^2/kT, kT=1[hard]
                   lambda=sigma*sigma/(8.0*(S11+mpfr_get_d(s1122,MPFR_RNDN))), dLambda=sigma*sigma*(fabs(dS11)+fabs(mpfr_get_d(dS1122,MPFR_RNDN)))/(8.0*fabs(S11+mpfr_get_d(s1122,MPFR_RNDN))*fabs(S11+mpfr_get_d(s1122,MPFR_RNDN))),

                   avB=4.0*lambda, dAvB=4.0*dLambda,
                   avMy=(avB-avB*avNu)/(1.0+avNu), dAvMy=fabs(dAvB*(1.0-avNu)/(1.0+avNu))+fabs(dAvNu*(-avB/(1.0+avNu)-(avB-avB*avNu)/(1.0+avNu)/(1.0+avNu))),
                   avE=4.0*avB*avMy/(avB+avMy), dAvE=4.0*(fabs(avB*avB*dAvMy)+fabs(dAvB*avMy*avMy))/fabs(avB+avMy)/fabs(avB+avMy),

                   //B=kT*Vp/<(V-Vp)^2>  OR  B=kT/(Vp*<(Exx+Eyy)^2>)
                   avBDirect1=mpfr_get_d(avVolume,MPFR_RNDN)/mpfr_get_d(eVSQ,MPFR_RNDN), dAvBDirect1=fabs(mpfr_get_d(dAvVolume,MPFR_RNDN)/mpfr_get_d(eVSQ,MPFR_RNDN))+fabs(mpfr_get_d(dEVSQ,MPFR_RNDN)*mpfr_get_d(avVolume,MPFR_RNDN)/mpfr_get_d(eVSQ,MPFR_RNDN)/mpfr_get_d(eVSQ,MPFR_RNDN)),
                   avBDirect2=1.0/(mpfr_get_d(avVolume,MPFR_RNDN)*mpfr_get_d(eXXPeYYSQ,MPFR_RNDN)), dAvBDirect2=fabs(mpfr_get_d(dAvVolume,MPFR_RNDN)/mpfr_get_d(eXXPeYYSQ,MPFR_RNDN)/mpfr_get_d(avVolume,MPFR_RNDN)/mpfr_get_d(avVolume,MPFR_RNDN))+fabs(mpfr_get_d(dEXXPeYYSQ,MPFR_RNDN)/mpfr_get_d(eXXPeYYSQ,MPFR_RNDN)/mpfr_get_d(eXXPeYYSQ,MPFR_RNDN)/mpfr_get_d(avVolume,MPFR_RNDN)),
                   //My=kT/(Vp*<(Exx-Eyy)^2>)  OR  My=kT/(4Vp*<Exy^2>)
                   avMyDirect1=1.0/(mpfr_get_d(avVolume,MPFR_RNDN)*mpfr_get_d(eXXMeYYSQ,MPFR_RNDN)), dAvMyDirect1=fabs(mpfr_get_d(dAvVolume,MPFR_RNDN)/mpfr_get_d(eXXMeYYSQ,MPFR_RNDN)/mpfr_get_d(avVolume,MPFR_RNDN)/mpfr_get_d(avVolume,MPFR_RNDN))+fabs(mpfr_get_d(dEXXMeYYSQ,MPFR_RNDN)/mpfr_get_d(eXXMeYYSQ,MPFR_RNDN)/mpfr_get_d(eXXMeYYSQ,MPFR_RNDN)/mpfr_get_d(avVolume,MPFR_RNDN)),
                   avMyDirect2=1.0/(4.0*mpfr_get_d(avVolume,MPFR_RNDN)*mpfr_get_d(e1212,MPFR_RNDN)), dAvMyDirect2=0.25*(fabs(mpfr_get_d(dAvVolume,MPFR_RNDN)/mpfr_get_d(e1212,MPFR_RNDN)/mpfr_get_d(avVolume,MPFR_RNDN)/mpfr_get_d(avVolume,MPFR_RNDN))+fabs(mpfr_get_d(dE1212,MPFR_RNDN)/mpfr_get_d(e1212,MPFR_RNDN)/mpfr_get_d(e1212,MPFR_RNDN)/mpfr_get_d(avVolume,MPFR_RNDN))),
                   //Nu=(B-My)/(B+My)
                   avNuDirect1=(avBDirect1-avMyDirect1)/(avBDirect1+avMyDirect1), dAvNuDirect1=2.0*(fabs(avBDirect1*dAvMyDirect1)+fabs(dAvBDirect1*avMyDirect1))/pow(avBDirect1+avMyDirect1,2),
                   avNuDirect2=(avBDirect2-avMyDirect2)/(avBDirect2+avMyDirect2), dAvNuDirect2=2.0*(fabs(avBDirect2*dAvMyDirect2)+fabs(dAvBDirect2*avMyDirect2))/pow(avBDirect2+avMyDirect2,2);
            mpfr_clear(buff1); mpfr_clear(buff2);
            printf("done\n");

            //wyznaczenie Probability Density Distribution Function
            /*printf("Determination of probability density distribution function... "); fflush(stdout);     //probDensDistFunMode - comment if not desired 4/7
            fileProbDensDistFun=fopen(bufferProbDensDistFun,"rt");
            double **avEdgeDistance=new double*[activeN]; for (int i=0;i<activeN;i++) avEdgeDistance[i]=new double[activeN];  //tablica tworzona dynamicznie zapisywana jest w HEAP a nie STACK (stack może zostać przepełniony i jest stackoverflow/segmentation fault   - pamiętać o delete[] tabel
            double licznik=0;
            for (int i=0;i<activeN;i++) for (int j=0;j<activeN;j++) avEdgeDistance[i][j]=0;
            char data[3][50]={"","",""}; int actIndex=0,character,dataType;
            while ((character=fgetc(fileProbDensDistFun))!=EOF) {
                if (character=='{') dataType=0;
                else if (character!=',' && character!='}') data[dataType][actIndex++]=character;
                else {
                    data[dataType++][actIndex++]=' '; actIndex=0;

                    if (dataType==3) {
                        int index[2]={(int)strtol(data[0],NULL,10),(int)strtol(data[1],NULL,10)};
                        if (index[0]==0 && index[1]==1) licznik++;  //i=0, j=1  identyfikuja pelen 'cykl' par (to pierwsza para jaka powinna zawsze byc zapisana wg metody updateNeighbours())
                        avEdgeDistance[index[0]][index[1]]+=strtod(data[2],NULL);
                    }
                    if (character=='}') fgetc(fileProbDensDistFun);
                }
            }
            fclose(fileProbDensDistFun);
            for (int i=0;i<activeN;i++) for (int j=0;j<activeN;j++) avEdgeDistance[i][j]/=licznik;
            printf("done\n");*/

            long timeEndMath=time(0);





/////////////////////////////////////////////// ZAPIS DANYCH DO PLIKU

            printf("Saving data to files... "); fflush(stdout);
            timeEq+=(timeEquilibration-timeStart); timeMe+=(timeEnd-timeEquilibration); timeMath+=(timeEndMath-timeEnd);

            fileResults = fopen(resultsFileName,"a");
            fileExcelResults = fopen(excelResultsFileName,"a");
            if (!onlyMath[0]) {
                fileConfigurations = fopen(bufferConfigurations,"w");
                fileConfigurationsList = fopen(configurationsListFileName,"a");
            }
            //fileProbDensDistFunResults = fopen(bufferProbDensDistFunResults,"w");     //probDensDistFunMode - comment if not desired 5/7

            fprintf(fileResults,"%ld\t%.17E\t%.17E\t%.17E\t%.17E\t%.17E\t%.17E\t%.17E\t%.17E\t%.17E\t%.17E\t%.17E\t%.17E\t%.17E\t%.17E\t%.17E\t%.17E\t%.17E\t%.17E\t%.17E\t%.17E\t%.17E\t%.17E\t%.17E\t%.17E\t%.17E\t%.17E\t%.17E\t%.17E\t%.17E\t%.17E\t%.17E\t%.17E\t%.17E\t%.17E\t%.17E\t%.17E\t%.17E\t%.17E\t%.17E\t%.17E\t%.17E\t%.17E\t%.17E\t%.17E\t%.17E\t%.17E\t%.17E\n",(cycle+(long)args[4]),pressureReduced,mpfr_get_d(avVolume,MPFR_RNDN),mpfr_get_d(dAvVolume,MPFR_RNDN),mpfr_get_d(avBoxMatrix[0],MPFR_RNDN),mpfr_get_d(dAvBoxMatrix[0],MPFR_RNDN),mpfr_get_d(avBoxMatrix[1],MPFR_RNDN),mpfr_get_d(dAvBoxMatrix[1],MPFR_RNDN),mpfr_get_d(avBoxMatrix[2],MPFR_RNDN),mpfr_get_d(dAvBoxMatrix[2],MPFR_RNDN),mpfr_get_d(avRho,MPFR_RNDN),mpfr_get_d(dAvRho,MPFR_RNDN),mpfr_get_d(avPacFrac,MPFR_RNDN),mpfr_get_d(dAvPacFrac,MPFR_RNDN),mpfr_get_d(s1111,MPFR_RNDN),mpfr_get_d(dS1111,MPFR_RNDN),mpfr_get_d(s1122,MPFR_RNDN),mpfr_get_d(dS1122,MPFR_RNDN),mpfr_get_d(s1212,MPFR_RNDN),mpfr_get_d(dS1212,MPFR_RNDN),mpfr_get_d(s2222,MPFR_RNDN),mpfr_get_d(dS2222,MPFR_RNDN),mpfr_get_d(s1112,MPFR_RNDN),mpfr_get_d(dS1112,MPFR_RNDN),mpfr_get_d(s1222,MPFR_RNDN),mpfr_get_d(dS1222,MPFR_RNDN),avNu,dAvNu,avNu2,dAvNu2,avB,dAvB,avMy,dAvMy,avE,dAvE,avNuDirect1,dAvNuDirect1,avNuDirect2,dAvNuDirect2,avBDirect1,dAvBDirect1,avBDirect2,dAvBDirect2,avMyDirect1,dAvMyDirect1,avMyDirect2,dAvMyDirect2);
            fprintf(fileExcelResults,"%.17E\t%.17E\t%.17E\t%.17E\t%.17E\t%.17E\t%.17E\n",pressureReduced,mpfr_get_d(avPacFrac,MPFR_RNDN),avNu,avNu2,avB,avMy,avE);

            if (!onlyMath[0]) {
                rho=N/volume; pacFrac=1.0/VcpPerParticle/rho;
                fprintf(fileConfigurations,"Rho: %.12E\tV/V_cp: %.12E\tPressureRed: %.12E\tRandStart: %.1f\tRandSteps: %.1f\tCycles: %ld\tEquilTime: %ldsec\tMeasuTime: %ldsec\n\n",rho,pacFrac,
                        pressureReduced,randomStartStep[0],randomStartStep[1],(cycle+(long)args[4]),(timeEquilibration-timeStart),(timeEnd-timeEquilibration));
                fprintf(fileConfigurations,"=====================\tConfigurations (data to load)\t=====================\n");
                fprintf(fileConfigurations,"%.1f %.1f %.17E %.17E %ld %.17E %.17E %.17E %.17E %.17E {",randomStartStep[0],randomStartStep[1],rho,pressureReduced,(cycle+(long)args[4]),boxMatrix[0][0],boxMatrix[1][1],boxMatrix[0][1],deltaR,deltaV);
                fprintf(fileConfigurationsList,"multimers[x_,y_,kI_]:={");
                for (int i=0;i<activeN;i++) {
                    fprintf(fileConfigurations,"%c[%.17E,%.17E,%.17E],",particles[i].type==0?'b':'w',particles[i].r[0],particles[i].r[1],particles[i].diameter);
                    fprintf(fileConfigurationsList,"%c[%.12E+x,%.12E+y,%.17E],",particles[i].type==0?'b':'w',particles[i].r[0],particles[i].r[1],particles[i].diameter);
                }
                fprintf(fileConfigurations,"{Opacity[0.2],Red,Polygon[{{0,0},{%.12E,%.12E},{%.12E,%.12E},{%.12E,%.12E}}]},{Opacity[0.2],Green,Disk[{%.12E,%.12E},%.12E]}}",boxMatrix[0][0],boxMatrix[1][0],boxMatrix[0][0]+boxMatrix[0][1],boxMatrix[1][0]+boxMatrix[1][1],boxMatrix[0][1],boxMatrix[1][1],particles[indexScanned].r[0],particles[indexScanned].r[1],neighRadius);
                fprintf(fileConfigurations,"\n==========================================\n\n\nboxMatrix[0][0]=%.12E, boxMatrix[1][1]=%.12E, boxMatrix[0][1]=boxMatrix[1][0]=%.12E",boxMatrix[0][0],boxMatrix[1][1],boxMatrix[0][1]);
                fprintf(fileConfigurationsList,"{Opacity[If[x==0 && y==0,0.4,0]],Red,Polygon[{{0,0},{%.12E,%.12E},{%.12E,%.12E},{%.12E,%.12E}}]},{Opacity[If[x==0 && y==0,0.4,0]],Green,Disk[{%.12E,%.12E},%.12E]}};\nconfigurationsList=Append[configurationsList,g[%.12E,multimers[%.12E,%.12E,kolorIndex=1],multimers[%.12E,%.12E,kolorIndex=1],multimers[%.12E,%.12E,kolorIndex=1],multimers[0,0,kolorIndex=1]]];\n",
                        boxMatrix[0][0],boxMatrix[1][0],boxMatrix[0][0]+boxMatrix[0][1],boxMatrix[1][0]+boxMatrix[1][1],boxMatrix[0][1],boxMatrix[1][1],particles[indexScanned].r[0],particles[indexScanned].r[1],neighRadius,pacFrac,boxMatrix[0][0],boxMatrix[1][0],boxMatrix[0][0]+boxMatrix[0][1],boxMatrix[1][0]+boxMatrix[1][1],boxMatrix[0][1],boxMatrix[1][1]);
            }

            /*for (int i=0;i<activeN;i++) {  //probDensDistFunMode - comment if not desired 6/7
                for (int j=0;j<activeN;j++) if (avEdgeDistance[i][j]!=0) fprintf(fileProbDensDistFunResults,"%d,%d,%.17E\n",i,j,avEdgeDistance[i][j]);
                delete [] avEdgeDistance[i];
            } delete [] avEdgeDistance;*/

            fclose(fileResults); fclose(fileExcelResults);
            if (!onlyMath[0]) {
                fclose(fileConfigurations); fclose(fileConfigurationsList);
            }
            //fclose(fileProbDensDistFunResults);  //probDensDistFunMode - comment if not desired 7/7
            printf("done\n\n");

            /*Clear of MPFR variables*/
            for (int i=0;i<10;i++) {
                mpfr_clear(avVolumeSet[i]);
                for (int j=0;j<3;j++) mpfr_clear(avBoxMatrixSet[i][j]);
                mpfr_clear(avRhoSet[i]);
                mpfr_clear(avPacFracSet[i]);
                mpfr_clear(e1111Set[i]); mpfr_clear(e1122Set[i]); mpfr_clear(e1212Set[i]);
                mpfr_clear(e2222Set[i]); mpfr_clear(e1112Set[i]); mpfr_clear(e1222Set[i]);
                mpfr_clear(eVSQSet[i]); mpfr_clear(eXXPeYYSQSet[i]); mpfr_clear(eXXMeYYSQSet[i]);
            }
            mpfr_clear(avVolume); mpfr_clear(avRho); mpfr_clear(avPacFrac);
            mpfr_clear(dAvVolume); mpfr_clear(dAvRho); mpfr_clear(dAvPacFrac);
            mpfr_clear(HxyHyx); mpfr_clear(HxxHyy); mpfr_clear(HxxHxy);
            mpfr_clear(Hxx2); mpfr_clear(Hyy2); mpfr_clear(mod0); mpfr_clear(mod1);
            mpfr_clear(e1111); mpfr_clear(e1122); mpfr_clear(e1212);
            mpfr_clear(e2222); mpfr_clear(e1112); mpfr_clear(e1222);
            mpfr_clear(eVSQ); mpfr_clear(eXXPeYYSQ); mpfr_clear(eXXMeYYSQ);
            mpfr_clear(dE1111); mpfr_clear(dE1122); mpfr_clear(dE1212);
            mpfr_clear(dE2222); mpfr_clear(dE1112); mpfr_clear(dE1222);
            mpfr_clear(dEVSQ); mpfr_clear(dEXXPeYYSQ); mpfr_clear(dEXXMeYYSQ);
            mpfr_clear(s1111); mpfr_clear(dS1111); mpfr_clear(s1122); mpfr_clear(dS1122);
            mpfr_clear(s1212); mpfr_clear(dS1212); mpfr_clear(s2222); mpfr_clear(dS2222);
            mpfr_clear(s1112); mpfr_clear(dS1112); mpfr_clear(s1222); mpfr_clear(dS1222);
            for (int j=0;j<3;j++) {
                mpfr_clear(avBoxMatrix[j]); mpfr_clear(dAvBoxMatrix[j]);
            }
        }




/////////////////////////////////////////////// PRZYGOTOWANIE DO KOLEJNEJ ITERACJI

        getNextArgument(arg,true);
        for (int i=0;i<2;i++) for (int j=0;j<2;j++) oldBoxMatrix[i][j]=boxMatrix[i][j];
        args[4]=0;
        loadedConfiguration=0;
        generatorStartPoint=0;
    }
    printf("\nTime for equilibrations: %ldsec, time for measurments: %ldsec, time for math: %ldsec.\n",timeEq,timeMe,timeMath);
    printf("\nSimulation has been completed.\n"); fflush(stdout);
}
