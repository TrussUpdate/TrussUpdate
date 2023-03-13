/***************
truss decomposition of probabilistic graphs
first step
generate several k and r
main cpp
****************/
#include <math.h>
#include <sys/time.h>
#include "common.h"
#include "myG.h"
#include "file.h"
#include "deTruss.h"
#include "calPek.h"

extern long g_lInitCnt;
extern long g_lUptCnt;
extern long g_lStoreTime;

/*****************
input:
        char *szFileName data file name
        long double ldMinP
description:
        global init
******************/
int init(char *szFileName, myG &oMainG)
{
    printf("read data\n");
    file_fillG(szFileName, oMainG);
    printf("start init\n");
    oMainG.init();
    printf("init done\n");
    return 0;
}

/*****************
input:
        Graph file name
        output data file name
        long double minimum possibility
        long double accuracy
        int kSpace
        dataset label
description:
        main function
******************/
int main(int argc, char *argv[])
{
	struct timeval tv;
	long lStartTime = 0;
	long lCurTime = 0;
	long lDelGTime = 0;
	long lSaveBaseTime = 0;
	double ldMinP = 0;
	double ldAccu = 0;
	myG oMainG;
	bool bStore = false;
	char * pcIndexPath = NULL;

    if (4 > argc)
    {
        printf("ERROR argc: %d\n",
           argc);
        ASSERT(0);
    }

    char * pcData = argv[1];
    sscanf(argv[2], "%lf", &ldMinP);
    sscanf(argv[3], "%lf", &ldAccu);
    if (4 < argc)
    {
        pcIndexPath = argv[4];
        bStore = true;
    }
    printf("get minimum p: %lf accuracy: %lf data: %s\n",
           ldMinP, ldAccu, pcData);
    ASSERT(( 0 <= ldMinP ) && (1 > ldMinP));
    ASSERT(( 0 < ldAccu ) && (1 > ldAccu));

    oMainG.m_dMinP = ldMinP;
    oMainG.m_dAccu = ldAccu;

    printf("start read file\n");
	gettimeofday(&tv, NULL);
	lStartTime = tv.tv_sec * 1000 + tv.tv_usec / 1000;
    init(pcData, oMainG);
	gettimeofday(&tv, NULL);
	lCurTime = tv.tv_sec * 1000 + tv.tv_usec / 1000;
	printf("read file Time=%ld ms\n", lCurTime - lStartTime);

    printf("start deTruss\n");
	lStartTime = lCurTime;
	/** first deTruss **/
	vector<vector<int> > vKG;
    deTruss_Detm(oMainG, vKG);
	gettimeofday(&tv, NULL);
	lCurTime = tv.tv_sec * 1000 + tv.tv_usec / 1000;
	long lDetrussTime = lCurTime - lStartTime;
	printf("deTruss Time=%ld ms\n", lDetrussTime);

    printf("start pTruss\n");
	lStartTime = lCurTime;
	/** second pTruss **/
	vector<vector<int> > vPKG;
    deTruss_ByP(oMainG, vPKG, ldMinP, ldAccu);
	gettimeofday(&tv, NULL);
	lCurTime = tv.tv_sec * 1000 + tv.tv_usec / 1000;
	long lPtrussTime = lCurTime - lStartTime;
	printf("pTruss Time=%ld ms\n", lPtrussTime);

    if (bStore)
    {
        gettimeofday(&tv, NULL);
        lStartTime = tv.tv_sec * 1000 + tv.tv_usec / 1000;
        file_saveBasic(oMainG, vKG, pcIndexPath);
        gettimeofday(&tv, NULL);
        lCurTime = tv.tv_sec * 1000 + tv.tv_usec / 1000;
        lSaveBaseTime += lCurTime;
    }
    vKG.clear();

    printf("start del\n");
	lStartTime = lCurTime;
	/** delete G **/
    deTruss_PT(oMainG, vPKG, ldMinP, ldAccu, bStore, pcIndexPath);
	gettimeofday(&tv, NULL);
	lCurTime = tv.tv_sec * 1000 + tv.tv_usec / 1000;
	lDelGTime = lCurTime - lStartTime;
	printf("peel all G Time=%ld ms\n", lDelGTime - g_lStoreTime);

	printf("Vector init: %ld update: %ld\n", g_lInitCnt, g_lUptCnt);

	printf("first step Time: %ld ms\n", lDetrussTime + lDelGTime - g_lStoreTime);
	/* store */
    /*if (4 < argc)
    {
        char * pcIndexPath = argv[4];
        file_saveIndex(oMainG, vKG, vPKG, vGPos, pcIndexPath);
    }*/
    return 0;
}
