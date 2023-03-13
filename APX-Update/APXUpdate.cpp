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
extern long g_lReadTime;
extern long g_lStoreTime;

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
	myG oMainG;
	vector<vector<int> > vKG;
	//vector<vector<int> > vPKG;
	//vector<vector<pair<int, double> > > vGPos;
	vector<TPST_INDEX_EDGE> vQuery;
	char * pcSavePath = NULL;
	bool bStore = false;

    if (3 > argc)
    {
        printf("ERROR argc: %d\n",
           argc);
        ASSERT(0);
    }

    char * pcData = argv[1];
    char * pcUptData = argv[2];
    if (3 < argc)
    {
        pcSavePath = argv[3];
        bStore = true;
    }

    printf("start read file\n");
	gettimeofday(&tv, NULL);
	lStartTime = tv.tv_sec * 1000 + tv.tv_usec / 1000;
    file_readBasic(oMainG, vKG, pcData);
    oMainG.init();
    file_readQuery(pcUptData, vQuery);
	gettimeofday(&tv, NULL);
	lCurTime = tv.tv_sec * 1000 + tv.tv_usec / 1000;
	printf("read file Time=%ld ms\n", lCurTime - lStartTime);

	printf("get index k_max: %d pK_max: %d miniP: %lf accu: %lf\n",
        oMainG.m_iMaxK, oMainG.m_iMaxPK, oMainG.m_dMinP, oMainG.m_dAccu);

    printf("start query batch: %d\n", vQuery.size());
	gettimeofday(&tv, NULL);
	lStartTime = tv.tv_sec * 1000 + tv.tv_usec / 1000;

    deTruss_update(oMainG, vQuery, pcData, bStore, pcSavePath);

	gettimeofday(&tv, NULL);
	lCurTime = tv.tv_sec * 1000 + tv.tv_usec / 1000;
	long lUpdateTime = lCurTime - lStartTime - g_lStoreTime - g_lReadTime;

	printf("Update Time: %ld ms\n", lUpdateTime);
	printf("Vector init: %ld update: %ld\n", g_lInitCnt, g_lUptCnt);

	/* store */
    /*if (3 < argc)
    {
        char * pcIndexPath = argv[3];
        deTruss_save(oMainG, vPKG, vGPos);
    }*/
	if (bStore)
    {
        file_saveBasic(oMainG, vKG, pcSavePath);
    }
    return 0;
}
