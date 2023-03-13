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
    /* only simple member variables are useful */
	myG oMainG;
	/* edges information */
    vector<TPST_INDEX_EDGE> vDataE;
    /* useless */
    vector<pair<int, int> > vPos;
	//TPST_INDEX_EDGE stNode = {0};
	bool bStoreG = false;

    if (6 > argc)
    {
        printf("ERROR argc: %d\n",
           argc);
        ASSERT(0);
    }

    char * pcData = argv[1];
    char * pcQueryPath = argv[2];
    int iSample = 0;
    float fDelta = 0;
    int iBatch = 0;
    sscanf(argv[3], "%d", &iSample);
    sscanf(argv[4], "%f", &fDelta);
    sscanf(argv[5], "%d", &iBatch);
	char * pcSavePath = NULL;
    if (6 < argc)
    {
        pcSavePath = argv[6];
        bStoreG = true;
    }

    file_readBlock(oMainG, vDataE, vPos, pcData);
    printf("INIT done\n");

    /* sample */
    vector<int> vEdges(vDataE.size());
    std::iota(vEdges.begin(), vEdges.end(), 1);
    srand (static_cast <unsigned> (time(0)));
    random_shuffle(vEdges.begin(), vEdges.end());
    printf("shuffle done\n");

    for (int iCurBat = 0; iCurBat < iBatch; ++iCurBat)
    {
        char caQuery[100] = {0};
        sprintf(caQuery, "%s/query_%d_%d_%d.txt", pcQueryPath, (int)(fDelta * 10), iSample, iCurBat + 1);
        printf("save file: %s\n", caQuery);
        vector<TPST_INDEX_EDGE> vQuery;
        int iIdStart = rand() % (vEdges.size() - iSample);
        for (int i = 0; i < iSample; ++i)
        {
            int iEid = vEdges[iIdStart + i];
            TPST_INDEX_EDGE stNode = vDataE[iEid];

            //printf("sample %d (%d, %d) %lf\n", iEid, stNode.x, stNode.y, stNode.p);

            float fStart = 0;
            float fEnd = 1;

            if (fDelta > 1)
            {
                /* [self, 1] */
                fStart = stNode.p;
            }
            else if (fDelta < -1)
            {
                /* [0, self] */
                fEnd = stNode.p;
            }
            else if (fDelta > 0)
            {
                /* [self, self + delta] */
                fStart = stNode.p;
                fEnd = min(stNode.p + fDelta, (float)1.0);
            }
            else
            {
                /* [self - delta, self] */
                fStart = max((float)0.0, stNode.p + fDelta);
                fEnd = stNode.p;
            }

            if (0 >= fEnd - fStart)
            {
                /* invalid */
                // printf("sample %d [%f, %f)\n", i, fStart, fEnd);
                --i;
                iIdStart = rand() % (vEdges.size() - iSample);
                continue;
            }

            float fOldP = stNode.p;

            stNode.p = fStart + static_cast <float> (rand()) / (static_cast <float> (RAND_MAX / (fEnd - fStart)));
            //printf("random (%d, %d) old p: %d p: %.6f\n", stNode.x, stNode.y, fOldP, stNode.p);

            vQuery.push_back(stNode);
        }

        file_saveQuery(caQuery, vQuery);
    }
    printf("sample done\n");


	/*if (bStoreG)
    {
        file_saveBasic(oMainG, vKG, pcSavePath);
    }*/
    return 0;
}
