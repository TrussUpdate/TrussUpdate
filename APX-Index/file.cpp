/***************
input file output file function
****************/

#include <unistd.h>
#include <fcntl.h>

#include "common.h"
#include "myG.h"
#include "file.h"

extern long double g_ldAccu;
int g_iOffset = 1;

/*****************
input:
        char *pcFileName
        long double ldMinP
        myG &mpG
description:
        read file and fill the list g_initG
******************/
int file_fillG(char *pcFileName, myG &mpG)
{
    FILE *fp = NULL;
    int iTempX = 0;
    int iTempY = 0;
    char caBuffer[ONE_LINE_BUFFER] = {0};
    long double fP = 0;

    fp = fopen(pcFileName, "rt");
    if (NULL == fp)
    {
        printf("error no file: %s\n", pcFileName);
        ASSERT(0);
    }
    while (fgets(caBuffer, ONE_LINE_BUFFER - 1, fp) > 0)
    {
        sscanf(caBuffer, "%d,%d,%Lf", &iTempX, &iTempY, &(fP));
        if ((fP < 0) || (fP > 1))
        {
            printf("error %d,%d,%Lf\n", iTempX, iTempY, fP);
            continue;
        }

        if (iTempX == iTempY)
        {
            continue;
        }

        mpG.add(iTempX + g_iOffset, iTempY + g_iOffset, fP);
    }
    fclose(fp);
    return 0;
}
/*****************
input:
        myG &mpG
        int k
        vector<int> vEdges #eid
        vector<pair<int, double> > vPos # <count, r>
        char *pcSavePath
description:
        save Index
******************/
int file_saveBasic(myG &obG, vector<vector<int> > &vKG, char *pcSavePath)
{
    char pcSaveFile[FILE_NAME_BUFFER] = {0};

    ASSERT(!vKG.empty());
    sprintf(pcSaveFile, "%s/basic.pt", pcSavePath);

    vector<int> vEdges;
    /* k, number */
    vector<pair<int, int> > vPos;
    vEdges.reserve(obG.m_iMaxEId);

    for (int iCurK = 2; iCurK < vKG.size(); ++iCurK)
    {
        for (int iEid : vKG[iCurK])
        {
            vEdges.push_back(iEid);
        }
        if (0 < vKG[iCurK].size())
        {
            vPos.push_back({iCurK, vKG[iCurK].size()});
        }
    }

    vector<TPST_INDEX_EDGE> vSaveE(vEdges.size());
    char caFirstLine[FIRST_LINE] = {0};

    sprintf(caFirstLine, "E: %d\nR: %d\nK: %d\nP: %d\nMP: %lf\nAC: %lf\n",
            vEdges.size(), vPos.size(), obG.m_iMaxPK, obG.m_iMaxPId, obG.m_dMinP, obG.m_dAccu);

    int uMaxBuffer = 100000000;

    TPST_E *pstNode = NULL;
    for (int i = 0; i < vEdges.size(); ++i)
    {
        pstNode = obG.findNode(vEdges[i]);
        if (vEdges[i] != pstNode->eid)
        {
            /* removed */
            continue;
        }
        vSaveE[i].x = pstNode->paXY.first - g_iOffset;
        vSaveE[i].y = pstNode->paXY.second - g_iOffset;
        vSaveE[i].p = pstNode->p;
    }

    /* write */
    int fd = open(pcSaveFile, O_CREAT|O_WRONLY|O_TRUNC, 0644);

    int ret = write(fd, caFirstLine, 100 * sizeof(char));
    ASSERT(100 == ret);

    char * pcData = (char *)&(vSaveE[0]);
    long long llCurPos = 0;
    long long llTotoalSize = vEdges.size() * sizeof(TPST_INDEX_EDGE);
    int uStep = (uMaxBuffer < llTotoalSize)?uMaxBuffer:llTotoalSize;
    while (llCurPos < llTotoalSize)
    {
        if (uStep <= uMaxBuffer)
        {
            ret = write(fd, pcData + llCurPos, uStep * sizeof(char));
            ASSERT(uStep == ret);
            llCurPos += uStep;

            uStep = uStep < (llTotoalSize - llCurPos)? uStep : (llTotoalSize - llCurPos);
        }
        else
        {
            uStep = uMaxBuffer;
            ASSERT(1 < uStep);
        }
    }
    /* second part */

    pcData = (char *)&(vPos[0]);
    llCurPos = 0;
    llTotoalSize = vPos.size() * sizeof(pair<int, int>);
    uStep = (uMaxBuffer < llTotoalSize)?uMaxBuffer:llTotoalSize;
    while (llCurPos < llTotoalSize)
    {
        if (uStep <= uMaxBuffer)
        {
            ret = write(fd, pcData + llCurPos, uStep * sizeof(char));
            ASSERT(uStep == ret);
            llCurPos += uStep;

            uStep = uStep < (llTotoalSize - llCurPos)? uStep : (llTotoalSize - llCurPos);
        }
        else
        {
            uStep = uMaxBuffer;
            ASSERT(1 < uStep);
        }
    }

    close(fd);
    pcData = NULL;

    return 0;
}
/*****************
input:
        myG &mpG
        int k
        vector<int> vEdges #eid
        vector<pair<int, double> > vPos # <count, r>
        char *pcSavePath
description:
        save Index
******************/
int file_saveIndex(myG &mpG, int k, vector<int> &vEdges, vector<pair<int, double> > &vPos, char *pcSavePath)
{
    char pcSaveFile[FILE_NAME_BUFFER] = {0};
    sprintf(pcSaveFile, "%s/%d.pt", pcSavePath, k);

    vector<TPST_INDEX_EDGE> vSaveE(vEdges.size());
    char caFirstLine[FIRST_LINE] = {0};

    ASSERT(k > 2);

    sprintf(caFirstLine, "k: %d\nE: %d\nR: %d\n", k, vEdges.size(), vPos.size());

    int uMaxBuffer = 100000000;

    TPST_E *pstNode = NULL;
    for (int i = 0; i < vEdges.size(); ++i)
    {
        pstNode = mpG.findNode(vEdges[i]);
        if (vEdges[i] != pstNode->eid)
        {
            /* removed */
            continue;
        }
        vSaveE[i].x = pstNode->paXY.first - g_iOffset;
        vSaveE[i].y = pstNode->paXY.second - g_iOffset;
    }

    /* write */
    int fd = open(pcSaveFile, O_CREAT|O_WRONLY|O_TRUNC, 0644);

    int ret = write(fd, caFirstLine, 100 * sizeof(char));
    ASSERT(100 == ret);

    char * pcData = (char *)&(vSaveE[0]);
    long long llCurPos = 0;
    long long llTotoalSize = vEdges.size() * sizeof(TPST_INDEX_EDGE);
    int uStep = (uMaxBuffer < llTotoalSize)?uMaxBuffer:llTotoalSize;
    while (llCurPos < llTotoalSize)
    {
        if (uStep <= uMaxBuffer)
        {
            ret = write(fd, pcData + llCurPos, uStep * sizeof(char));
            ASSERT(uStep == ret);
            llCurPos += uStep;

            uStep = uStep < (llTotoalSize - llCurPos)? uStep : (llTotoalSize - llCurPos);
        }
        else
        {
            uStep = uMaxBuffer;
            ASSERT(1 < uStep);
        }
    }
    /* second part */

    pcData = (char *)&(vPos[0]);
    llCurPos = 0;
    llTotoalSize = vPos.size() * sizeof(pair<int, double>);
    uStep = (uMaxBuffer < llTotoalSize)?uMaxBuffer:llTotoalSize;
    while (llCurPos < llTotoalSize)
    {
        if (uStep <= uMaxBuffer)
        {
            ret = write(fd, pcData + llCurPos, uStep * sizeof(char));
            ASSERT(uStep == ret);
            llCurPos += uStep;

            uStep = uStep < (llTotoalSize - llCurPos)? uStep : (llTotoalSize - llCurPos);
        }
        else
        {
            uStep = uMaxBuffer;
            ASSERT(1 < uStep);
        }
    }

    close(fd);
    pcData = NULL;

    return 0;
}
#if 0
/*****************
input:
        myG &mpG
        int k
        vector<int> vEdges #eid
        vector<pair<int, double> > vPos # <count, r>
        char *pcSaveFile
description:
        save Index
******************/
int file_saveIndex(myG &mpG, vector<vector<int> > &vKG, vector<vector<int> > &vPKG, vector<vector<pair<int, double> > > &vRetPos, char *pcSavePath)
{
    char caFileName[FILE_NAME_BUFFER] = {0};

    ASSERT(!vKG.empty());
    sprintf(caFileName, "%s/basic.pt", pcSavePath);
    fileSaveBasic(mpG, vKG, caFileName);

    for (int iCurK = 3; iCurK < vPKG.size(); ++iCurK)
    {
        sprintf(caFileName, "%s/%d.pt", pcSavePath, iCurK);
        fileSaveIndex(mpG, iCurK, vPKG[iCurK], vRetPos[iCurK], caFileName);
    }
    return 0;
}
#endif
/*****************
input:
        myG &mpG
        int k
        vector<int> vEdges #eid
        vector<pair<int, double> > vPos # <count, r>
        char *pcSaveFile
description:
        save Index
******************/
int fileRead2T(myG &mpG, vector<int> &vEdges, vector<pair<int, double> > &vPos, char *pcFile)
{
    vEdges.clear();
    vPos.clear();
    char caFirstLine[FIRST_LINE] = {0};
    /* read */
    int fd = open(pcFile, O_RDONLY);
    ASSERT(-1 != fd);
    int ret = read(fd, caFirstLine, FIRST_LINE * sizeof(char));
    //printf("read file %s\n", caFirstLine);
    ASSERT(FIRST_LINE == ret);
    int iCurK, iGE, iAllR, iMaxK, iMaxP;
    ret = sscanf(caFirstLine, "k: %d\nE: %dR: %d\nK: %d\nP: %d\n", &iCurK, &iGE, &iAllR, &iMaxK, &iMaxP);
    ASSERT(5 == ret);
    ASSERT(2 == iCurK);

    //mpG.m_iMaxEId = iGE;
    //mpG.m_iMaxPId = iMaxP;
    mpG.m_iMaxK = iMaxK;
    vector<TPST_INDEX_EDGE> vDataE(iGE);
    vPos.resize(iAllR);

    int uMaxBuffer = 100000000;

    long long llCurPos = 0;
    long long llTotoalSize = iGE * sizeof(TPST_INDEX_EDGE);
    uint32_t uStep = (uMaxBuffer < llTotoalSize)?uMaxBuffer:llTotoalSize;
    ASSERT(1 < uStep);
    char *pcData = (char *)(&(vDataE[0]));
    while (llCurPos < llTotoalSize)
    {
        ret = read(fd, pcData + llCurPos, uStep * sizeof(char));
        ASSERT(uStep == ret);
        llCurPos += uStep;

        uStep = uStep < (llTotoalSize - llCurPos)? uStep : (llTotoalSize - llCurPos);
    }
    /* second part */
    llCurPos = 0;
    llTotoalSize = iAllR * sizeof(pair<int, double>);
    uStep = (uMaxBuffer < llTotoalSize)?uMaxBuffer:llTotoalSize;
    ASSERT(1 < uStep);
    pcData = (char *)(&(vPos[0]));
    while (llCurPos < llTotoalSize)
    {
        ret = read(fd, pcData + llCurPos, uStep * sizeof(char));
        ASSERT(uStep == ret);
        llCurPos += uStep;

        uStep = uStep < (llTotoalSize - llCurPos)? uStep : (llTotoalSize - llCurPos);
    }
    close(fd);

    /* save data to obG */
    TPST_E *pstNode = NULL;
    int i = 0;
    for (auto atPos : vPos)
    {
        for (int j = 0; j < atPos.first; ++j)
        {
            int x = vDataE[i + j].x + g_iOffset;
            int y = vDataE[i + j].y + g_iOffset;
            mpG.add(x, y, atPos.second);
            vEdges.push_back(mpG.m_iMaxEId);
        }
        i += atPos.first;
    }
    mpG.init();

    return 0;
}
/*****************
input:
        myG &mpG
        int k
        vector<int> vEdges #eid
        vector<pair<int, double> > vPos # <count, r>
        char *pcSaveFile
description:
        save Index
******************/
int fileReadIndex(myG &mpG, vector<int> &vEdges, vector<pair<int, double> > &vPos, char *pcFile)
{
    vEdges.clear();
    vPos.clear();
    char caFirstLine[FIRST_LINE] = {0};
    /* read */
    int fd = open(pcFile, O_RDONLY);
    ASSERT(-1 != fd);
    int ret = read(fd, caFirstLine, FIRST_LINE * sizeof(char));
    //printf("read file %s\n", caFirstLine);
    ASSERT(FIRST_LINE == ret);
    int iCurK, iGE, iAllR;
    ret = sscanf(caFirstLine, "k: %d\nE: %d\nR: %d\n", &iCurK, &iGE, &iAllR);
    ASSERT(3 == ret);

    vector<TPST_INDEX_EDGE> vDataE(iGE);
    vPos.resize(iAllR);

    int uMaxBuffer = 100000000;

    long long llCurPos = 0;
    long long llTotoalSize = iGE * sizeof(TPST_INDEX_EDGE);
    uint32_t uStep = (uMaxBuffer < llTotoalSize)?uMaxBuffer:llTotoalSize;
    ASSERT(1 < uStep);
    char *pcData = (char *)(&(vDataE[0]));
    while (llCurPos < llTotoalSize)
    {
        ret = read(fd, pcData + llCurPos, uStep * sizeof(char));
        ASSERT(uStep == ret);
        llCurPos += uStep;

        uStep = uStep < (llTotoalSize - llCurPos)? uStep : (llTotoalSize - llCurPos);
    }
    /* second part */
    llCurPos = 0;
    llTotoalSize = iAllR * sizeof(pair<int, double>);
    uStep = (uMaxBuffer < llTotoalSize)?uMaxBuffer:llTotoalSize;
    ASSERT(1 < uStep);
    pcData = (char *)(&(vPos[0]));
    while (llCurPos < llTotoalSize)
    {
        ret = read(fd, pcData + llCurPos, uStep * sizeof(char));
        ASSERT(uStep == ret);
        llCurPos += uStep;

        uStep = uStep < (llTotoalSize - llCurPos)? uStep : (llTotoalSize - llCurPos);
    }
    close(fd);

    /* save data to vector */
    TPST_E *pstNode = NULL;
    for (int i = 0; i < iGE; ++i)
    {
        int x = vDataE[i].x + g_iOffset;
        int y = vDataE[i].y + g_iOffset;
        pstNode = mpG.findNode(x, y);
        ASSERT(NULL != pstNode);
        vEdges.push_back(pstNode->eid);
    }

    return 0;
}
/*****************
input:
        myG &mpG
        int k
        vector<int> vEdges #eid
        vector<pair<int, double> > vPos # <count, r>
        char *pcSaveFile
description:
        save Index
******************/
int file_printIndex(char *pcFile)
{
    char caFirstLine[FIRST_LINE] = {0};
    /* read */
    int fd = open(pcFile, O_RDONLY);
    ASSERT(-1 != fd);
    int ret = read(fd, caFirstLine, FIRST_LINE * sizeof(char));
    //printf("read file %s\n", caFirstLine);
    ASSERT(FIRST_LINE == ret);
    int iCurK, iGE, iAllR;
    ret = sscanf(caFirstLine, "k: %d\nE: %d\nR: %d\n", &iCurK, &iGE, &iAllR);
    ASSERT(3 == ret);

    printf("%s", caFirstLine);

    vector<TPST_INDEX_EDGE> vDataE(iGE);
    vector<pair<int, double> > vPos(iAllR);

    int uMaxBuffer = 100000000;

    long long llCurPos = 0;
    long long llTotoalSize = iGE * sizeof(TPST_INDEX_EDGE);
    //uint32_t uStep = min(uMaxBuffer, llTotoalSize);
    uint32_t uStep = (uMaxBuffer < llTotoalSize)?uMaxBuffer:llTotoalSize;
    ASSERT(1 < uStep);
    char *pcData = (char *)(&(vDataE[0]));
    while (llCurPos < llTotoalSize)
    {
        ret = read(fd, pcData + llCurPos, uStep * sizeof(char));
        ASSERT(uStep == ret);
        llCurPos += uStep;

        uStep = uStep < (llTotoalSize - llCurPos)? uStep : (llTotoalSize - llCurPos);
    }
    /* second part */
    llCurPos = 0;
    llTotoalSize = iAllR * sizeof(pair<int, double>);
    uStep = (uMaxBuffer < llTotoalSize)?uMaxBuffer:llTotoalSize;
    ASSERT(1 < uStep);
    pcData = (char *)(&(vPos[0]));
    while (llCurPos < llTotoalSize)
    {
        ret = read(fd, pcData + llCurPos, uStep * sizeof(char));
        ASSERT(uStep == ret);
        llCurPos += uStep;

        uStep = uStep < (llTotoalSize - llCurPos)? uStep : (llTotoalSize - llCurPos);
    }
    close(fd);

    /* print */
    for (int i = 0; i < iGE; ++i)
    {
        int x = vDataE[i].x;
        int y = vDataE[i].y;
        printf("%d %d\n", x, y);
    }
    for (int i = 0; i < iAllR; ++i)
    {
        printf("%d %.3f\n", vPos[i].first, (float)(vPos[i].second));
    }

    return 0;
}
/*****************
input:
        myG &mpG
        int k
        vector<int> vEdges #eid
        vector<pair<int, double> > vPos # <count, r>
        char *pcSaveFile
description:
        save Index
******************/
int file_readIndex(myG &mpG, vector<vector<int> > &vKG, vector<vector<pair<int, double> > > &vRetPos, char *pcSavePath)
{
    char caFileName[FILE_NAME_BUFFER] = {0};

    int iCurK = 2;
    sprintf(caFileName, "%s/%d.pt", pcSavePath, iCurK);
    vKG.resize(3);
    vRetPos.resize(3);
    fileRead2T(mpG, vKG[iCurK], vRetPos[iCurK], pcSavePath);
    vKG.resize(mpG.m_iMaxK + 1);
    vRetPos.resize(mpG.m_iMaxK + 1);

    for (iCurK = 3; iCurK < mpG.m_iMaxK; ++iCurK)
    {
        sprintf(caFileName, "%s/%d.pt", pcSavePath, iCurK);
        fileReadIndex(mpG, vKG[iCurK], vRetPos[iCurK], pcSavePath);
    }
    return 0;
}
