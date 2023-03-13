/***************
input file output file function
****************/

#include <unistd.h>
#include <fcntl.h>

#include "common.h"
#include "myG.h"
#include "file.h"
#include "deTruss.h"

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
        char *szFileName
        map<int, vector<TPST_INDEX_EDGE> > &mpQuery
description:
        read file and fill the list g_initG
******************/
int file_readQuery(char *szFileName, vector<TPST_INDEX_EDGE> &query)
{
    FILE *fp = NULL;
    //int iTimeStamp = 0;
    int iTempX = 0;
    int iTempY = 0;
    char caBuffer[ONE_LINE_BUFFER] = {0};
    double fP = 0;
    TPST_INDEX_EDGE stNode = {0};

    fp = fopen(szFileName, "rt");
    if (NULL == fp)
    {
        printf("error no file: %s\n", szFileName);
        ASSERT(0);
    }
    while (fgets(caBuffer, ONE_LINE_BUFFER - 1, fp) > 0)
    {
        sscanf(caBuffer, "%d,%d,%lf", &iTempX, &iTempY, &(fP));
        if ((fP < 0) || (fP > 1))
        {
            printf("error %d,%d,%lf\n", iTempX, iTempY, fP);
            continue;
        }

        if (iTempX == iTempY)
        {
            continue;
        }

        stNode.x = iTempX + g_iOffset;
        stNode.y = iTempY + g_iOffset;
        stNode.p = fP;
        //mpQuery[iTimeStamp].push_back(stNode);
        query.push_back(stNode);
    }
    fclose(fp);
    return 0;
}
/*****************
input:
        char *szFileName
        map<int, vector<TPST_INDEX_EDGE> > &mpQuery
description:
        read file and fill the list g_initG
******************/
int file_saveQuery(char *szFileName, vector<TPST_INDEX_EDGE> &query)
{
    FILE *fp = NULL;
    int iTempX = 0;
    int iTempY = 0;
    char caBuffer[ONE_LINE_BUFFER] = {0};
    float fP = 0;

    fp = fopen(szFileName, "wt");
    if (NULL == fp)
    {
        printf("error no file: %s\n", szFileName);
        ASSERT(0);
    }
    for (auto atNode : query)
    {
        iTempX = atNode.x;
        iTempY = atNode.y;
        fP = atNode.p;

        fprintf(fp, "%d,%d,%.6f\n", iTempX, iTempY, fP);
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
        char *pcSaveFile
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
        myG &mpG,
        vector<TPST_INDEX_EDGE> &vDataE,
        vector<pair<int, int> > &vPos,  // <k, number>
        char *pcSavePath
description:
        save Index
******************/
int file_readBlock(myG &mpG, vector<TPST_INDEX_EDGE> &vDataE, vector<pair<int, int> > &vPos, char *pcSavePath)
{
    char pcFile[FILE_NAME_BUFFER] = {0};

    sprintf(pcFile, "%s/basic.pt", pcSavePath);

    char caFirstLine[FIRST_LINE] = {0};
    /* read */
    int fd = open(pcFile, O_RDONLY);
    ASSERT(-1 != fd);
    int ret = read(fd, caFirstLine, FIRST_LINE * sizeof(char));
    //printf("read file %s first size: %d\n", pcFile, ret);
    ASSERT(FIRST_LINE == ret);
    int iGE, iPos, iMaxK, iMaxP;
    double dMinP, dAccu;
    ret = sscanf(caFirstLine, "E: %d\nR: %d\nK: %d\nP: %d\nMP: %lf\nAC: %lf\n",
                 &iGE, &iPos, &iMaxK, &iMaxP, &dMinP, &dAccu);
    ASSERT(6 == ret);

    mpG.m_iMaxPK = iMaxK;
    mpG.m_dMinP = dMinP;
    mpG.m_dAccu = dAccu;
    vDataE.resize(iGE);
    vPos.resize(iPos);

    int uMaxBuffer = 100000000;

    long long llCurPos = 0;
    long long llTotoalSize = iGE * sizeof(TPST_INDEX_EDGE);
    uint32_t uStep = (uMaxBuffer < llTotoalSize)?uMaxBuffer:llTotoalSize;
    ASSERT(1 < uStep);
    char *pcData = (char *)(&(vDataE[0]));
    while (llCurPos < llTotoalSize)
    {
        ret = read(fd, pcData + llCurPos, uStep * sizeof(char));
        //printf("FILE block want: %d get: %d\n", uStep, ret);
        ASSERT(uStep == ret);
        llCurPos += uStep;

        uStep = uStep < (llTotoalSize - llCurPos)? uStep : (llTotoalSize - llCurPos);
    }
    /* second part */
    llCurPos = 0;
    llTotoalSize = iPos * sizeof(pair<int, int>);
    uStep = (uMaxBuffer < llTotoalSize)?uMaxBuffer:llTotoalSize;
    ASSERT(1 < uStep);
    pcData = (char *)(&(vPos[0]));
    while (llCurPos < llTotoalSize)
    {
        ret = read(fd, pcData + llCurPos, uStep * sizeof(char));
        //printf("FILE block want: %d get: %d\n", uStep, ret);
        ASSERT(uStep == ret);
        llCurPos += uStep;

        uStep = uStep < (llTotoalSize - llCurPos)? uStep : (llTotoalSize - llCurPos);
    }
    close(fd);

    ASSERT(!vPos.empty());
    mpG.m_iMaxK = vPos.back().first;

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
int file_readBasic(myG &mpG, vector<vector<int> > &vKG, char *pcSavePath)
{
    vector<int> vEdges;
    vector<TPST_INDEX_EDGE> vDataE;
    /* k, number */
    vector<pair<int, int> > vPos;

    file_readBlock(mpG, vDataE, vPos, pcSavePath);

    vEdges.reserve(vDataE.size());
    /* save data to obG */
    for (int j = 0; j < vDataE.size(); ++j)
    {
        int x = vDataE[j].x + g_iOffset;
        int y = vDataE[j].y + g_iOffset;
        float p = vDataE[j].p;
        mpG.add(x, y, p);
        vEdges.push_back(mpG.m_iMaxEId);
        //printf("FILE add %d (%d, %d)\n", mpG.m_iMaxEId, x, y);
    }
    //mpG.init();

    vKG.resize(mpG.m_iMaxK + 1);
    TPST_E *pstNode = NULL;
    int i = 0;
    for (auto atPos : vPos)
    {
        int k = atPos.first;
        int cnt = atPos.second;
        //printf("FILE k: %d cnt: %d\n", k, cnt);
        for (int j = 0; j < cnt; ++j)
        {
            int eid = vEdges[i + j];
            pstNode = mpG.findNode(eid);
            pstNode->iTrussness = k;
            vKG[k].push_back(eid);
            //printf("FILE init %d k: %d\n", eid, pstNode->iTrussness);
        }
        i += cnt;
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
int file_readIndex(myG &mpG, int iCurK, vector<int> &vEdges, vector<pair<int, double> > &vPos, char *pcSavePath)
{
    char pcFile[FILE_NAME_BUFFER] = {0};
    sprintf(pcFile, "%s/%d.pt", pcSavePath, iCurK);

    vEdges.clear();
    vPos.clear();
    char caFirstLine[FIRST_LINE] = {0};
    /* read */
    int fd = open(pcFile, O_RDONLY);
    if (-1 == fd)
    {
        /* no such (k, \gamma )-truss */
        return 0;
    }
    int ret = read(fd, caFirstLine, FIRST_LINE * sizeof(char));
    //printf("read file %s\n", caFirstLine);
    ASSERT(FIRST_LINE == ret);
    int iGetK, iGE, iAllR;
    ret = sscanf(caFirstLine, "k: %d\nE: %d\nR: %d\n", &iGetK, &iGE, &iAllR);
    ASSERT(3 == ret);
    ASSERT(iCurK == iGetK);

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
    vEdges.reserve(iGE);
    TPST_E *pstNode = NULL;
    for (int i = 0; i < iGE; ++i)
    {
        int x = vDataE[i].x + g_iOffset;
        int y = vDataE[i].y + g_iOffset;
        //printf("File want: (%d, %d)\n", x, y);
        pstNode = mpG.findNode(x, y);
        ASSERT(NULL != pstNode);
        vEdges.push_back(pstNode->eid);
    }
    int i = 0;
    for (auto atPos : vPos)
    {
        int cnt = atPos.first;
        double dGamma = atPos.second;
        int iCurR = DE_TRUSS_P2R(dGamma, mpG.m_dMinP, mpG.m_dAccu);
        for (int j = 0; j < cnt; ++j)
        {
            int eid = vEdges[i + j];
            pstNode = mpG.findNode(eid);
            if (iCurK >= pstNode->vSky.size())
            {
                pstNode->vSky.resize(iCurK + 1);
            }
            pstNode->vSky[iCurK] = iCurR;
        }
        i += cnt;
    }

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
int file_readIndex(myG &mpG, vector<vector<int> > &vKG, vector<vector<int> > &vPKG, vector<vector<pair<int, double> > > &vRetPos, char *pcSavePath)
{
    char caFileName[FILE_NAME_BUFFER] = {0};

    sprintf(caFileName, "%s/basic.pt", pcSavePath);
    fileReadBasic(mpG, vKG, caFileName);
    vPKG.resize(mpG.m_iMaxPK + 1);
    vRetPos.resize(mpG.m_iMaxPK + 1);
    for (int iCurK = mpG.m_iMaxPK; iCurK > 2; --iCurK)
    {
        sprintf(caFileName, "%s/%d.pt", pcSavePath, iCurK);
        fileReadIndex(mpG, vPKG[iCurK], vRetPos[iCurK], caFileName);
    }
    return 0;
}
#endif
