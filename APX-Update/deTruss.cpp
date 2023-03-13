/***************
truss decomposition
****************/
#include <sys/time.h>
#include <queue>

#include "common.h"
#include "myG.h"
#include "deTruss.h"
#include "calPek.h"
#include "file.h"

vector<int> g_vLfE;
vector<int> g_vRtE;
vector<double> g_vLfP;
vector<double> g_vRtP;
bool g_bLock;

long g_lReadTime;
long g_lStoreTime;
/*****************
input:
        vector <long double> vfDtE
        long double ldP
        long double ldAccu
output:
        real support number
description:
        set key
******************/
int getR(vector <long double> &vfDtE, double ldP, double ldMiniP, double ldAccu, int iDesK)
{
	long double sum = 0;

    if(vfDtE.size() - 1 + 2 < iDesK)
    {
        return 0;
    }

	if (iDesK - 2 < vfDtE.size() / 2)
    {
        /* small k */
        for (int i = 0; i < iDesK - 2; ++i)
        {
            sum += vfDtE[i];
        }
        sum = 1 - sum;
    }
    else
    {
        for (int i = vfDtE.size() - 1; i >= iDesK - 2; --i)
        {
            sum += vfDtE[i];
        }
    }

    sum *= ldP;

    /*if (sum < ldMiniP)
    {
        return 0;
    }*/


    // DEBUG
    /*for (auto atR : vfDtE)
    {
        printf("%Lf ", atR);
    }
    printf("\nGET_R p: %f real: %Lf number: %d\n", ldP, sum, (int)(sum / ldAccu));
    printf("GET_R ldMiniP: %lf accu: %lf size: %d k: %d\n", ldMiniP, ldAccu, vfDtE.size(), iDesK);*/
    //ASSERT(0);

    return DE_TRUSS_P2R(sum, ldMiniP, ldAccu);
}

#if 0
/*****************
input:
        vector <long double> vfDtE
        long double ldP
        long double ldAccu
output:
        real support number
description:
        set key
******************/
int getSup(vector <long double> &vfDtE, double ldP, double ldMiniP)
{
	long double sum = 0;
	int iSup = 0;

	ASSERT(ldP >= ldMiniP);

    for (iSup = vfDtE.size() - 1; iSup > -1; --iSup)
    {
        sum += vfDtE[iSup];
        if (sum * ldP >= ldMiniP)
        {
            break;
        }
    }

    if (iSup < 0)
    {
        return 0;
    }

    return iSup;
}

/*****************
input:
        vector <long double> vfDtE
        long double ldP
        long double ldAccu
output:
        real support number
description:
        set key
******************/
double recoverR(int iCurR, double ldAccu)
{
    if (0 == iCurR)
    {
        return 0;
    }
    return (iCurR * ldAccu);
}
/*****************
input:
        vector <long double> vfDtE
        long double ldP
        long double ldAccu
output:
        real support number
description:
        set key
******************/
double PtoR(myG &obG, vector<int> &vEdges, vector<pair<int, double> > &vPos)
{
    //printf("PtoR start max eid: %d\n", obG.m_iMaxEId);
    vEdges.clear();
    vPos.clear();
    for (int iEid = 1; iEid <= obG.m_iMaxEId; ++iEid)
    {
        if (iEid != obG.findNode(iEid)->eid)
        {
            /* removed */
            continue;
        }
        vEdges.push_back(iEid);
    }
    //printf("PtoR edges: %d\n", vEdges.size());
    sort(vEdges.begin(), vEdges.end(), [&obG](const auto e1, const auto e2) {
         return (obG.findNode(e1)->p < obG.findNode(e2)->p);
      });

    /* debug */
    /*for (int iEid : vEdges)
    {
        printf("(%d, %.3f) ", iEid, (float)obG.findNode(iEid)->p);
    }
    printf("PtoR edges: %d sort done\n", vEdges.size());*/

    /* count */
    int iCnt = 0;
    double ldCurV = 0;
    for (auto eid : vEdges)
    {
        if (ldCurV != obG.findNode(eid)->p)
        {
            if (0 < iCnt)
            {
                vPos.push_back({iCnt, ldCurV});
                iCnt = 0;
            }
            ldCurV = obG.findNode(eid)->p;
        }
        ++iCnt;
    }
    // last
    if (0 < iCnt)
    {
        vPos.push_back({iCnt, ldCurV});
    }
    return 0;
}

/*****************
input:
        myG &mapInitG
        LIST_DECOMP_G &lstDeG
        long double ldMiniP
output:
        LIST_DECOMP_G &lstDeG
description:
        truss decomposition
******************/
int deTruss_ByP(myG &obG, vector<vector<int> > &vKG, double ldMiniP, double ldAccu)
{
    vector<int> vERaw;
    /* init */
    vERaw.reserve(obG.m_iMaxEId);
    for (int iEid = 1; iEid <= obG.m_iMaxEId; ++iEid)
    {
        if (obG.findNode(iEid)->eid == iEid)
        {
            if ((obG.findNode(iEid)->p >= ldMiniP) && (2 < obG.findNode(iEid)->iTrussness))
            {
                vERaw.push_back(iEid);
            }
        }
    }
    printf("DETRUSS_P valid graph size: %d\n", vERaw.size());
    /* init vector */
    calPek_part(obG, vERaw);
    printf("DETRUSS_P init pek done\n");

    /* init support */
    int iMaxSup = 0;
    for (int iEid : vERaw)
    {
        TPST_E *pstNode = obG.findNode(iEid);
        pstNode->iSup = getSup(pstNode->vldPek, pstNode->p, ldMiniP);
        iMaxSup = max(iMaxSup, pstNode->iSup);
    }
    printf("DETRUSS_P max support: %d\n", iMaxSup);

    /* <eid, sup> */
    auto cmp = [](pair<int, int> left, pair<int, int> right) { return left.second > right.second; };
    priority_queue<pair<int, int>, vector<pair<int, int> >, decltype(cmp)> prQ(cmp);
    for (int iEid : vERaw) prQ.push({iEid, obG.findNode(iEid)->iSup});
    printf("DETRUSS_P init Q done\n");

    vector<int> vKT(obG.m_iMaxEId + 1, 0);
    vector<bool> removed(obG.m_iMaxEId + 1, false);
    vector<bool> queueFlag(obG.m_iMaxEId + 1, false);
    int iCurSup = 0;
    int iRmCnt = 0;
    vector<int> vWait(obG.m_iMaxEId + 1);
    // 2.2.1. process the edges layer by layer
    while (iRmCnt < vERaw.size())
    {
        vWait.clear();
        if (prQ.empty())
        {
            break;
        }
        iCurSup = prQ.top().second;
        ASSERT(iCurSup >= 0);
        while (!prQ.empty())
        {
            int iEid = prQ.top().first;

            if (prQ.top().second <= iCurSup)
            {
                if ((!removed[iEid]) && (!queueFlag[iEid]))
                  {
                      vWait.push_back(iEid);
                      queueFlag[iEid] = true;
                      vKT[iEid] = iCurSup + 2;
                      //ASSERT(iCurSup == obG.findNode(iEid)->iSup);
                  }
                /*printf("DETRUSS_P rm cnt: %d pop eid: %d support: %d\n",
                       iRmCnt, prQ.top().first, prQ.top().second);*/
                prQ.pop();
            }
            else
            {
                /* not now */
                break;
            }
        }

        if (vWait.empty())
        {
            continue;
        }

        /*printf("DETRUSS_P support: %d edges: %d\n",
               iCurSup, vWait.size());*/
        for (auto eid : vWait)
        {
            ASSERT(!obG.findNode(eid)->bRm);
            removed[eid] = true;
            queueFlag[eid] = false;
            ++iRmCnt;
            //printf("start find triangles\n");
            // find triangles containing the edge with ID eid
            g_vLfE.clear();
            g_vRtE.clear();
            obG.findNotRmNeb(obG.findNode(eid)->paXY.first, obG.findNode(eid)->paXY.second, g_vLfE, g_vRtE);

            /*printf("DETRUSS_P rm edge: %d support: %d neighbors: %d\n",
                   eid, iCurSup, g_vLfE.size());*/

            //printf("start rm triangles\n");
            for (int i = 0; i < g_vLfE.size(); ++i)
            {
                const int e1 = g_vLfE[i];
                const int e2 = g_vRtE[i];
                if (removed[e1] || removed[e2]) continue;

                if (obG.findNode(e1)->iSup > iCurSup)
                {
                    /*printf("DETRUSS_P rm begin eid: %d support: %d neighbors: %d\n",
                        e1, obG.findNode(e1)->iSup, obG.findNode(e1)->vldPek.size() - 1);*/
                    calPek_rm(obG, e1, eid, e2, obG.findNode(e1)->vldPek);
                    int iOldSup = obG.findNode(e1)->iSup;
                    obG.findNode(e1)->iSup = getSup(obG.findNode(e1)->vldPek, obG.findNode(e1)->p, ldMiniP);
                    if (obG.findNode(e1)->iSup < iOldSup)
                    {
                        //vKT[e1] = iCurSup + 2;
                        prQ.push({e1, max(obG.findNode(e1)->iSup, iCurSup)});
                    }
                    /*printf("DETRUSS_P rm end eid: %d support: %d neighbors: %d\n",
                        e1, obG.findNode(e1)->iSup, obG.findNode(e1)->vldPek.size() - 1);*/
                }
                if (obG.findNode(e2)->iSup > iCurSup)
                {
                    /*printf("DETRUSS_P rm begin eid: %d support: %d neighbors: %d\n",
                        e2, obG.findNode(e2)->iSup, obG.findNode(e2)->vldPek.size() - 1);*/
                    calPek_rm(obG, e2, eid, e1, obG.findNode(e2)->vldPek);

                    int iOldSup = obG.findNode(e2)->iSup;
                    obG.findNode(e2)->iSup = getSup(obG.findNode(e2)->vldPek, obG.findNode(e2)->p, ldMiniP);
                    if (obG.findNode(e2)->iSup < iOldSup)
                    {
                        //vKT[e2] = iCurSup + 2;
                        prQ.push({e2, max(obG.findNode(e2)->iSup, iCurSup)});
                    }
                    /*printf("DETRUSS_P rm end eid: %d support: %d neighbors: %d\n",
                        e2, obG.findNode(e2)->iSup, obG.findNode(e2)->vldPek.size() - 1);*/
                }
            }
            obG.findNode(eid)->bRm = true;
            //vEDone.push_back(eid);
        }

    }

    obG.m_iMaxPK = iCurSup + 2;

    printf("DETRUSS_P valid k max: %d\n", obG.m_iMaxPK);
    vKG.clear();
    vKG.resize(obG.m_iMaxPK + 1);
    for (int iEid : vERaw)
    {
        ASSERT(vKT[iEid] <= obG.m_iMaxPK);
        vKG[vKT[iEid]].push_back(iEid);
    }

    /* debug */
    /*for (int iCurK = 0; iCurK < vKG.size(); ++iCurK)
    {
        printf("DETRUSS_P %d-class: %d\n", iCurK, vKG[iCurK].size());
    }*/
    /*for (int iEid : vKG[5])
    {
        printf("DETRUSS_P 5-class edge: %d\n", iEid);
    }*/
    return 0;
}

/*****************
input:
        myG &obG
        int iDesK
        double ldMiniP
        double ldAccu
        vector<int> &vEdges
        vector<pair<int, int> > &vPos
output:
        LIST_DECOMP_G &lstDeG
description:
        vEdges $\cap vUpE = $\emptyset
******************/
int deTruss_topL(myG &obG, int iDesK, double ldMiniP, double ldAccu, vector<int> &vEdges, vector<pair<int, int> > &vPos)
{
    int iUpR = MAX_R(ldAccu);
    vPos.clear();
    vector<int> vERaw;
    vector<int> vEDone;
    /* init */
    for (int iEid = 1; iEid <= obG.m_iMaxEId; ++iEid)
    {
        obG.findNode(iEid)->bRm = true;
    }
    for (int iEid : vEdges)
    {
        if (obG.findNode(iEid)->p >= ldMiniP)
        {
            obG.findNode(iEid)->bRm = false;
            vERaw.push_back(iEid);
        }
        else
        {
            vEDone.push_back(iEid);
        }
        obG.findNode(iEid)->iColR = -1;
    }

    /* init vector */
    calPek_part(obG, vERaw);

    //printf("DE_TRUSS first raw: %d\n", vERaw.size());

    vector<int> vER(obG.m_iMaxEId + 1, 0);
    auto cmp = [](pair<int, int> &left, pair<int, int> &right) { return left.second > right.second; };
    priority_queue<int, vector<pair<int, int> >, decltype(cmp)> prQ(cmp);
    /* init R */
    for (int iEid : vERaw)
    {
        vER[iEid] = getR(obG.findNode(iEid)->vldPek, obG.findNode(iEid)->p, ldMiniP, ldAccu, iDesK);
        //printf("DE_TRUSS edge: %d R: %d size: %d\n", iEid, vER[iEid], obG.findNode(iEid)->vldPek.size());
        prQ.push({iEid, vER[iEid]});
    }

    vector<int> vPT(obG.m_iMaxEId + 1, -1);
    vector<bool> removed(obG.m_iMaxEId + 1, false);
    int iCurR = 0;
    int iRCnt = 0;
    int iRmCnt = 0;

    if (0 < vEDone.size())
    {
        iRCnt = vEDone.size();
        iRmCnt = iRCnt;
    }

    vector<int> vWait(obG.m_iMaxEId + 1);
    // 2.2.1. process the edges layer by layer
    while (iRmCnt < vERaw.size())
    {
        if (iCurR >= iUpR)
        {
            /* end */
            break;
        }
        vWait.clear();
        while (!prQ.empty())
        {
            int iEid = prQ.top().first;
            if (vER[iEid] > iCurR)
            {
                /* end */
                break;
            }
            if ((!removed[iEid]) && (iCurR != vPT[iEid]))
              {
                  vWait.push_back(iEid);
                  vPT[iEid] = iCurR;
                  if (iCurR != obG.findNode(iEid)->iColR)
                  {
                    obG.findNode(iEid)->iColR = iCurR;
                    obG.findNode(iEid)->vldResPek = obG.findNode(iEid)->vldPek;
                    //printf("FIRST debug 1 save %d size: %d\n", iEid, obG.findNode(iEid)->vldResPek.size());
                  }
                ++iRCnt;
              }
            prQ.pop();
        }

        if (vWait.empty())
        {
          if (0 < iRCnt)
          {
            vPos.push_back({iRCnt, iCurR});
            iRCnt = 0;
          }
          if (prQ.empty())
          {
              break;
          }
          iCurR = vER[prQ.top().first];
          //printf("DE_TRUSS first get new R: %d\n", iCurR);
          continue;
        }

        for (auto eid : vWait)
        {
            removed[eid] = true;
            ++iRmCnt;
            //printf("start find triangles\n");
            // find triangles containing the edge with ID eid
            g_vLfE.clear();
            g_vRtE.clear();
            obG.findNotRmNeb(obG.findNode(eid)->paXY.first, obG.findNode(eid)->paXY.second, g_vLfE, g_vRtE);

            //printf("start rm triangles\n");
            for (int i = 0; i < g_vLfE.size(); ++i)
            {
                const int e1 = g_vLfE[i];
                const int e2 = g_vRtE[i];
                if (removed[e1] || removed[e2]) continue;

                if (iCurR != obG.findNode(e1)->iColR)
                {
                    obG.findNode(e1)->iColR = iCurR;
                    obG.findNode(e1)->vldResPek = obG.findNode(e1)->vldPek;

                    //printf("FIRST debug 2 save %d size: %d\n", e1, obG.findNode(e1)->vldResPek.size());
                }
                if (vER[e1] > iCurR)
                {
                    calPek_rm(obG, e1, eid, e2, obG.findNode(e1)->vldPek);
                    int iTpR = getR(obG.findNode(e1)->vldPek, obG.findNode(e1)->p, ldMiniP, ldAccu, iDesK);
                    if (iTpR < vER[e1])
                    {
                        vER[e1] = iTpR;
                        prQ.push({e1, iTpR});
                    }
                }
                if (iCurR != obG.findNode(e2)->iColR)
                {
                    obG.findNode(e2)->iColR = iCurR;
                    obG.findNode(e2)->vldResPek = obG.findNode(e2)->vldPek;
                    //printf("FIRST debug 3 save %d size: %d\n", e2, obG.findNode(e2)->vldResPek.size());
                }
                if (vER[e2] > iCurR)
                {
                    calPek_rm(obG, e2, eid, e1, obG.findNode(e2)->vldPek);

                    int iTpR = getR(obG.findNode(e2)->vldPek, obG.findNode(e2)->p, ldMiniP, ldAccu, iDesK);
                    if (iTpR < vER[e2])
                    {
                        vER[e2] = iTpR;
                        prQ.push({e2, iTpR});
                    }
                }
            }
            obG.findNode(eid)->bRm = true;
            vEDone.push_back(eid);
        }

    }

    /* last */
    if (0 < iRCnt)
    {
        vPos.push_back({iRCnt, iCurR});
        iRCnt = 0;
    }
    if (iCurR >= iUpR)
    {
        iRCnt = 0;
        for (int iEid : vERaw)
        {
            if (!(removed[iEid]))
            {
                ++iRCnt;
                vEDone.push_back(iEid);
                obG.findNode(iEid)->vldResPek = obG.findNode(iEid)->vldPek;
                obG.findNode(iEid)->iColR = iUpR;
                //printf("FIRST debug 4 save %d size: %d\n", iEid, obG.findNode(iEid)->vldResPek.size());
            }
        }
        if (0 < iRCnt)
        {
            vPos.push_back({iRCnt, iUpR});
        }
    }
    vEdges.swap(vEDone);
    for (int iEid : vEdges)
    {
        obG.findNode(iEid)->bUp = true;
    }
    //printf("DE_TRUSS first done\n");
    //ASSERT(0);
    return 0;
}
/*****************
input:
        myG &mapInitG
        LIST_DECOMP_G &lstDeG
        long double ldMiniP
output:
        vEdges $\subseteq vEdges $\cup vAffE
description:
        vEdges $\cap vUpE = $\emptyset
        vAffE $\cap vUpE = $\emptyset
        vEdges $\cap vAffE = $\emptyset
        vEdges in k-truss
        vAffE in (k+1)-truss but smaller than iUpR
        vUpE in (k+1)-truss and no smaller than iUpR
******************/
int deTruss_firstL(myG &obG, int iDesK, double ldMiniP, double ldAccu, vector<int> &vEdges, vector<int> &vAffE, vector<pair<int, int> > &vPos, vector<int> &vUpE, int iUpR)
{
    //printf("DE_TRUSS get upper R: %d affected: %d\n", iUpR, vAffE.size());
    vector<int> vERaw;
    vector<int> vENew;
    vector<int> vValidAffE;
    vector<int> vEDone;

    vERaw.reserve(vEdges.size() + vAffE.size());

    for (int iEid : vEdges)
    {
        if (obG.findNode(iEid)->p >= ldMiniP)
        {
            obG.findNode(iEid)->bRm = false;
            vERaw.push_back(iEid);
        }
        else
        {
            vEDone.push_back(iEid);
        }
        obG.findNode(iEid)->iColR = -1;
        obG.findNode(iEid)->iPT = -1;
    }
    vENew = vERaw;

    for (int iEid : vAffE)
    {
        if (obG.findNode(iEid)->p >= ldMiniP)
        {
            obG.findNode(iEid)->bRm = false;
            vValidAffE.push_back(iEid);

            obG.findNode(iEid)->vldPek.clear();
            obG.findNode(iEid)->vldPek.swap(obG.findNode(iEid)->vldResPek);
        }
        else
        {
            vEDone.push_back(iEid);
        }
        /*printf("DE_TRUSS FIRST eid: %d iColR: %d iPT: %d pek: %d\n",
               iEid, obG.findNode(iEid)->iColR, obG.findNode(iEid)->iPT,
               obG.findNode(iEid)->vldPek.size());*/
        obG.findNode(iEid)->iColR = -1;
        obG.findNode(iEid)->iPT = -1;
    }
    vERaw.insert(vERaw.end(), vValidAffE.begin(), vValidAffE.end());
    //printf("DE_TRUSS raw: %d\n", vERaw.size());

    // int upper layer
    for (int iEid : vUpE)
    {
        obG.findNode(iEid)->bUp = true;
        obG.findNode(iEid)->bRm = false;
    }
    for (int iEid : vENew)
    {
        obG.findNode(iEid)->bNew = true;
    }
    /* init vector */
    calPek_part(obG, vENew);
    for (int iEid : vValidAffE)
    {
        obG.findNode(iEid)->bUp = false;
        g_vLfE.clear();
        g_vRtE.clear();
        obG.findNotRmNeb(obG.findNode(iEid)->paXY.first, obG.findNode(iEid)->paXY.second, g_vLfE, g_vRtE);
        for (int i = 0; i < g_vLfE.size(); ++i)
        {
            const int e1 = g_vLfE[i];
            const int e2 = g_vRtE[i];
            if (obG.findNode(e1)->bNew || obG.findNode(e2)->bNew)
            {
                /* new triangle */
                int iRes = calPek_add(obG, iEid, e1, e2, obG.findNode(iEid)->vldPek);
                if (0 != iRes)
                {
                    /* recalculate */
                    break;
                }
            }
        }
        if (g_vLfE.size() + 1 != obG.findNode(iEid)->vldPek.size())
        {
            /* ERROR */
            printf("ERROR edge: %d p: %.3f R: %d left size: %d pek: %d\n",
                   iEid, obG.findNode(iEid)->p, obG.findNode(iEid)->iPT,
                   g_vLfE.size(), obG.findNode(iEid)->vldPek.size());
            for (int i = 0; i < g_vLfE.size(); ++i)
            {
                const int e1 = g_vLfE[i];
                const int e2 = g_vRtE[i];
                if (obG.findNode(e1)->bNew || obG.findNode(e2)->bNew)
                {
                    printf("new: %lf %lf\n", obG.findNode(e1)->p, obG.findNode(e2)->p);
                }
                else
                {
                    printf("old: %lf %lf\n", obG.findNode(e1)->p, obG.findNode(e2)->p);
                }
            }
            for (auto atR : obG.findNode(iEid)->vldPek)
            {
                printf("%Lf\n", atR);
            }
            ASSERT(0);
        }
    }

    //printf("DE_TRUSS init R size: %d\n", vERaw.size());
    int iCurR = 0;
    int iRCnt = 0;
    int iRmCnt = 0;
    if (0 < vEDone.size())
    {
        iRCnt = vEDone.size();
        //iRmCnt = iRCnt;
    }
    //printf("DE_TRUSS init R done\n");
    auto cmp = [](pair<int, int> &left, pair<int, int> &right) { return left.second > right.second; };
    priority_queue<int, vector<pair<int, int> >, decltype(cmp)> prQ(cmp);
    /* init R */
    for (int iEid : vERaw)
    {
        obG.findNode(iEid)->iSelfR = getR(obG.findNode(iEid)->vldPek, obG.findNode(iEid)->p, ldMiniP, ldAccu, iDesK);
        prQ.push({iEid, obG.findNode(iEid)->iSelfR});
    }

    //vector<int> vPT(obG.m_iMaxEId + 1, 0);
    //vector<bool> removed(obG.m_iMaxEId + 1, false);
    vector<int> vWait(obG.m_iMaxEId + 1);
    // 2.2.1. process the edges layer by layer

    //int iDebugFirstEid = 0;
    //int iDebugR = 0;
    while (iRmCnt < vERaw.size())
    {
        //printf("FIRST curR: %d Q: %d\n", iCurR, prQ.size());
        if (iCurR >= iUpR)
        {
            /* end */
            break;
        }
        vWait.clear();
        while (!prQ.empty())
        {
            int iEid = prQ.top().first;
            if (obG.findNode(iEid)->iSelfR > iCurR)
            {
                /* end */
                //printf("FIRST get larger R: %d\n", obG.findNode(iEid)->iSelfR);
                break;
            }
            if ((!obG.findNode(iEid)->bRm) && (iCurR != obG.findNode(iEid)->iPT))
              {
                  /*printf("FIRST get edge: %d from Q selfR: %d curR: %d\n",
                         iEid, obG.findNode(iEid)->iSelfR, iCurR);*/

                  /* debug */
                  /*if (0 == iDebugFirstEid)
                  {
                      iDebugFirstEid = iEid;
                      iDebugR = obG.findNode(iEid)->iSelfR;
                  }*/

                  vWait.push_back(iEid);
                  obG.findNode(iEid)->iPT = iCurR;
                  if (iCurR != obG.findNode(iEid)->iColR)
                  {
                    obG.findNode(iEid)->iColR = iCurR;
                    obG.findNode(iEid)->vldResPek = obG.findNode(iEid)->vldPek;
                  }
                ++iRCnt;
              }
            prQ.pop();
        }

        //printf("FIRST raw R: %d size: %d\n", iCurR, vWait.size());
        if (vWait.empty())
        {
          if (0 < iRCnt)
          {
            vPos.push_back({iRCnt, iCurR});
            iRCnt = 0;
          }
          if (prQ.empty())
          {
              break;
          }
          /*printf("DE_TRUSS first get eid: %d old R: %d new R: %d\n",
                 prQ.top().first, iCurR, obG.findNode(prQ.top().first)->iSelfR);*/
          iCurR = obG.findNode(prQ.top().first)->iSelfR;
          ASSERT(!(obG.findNode(prQ.top().first)->bRm));
          continue;
        }

        //printf("FIRST R: %d size: %d\n", iCurR, vWait.size());

        for (auto eid : vWait)
        {
            ++iRmCnt;
            //printf("start find triangles\n");
            // find triangles containing the edge with ID eid
            g_vLfE.clear();
            g_vRtE.clear();
            obG.findNotRmNeb(obG.findNode(eid)->paXY.first, obG.findNode(eid)->paXY.second, g_vLfE, g_vRtE);

            //printf("start rm triangles left: %d right: %d\n", obG.m_vLfE.size(), obG.m_vRtE.size());
            ASSERT(g_vLfE.size() == g_vRtE.size());
            for (int i = 0; i < g_vLfE.size(); ++i)
            {
                const int e1 = g_vLfE[i];
                const int e2 = g_vRtE[i];
                if (obG.findNode(e1)->bRm || obG.findNode(e2)->bRm) continue;

                if (!(obG.findNode(e1)->bUp) && (obG.findNode(e1)->iSelfR > iCurR))
                {
                    if (iCurR != obG.findNode(e1)->iColR)
                    {
                        obG.findNode(e1)->iColR = iCurR;
                        obG.findNode(e1)->vldResPek = obG.findNode(e1)->vldPek;
                    }
                    calPek_rm(obG, e1, eid, e2, obG.findNode(e1)->vldPek);
                    /*printf("FIRST rm Pek of edge: %d (%d, %d) neighbor: %d (%d, %d) %d (%d, %d) 1\n",
                           e1, obG.findNode(e1)->paXY.first, obG.findNode(e1)->paXY.second,
                           eid, obG.findNode(eid)->paXY.first, obG.findNode(eid)->paXY.second,
                           e2, obG.findNode(e2)->paXY.first, obG.findNode(e2)->paXY.second);*/
                    /* debug */
                    /*if ((46695 == e1) && (50490 == eid))
                    {
                        for (int iTpId = 0; iTpId < g_vLfE.size(); ++iTpId)
                        {
                            int iTpLeft = g_vLfE[iTpId];
                            int iTpRight = g_vRtE[iTpId];
                            printf("%d (%d, %d) %d (%d, %d)\n",
                                   iTpLeft, obG.findNode(iTpLeft)->paXY.first, obG.findNode(iTpLeft)->paXY.second,
                                   iTpRight, obG.findNode(iTpRight)->paXY.first, obG.findNode(iTpRight)->paXY.second);
                        }

                        vector<int> vTpLfE;
                        vector<int> vTpRtE;
                        obG.findNotRmNeb(obG.findNode(eid)->paXY.first, obG.findNode(eid)->paXY.second, vTpLfE, vTpRtE);
                        for (int iTpId = 0; iTpId < vTpLfE.size(); ++iTpId)
                        {
                            int iTpLeft = vTpLfE[iTpId];
                            int iTpRight = vTpRtE[iTpId];
                            printf("DEBUG %d (%d, %d) %d (%d, %d)\n",
                                   iTpLeft, obG.findNode(iTpLeft)->paXY.first, obG.findNode(iTpLeft)->paXY.second,
                                   iTpRight, obG.findNode(iTpRight)->paXY.first, obG.findNode(iTpRight)->paXY.second);
                        }
                    }*/
                    int iTpR = getR(obG.findNode(e1)->vldPek, obG.findNode(e1)->p, ldMiniP, ldAccu, iDesK);
                    if (iTpR < obG.findNode(e1)->iSelfR)
                    {
                        obG.findNode(e1)->iSelfR = iTpR;
                        prQ.push({e1, iTpR});
                        //printf("FIRST push edge: %d into Q selfR: %d curR: %d\n", e1, iTpR, iCurR);
                    }
                }
                if (!(obG.findNode(e2)->bUp) && (obG.findNode(e2)->iSelfR > iCurR))
                {
                    if (iCurR != obG.findNode(e2)->iColR)
                    {
                        obG.findNode(e2)->iColR = iCurR;
                        obG.findNode(e2)->vldResPek = obG.findNode(e2)->vldPek;
                    }

                    calPek_rm(obG, e2, eid, e1, obG.findNode(e2)->vldPek);
                    /*printf("FIRST rm Pek of edge: %d (%d, %d) neighbor: %d (%d, %d) %d (%d, %d) 2\n",
                           e2, obG.findNode(e2)->paXY.first, obG.findNode(e2)->paXY.second,
                           eid, obG.findNode(eid)->paXY.first, obG.findNode(eid)->paXY.second,
                           e1, obG.findNode(e1)->paXY.first, obG.findNode(e1)->paXY.second);*/

                    int iTpR = getR(obG.findNode(e2)->vldPek, obG.findNode(e2)->p, ldMiniP, ldAccu, iDesK);
                    if (iTpR < obG.findNode(e2)->iSelfR)
                    {
                        obG.findNode(e2)->iSelfR = iTpR;
                        prQ.push({e2, iTpR});
                        /*printf("FIRST push edge: %d into Q selfR: %d curR: %d\n",
                               e2, iTpR, iCurR);*/
                    }
                }
            }

            //printf("rm triangles done\n");
            obG.findNode(eid)->bRm = true;
            vEDone.push_back(eid);
        }

    }

    /* last */
    if (0 < iRCnt)
    {
        vPos.push_back({iRCnt, iCurR});

        /* debug */
        /*vector<int> vSubE;
        vSubE.insert(vSubE.end(), vEDone.end() - iRCnt, vEDone.end());
        pair<int, int> paMinXY;
        int iDiff = deTruss_verify(obG, vSubE, iCurR, iDesK, ldMiniP, ldAccu, paMinXY);
        if(abs(iDiff) > 1)
        {
          printf("ERROR low accuracy k: %d |H|: %d estimate R: %d difference: %d\n",
                 iDesK, vSubE.size(), iCurR, iDiff);

          TPST_E *pstNode = NULL;
          for (int iEid : vERaw)
          {
              pstNode = obG.findNode(iEid);
              pstNode->bRm = false;
          }
          vector<int> vDebugE;
          pstNode = obG.findNode(paMinXY.first, paMinXY.second);
          vDebugE.push_back(pstNode->eid);
          calPek_part(obG, vDebugE);
          int iAfterR = getR(pstNode->vldPek, pstNode->p, ldMiniP, ldAccu, iDesK);
          printf("ERROR minimum edge: %d reR: %d new: %d rm: %d up: %d\n",
                 pstNode->eid, iAfterR, pstNode->bNew, pstNode->bRm, pstNode->bUp);

          vDebugE.clear();
          vDebugE.push_back(iDebugFirstEid);
          calPek_part(obG, vDebugE);
          pstNode = obG.findNode(iDebugFirstEid);
          iAfterR = getR(pstNode->vldPek, pstNode->p, ldMiniP, ldAccu, iDesK);
          printf("ERROR first edge: %d R: %d recalculate R: %d new: %d rm: %d up: %d\n",
                 iDebugFirstEid, iDebugR, iAfterR, pstNode->bNew, pstNode->bRm, pstNode->bUp);
          ASSERT(0);
        }*/

        iRCnt = 0;
    }

    //printf("DE_TRUSS first last R: %d done: %d\n", iCurR, vEDone.size());
    if (iCurR >= iUpR)
    {
        iRCnt = 0;
        for (int iEid : vERaw)
        {
            if (!(obG.findNode(iEid)->bRm))
            {
                ++iRCnt;
                vEDone.push_back(iEid);
                obG.findNode(iEid)->vldResPek = obG.findNode(iEid)->vldPek;
                obG.findNode(iEid)->iColR = iUpR;
            }
        }
        if (0 < iRCnt)
        {
            vPos.push_back({iRCnt, iUpR});
        }
    }
    else
    {
        /* no remained edges */
        vPos.push_back({0, iUpR});
    }
    /* restore */
    for (int iEid : vENew)
    {
        obG.findNode(iEid)->bNew = false;
    }
    vEdges.swap(vEDone);
    //printf("DE_TRUSS first done\n");
    return 0;
}

/*****************
input:
        myG &mapInitG
        LIST_DECOMP_G &lstDeG
        long double ldMiniP
output:
        vEdges $\subseteq vEdges $\cup vAffE
description:
        vEdges $\cap vUpE = $\emptyset
        vAffE $\cap vUpE = $\emptyset
        vEdges $\cap vAffE = $\emptyset
        vEdges in k-truss
        vAffE in (k+1)-truss but smaller than iUpR
        vUpE in (k+1)-truss and no smaller than iUpR
******************/
int deTruss_middleL(myG &obG, int iDesK, double ldMiniP, double ldAccu, vector<int> &vEdges, vector<int> &vAffE, vector<pair<int, int> > &vPos, vector<int> &vUpE, int iUpR)
{
    //printf("DE_TRUSS get upper R: %d affected: %d\n", iUpR, vAffE.size());
    vector<int> vERaw = vEdges;
    vector<int> vEDone;
    vERaw.insert(vERaw.end(), vAffE.begin(), vAffE.end());
    //printf("DE_TRUSS raw: %d\n", vERaw.size());
    for (int iEid : vERaw)
    {
        /*printf("DE_TRUSS SECOND eid: %d selfR: %d iColR: %d iPT: %d pek: %d\n",
               iEid, obG.findNode(iEid)->iSelfR, obG.findNode(iEid)->iColR, obG.findNode(iEid)->iPT,
               obG.findNode(iEid)->vldResPek.size());*/
        obG.findNode(iEid)->iColR = -1;
        obG.findNode(iEid)->iPT = -1;
        ASSERT(!obG.findNode(iEid)->bRm);
    }
    for (int iEid : vAffE)
    {
        obG.findNode(iEid)->vldPek.clear();
        obG.findNode(iEid)->vldPek.swap(obG.findNode(iEid)->vldResPek);
        obG.findNode(iEid)->bUp = false;
        ASSERT(!obG.findNode(iEid)->bRm);
        /*printf("DE_TRUSS SECOND successor eid: %d pek: %d\n",
               iEid, obG.findNode(iEid)->vldPek.size());*/
    }
    for (int iEid : vEdges)
    {
        obG.findNode(iEid)->bNew = true;
    }
    /* init vector */
    //calPek_part(obG, vEdges);
    for (int iEid : vAffE)
    {
        g_vLfE.clear();
        g_vRtE.clear();
        obG.findNotRmNeb(obG.findNode(iEid)->paXY.first, obG.findNode(iEid)->paXY.second, g_vLfE, g_vRtE);
        int iOldPek = obG.findNode(iEid)->vldPek.size();
        for (int i = 0; i < g_vLfE.size(); ++i)
        {
            const int e1 = g_vLfE[i];
            const int e2 = g_vRtE[i];
            if (obG.findNode(e1)->bNew || obG.findNode(e2)->bNew)
            {
                /* new triangle */
                int iRes = calPek_add(obG, iEid, e1, e2, obG.findNode(iEid)->vldPek);
                if (0 != iRes)
                {
                    /* recalculate */
                    break;
                }
            }
        }
        if (g_vLfE.size() + 1 != obG.findNode(iEid)->vldPek.size())
        {
            /* ERROR */
            printf("ERROR edge: %d p: %.3f R: %d left size: %d pek: %d old pek: %d upR: %d\n",
                   iEid, obG.findNode(iEid)->p, obG.findNode(iEid)->iPT,
                   g_vLfE.size(), obG.findNode(iEid)->vldPek.size(), iOldPek, iUpR);
            for (int i = 0; i < g_vLfE.size(); ++i)
            {
                const int e1 = g_vLfE[i];
                const int e2 = g_vRtE[i];
                if (obG.findNode(e1)->bNew || obG.findNode(e2)->bNew)
                {
                    printf("new: %lf %lf\n", obG.findNode(e1)->p, obG.findNode(e2)->p);
                }
                else
                {
                    printf("old: %lf %lf\n", obG.findNode(e1)->p, obG.findNode(e2)->p);
                }
            }
            for (auto atR : obG.findNode(iEid)->vldPek)
            {
                printf("%Lf\n", atR);
            }
            g_vLfE.clear();
            g_vRtE.clear();
            obG.findNeb(obG.findNode(iEid)->paXY.first, obG.findNode(iEid)->paXY.second, g_vLfE, g_vRtE);
            for (int i = 0; i < g_vLfE.size(); ++i)
            {
                const int e1 = g_vLfE[i];
                const int e2 = g_vRtE[i];
                printf("edge: %d %d new: %d %d rm: %d %d up: %d %d PT: %d %d oldPek: %d %d\n",
                       e1, e2,
                       obG.findNode(e1)->bNew, obG.findNode(e2)->bNew,
                       obG.findNode(e1)->bRm, obG.findNode(e2)->bRm,
                       obG.findNode(e1)->bUp, obG.findNode(e2)->bUp,
                       obG.findNode(e1)->iPT, obG.findNode(e2)->iPT,
                       obG.findNode(e1)->vldResPek.size(), obG.findNode(e2)->vldResPek.size());
            }
            ASSERT(0);
        }
    }

    //printf("DE_TRUSS init R size: %d\n", vERaw.size());
    //printf("DE_TRUSS init R done\n");
    auto cmp = [](pair<int, int> &left, pair<int, int> &right) { return left.second > right.second; };
    priority_queue<int, vector<pair<int, int> >, decltype(cmp)> prQ(cmp);
    /* init R */
    for (int iEid : vERaw)
    {
        obG.findNode(iEid)->iSelfR = getR(obG.findNode(iEid)->vldPek, obG.findNode(iEid)->p, ldMiniP, ldAccu, iDesK);
        prQ.push({iEid, obG.findNode(iEid)->iSelfR});
    }

    //vector<int> vPT(obG.m_iMaxEId + 1, 0);
    //vector<bool> removed(obG.m_iMaxEId + 1, false);
    int iCurR = 0;
    int iRCnt = 0;
    int iRmCnt = 0;
    vector<int> vWait(obG.m_iMaxEId + 1);
    // 2.2.1. process the edges layer by layer
    while (iRmCnt < vERaw.size())
    {
        if (iCurR >= iUpR)
        {
            /* end */
            break;
        }
        vWait.clear();
        while (!prQ.empty())
        {
            int iEid = prQ.top().first;
            if (obG.findNode(iEid)->iSelfR > iCurR)
            {
                /* end */
                break;
            }
            if ((!obG.findNode(iEid)->bRm) && (iCurR != obG.findNode(iEid)->iPT))
              {
                  vWait.push_back(iEid);
                  obG.findNode(iEid)->iPT = iCurR;
                  if (iCurR != obG.findNode(iEid)->iColR)
                  {
                    obG.findNode(iEid)->iColR = iCurR;
                    obG.findNode(iEid)->vldResPek = obG.findNode(iEid)->vldPek;
                  }
                ++iRCnt;
              }
            prQ.pop();
        }

        if (vWait.empty())
        {
            if (0 < iRCnt)
            {
                vPos.push_back({iRCnt, iCurR});
                iRCnt = 0;
            }
            if (prQ.empty())
            {
                break;
            }
            iCurR = obG.findNode(prQ.top().first)->iSelfR;
            //printf("DE_TRUSS second get new R: %d\n", iCurR);

            /* debug */
            /*if (iCurR >= iUpR)
            {
                continue;
            }
            vector<int> vSubE = vUpE;
            for (int iEid : vERaw)
            {
                if (!obG.findNode(iEid)->bRm)
                {
                    vSubE.push_back(iEid);
                }
            }
            pair<int, int> paMinXY;
            int iDiff = deTruss_verify(obG, vSubE, iCurR, iDesK, ldMiniP, ldAccu, paMinXY);
            if(abs(iDiff) > 1)
            {
                printf("ERROR low accuracy k: %d |H|: %d upR: %d estimate R: %d difference: %d\n",
                    iDesK, vSubE.size(), iUpR, iCurR, iDiff);
                TPST_E *pstNode = obG.findNode(paMinXY.first, paMinXY.second);
                vector<int> vLfE;
                vector<int> vRtE;
                vector <long double> vfDtE;
                obG.findNotRmNeb(pstNode->paXY.first, pstNode->paXY.second, vLfE, vRtE);
                cal_OneEdgePek(obG, vLfE, vRtE, vfDtE);
                int iDebugR = getR(vfDtE, pstNode->p, ldMiniP, ldAccu, iDesK);
                printf("ERROR minimum edge: %d selfR: %d reR: %d new: %d rm: %d up: %d neighbors: %d\n",
                    pstNode->eid, pstNode->iSelfR, iDebugR, pstNode->bNew, pstNode->bRm, pstNode->bUp,
                    vLfE.size());

                for (int i = 0; i < vLfE.size(); ++i)
                {
                    printf("SECOND %d(%d, %d) %d(%d, %d) new: %d %d rm: %d %d up: %d %d R: %d %d p: %f %f t: %d %d\n",
                           vLfE[i], obG.findNode(vLfE[i])->paXY.first, obG.findNode(vLfE[i])->paXY.second,
                           vRtE[i], obG.findNode(vRtE[i])->paXY.first, obG.findNode(vRtE[i])->paXY.second,
                           obG.findNode(vLfE[i])->bNew, obG.findNode(vRtE[i])->bNew,
                           obG.findNode(vLfE[i])->bRm, obG.findNode(vRtE[i])->bRm,
                           obG.findNode(vLfE[i])->bUp, obG.findNode(vRtE[i])->bUp,
                           obG.findNode(vLfE[i])->iSelfR, obG.findNode(vRtE[i])->iSelfR,
                           obG.findNode(vLfE[i])->p, obG.findNode(vRtE[i])->p,
                           obG.findNode(vLfE[i])->iTrussness, obG.findNode(vRtE[i])->iTrussness);
                }
                for (int iEid : vAffE)
                {
                    printf("ERROR affected eid: %d\n",
                        iEid);
                }
                for (int iEid : vEDone)
                {
                    printf("ERROR done eid: %d\n",
                        iEid);
                }
                ASSERT(0);
            }*/
            continue;
        }

        for (auto eid : vWait)
        {
            ++iRmCnt;
            //printf("start find triangles\n");
            // find triangles containing the edge with ID eid
            g_vLfE.clear();
            g_vRtE.clear();
            obG.findNotRmNeb(obG.findNode(eid)->paXY.first, obG.findNode(eid)->paXY.second, g_vLfE, g_vRtE);

            //printf("start rm triangles\n");
            for (int i = 0; i < g_vLfE.size(); ++i)
            {
                const int e1 = g_vLfE[i];
                const int e2 = g_vRtE[i];
                if (obG.findNode(e1)->bRm || obG.findNode(e2)->bRm) continue;

                if (!(obG.findNode(e1)->bUp) && (obG.findNode(e1)->iSelfR > iCurR))
                {
                    if (iCurR != obG.findNode(e1)->iColR)
                    {
                        obG.findNode(e1)->iColR = iCurR;
                        obG.findNode(e1)->vldResPek = obG.findNode(e1)->vldPek;
                    }
                    calPek_rm(obG, e1, eid, e2, obG.findNode(e1)->vldPek);
                    int iTpR = getR(obG.findNode(e1)->vldPek, obG.findNode(e1)->p, ldMiniP, ldAccu, iDesK);
                    if (iTpR < obG.findNode(e1)->iSelfR)
                    {
                        obG.findNode(e1)->iSelfR = iTpR;
                        prQ.push({e1, iTpR});
                    }
                }
                if (!(obG.findNode(e2)->bUp) && (obG.findNode(e2)->iSelfR > iCurR))
                {
                    if (iCurR != obG.findNode(e2)->iColR)
                    {
                        obG.findNode(e2)->iColR = iCurR;
                        obG.findNode(e2)->vldResPek = obG.findNode(e2)->vldPek;
                    }

                    calPek_rm(obG, e2, eid, e1, obG.findNode(e2)->vldPek);

                    int iTpR = getR(obG.findNode(e2)->vldPek, obG.findNode(e2)->p, ldMiniP, ldAccu, iDesK);
                    if (iTpR < obG.findNode(e2)->iSelfR)
                    {
                        obG.findNode(e2)->iSelfR = iTpR;
                        prQ.push({e2, iTpR});
                    }
                }
            }
            obG.findNode(eid)->bRm = true;
            vEDone.push_back(eid);
        }

    }

    /* last */
    if (0 < iRCnt)
    {
        vPos.push_back({iRCnt, iCurR});

        /* debug */
        /*vector<int> vSubE;
        vSubE.insert(vSubE.end(), vEDone.end() - iRCnt, vEDone.end());
        pair<int, int> paMinXY;
        int iDiff = deTruss_verify(obG, vSubE, iCurR, iDesK, ldMiniP, ldAccu, paMinXY);
        if(abs(iDiff) > 1)
        {
          printf("ERROR low accuracy k: %d estimate R: %d difference: %d\n",
                 iDesK, iCurR, iDiff);
          TPST_E *pstNode = obG.findNode(paMinXY.first, paMinXY.second);
          printf("ERROR minimum edge: %d new: %d rm: %d up: %d\n",
                 pstNode->eid, pstNode->bNew, pstNode->bRm, pstNode->bUp);
          ASSERT(0);
        }*/

        iRCnt = 0;
    }

    //printf("DE_TRUSS second last R: %d done: %d\n", iCurR, vEDone.size());
    if (iCurR >= iUpR)
    {
        iRCnt = 0;
        for (int iEid : vERaw)
        {
            if (!(obG.findNode(iEid)->bRm))
            {
                ++iRCnt;
                vEDone.push_back(iEid);
                obG.findNode(iEid)->vldResPek = obG.findNode(iEid)->vldPek;
                obG.findNode(iEid)->iColR = iUpR;
            }
        }
        if (0 < iRCnt)
        {
            vPos.push_back({iRCnt, iUpR});
        }
    }
    else
    {
        /* no remained edges */
        vPos.push_back({0, iUpR});
    }
    /* restore */
    for (int iEid : vEdges)
    {
        obG.findNode(iEid)->bNew = false;
    }
    vEdges.swap(vEDone);
    //printf("DE_TRUSS second done raw size: %d\n", vERaw.size());
    return 0;
}

/*****************
input:
        myG &mapInitG
        LIST_DECOMP_G &lstDeG
        long double ldMiniP
output:
        LIST_DECOMP_G &lstDeG
description:
        truss decomposition
******************/
int deTruss_PT(myG &obG, vector<vector<int> > &vKG, vector<vector<pair<int, double> > > &vRetPos, double ldMiniP, double ldAccu)
{
    int iKTrussSize = 0;
    vector<int> vUpG;
    vector<vector<pair<int, int> > > vGPos;
    int iCurK = vKG.size() - 1;
    vGPos.resize(vKG.size());
    vRetPos.resize(vKG.size());
    printf("=============PEEL top k: %d k-class: %d\n", iCurK, vKG[iCurK].size());
    iKTrussSize = vKG[iCurK].size();
    if (2 >= iCurK)
    {
        PtoR(obG, vKG[iCurK], vRetPos[iCurK]);
        return 0;
    }
    /* top layer */
    ASSERT(!vKG[iCurK].empty());
    deTruss_topL(obG, iCurK, ldMiniP, ldAccu, vKG[iCurK], vGPos[iCurK]);
    --iCurK;
    /* next layers */
    for (; iCurK > 0; --iCurK)
    {
        if (2 >= iCurK)
        {
            PtoR(obG, vKG[iCurK], vRetPos[iCurK]);
            break;
        }
        /* debug */
        //verify(obG, vKG[iCurK+1], vGPos[iCurK+1], iCurK+1, ldMiniP, ldAccu);
        printf("=============PEEL current k: %d k-class: %d\n", iCurK, vKG[iCurK].size());
        iKTrussSize += vKG[iCurK].size();
        /* debug */
        /*for (auto atPos : vGPos[iCurK+1])
        {
            printf("PEEL upper R: (%d, %d)\n", atPos.first, atPos.second);
        }*/
        /* first layer */
        vector<pair<int, int> >::iterator itCircle = vGPos[iCurK+1].begin();

        //printf("PEEL upper edges: %d upper groups: %d\n", vKG[iCurK+1].size(), vGPos[iCurK+1].size());

        int iUpFirstCnt = itCircle->first;
        int iUpR = itCircle->second;
        /* offset in (k+1)-truss */
        int iUpCnt = 0;
        vector <int> vUpE;
        vector <int> vAffE;

        //printf("PEEL get top R: %d\n", iUpR);
        if (0 == iUpR)
        {
            vAffE.insert(vAffE.end(), vKG[iCurK+1].begin(), vKG[iCurK+1].begin() + iUpFirstCnt);
            //printf("PEEL vAffEs size: %d\n", vAffE.size());
            //printf("PEEL vKG: %d cnt: %d\n", vKG[iCurK+1].size(), iUpFirstCnt);
            vUpE.insert(vUpE.end(), vKG[iCurK+1].begin() + iUpFirstCnt, vKG[iCurK+1].end());
            //printf("PEEL vUpE size: %d\n", vUpE.size());
            iUpCnt = iUpFirstCnt;
            ++itCircle;
            if (itCircle != vGPos[iCurK+1].end())
            {
                iUpR = itCircle->second;
            }
            else
            {
                iUpR = MAX_R(ldAccu);
            }
        }
        else
        {
            vUpE = vKG[iCurK+1];
        }

        vector<int> vEdges;
        vEdges.swap(vKG[iCurK]);
        vKG[iCurK].clear();
        vGPos[iCurK].clear();
        //printf("PEEL first top R: %d edges: %d\n", iUpR, vEdges.size() + vAffE.size() + vUpE.size());
        //deTruss_firstL(obG, iCurK, ldMiniP, ldAccu, vEdges, vGPos[iCurK], vUpE, iUpR);
        deTruss_firstL(obG, iCurK, ldMiniP, ldAccu, vEdges, vAffE, vGPos[iCurK], vUpE, iUpR);
        ASSERT(iKTrussSize == vEdges.size() + vUpE.size());

        //printf("PEEL first done edges: %d up edges: %d upR: %d\n", vEdges.size(), vUpE.size(), iUpR);

        for (; itCircle != vGPos[iCurK+1].end();)
        {
            int iCurCnt = vGPos[iCurK].back().first;
            int iCurR = vGPos[iCurK].back().second;
            vGPos[iCurK].pop_back();
            int iAffCnt = itCircle->first;
            //printf("PEEL current: %d done: %d, last R: %d\n", vEdges.size(), vEdges.size() - iCurCnt, iCurR);
            ++itCircle;
            if (itCircle != vGPos[iCurK+1].end())
            {
                iUpR = itCircle->second;
            }
            else
            {
                iUpR = MAX_R(ldAccu);
            }
            //printf("PEEL top R: %d\n", iUpR);
            vKG[iCurK].insert(vKG[iCurK].end(), vEdges.begin(), vEdges.end() - iCurCnt);
            vEdges.erase(vEdges.begin(), vEdges.end() - iCurCnt);
            vAffE.clear();
            vAffE.insert(vAffE.end(), vKG[iCurK+1].begin() + iUpCnt, vKG[iCurK+1].begin() + iUpCnt + iAffCnt);
            //printf("PEEL affE: %d\n", vAffE.size());
            iUpCnt += iAffCnt;
            vUpE.clear();
            if (iUpCnt < vKG[iCurK+1].size())
            {
                vUpE.insert(vUpE.end(), vKG[iCurK+1].begin() + iUpCnt, vKG[iCurK+1].end());
            }
            //printf("PEEL used: %d upE: %d\n", iUpCnt, vUpE.size());
            deTruss_middleL(obG, iCurK, ldMiniP, ldAccu, vEdges, vAffE, vGPos[iCurK], vUpE, iUpR);
            //printf("PEEL second done edges: %d up edges: %d upR: %d\n", vEdges.size(), vUpE.size(), iUpR);
        }

        /* last */
        vKG[iCurK].insert(vKG[iCurK].end(), vEdges.begin(), vEdges.end());
        if (iKTrussSize != vKG[iCurK].size())
        {
            printf("ERROR k-truss: %d get %d\n", iKTrussSize, vKG[iCurK].size());
            ASSERT(0);
        }

        int iLastCnt = vGPos[iCurK].back().first;
        if (0 == iLastCnt)
        {
            vGPos[iCurK].pop_back();
        }
        //printf("PEEL k-truss: %d\n", vKG[iCurK].size());
    }

    //printf("PEEL save R size: %d\n", vGPos.size());
    /* save */
    vRetPos.resize(vGPos.size());
    for (iCurK = 3; iCurK < vGPos.size(); ++iCurK)
    {
        for (auto atPos : vGPos[iCurK])
        {
            int iCnt = atPos.first;
            double ldP = recoverR(atPos.second, ldAccu);
            vRetPos[iCurK].push_back({atPos.first, ldP});
            /*printf("PEEL k: %d (%d, %d) (%d, %.3f)\n", iCurK,
                   atPos.first, atPos.second,
                   atPos.first, (float)ldP);*/
        }
    }
    //printf("PEEL save done R size: %d\n", vRetPos.size());
    return 0;
}

/*****************
input:
        myG &mapInitG
        LIST_DECOMP_G &lstDeG
output:
        LIST_DECOMP_G &lstDeG
description:
        truss decomposition
******************/
int deTruss_Detm(myG &obG, vector<vector<int> > &vKG)
{
    uint32_t n_ = obG.m_iMaxPId;
    uint32_t m_ = obG.m_iMaxEId;
  // truss decomposition
  // 1. compute the support of each edge by triangle listing
  // 1.1. define a total order over the vertices
  const auto pred = [&obG](const uint32_t v1, const uint32_t v2) {
    const size_t deg1 = obG.m_Adj[v1].size();
    const size_t deg2 = obG.m_Adj[v2].size();
    if (deg1 != deg2) return deg1 > deg2;
    else return v1 > v2;
  };
  // 1.2. sort the vertices in non-ascending order of degree
  //printf("DE_TRUSS start sort nodes\n");
  std::vector<uint32_t> verts(n_ + 1);
  std::iota(verts.begin(), verts.end(), 0);
  std::sort(verts.begin(), verts.end(), pred);
  // 1.3. call the "forward" algorithm to list triangles
  //printf("DE_TRUSS start count support\n");
  std::vector<uint32_t> sup(m_ + 1, 0);
  std::vector<std::vector<TPST_ADJ>> A(n_ + 1);
  for (const uint32_t v : verts) {
    for (const auto ae : obG.m_Adj[v]) {
      const uint32_t u = ae.pid;
      const uint32_t e = ae.eid;
      if (!pred(v, u)) continue;
      size_t pv = 0, pu = 0;
      while (pv < A[v].size() && pu < A[u].size()) {
        if (A[v][pv].pid == A[u][pu].pid) {
          ++sup[A[v][pv].eid]; ++sup[A[u][pu].eid];
          ++sup[e];
          ++pv; ++pu;
        } else if (pred(A[v][pv].pid, A[u][pu].pid)) {
          ++pv;
        } else {
          ++pu;
        }
      }
      A[u].push_back({v, e});
    }
  }
  //decltype(A)().swap(A);
  decltype(verts)().swap(verts);

  // 2. decomposition
  // 2.1. sort the edges according to their supports

  //printf("DE_TRUSS start fill bins\n");
  const uint32_t max_sup = *std::max_element(sup.cbegin(), sup.cend());
  std::vector<std::vector<uint32_t>> bin(max_sup + 1);
  for (uint32_t eid = 1; eid <= m_; ++eid)
  {
      if (obG.findNode(eid)->bRm)
      {
          /* removed */
          continue;
      }
    bin[sup[eid]].push_back(eid);
    /* debug */
    g_vLfE.clear();
    g_vRtE.clear();
    obG.findNeb(obG.findNode(eid)->paXY.first, obG.findNode(eid)->paXY.second, g_vLfE, g_vRtE);
    if (g_vLfE.size() != sup[eid])
    {
        printf("ERROR eid: %d bool: %d sup: %d neighbor: %d\n",
               eid, obG.findNode(eid)->bRm, sup[eid], g_vLfE.size());

        printf("node %d: ", obG.findNode(eid)->paXY.first);
        for (auto atNode : obG.m_Adj[obG.findNode(eid)->paXY.first])
        {
            printf("%d(%d) ", atNode.pid, atNode.eid);
        }
        printf("\n A node %d: ", obG.findNode(eid)->paXY.first);
        for (auto atNode : A[obG.findNode(eid)->paXY.first])
        {
            printf("%d ", atNode.pid);
        }
        printf("\n node %d: ", obG.findNode(eid)->paXY.second);
        for (auto atNode : obG.m_Adj[obG.findNode(eid)->paXY.second])
        {
            printf("%d(%d) ", atNode.pid, atNode.eid);
        }
        printf("\n A node %d: ", obG.findNode(eid)->paXY.second);
        for (auto atNode : A[obG.findNode(eid)->paXY.second])
        {
            printf("%d ", atNode.pid);
        }
        printf("\n");
        ASSERT(0);
    }
  }

  //printf("start peeling max sup: %d\n", max_sup);
  std::vector<uint32_t> vET(m_ + 1, 0);
  // 2.2. peeling
  std::vector<bool> removed(m_ + 1, false);
  uint32_t uiK = 2;
  uint32_t uiMaxK = 2;
  uint32_t uiRmCnt = 0;
  std::vector<uint32_t> vWait(m_ + 1);
  // 2.2.1. process the edges layer by layer
    while (uiK <= max_sup + 2)
    {
        //printf("DE_TRUSS rm cnt: %d total: %d\n", uiRmCnt, m_);
        vWait.clear();
        for (uint32_t i = 0; i <= uiK - 2; ++i)
        {
          for (auto eid : bin[i])
          {
                /* avoid repetition */
              if ((!removed[eid]) && (uiK != vET[eid]))
              {
                  vWait.push_back(eid);
                  vET[eid] = uiK;
              }
          }
          bin[i].clear();
        }

        if (vWait.empty())
        {
          ++uiK;
          continue;
        }

        //printf("DE_TRUSS peel %d edges in %d-truss\n", vWait.size(), uiK);
        uiMaxK = max(uiMaxK, uiK);

        for (auto eid : vWait)
        {
            ASSERT(eid < removed.size());
            removed[eid] = true;
            ++uiRmCnt;
            //printf("start find triangles\n");
            // find triangles containing the edge with ID eid
            std::vector<std::pair<uint32_t, uint32_t>> tris;
            {
              const uint32_t v1 = obG.findNode(eid)->paXY.first;
              const uint32_t v2 = obG.findNode(eid)->paXY.second;
              size_t p1 = 0, p2 = 0;
                ASSERT(v1 < obG.m_Adj.size());
                ASSERT(v2 < obG.m_Adj.size());
                ASSERT(p1 < obG.m_Adj[v1].size());
                ASSERT(p2 < obG.m_Adj[v2].size());
              while (p1 < obG.m_Adj[v1].size() && p2 < obG.m_Adj[v2].size()) {
                if (obG.m_Adj[v1][p1].pid == obG.m_Adj[v2][p2].pid) {
                  tris.push_back({obG.m_Adj[v1][p1].eid, obG.m_Adj[v2][p2].eid});
                  ++p1; ++p2;
                } else if (obG.m_Adj[v1][p1].pid < obG.m_Adj[v2][p2].pid) {
                  ++p1;
                } else {
                  ++p2;
                }
              }
            }

            //printf("start rm triangles\n");
            for (const auto tri : tris)
            {
                const uint32_t e1 = tri.first;
                const uint32_t e2 = tri.second;
                ASSERT(e1 < removed.size());
                ASSERT(e2 < removed.size());
                if (removed[e1] || removed[e2]) continue;

                /*printf("DEBUG self: (%d, %d) get (%d, %d) (%d, %d) Pos: %d %d %d\n",
                       obG.findNode(eid)->paXY.first,
                       obG.findNode(eid)->paXY.second,
                       obG.findNode(e1)->paXY.first,
                       obG.findNode(e1)->paXY.second,
                       obG.findNode(e2)->paXY.first,
                       obG.findNode(e2)->paXY.second,
                       obG.m_Adj[obG.findNode(eid)->paXY.first].size(),
                       obG.m_Adj[obG.findNode(eid)->paXY.second].size(),
                       obG.m_Adj[obG.findNode(e1)->paXY.second].size());*/

                for (const uint32_t e : {e1, e2})
                {
                    ASSERT(e < sup.size());
                    --sup[e];
                    ASSERT_MSG(sup[e] < bin.size(), "get support: " << sup[e]);
                    bin[sup[e]].push_back(e);
                }
            }
        }
        //printf("DE_TRUSS peel end\n");

    }
    obG.m_iMaxK = uiMaxK;
    vKG.resize(obG.m_iMaxK + 1);
    //printf("Detruss max k: %d\n", obG.m_iMaxK);
    for (uint32_t eid = 1; eid <= m_; ++eid)
    {
        if (obG.findNode(eid)->bRm)
        {
          /* removed */
          continue;
        }
        vKG[vET[eid]].push_back(eid);
        obG.findNode(eid)->iTrussness = vET[eid];
    }
    return 0;
}

#endif
/*****************
input:
        myG &mapInitG
        LIST_DECOMP_G &lstDeG
        long double ldMiniP
output:
        LIST_DECOMP_G &lstDeG
description:
        truss decomposition
******************/
int verify(myG &obG, vector<int> &vKG, vector<pair<int, int> > &vPos, int iCurK, double ldMiniP, double ldAccu)
{
    vector<int>::iterator itBegin = vKG.begin();
    vector<int> vSubE;
    pair<int, int> paXY;
    for (auto atPos : vPos)
    {
        int iOffset = atPos.first;
        int iEsR = atPos.second;
        vSubE.clear();
        vSubE.insert(vSubE.end(), itBegin, vKG.end());
        int iDiff = deTruss_verify(obG, vSubE, iEsR, iCurK, ldMiniP, ldAccu, paXY);
        itBegin += iOffset;

        printf("VERIFY k: %d calculated R: %d difference: %d\n", iCurK, iEsR, iDiff);
        if (abs(iDiff) > 1)
        {
            TPST_E *pstNode = obG.findNode(paXY.first, paXY.second);
            printf("ERROR %d-truss edge: %d (%d, %d) self R: %d\nNeighbors:\n",
                   iCurK, pstNode->eid, paXY.first, paXY.second, pstNode->vSky[iCurK]);
            ASSERT(!g_bLock);
            g_bLock = true;
            g_vLfE.clear();
            g_vRtE.clear();
            obG.findNeb(paXY.first, paXY.second, g_vLfE, g_vRtE);
            int iNeighbors = 0;
            for (int iPos = 0; iPos < g_vLfE.size(); ++iPos)
            {
                TPST_E *pstLf = obG.findNode(g_vLfE[iPos]);
                TPST_E *pstRt = obG.findNode(g_vRtE[iPos]);
                printf("%d (%d, %d) %d (%d, %d)\n",
                       pstLf->eid, pstLf->paXY.first, pstLf->paXY.second,
                       pstRt->eid, pstRt->paXY.first, pstRt->paXY.second);
                if ((iCurK >= pstLf->vSky.size()) || (iCurK >= pstRt->vSky.size()))
                {
                    /* not in it */
                    continue;
                }
                printf("gamma: %d %d\n", pstLf->vSky[iCurK], pstRt->vSky[iCurK]);
                int iMinR = min(pstLf->vSky[iCurK], pstRt->vSky[iCurK]);
                if (iMinR >= pstNode->vSky[iCurK])
                {
                    ++iNeighbors;
                }
            }
            g_bLock = false;
            printf("neighbors: %d\n", iNeighbors);
            ASSERT(0);
        }
    }
    return 0;
}

/*****************
input:
        myG &mapInitG
        LIST_DECOMP_G &lstDeG
        long double ldMiniP
output:
        LIST_DECOMP_G &lstDeG
description:
        truss decomposition
******************/
int deTruss_verify(myG &obG, vector<int> &vSubG, int iResR, int iCurK, double ldMiniP, double ldAccu, pair<int, int> &paXY)
{
	myG oSubG;

	for(int iEid : vSubG)
    {
        oSubG.add(obG.findNode(iEid)->paXY.first, obG.findNode(iEid)->paXY.second, obG.findNode(iEid)->p);
    }

	oSubG.init();

	vector<int> vNewE;
	for(int i = 1; i <= oSubG.m_iMaxEId; ++i)
    {
        vNewE.push_back(i);
    }
	/* calculate pek */
    calPek_part(oSubG, vNewE);

    /* init R */
    int iMinR = MAX_R(ldAccu);
    for (int iEid : vNewE)
    {
        oSubG.findNode(iEid)->iSelfR = getR(oSubG.findNode(iEid)->vldPek, oSubG.findNode(iEid)->p, ldMiniP, ldAccu, iCurK);
        //iMinR = min(iMinR, oSubG.findNode(iEid)->iSelfR);
        if (iMinR > oSubG.findNode(iEid)->iSelfR)
        {
            iMinR = oSubG.findNode(iEid)->iSelfR;
            paXY = oSubG.findNode(iEid)->paXY;
        }
    }

    //printf("VERIFY k: %d calculated R: %d real R: %d difference: %d\n", iCurK, iResR, iMinR, iMinR - iResR);

    if (abs(iMinR - iResR) > 1)
    {
        vector<int> vLfE;
        vector<int> vRtE;
        oSubG.findNotRmNeb(paXY.first, paXY.second, vLfE, vRtE);
        printf("SUBGRAPH neighbors: %d\n", vLfE.size());
        for (int i = 0; i < vLfE.size(); ++i)
        {
            printf("SUBGRAPH %d(%d, %d) %d(%d, %d)\n",
                   vLfE[i], oSubG.findNode(vLfE[i])->paXY.first, oSubG.findNode(vLfE[i])->paXY.second,
                   vRtE[i], oSubG.findNode(vRtE[i])->paXY.first, oSubG.findNode(vRtE[i])->paXY.second);
        }
    }

    return iMinR - iResR;
}

/*****************
input:
        myG &mapInitG
        LIST_DECOMP_G &lstDeG
        long double ldMiniP
output:
        LIST_DECOMP_G &lstDeG
description:
        truss decomposition
******************/
int deTruss_rePek(myG &obG, int iEid, int iCurK, bool bRmValid)
{
    TPST_E *pstNode = obG.findNode(iEid);
    int iDesSize = pstNode->vldPek.size();
    int iSup = pstNode->vLfE.size();

    ASSERT(!g_bLock);
    g_bLock = true;
    g_vLfE.clear();
    g_vRtE.clear();
    for (int i = 0; i < iSup; ++i)
    {
        TPST_E *pstLf = obG.findNode(pstNode->vLfE[i]);
        TPST_E *pstRt = obG.findNode(pstNode->vRtE[i]);

        if (bRmValid)
        {
            if (pstLf->bRm || pstRt->bRm)
            {
                continue;
            }
        }
        else
        {
            if (pstLf->bRm && pstLf->vSky[iCurK] < pstNode->vSky[iCurK])
            {
                continue;
            }
            if (pstRt->bRm && pstRt->vSky[iCurK] < pstNode->vSky[iCurK])
            {
                continue;
            }
        }
        if ((!pstLf->bLc) && (pstLf->vSky[iCurK] < pstNode->vSky[iCurK]))
        {
            continue;
        }
        if ((!pstRt->bLc) && (pstRt->vSky[iCurK] < pstNode->vSky[iCurK]))
        {
            continue;
        }
         g_vLfE.push_back(pstLf->eid);
         g_vRtE.push_back(pstRt->eid);

    }
    pstNode->vldPek.clear();
    pstNode->vldPek.resize(g_vLfE.size() + 1);
    cal_OneEdgePek(obG, g_vLfE, g_vRtE, pstNode->vldPek);

    g_bLock = false;

    if (iDesSize != pstNode->vldPek.size())
    {
#if 0
        map <int, int> mpTemp = pstNode->mpDebugPek;
        printf("ERROR CAL eid: %d destinate size: %d check size: %d real: %d neighbors: %d\n",
               pstNode->eid, iDesSize, pstNode->vldPek.size(), pstNode->mpDebugPek.size() / 2, iSup);
        for (int i = 0; i < iSup; ++i)
        {
            TPST_E *pstLf = obG.findNode(pstNode->vLfE[i]);
            TPST_E *pstRt = obG.findNode(pstNode->vRtE[i]);
            printf("CAL current R: %d gamma: %d neighbors %d %d p: %f %f rm: %d %d R: %d %d self R: %d %d local: %d %d\n",
                   pstNode->iSelfR, pstNode->vSky[iCurK],
                   pstLf->eid, pstRt->eid,
                   pstLf->p, pstRt->p,
                   pstLf->bRm, pstRt->bRm,
                   pstLf->vSky[iCurK], pstRt->vSky[iCurK],
                   pstLf->iSelfR, pstRt->iSelfR,
                   pstLf->bLc, pstRt->bLc);
            if (bRmValid)
            {
                if (pstLf->bRm || pstRt->bRm)
                {
                    continue;
                }
            }
            else
            {
                if (pstLf->bRm && pstLf->vSky[iCurK] < pstNode->vSky[iCurK])
                {
                    continue;
                }
                if (pstRt->bRm && pstRt->vSky[iCurK] < pstNode->vSky[iCurK])
                {
                    continue;
                }
            }
            if ((!pstLf->bLc) && (pstLf->vSky[iCurK] < pstNode->vSky[iCurK]))
            {
                continue;
            }
            if ((!pstRt->bLc) && (pstRt->vSky[iCurK] < pstNode->vSky[iCurK]))
            {
                continue;
            }
            printf("True\n");
            if (mpTemp.find(pstLf->eid) == mpTemp.end())
            {
                printf("Not in it: %d %d\n", pstLf->eid, pstRt->eid);
            }
            else
            {
                mpTemp.erase(pstLf->eid);
                mpTemp.erase(pstRt->eid);
            }

        }
        for (auto atN : mpTemp)
        {
            if (atN.first < atN.second)
            {
                printf("Remained: %d %d\n", atN.first, atN.second);
            }
        }
#endif
        ASSERT(0);
    }
    ASSERT(cal_check(pstNode->vldPek));

    return 0;
}
/*****************
input:
        myG &mapInitG
        LIST_DECOMP_G &lstDeG
        long double ldMiniP
output:
        LIST_DECOMP_G &lstDeG
description:
        truss decomposition
******************/
int init(myG &obG, int iEid, int iCurK, int iDesR)
{
    TPST_E *pstNode = obG.findNode(iEid);
    pstNode->vldPek.clear();
    //pstNode->mpDebugPek.clear();

    ASSERT(!g_bLock);
    g_bLock = true;
    g_vLfE.clear();
    g_vRtE.clear();

    pstNode->vLfE.clear();
    pstNode->vRtE.clear();

    vector<int> vLfE;
    vector<int> vRtE;

    obG.findNeb(pstNode->paXY.first, pstNode->paXY.second, g_vLfE, g_vRtE);

    for (int iPos = 0; iPos < g_vLfE.size(); ++iPos)
    {
        int iLfE = g_vLfE[iPos];
        int iRtE = g_vRtE[iPos];
        TPST_E *pstLf = obG.findNode(iLfE);
        TPST_E *pstRt = obG.findNode(iRtE);

        /* printf("INIT k: %d, iDesR: %d find %d %d pek: %d %d p: %lf %lf by %d\n",
               iCurK, iDesR, iLfE, iRtE, pstLf->vSky.size(), pstRt->vSky.size(),
               pstLf->p, pstRt->p, pstNode->eid); */
        if ((pstLf->p < obG.m_dMinP) || (pstRt->p < obG.m_dMinP))
        {
            /* not in it */
            continue;
        }
#if 0
        if (0 < iDesR)
        {
            if ((iCurK >= pstLf->vSky.size()) || (iCurK >= pstRt->vSky.size()))
            {
                /* not in it */
                continue;
            }
        }
        else
        {
            /* bottom layer */
            if ((iCurK > pstLf->iTrussness) || (iCurK > pstRt->iTrussness))
            {
                /* not in it */
                continue;
            }
            if (iCurK >= pstLf->vSky.size())
            {
                int iOldSize = pstLf->vSky.size();
                pstLf->vSky.resize(iCurK + 1);
                for (int iTpK = iOldSize; iTpK <= iCurK; ++iTpK)
                {
                    /* gamma = 0 */
                    pstLf->vSky[iTpK] = 0;
                }
            }
            if (iCurK >= pstRt->vSky.size())
            {
                int iOldSize = pstRt->vSky.size();
                pstRt->vSky.resize(iCurK + 1);
                for (int iTpK = iOldSize; iTpK <= iCurK; ++iTpK)
                {
                    /* gamma = 0 */
                    pstRt->vSky[iTpK] = 0;
                }
            }
        }
#endif

        if ((iCurK > pstLf->iTrussness) || (iCurK > pstRt->iTrussness))
        {
            /* not in it */
            continue;
        }
        if (iCurK >= pstLf->vSky.size())
        {
            int iOldSize = pstLf->vSky.size();
            pstLf->vSky.resize(iCurK + 1);
            for (int iTpK = iOldSize; iTpK <= iCurK; ++iTpK)
            {
                /* gamma = 0 */
                pstLf->vSky[iTpK] = 0;
            }
        }
        if (iCurK >= pstRt->vSky.size())
        {
            int iOldSize = pstRt->vSky.size();
            pstRt->vSky.resize(iCurK + 1);
            for (int iTpK = iOldSize; iTpK <= iCurK; ++iTpK)
            {
                /* gamma = 0 */
                pstRt->vSky[iTpK] = 0;
            }
        }
        pstNode->vLfE.push_back(iLfE);
        pstNode->vRtE.push_back(iRtE);
        if ((iDesR > pstLf->vSky[iCurK]) || (iDesR > pstRt->vSky[iCurK]))
        {
            continue;
        }
        //printf("INIT %d %d by %d\n", iLfE, iRtE, pstNode->eid);
        vLfE.push_back(iLfE);
        vRtE.push_back(iRtE);
        //pstNode->mpDebugPek[iLfE] = iRtE;
        //pstNode->mpDebugPek[iRtE] = iLfE;
    }

    g_bLock = false;
    cal_OneEdgePek(obG, vLfE, vRtE, pstNode->vldPek);

    pstNode->vldOldPek = pstNode->vldPek;
    //pstNode->mpDebugOldPek = pstNode->mpDebugPek;

    pstNode->bInit = true;
    obG.m_vInitE.push_back(pstNode->eid);
    /*printf("INIT k: %d eid: %d pek: %d real: %d desR: %d neighbors: %d\n",
           iCurK, pstNode->eid, pstNode->vldPek.size(),
           pstNode->mpDebugPek.size() / 2, iDesR, pstNode->vLfE.size());*/

    //ASSERT(pstNode->vldPek.size() - 1 == pstNode->mpDebugPek.size() / 2);

    /*if ((9761 == iEid) || (9411 == iEid) || (7813 == iEid))
    {
        for (int iPos = 0; iPos < pstNode->vLfE.size(); ++iPos)
        {
            TPST_E *pstLf = obG.findNode(pstNode->vLfE[iPos]);
            TPST_E *pstRt = obG.findNode(pstNode->vRtE[iPos]);
            printf("INIT neighbors: %d (%d, %d) %d (%d, %d) R: %d %d by %d\n",
                   pstLf->eid, pstLf->paXY.first, pstLf->paXY.second,
                   pstRt->eid, pstRt->paXY.first, pstRt->paXY.second,
                   pstLf->vSky[iCurK], pstRt->vSky[iCurK], iEid);
        }
        printf("INIT pek: %d neighbors: %d by %d\n",
               pstNode->vldPek.size(), pstNode->vLfE.size(), iEid);
        for (int iPos = 0; iPos < g_vLfE.size(); ++iPos)
        {
            TPST_E *pstLf = obG.findNode(g_vLfE[iPos]);
            TPST_E *pstRt = obG.findNode(g_vRtE[iPos]);
            printf("INIT neighbors: %d (%d, %d) %d (%d, %d) p: %f %f t: %d %d R: %d %d sky: %d %d init: %d %d queue: %d %d by %d\n",
                   pstLf->eid, pstLf->paXY.first, pstLf->paXY.second,
                   pstRt->eid, pstRt->paXY.first, pstRt->paXY.second,
                   pstLf->p, pstRt->p,
                   pstLf->iTrussness, pstRt->iTrussness,
                   pstLf->vSky[iCurK], pstRt->vSky[iCurK],
                   pstLf->vSky.size(), pstRt->vSky.size(),
                   pstLf->bInit, pstRt->bInit,
                   pstLf->bInQ, pstRt->bInQ,
                   iEid);
        }
        printf("INIT original neighbors: %d by %d\n",
               g_vLfE.size(), iEid);
    }*/
    return 0;
}
/*****************
input:
        myG &mapInitG
        LIST_DECOMP_G &lstDeG
        long double ldMiniP
output:
        LIST_DECOMP_G &lstDeG
description:
        truss decomposition
******************/
int insOneHop(myG &obG, int iCurK, const vector<pair<int, double> > &vNewIncE, vector<int> &vSeed)
{
    vector<int> vRawE;

    for (auto atN : vNewIncE)
    {
        TPST_E *pstNode = obG.findNode(atN.first);

        if (iCurK >= pstNode->vSky.size())
        {
            if (pstNode->iTrussness >= iCurK)
            {
                int iOldSize = pstNode->vSky.size();
                pstNode->vSky.resize(iCurK + 1);
                for (int iTpK = iOldSize; iTpK <= iCurK; ++iTpK)
                {
                    /* gamma = 0 */
                    pstNode->vSky[iTpK] = 0;
                }
            }
            else
            {
                /* not in it */
                continue;
            }
        }
        vRawE.push_back(pstNode->eid);
    }
    /* sort by R */
    sort(vRawE.begin(), vRawE.end(),
        [&obG, &iCurK](const int e1, const int e2) {
         if (obG.findNode(e1)->vSky[iCurK] == obG.findNode(e2)->vSky[iCurK]) return e1 < e2;
         else return obG.findNode(e1)->vSky[iCurK] < obG.findNode(e2)->vSky[iCurK];
        });

    int iPrevR = 0;
    for (int iEid : vRawE)
    {
        TPST_E *pstNode = obG.findNode(iEid);

        int iOldR = pstNode->vSky[iCurK];
        if (iOldR < iPrevR)
        {
            printf("ERROR sort R: %d %d\n", iPrevR, iOldR);
            ASSERT(0);
        }
        iPrevR = iOldR;
        if (!pstNode->bInit)
        {
            init(obG, pstNode->eid, iCurK, iOldR);
        }
        int iNewR = getR(pstNode->vldPek, pstNode->p, obG.m_dMinP, obG.m_dAccu, iCurK);
        /*printf("INS_ONE_HOP root eid: %d (%d, %d) old R: %d new R: %d\n",
               pstNode->eid, pstNode->paXY.first, pstNode->paXY.second, iOldR, iNewR);*/

        if (iNewR > iOldR)
        {
            vSeed.push_back(pstNode->eid);
        }
        /* neighbors */
        for (int iPos = 0; iPos < pstNode->vLfE.size(); ++iPos)
        {
            int iLfE = pstNode->vLfE[iPos];
            int iRtE = pstNode->vRtE[iPos];
            TPST_E *pstLf = obG.findNode(iLfE);
            TPST_E *pstRt = obG.findNode(iRtE);

            /*printf("INS_ONE_HOP neighbors eid: %d %d (%d, %d) (%d, %d)\n",
                   pstLf->eid, pstRt->eid,
                   pstLf->paXY.first, pstLf->paXY.second,
                   pstRt->paXY.first, pstRt->paXY.second);*/

            iOldR = pstLf->vSky[iCurK];
            if (!pstLf->bInit)
            {
                init(obG, pstLf->eid, iCurK, iOldR);
            }
            iNewR = getR(pstLf->vldPek, pstLf->p, obG.m_dMinP, obG.m_dAccu, iCurK);
            //printf("INS_ONE_HOP left old R: %d new R: %d\n", iOldR, iNewR);
            if (iNewR > iOldR)
            {
                vSeed.push_back(pstLf->eid);
            }

            iOldR = pstRt->vSky[iCurK];
            if (!pstRt->bInit)
            {
                init(obG, pstRt->eid, iCurK, iOldR);
            }
            iNewR = getR(pstRt->vldPek, pstRt->p, obG.m_dMinP, obG.m_dAccu, iCurK);
            //printf("INS_ONE_HOP right old R: %d new R: %d\n", iOldR, iNewR);
            if (iNewR > iOldR)
            {
                vSeed.push_back(pstRt->eid);
            }
        }
    }
    return 0;
}

/*****************
input:
        myG &mapInitG
        LIST_DECOMP_G &lstDeG
        long double ldMiniP
output:
        LIST_DECOMP_G &lstDeG
description:
        truss decomposition
******************/
int decOneHop(myG &obG, int iCurK, const vector<pair<int, double> > &vNewDecE, vector<int> &vSeed)
{
    vector<int> vRawE;

    for (auto atN : vNewDecE)
    {
        TPST_E *pstNode = obG.findNode(atN.first);

        if (iCurK >= pstNode->vSky.size())
        {
            /* not in it */
            continue;
        }
        if (atN.second < obG.m_dMinP)
        {
            /* cannot change */
            continue;
        }
        if (pstNode->vSky[iCurK] <= 0)
        {
            /* cannot change */
            continue;
        }
        vRawE.push_back(pstNode->eid);
    }
    /* sort by R, decreasing order */
    sort(vRawE.begin(), vRawE.end(),
        [&obG, &iCurK](const int e1, const int e2) {
         if (obG.findNode(e1)->vSky[iCurK] == obG.findNode(e2)->vSky[iCurK]) return e1 > e2;
         else return obG.findNode(e1)->vSky[iCurK] > obG.findNode(e2)->vSky[iCurK];
        });
    for (int iEid : vRawE)
    {
        TPST_E *pstNode = obG.findNode(iEid);

        int iOldR = pstNode->vSky[iCurK];
        if (!pstNode->bInit)
        {
            init(obG, pstNode->eid, iCurK, iOldR);
        }
        int iNewR = getR(pstNode->vldPek, pstNode->p, obG.m_dMinP, obG.m_dAccu, iCurK);
        /*printf("INS_ONE_HOP root eid: %d (%d, %d) old R: %d new R: %d\n",
               pstNode->eid, pstNode->paXY.first, pstNode->paXY.second, iOldR, iNewR);*/

        if (iNewR < iOldR)
        {
            vSeed.push_back(pstNode->eid);
        }
        /* neighbors */
        for (int iPos = 0; iPos < pstNode->vLfE.size(); ++iPos)
        {
            int iLfE = pstNode->vLfE[iPos];
            int iRtE = pstNode->vRtE[iPos];
            TPST_E *pstLf = obG.findNode(iLfE);
            TPST_E *pstRt = obG.findNode(iRtE);

            /*printf("INS_ONE_HOP neighbors eid: %d %d (%d, %d) (%d, %d)\n",
                   pstLf->eid, pstRt->eid,
                   pstLf->paXY.first, pstLf->paXY.second,
                   pstRt->paXY.first, pstRt->paXY.second);*/

            iOldR = pstLf->vSky[iCurK];
            if (!pstLf->bInit)
            {
                init(obG, pstLf->eid, iCurK, iOldR);
            }
            iNewR = getR(pstLf->vldPek, pstLf->p, obG.m_dMinP, obG.m_dAccu, iCurK);
            //printf("INS_ONE_HOP left old R: %d new R: %d\n", iOldR, iNewR);
            if (iNewR < iOldR)
            {
                vSeed.push_back(pstLf->eid);
            }

            iOldR = pstRt->vSky[iCurK];
            if (!pstRt->bInit)
            {
                init(obG, pstRt->eid, iCurK, iOldR);
            }
            iNewR = getR(pstRt->vldPek, pstRt->p, obG.m_dMinP, obG.m_dAccu, iCurK);
            //printf("INS_ONE_HOP right old R: %d new R: %d\n", iOldR, iNewR);
            if (iNewR < iOldR)
            {
                vSeed.push_back(pstRt->eid);
            }
        }
    }
    return 0;
}
/*****************
input:
        myG &mapInitG
        LIST_DECOMP_G &lstDeG
        long double ldMiniP
output:
        LIST_DECOMP_G &lstDeG
description:
        truss decomposition
******************/
int expend(myG &obG, int iCurK, int iCurR, vector<int> &vSeed, int *piLB, int *piUB)
{
    vector<int> vWait = vSeed;
    for (int iEid : vWait)
    {
        TPST_E *pstNode = obG.findNode(iEid);
        pstNode->bBFS = true;
        //printf("EXPEND seed: %d init: %d\n", iEid, pstNode->bInit);
        ASSERT(pstNode->vSky.size() > iCurK);
        ASSERT(pstNode->p >= obG.m_dMinP);
    }

    *piLB = 0;
    *piUB = MAX_R(obG.m_dAccu);

    while (!vWait.empty())
    {
        int iEid = vWait.back();
        vWait.pop_back();
        TPST_E *pstNode = obG.findNode(iEid);
        //printf("EXPEND pop: %d init: %d\n", iEid, pstNode->bInit);

        /*g_vLfE.clear();
        g_vRtE.clear();

        obG.findNeb(pstNode->paXY.first, pstNode->paXY.second, g_vLfE, g_vRtE);*/

        ASSERT(pstNode->bInit);

        for (int iPos = 0; iPos < pstNode->vLfE.size(); ++iPos)
        {
            int iLfE = pstNode->vLfE[iPos];
            int iRtE = pstNode->vRtE[iPos];
            TPST_E *pstLf = obG.findNode(iLfE);
            TPST_E *pstRt = obG.findNode(iRtE);

            if ((iCurK >= pstLf->vSky.size()) || (iCurK >= pstRt->vSky.size()))
            {
                /* not in it */
                continue;
            }
            if ((pstLf->p < obG.m_dMinP) || (pstRt->p < obG.m_dMinP))
            {
                continue;
            }
            /*printf("EXPEND current eid: %d R: %d neighbors: %d %d R: %d %d BFS flag: %d %d\n",
                   iEid, iCurR, iLfE, iRtE,
                   pstLf->vSky[iCurK], pstRt->vSky[iCurK],
                   pstLf->bBFS, pstRt->bBFS);*/
            int iMinR = min(pstLf->vSky[iCurK], pstRt->vSky[iCurK]);

            if (iMinR == iCurR)
            {
                /* extend */
                if ((pstLf->vSky[iCurK] == iCurR) && (!pstLf->bBFS))
                {
                    pstLf->bBFS = true;
                    vWait.push_back(iLfE);
                    vSeed.push_back(iLfE);
                    if (!pstLf->bInit)
                    {
                        init(obG, pstLf->eid, iCurK, iCurR);
                    }
                }
                if ((pstRt->vSky[iCurK] == iCurR) && (!pstRt->bBFS))
                {
                    pstRt->bBFS = true;
                    vWait.push_back(iRtE);
                    vSeed.push_back(iRtE);
                    if (!pstRt->bInit)
                    {
                        init(obG, pstRt->eid, iCurK, iCurR);
                    }
                }

                /* upper bound */
                if (pstLf->vSky[iCurK] > pstRt->vSky[iCurK])
                {
                    *piUB = min(*piUB, pstLf->vSky[iCurK]);
                }
                else if (pstRt->vSky[iCurK] > pstLf->vSky[iCurK])
                {
                    *piUB = min(*piUB, pstRt->vSky[iCurK]);
                }
            }
            else if (iMinR > iCurR)
            {
                /* upper bound */
                *piUB = min(*piUB, iMinR);
            }
            else
            {
                /* lower bound */
                *piLB = max(*piLB, iMinR);
            }
        }

    }

    /* restore */
    for (int iEid : vSeed)
    {
        TPST_E *pstNode = obG.findNode(iEid);
        pstNode->bBFS = false;
    }
    return 0;
}


/*****************
input:
        myG &mapInitG
        LIST_DECOMP_G &lstDeG
        long double ldMiniP
output:
        LIST_DECOMP_G &lstDeG
description:
        truss decomposition
******************/
int changeR(myG &obG, int iCurK, int iEid, int iOldR, int iNewR)
{
    TPST_E *pstNode = obG.findNode(iEid);

    map<int, int> mpRecal;
    //int iOldR = pstNode->vSky[iCurK];

    pstNode->vSky[iCurK] = iNewR;
    /*printf("CHANGE_R save %d (%d, %d) old R: %d new R: %d pek: %d old pek: %d real: %d old real: %d\n",
           iEid, pstNode->paXY.first, pstNode->paXY.second, iOldR, iNewR,
           pstNode->vldPek.size(), pstNode->vldOldPek.size(),
           pstNode->mpDebugPek.size() / 2, pstNode->mpDebugOldPek.size() / 2);*/

    if (iNewR > iOldR)
    {
        /* restore pek */
        /* printf("CHANGE_R restore %d neighbors real: %d to %d\n",
               pstNode->eid, pstNode->mpDebugPek.size() / 2, pstNode->mpDebugOldPek.size() / 2); */
        pstNode->vldPek = pstNode->vldOldPek;
        //pstNode->mpDebugPek = pstNode->mpDebugOldPek;
        for (int i = 0; i < pstNode->vLfE.size(); ++i)
        {
            const int e1 = pstNode->vLfE[i];
            const int e2 = pstNode->vRtE[i];
            /*printf("CHANGE_R DEBUG neighbors: %d %d R: %d %d rm: %d %d by %d real: %d\n", e1, e2,
                   obG.findNode(e1)->vSky[iCurK], obG.findNode(e2)->vSky[iCurK],
                   obG.findNode(e1)->bRm, obG.findNode(e2)->bRm,
                   iEid, pstNode->mpDebugPek.size() / 2);*/
            int iMinR = min(obG.findNode(e1)->vSky[iCurK], obG.findNode(e2)->vSky[iCurK]);
            //if ((iMinR < iOldR) || (iMinR > iNewR)) continue;
            if (iMinR > iNewR) continue;
            int iBeforeMinR = iOldR;
            if (!obG.findNode(e1)->bLc)
            {
                iBeforeMinR = min(iBeforeMinR, obG.findNode(e1)->vSky[iCurK]);
            }
            if (!obG.findNode(e2)->bLc)
            {
                iBeforeMinR = min(iBeforeMinR, obG.findNode(e2)->vSky[iCurK]);
            }
            if (iBeforeMinR < iOldR) continue;

            /* self */
            if ((obG.findNode(e1)->bLc) && (obG.findNode(e2)->bLc))
            {
                if (obG.findNode(e1)->bRm && obG.findNode(e2)->bRm)
                {
                    if (iMinR < iNewR)
                    {
                        /* peel */
                        //pstNode->mpDebugPek.erase(e1);
                        //pstNode->mpDebugPek.erase(e2);
                        bool bRet = calPek_rm(obG, iEid, e1, e2, pstNode->vldPek, iCurK);
                        if (!bRet)
                        {
                            mpRecal[iEid];
                        }
                        //printf("CHANGE_R rm %d pek: %d\n", pstNode->eid, pstNode->vldPek.size());
                    }
                }
                else if (obG.findNode(e1)->bRm)
                {
                    if (obG.findNode(e1)->vSky[iCurK] < iNewR)
                    {
                        /* peel */
                        //pstNode->mpDebugPek.erase(e1);
                        //pstNode->mpDebugPek.erase(e2);
                        bool bRet = calPek_rm(obG, iEid, e1, e2, pstNode->vldPek, iCurK);
                        if (!bRet)
                        {
                            mpRecal[iEid];
                        }
                        //printf("CHANGE_R rm %d pek: %d\n", pstNode->eid, pstNode->vldPek.size());
                    }
                }
                else if (obG.findNode(e2)->bRm)
                {
                    if (obG.findNode(e2)->vSky[iCurK] < iNewR)
                    {
                        /* peel */
                        //pstNode->mpDebugPek.erase(e1);
                        //pstNode->mpDebugPek.erase(e2);
                        bool bRet = calPek_rm(obG, iEid, e1, e2, pstNode->vldPek, iCurK);
                        if (!bRet)
                        {
                            mpRecal[iEid];
                        }
                        //printf("CHANGE_R rm %d pek: %d\n", pstNode->eid, pstNode->vldPek.size());
                    }
                }
            }
            else if (obG.findNode(e1)->bLc)
            {
                if (obG.findNode(e1)->bRm)
                {
                    if (iMinR < iNewR)
                    {
                        /* peel */
                        //pstNode->mpDebugPek.erase(e1);
                        //pstNode->mpDebugPek.erase(e2);
                        bool bRet = calPek_rm(obG, iEid, e1, e2, pstNode->vldPek, iCurK);
                        if (!bRet)
                        {
                            mpRecal[iEid];
                        }
                        //printf("CHANGE_R rm %d pek: %d\n", pstNode->eid, pstNode->vldPek.size());
                    }
                }
                else
                {
                    /* neighbors */
                    if (obG.findNode(e2)->bInit && (iNewR >= obG.findNode(e2)->vSky[iCurK]))
                    {
                        if (iEid < e1)
                        {
                            /* add neighbors */
                            //obG.findNode(e2)->mpDebugPek[iEid] = e1;
                            //obG.findNode(e2)->mpDebugPek[e1] = iEid;
                            bool bRet = calPek_add(obG, e2, iEid, e1, obG.findNode(e2)->vldPek, iCurK);
                            if (!bRet)
                            {
                                mpRecal[e2];
                            }
                            obG.findNode(e2)->vldOldPek = obG.findNode(e2)->vldPek;
                            //obG.findNode(e2)->mpDebugOldPek = obG.findNode(e2)->mpDebugPek;

                            /*printf("CHANGE_R %d add neighbors: %d real: %d\n", e2,
                                   obG.findNode(e2)->vldPek.size() - 1,
                                   obG.findNode(e2)->mpDebugPek.size() / 2);*/
                        }
                    }
                }
            }
            else if (obG.findNode(e2)->bLc)
            {
                if (obG.findNode(e2)->bRm)
                {
                    if (iMinR < iNewR)
                    {
                        /* peel */
                        //pstNode->mpDebugPek.erase(e1);
                        //pstNode->mpDebugPek.erase(e2);
                        bool bRet = calPek_rm(obG, iEid, e1, e2, pstNode->vldPek, iCurK);
                        if (!bRet)
                        {
                            mpRecal[iEid];
                        }
                        //printf("CHANGE_R rm %d pek: %d\n", pstNode->eid, pstNode->vldPek.size());
                    }
                }
                else
                {
                    /* neighbors */
                    if (obG.findNode(e1)->bInit && (iNewR >= obG.findNode(e1)->vSky[iCurK]))
                    {
                        if (iEid < e2)
                        {
                            /* add neighbors */
                            //obG.findNode(e1)->mpDebugPek[iEid] = e2;
                            //obG.findNode(e1)->mpDebugPek[e2] = iEid;
                            bool bRet = calPek_add(obG, e1, iEid, e2, obG.findNode(e1)->vldPek, iCurK);
                            if (!bRet)
                            {
                                mpRecal[e1];
                            }

                            obG.findNode(e1)->vldOldPek = obG.findNode(e1)->vldPek;
                            //obG.findNode(e1)->mpDebugOldPek = obG.findNode(e1)->mpDebugPek;

                            /*printf("CHANGE_R %d add neighbors: %d real: %d\n", e1,
                                   obG.findNode(e1)->vldPek.size() - 1,
                                   obG.findNode(e1)->mpDebugPek.size() / 2);*/
                        }
                    }
                }
            }
            else
            {
                /* neighbors */
                if (iNewR >= iMinR)
                {
                    if (obG.findNode(e1)->bInit && (obG.findNode(e1)->vSky[iCurK] == iMinR))
                    {
                        /* add neighbors */
                        //obG.findNode(e1)->mpDebugPek[iEid] = e2;
                        //obG.findNode(e1)->mpDebugPek[e2] = iEid;
                        bool bRet = calPek_add(obG, e1, iEid, e2, obG.findNode(e1)->vldPek, iCurK);
                        if (!bRet)
                        {
                            mpRecal[e1];
                        }

                        obG.findNode(e1)->vldOldPek = obG.findNode(e1)->vldPek;

                        //obG.findNode(e1)->mpDebugOldPek = obG.findNode(e1)->mpDebugPek;

                        /*printf("CHANGE_R %d add neighbors: %d real: %d\n", e1,
                               obG.findNode(e1)->vldPek.size() - 1,
                               obG.findNode(e1)->mpDebugPek.size() / 2);*/
                    }
                    if (obG.findNode(e2)->bInit && (obG.findNode(e2)->vSky[iCurK] == iMinR))
                    {
                        /* add neighbors */
                        //obG.findNode(e2)->mpDebugPek[iEid] = e1;
                        //obG.findNode(e2)->mpDebugPek[e1] = iEid;
                        bool bRet = calPek_add(obG, e2, iEid, e1, obG.findNode(e2)->vldPek, iCurK);
                        if (!bRet)
                        {
                            mpRecal[e2];
                        }
                        obG.findNode(e2)->vldOldPek = obG.findNode(e2)->vldPek;

                        //obG.findNode(e2)->mpDebugOldPek = obG.findNode(e2)->mpDebugPek;

                        /*printf("CHANGE_R %d add neighbors: %d real: %d\n", e2,
                               obG.findNode(e2)->vldPek.size() - 1,
                               obG.findNode(e2)->mpDebugPek.size() / 2);*/
                    }
                }
            }
#if 0
            if (obG.findNode(e1)->bInit && (!obG.findNode(e1)->bLc) &&
                (obG.findNode(e1)->vSky[iCurK] == iMinR))
            {
                calPek_add(obG, e1, iEid, e2, obG.findNode(e1)->vldPek);
                //printf("CHANGE_R add %d pek: %d by %d\n", e1, obG.findNode(e1)->vldPek.size(), iEid);
            }

            if (obG.findNode(e2)->bInit && (!obG.findNode(e2)->bLc) &&
                (obG.findNode(e2)->vSky[iCurK] == iMinR))
            {
                calPek_add(obG, e2, iEid, e1, obG.findNode(e2)->vldPek);
                //printf("CHANGE_R add %d pek: %d by %d\n", e2, obG.findNode(e2)->vldPek.size(), iEid);
            }
#endif
        }
    }
    else if (iNewR < iOldR)
    {
        //pstNode->bInit = false;

        /* restore pek */
        /*printf("CHANGE_R restore %d neighbors real: %d to %d R: %d %d rm: %d %d by %d real: %d\n",
               pstNode->eid, pstNode->mpDebugPek.size() / 2, pstNode->mpDebugOldPek.size() / 2);*/
        //pstNode->mpDebugPek = pstNode->mpDebugOldPek;
        pstNode->vldPek = pstNode->vldOldPek;
        /* increase self pek */
        for (int i = 0; i < pstNode->vLfE.size(); ++i)
        {
            const int e1 = pstNode->vLfE[i];
            const int e2 = pstNode->vRtE[i];
            /*printf("CHANGE_R DEBUG neighbors: %d %d R: %d %d rm: %d %d by %d\n", e1, e2,
                   obG.findNode(e1)->vSky[iCurK], obG.findNode(e2)->vSky[iCurK],
                   obG.findNode(e1)->bRm, obG.findNode(e2)->bRm,
                   iEid);*/
            /* not local edges */
            int iOldMinR = iOldR;
            if (!obG.findNode(e1)->bLc)
            {
                iOldMinR = min(obG.findNode(e1)->vSky[iCurK], iOldMinR);
            }
            if (!obG.findNode(e2)->bLc)
            {
                iOldMinR = min(obG.findNode(e2)->vSky[iCurK], iOldMinR);
            }
            //int iMinR = min(obG.findNode(e1)->vSky[iCurK], obG.findNode(e2)->vSky[iCurK]);
            if ((iOldMinR < iOldR) && (iOldMinR >= iNewR))
            {
                //pstNode->mpDebugPek[e2] = e1;
                //pstNode->mpDebugPek[e1] = e2;
                bool bRet = calPek_add(obG, iEid, e1, e2, pstNode->vldPek, iCurK);
                if (!bRet)
                {
                    mpRecal[iEid];
                }
                //printf("CHANGE_R add %d pek: %d\n", pstNode->eid, pstNode->vldPek.size());
            }

            /* maybe some peers degrade */
            int iNewMinR = min(obG.findNode(e1)->vSky[iCurK], obG.findNode(e2)->vSky[iCurK]);

            if ((iOldMinR >= iOldR) && (iNewMinR < iNewR))
            {
                //pstNode->mpDebugPek.erase(e1);
                //pstNode->mpDebugPek.erase(e2);
                bool bRet = calPek_rm(obG, iEid, e1, e2, pstNode->vldPek, iCurK);
                if (!bRet)
                {
                    mpRecal[iEid];
                }
                //printf("CHANGE_R rm %d pek: %d\n", pstNode->eid, pstNode->vldPek.size());
            }
        }
    }
    else
    {
        /* maybe some peers degrade */
        /* restore pek */
        /*printf("CHANGE_R restore %d neighbors real: %d to %d R: %d %d rm: %d %d by %d real: %d\n",
               pstNode->eid, pstNode->mpDebugPek.size() / 2, pstNode->mpDebugOldPek.size() / 2);*/
        //pstNode->mpDebugPek = pstNode->mpDebugOldPek;
        pstNode->vldPek = pstNode->vldOldPek;
        /* decrease self pek */
        for (int i = 0; i < pstNode->vLfE.size(); ++i)
        {
            const int e1 = pstNode->vLfE[i];
            const int e2 = pstNode->vRtE[i];
            /*printf("CHANGE_R DEBUG neighbors: %d %d R: %d %d rm: %d %d by %d\n", e1, e2,
                   obG.findNode(e1)->vSky[iCurK], obG.findNode(e2)->vSky[iCurK],
                   obG.findNode(e1)->bRm, obG.findNode(e2)->bRm,
                   iEid);*/

            int iOldMinR = iOldR;
            if (!obG.findNode(e1)->bLc)
            {
                iOldMinR = min(obG.findNode(e1)->vSky[iCurK], iOldMinR);
            }
            if (!obG.findNode(e2)->bLc)
            {
                iOldMinR = min(obG.findNode(e2)->vSky[iCurK], iOldMinR);
            }
            int iNewMinR = min(obG.findNode(e1)->vSky[iCurK], obG.findNode(e2)->vSky[iCurK]);

            //int iMinR = min(obG.findNode(e1)->vSky[iCurK], obG.findNode(e2)->vSky[iCurK]);
            if ((iOldMinR >= iOldR) && (iNewMinR < iOldR))
            {
                //pstNode->mpDebugPek.erase(e1);
                //pstNode->mpDebugPek.erase(e2);
                bool bRet = calPek_rm(obG, iEid, e1, e2, pstNode->vldPek, iCurK);
                if (!bRet)
                {
                    mpRecal[iEid];
                }
                //printf("CHANGE_R rm %d pek: %d\n", pstNode->eid, pstNode->vldPek.size());
            }
        }
    }

    for (auto atN : mpRecal)
    {
        if (obG.findNode(atN.first)->bLc && (!obG.findNode(atN.first)->bRm))
        {
            deTruss_rePek(obG, atN.first, iCurK, true);
        }
        else
        {
            deTruss_rePek(obG, atN.first, iCurK, false);
        }
    }

    pstNode->vldOldPek = pstNode->vldPek;
    //pstNode->mpDebugOldPek = pstNode->mpDebugPek;

    /* debug */
#if 0
    int iDebugValidCnt = 0;
    bool bError = false;
    for (int i = 0; i < pstNode->vLfE.size(); ++i)
    {
        const int e1 = pstNode->vLfE[i];
        const int e2 = pstNode->vRtE[i];
        /*printf("CHANGE_R DEBUG neighbors: %d %d R: %d %d rm: %d %d by %d\n", e1, e2,
               obG.findNode(e1)->vSky[iCurK], obG.findNode(e2)->vSky[iCurK],
               obG.findNode(e1)->bRm, obG.findNode(e2)->bRm,
               iEid);*/

        /*if (obG.findNode(e1)->bInit &&
            (obG.findNode(e1)->vldPek.size() - 1 != obG.findNode(e1)->mpDebugPek.size() / 2))
        {
            bError = true;
        }
        if (obG.findNode(e2)->bInit &&
            (obG.findNode(e2)->vldPek.size() - 1 != obG.findNode(e2)->mpDebugPek.size() / 2))
        {
            bError = true;
        }*/

        int iNewMinR = iNewR;
#if 0
        if (obG.findNode(e1)->bRm)
        {
            d
        }
        else
        {
            if (obG.findNode(e1)->bLc)
            {
                /* infinity */
            }
            else
            {

            }
        }
#endif
        if ((!obG.findNode(e1)->bRm) && (obG.findNode(e1)->bLc))
        {
            /* infinity */
        }
        else
        {
            iNewMinR = min(iNewMinR, obG.findNode(e1)->vSky[iCurK]);
        }

        if ((!obG.findNode(e2)->bRm) && (obG.findNode(e2)->bLc))
        {
            /* infinity */
        }
        else
        {
            iNewMinR = min(iNewMinR, obG.findNode(e2)->vSky[iCurK]);
        }

        //int iMinR = min(obG.findNode(e1)->vSky[iCurK], obG.findNode(e2)->vSky[iCurK]);
        if (iNewMinR >= iNewR)
        {
            ++iDebugValidCnt;
        }
    }
    if ((iDebugValidCnt != pstNode->vldPek.size() - 1) || (pstNode->vldPek.size() - 1 != pstNode->mpDebugPek.size() / 2))
    {
        bError = true;
    }
    if (bError)
    {
        map <int, int> mpTemp = pstNode->mpDebugPek;
        for (int i = 0; i < pstNode->vLfE.size(); ++i)
        {
            const int e1 = pstNode->vLfE[i];
            const int e2 = pstNode->vRtE[i];
            printf("CHANGE_R DEBUG neighbors: %d %d R: %d %d rm: %d %d local: %d %d by %d\n", e1, e2,
                   obG.findNode(e1)->vSky[iCurK], obG.findNode(e2)->vSky[iCurK],
                   obG.findNode(e1)->bRm, obG.findNode(e2)->bRm,
                   obG.findNode(e1)->bLc, obG.findNode(e2)->bLc,
                   iEid);
            int iNewMinR = iNewR;

            if ((!obG.findNode(e1)->bRm) && (obG.findNode(e1)->bLc))
            {
                /* infinity */
            }
            else
            {
                iNewMinR = min(iNewMinR, obG.findNode(e1)->vSky[iCurK]);
            }

            if ((!obG.findNode(e2)->bRm) && (obG.findNode(e2)->bLc))
            {
                /* infinity */
            }
            else
            {
                iNewMinR = min(iNewMinR, obG.findNode(e2)->vSky[iCurK]);
            }
            if (iNewMinR >= iNewR)
            {
                printf("True\n");

                if (mpTemp.find(e1) == mpTemp.end())
                {
                    printf("Not in it: %d %d\n", e1, e2);
                }
                else
                {
                    mpTemp.erase(e1);
                    mpTemp.erase(e2);
                }
            }
        }
        for (auto atN : mpTemp)
        {
            if (atN.first < atN.second)
            {
                printf("Remained: %d %d\n", atN.first, atN.second);
            }
        }
        printf("CHANGE_R ERROR cnt: %d get: %d real: %d oldR: %d newR: %d\n",
               iDebugValidCnt, pstNode->vldPek.size() - 1, pstNode->mpDebugPek.size() / 2,
               iOldR, iNewR);
        ASSERT(0);
    }
#endif
    return 0;
}
/*****************
input:
        myG &mapInitG
        LIST_DECOMP_G &lstDeG
        long double ldMiniP
output:
        LIST_DECOMP_G &lstDeG
description:
        truss decomposition
******************/
int peel(myG &obG, int iCurK, vector<int> &vSeed, int iLB, int iUB, vector<int> &vLRemain, vector<int> &vURemain)
{
    int iUpR = iUB;
    vector<int> vERaw;
    /* init */
    for (int iEid : vSeed)
    {
        TPST_E *pstNode = obG.findNode(iEid);
        ASSERT(pstNode->p >= obG.m_dMinP);

        vERaw.push_back(iEid);
    }

    auto cmp = [](pair<int, int> &left, pair<int, int> &right) { return left.second > right.second; };
    priority_queue<int, vector<pair<int, int> >, decltype(cmp)> prQ(cmp);

    int iDesR = 0;
    if (!vERaw.empty())
    {
        iDesR = obG.findNode(vERaw.front())->vSky[iCurK];
        //printf("PEEL first R: %d eid: %d\n", prQ.top().second, prQ.top().first);
    }
    /* init vector */
    for (int iEid : vERaw)
    {
        TPST_E *pstNode = obG.findNode(iEid);
        pstNode->bLc = true;
        /*if (!pstNode->bInit)
        {
            init(obG, pstNode->eid, iCurK, pstNode->vSky[iCurK]);
        }*/
        ASSERT(pstNode->bInit);
        pstNode->iSelfR = getR(pstNode->vldPek, pstNode->p, obG.m_dMinP, obG.m_dAccu, iCurK);
        prQ.push({iEid, pstNode->iSelfR});

        /* temporarily set max R */
        pstNode->vSky[iCurK] = iUB;
        /* printf("PEEL add %d R: %d pek: %d neighbors: %d real: %d LB: %d UB: %d\n",
               iEid, pstNode->iSelfR, pstNode->vldPek.size(), pstNode->vLfE.size(),
               pstNode->mpDebugPek.size() / 2, iLB, iUB); */
#if 0
        if (96132 == iEid)
        {
            /*for (int iPos = 0; iPos < pstNode->vLfE.size(); ++iPos)
            {
                TPST_E *pstLf = obG.findNode(pstNode->vLfE[iPos]);
                TPST_E *pstRt = obG.findNode(pstNode->vRtE[iPos]);
                printf("PEEL neighbors: %d (%d, %d) %d (%d, %d) R: %d %d\n",
                       pstLf->eid, pstLf->paXY.first, pstLf->paXY.second,
                       pstRt->eid, pstRt->paXY.first, pstRt->paXY.second,
                       pstLf->vSky[iCurK], pstRt->vSky[iCurK]);
            }*/
            printf("PEEL pek: %d neighbors: %d\n",
                   pstNode->vldPek.size(), pstNode->vLfE.size());
            /*for (auto atP : pstNode->vldPek)
            {
                printf("PEEL pek: %Lf\n", atP);
            }*/
        }
#endif
    }

    //printf("DE_TRUSS first raw: %d LB: %d UB: %d\n", vERaw.size(), iLB, iUB);

    //vector<int> vPT(obG.m_iMaxEId + 1, -1);
    //vector<bool> removed(obG.m_iMaxEId + 1, false);
    int iCurR = 0;
    int iRmCnt = 0;

    vector<int> vWait;
    // 2.2.1. process the edges layer by layer
    while (iRmCnt < vERaw.size())
    {
        if (iCurR >= iUpR)
        {
            /* end */
            break;
        }
        vWait.clear();
        while (!prQ.empty())
        {
            int iEid = prQ.top().first;
            if (prQ.top().second > iCurR)
            {
                /* end */
                break;
            }
            TPST_E *pstNode = obG.findNode(iEid);
            if (!pstNode->bRm)
            {
                vWait.push_back(iEid);

                if (iCurR < iLB)
                {
                    pstNode->iPT = iLB;
                    vLRemain.push_back(iEid);
                }
                else
                {
                    pstNode->iPT = iCurR;
                }
            }
            prQ.pop();
        }

        if (vWait.empty())
        {
          if (prQ.empty())
          {
              break;
          }
          iCurR = prQ.top().second;
          //printf("PEEL get new R: %d\n", iCurR);

          /* debug */
          /*if ((7 == iCurK) && (185 <= iCurR))
          {
            vector<int> vLfE;
            vector<int> vRtE;
            TPST_E *pstNode = obG.findNode(94066);
            int iOldSize = pstNode->vldPek.size();

            //printf("start rm triangles\n");
            for (int i = 0; i < pstNode->vLfE.size(); ++i)
            {
                const int e1 = pstNode->vLfE[i];
                const int e2 = pstNode->vRtE[i];
                printf("DEBUG current R: %d LB: %d neighbors %d %d rm: %d %d R: %d %d local: %d %d\n",
                       iCurR, iLB, e1, e2, obG.findNode(e1)->bRm, obG.findNode(e2)->bRm,
                       obG.findNode(e1)->vSky[iCurK], obG.findNode(e2)->vSky[iCurK],
                       obG.findNode(e1)->bLc, obG.findNode(e2)->bLc);

                if (obG.findNode(e1)->bRm || obG.findNode(e2)->bRm) continue;
                int iMinR = min(obG.findNode(e1)->vSky[iCurK], obG.findNode(e2)->vSky[iCurK]);
                if (iMinR <= iLB) continue;
                vLfE.push_back(e1);
                vRtE.push_back(e2);
                printf("True\n");
            }
            cal_OneEdgePek(obG, vLfE, vRtE, pstNode->vldPek);
            int iTpR = getR(pstNode->vldPek, pstNode->p, obG.m_dMinP, obG.m_dAccu, iCurK);
            printf("PEEL DEBUG %d R: %d old R: %d pek: %d old pek: %d neighbors: %d rm: %d local: %d\n",
                   pstNode->eid, iTpR, pstNode->iSelfR, pstNode->vldPek.size(), iOldSize,
                   pstNode->vLfE.size(), pstNode->bRm, pstNode->bLc);
              ASSERT(0);
          }*/
          continue;
        }

        for (auto eid : vWait)
        {
            ++iRmCnt;
            //printf("PEEL start %d find triangles\n", eid);
            TPST_E *pstNode = obG.findNode(eid);
            ASSERT(!pstNode->bRm);
            pstNode->bRm = true;

            //printf("start rm triangles\n");
            for (int i = 0; i < pstNode->vLfE.size(); ++i)
            {
                const int e1 = pstNode->vLfE[i];
                const int e2 = pstNode->vRtE[i];
                /*printf("PEEL current R: %d neighbors %d %d rm: %d %d R: %d %d self: %d %d local: %d %d LB: %d UB: %d K: %d by %d\n",
                       iCurR, e1, e2, obG.findNode(e1)->bRm, obG.findNode(e2)->bRm,
                       obG.findNode(e1)->vSky[iCurK], obG.findNode(e2)->vSky[iCurK],
                       obG.findNode(e1)->iSelfR, obG.findNode(e2)->iSelfR,
                       obG.findNode(e1)->bLc, obG.findNode(e2)->bLc, iLB, iUB, iCurK, eid);*/
                if (obG.findNode(e1)->bRm || obG.findNode(e2)->bRm) continue;
                int iMinR = min(obG.findNode(e1)->vSky[iCurK], obG.findNode(e2)->vSky[iCurK]);
                if ((iMinR <= iLB) && (0 < iDesR)) continue;
                /*printf("PEEL minimum R: %d neighbors %d %d rm: %d %d R: %d %d local: %d %d\n",
                       iMinR, e1, e2, obG.findNode(e1)->bRm, obG.findNode(e2)->bRm,
                       obG.findNode(e1)->iSelfR, obG.findNode(e2)->iSelfR,
                       obG.findNode(e1)->bLc, obG.findNode(e2)->bLc);*/

                if ((obG.findNode(e1)->bLc) && (obG.findNode(e1)->iSelfR > iCurR))
                {
                    /*printf("PEEL decrease before %d R: %d pek: %d real: %d by %d\n",
                           e1, obG.findNode(e1)->iSelfR, obG.findNode(e1)->vldPek.size(),
                           obG.findNode(e1)->mpDebugPek.size() / 2, eid);*/
                    //obG.findNode(e1)->mpDebugPek.erase(eid);
                    //obG.findNode(e1)->mpDebugPek.erase(e2);
                    bool bRet = calPek_rm(obG, e1, eid, e2, obG.findNode(e1)->vldPek, iCurK);
                    if (!bRet)
                    {
                        deTruss_rePek(obG, e1, iCurK, true);
                    }
                    int iTpR = getR(obG.findNode(e1)->vldPek, obG.findNode(e1)->p, obG.m_dMinP, obG.m_dAccu, iCurK);
                    /*printf("PEEL decrease %d R: %d old R: %d pek: %d by %d\n",
                           e1, iTpR, obG.findNode(e1)->iSelfR, obG.findNode(e1)->vldPek.size(), eid);*/
                    if (iTpR < obG.findNode(e1)->iSelfR)
                    {
                        obG.findNode(e1)->iSelfR = iTpR;
                        prQ.push({e1, iTpR});
                    }
                }

                if ((obG.findNode(e2)->bLc) && (obG.findNode(e2)->iSelfR > iCurR))
                {
                    /*printf("PEEL decrease before %d R: %d pek: %d real: %d by %d\n",
                           e2, obG.findNode(e2)->iSelfR, obG.findNode(e2)->vldPek.size(),
                           obG.findNode(e2)->mpDebugPek.size() / 2, eid);*/
                    //obG.findNode(e2)->mpDebugPek.erase(eid);
                    //obG.findNode(e2)->mpDebugPek.erase(e1);
                    bool bRet = calPek_rm(obG, e2, eid, e1, obG.findNode(e2)->vldPek, iCurK);
                    if (!bRet)
                    {
                        deTruss_rePek(obG, e2, iCurK, true);
                    }
                    int iTpR = getR(obG.findNode(e2)->vldPek, obG.findNode(e2)->p, obG.m_dMinP, obG.m_dAccu, iCurK);
                    /*printf("PEEL decrease %d R: %d old R: %d pek: %d by %d\n",
                           e2, iTpR, obG.findNode(e2)->iSelfR, obG.findNode(e2)->vldPek.size(), eid);*/
                    if (iTpR < obG.findNode(e2)->iSelfR)
                    {
                        obG.findNode(e2)->iSelfR = iTpR;
                        prQ.push({e2, iTpR});
                    }
                }
            }
            /* store */
            changeR(obG, iCurK, eid, iDesR, pstNode->iPT);
        }

    }

    /* last */
    if ((iCurR >= iUpR) && (iRmCnt < vERaw.size()))
    {
        for (int iEid : vERaw)
        {
            if (!(obG.findNode(iEid)->bRm))
            {
                TPST_E *pstNode = obG.findNode(iEid);
                pstNode->iPT = iUpR;
                vURemain.push_back(iEid);
                /* store */
                changeR(obG, iCurK, iEid, iDesR, pstNode->iPT);
            }
        }
    }

    /* restore */
    for (int iEid : vERaw)
    {
        obG.findNode(iEid)->bRm = false;
        obG.findNode(iEid)->bLc = false;
    }
    //printf("DE_TRUSS first done\n");
    //ASSERT(0);
    return iCurR;
}
/*****************
input:
        myG &mapInitG
        LIST_DECOMP_G &lstDeG
        long double ldMiniP
output:
        LIST_DECOMP_G &lstDeG
description:
        truss decomposition
******************/
int insMaintain(myG &obG, int iCurK, vector<int> &vSeed)
{
    /* eid, R */
    auto cmp = [](pair<int, int> &left, pair<int, int> &right) { return left.second > right.second; };
    priority_queue<pair<int, int>, vector<pair<int, int> >, decltype(cmp)> prQ(cmp);

    /* add Q */
    for (int iEid : vSeed)
    {
        TPST_E *pstNode = obG.findNode(iEid);
        if (!pstNode->bInQ)
        {
            prQ.push({iEid, pstNode->vSky[iCurK]});
            pstNode->bInQ = true;
        }
    }

    int iCurR = 0;
    int iLastR = 0;
    vector<int> vWait;
    while (!prQ.empty())
    {
        iCurR = prQ.top().second;
        vWait.clear();
        while (!prQ.empty())
        {
            int iEid = prQ.top().first;
            int iNextR = prQ.top().second;
            if (iNextR > iCurR)
            {
                /* next turn */
                break;
            }
            prQ.pop();
            obG.findNode(iEid)->bInQ = false;
            vWait.push_back(iEid);
        }

        int iLB = 0;
        int iUB = 0;
        //printf("MAINTAIN start expend R: %d seed: %d\n", iCurR, vWait.size());
        iLastR = expend(obG, iCurK, iCurR, vWait, &iLB, &iUB);
        //printf("MAINTAIN expend R: %d size: %d LB: %d UB: %d\n", iCurR, vWait.size(), iLB, iUB);
        vector<int> vLRemain;
        vector<int> vURemain;
        peel(obG, iCurK, vWait, iLB, iUB, vLRemain, vURemain);
        //printf("MAINTAIN get remain L: %d U: %d\n", vLRemain.size(), vURemain.size());
        if (!vLRemain.empty())
        {
            for (int iEid : vLRemain)
            {
                TPST_E *pstNode = obG.findNode(iEid);
                init(obG, pstNode->eid, iCurK, iLB);
                pstNode->iSelfR = getR(pstNode->vldPek, pstNode->p, obG.m_dMinP, obG.m_dAccu, iCurK);
                printf("MAINTAIN recalculate %d R: %d base R: %d\n", iEid, pstNode->iSelfR, iLB);

                init(obG, pstNode->eid, iCurK, iCurR);
                pstNode->iSelfR = getR(pstNode->vldPek, pstNode->p, obG.m_dMinP, obG.m_dAccu, iCurK);
                printf("MAINTAIN recalculate %d R: %d base R: %d\n", iEid, pstNode->iSelfR, iCurR);

                for (int iPos = 0; iPos < pstNode->vLfE.size(); ++iPos)
                {
                    TPST_E *pstLf = obG.findNode(pstNode->vLfE[iPos]);
                    TPST_E *pstRt = obG.findNode(pstNode->vRtE[iPos]);
                    printf("PEEL neighbors: %d (%d, %d) %d (%d, %d) R: %d %d\n",
                           pstLf->eid, pstLf->paXY.first, pstLf->paXY.second,
                           pstRt->eid, pstRt->paXY.first, pstRt->paXY.second,
                           pstLf->vSky[iCurK], pstRt->vSky[iCurK]);
                }
                printf("MAINTAIN pek: %d neighbors: %d\n",
                       pstNode->vldPek.size(), pstNode->vLfE.size());
                for (auto atP : pstNode->vldPek)
                {
                    printf("MAINTAIN pek: %Lf\n", atP);
                }
            }
            ASSERT(0);
        }
        if (!vURemain.empty())
        {
            for (int iEid : vURemain)
            {
                TPST_E *pstNode = obG.findNode(iEid);
                if (pstNode->bInQ)
                {
                    continue;
                }
                prQ.push({iEid, pstNode->vSky[iCurK]});
                pstNode->bInQ = true;
            }
        }
    }
    return max(iCurR, iLastR);
}
/*****************
input:
        myG &mapInitG
        LIST_DECOMP_G &lstDeG
        long double ldMiniP
output:
        LIST_DECOMP_G &lstDeG
description:
        truss decomposition
******************/
int decMaintain(myG &obG, int iCurK, vector<int> &vSeed)
{
    /* eid, R, decreasing order */
    auto cmp = [](pair<int, int> &left, pair<int, int> &right) { return left.second < right.second; };
    priority_queue<pair<int, int>, vector<pair<int, int> >, decltype(cmp)> prQ(cmp);

    /* add Q */
    for (int iEid : vSeed)
    {
        TPST_E *pstNode = obG.findNode(iEid);
        if (!pstNode->bInQ)
        {
            prQ.push({iEid, pstNode->vSky[iCurK]});
            pstNode->bInQ = true;
        }
    }

    int iCurR = 0;
    vector<int> vWait;
    while (!prQ.empty())
    {
        iCurR = prQ.top().second;
        vWait.clear();
        while (!prQ.empty())
        {
            int iEid = prQ.top().first;
            int iNextR = prQ.top().second;
            if (iNextR != iCurR)
            {
                /* next turn */
                break;
            }
            prQ.pop();
            obG.findNode(iEid)->bInQ = false;
            vWait.push_back(iEid);
        }

        int iLB = 0;
        int iUB = 0;
        //printf("MAINTAIN start expend R: %d\n", iCurR);
        expend(obG, iCurK, iCurR, vWait, &iLB, &iUB);
        //printf("MAINTAIN expend R: %d size: %d LB: %d UB: %d\n", iCurR, vWait.size(), iLB, iUB);
        vector<int> vLRemain;
        vector<int> vURemain;
        int iLastR = peel(obG, iCurK, vWait, iLB, iUB, vLRemain, vURemain);
        //printf("MAINTAIN get remain L: %d U: %d\n", vLRemain.size(), vURemain.size());
        if (!vURemain.empty())
        {
            /*printf("MAINTAIN expend R: %d size: %d LB: %d UB: %d last R: %d\n",
                   iCurR, vWait.size(), iLB, iUB, iLastR);*/
            /*for (int iEid : vURemain)
            {
                TPST_E *pstNode = obG.findNode(iEid);
                init(obG, pstNode->eid, iCurK, iLB);
                pstNode->iSelfR = getR(pstNode->vldPek, pstNode->p, obG.m_dMinP, obG.m_dAccu, iCurK);
                printf("MAINTAIN recalculate %d R: %d base R: %d\n", iEid, pstNode->iSelfR, iLB);

                init(obG, pstNode->eid, iCurK, iCurR);
                pstNode->iSelfR = getR(pstNode->vldPek, pstNode->p, obG.m_dMinP, obG.m_dAccu, iCurK);
                printf("MAINTAIN recalculate %d R: %d base R: %d\n", iEid, pstNode->iSelfR, iCurR);

                for (int iPos = 0; iPos < pstNode->vLfE.size(); ++iPos)
                {
                    TPST_E *pstLf = obG.findNode(pstNode->vLfE[iPos]);
                    TPST_E *pstRt = obG.findNode(pstNode->vRtE[iPos]);
                    printf("PEEL neighbors: %d (%d, %d) %d (%d, %d) R: %d %d\n",
                           pstLf->eid, pstLf->paXY.first, pstLf->paXY.second,
                           pstRt->eid, pstRt->paXY.first, pstRt->paXY.second,
                           pstLf->vSky[iCurK], pstRt->vSky[iCurK]);
                }
                printf("MAINTAIN pek: %d neighbors: %d\n",
                       pstNode->vldPek.size(), pstNode->vLfE.size());
                for (auto atP : pstNode->vldPek)
                {
                    printf("MAINTAIN pek: %Lf\n", atP);
                }
            }*/
            //ASSERT(0);
        }
        if (!vLRemain.empty())
        {
            for (int iEid : vLRemain)
            {
                TPST_E *pstNode = obG.findNode(iEid);
                if (pstNode->bInQ)
                {
                    continue;
                }
                prQ.push({iEid, pstNode->vSky[iCurK]});
                pstNode->bInQ = true;
            }
        }
    }
    return iCurR;
}

/*****************
input:
        myG &obG
        vector<vector<int> > &vKG   // input
        vector<vector<pair<int, double> > > &vPos   // output
output:
        LIST_DECOMP_G &lstDeG
description:
        truss decomposition
******************/
int save(myG &obG, int iCurK, vector<int> &vKG, vector<pair<int, double> > &vPos)
{
    /* debug */
    vector<pair<int, int> > vDebugPos;

    vPos.clear();
    printf("=============SAVE k: %d k-truss: %d\n", iCurK, vKG.size());
    sort(vKG.begin(), vKG.end(),
        [&iCurK, &obG](const int e1, const int e2) {
         int iR1 = obG.findNode(e1)->vSky[iCurK];
         int iR2 = obG.findNode(e2)->vSky[iCurK];
         if (iR1 == iR2) return e1 < e2;
         else return iR1 < iR2;
        });
    //printf("SAVE sort done\n");
    int iCnt = 0;
    int iCurR = 0;
    int iTotalCnt = 0;
    for (int iEid : vKG)
    {
        ASSERT(iCurK < obG.findNode(iEid)->vSky.size());
        int iGetR = obG.findNode(iEid)->vSky[iCurK];
        ASSERT(iGetR >= iCurR);
        if (iGetR > iCurR)
        {
            /* new R */
            if (0 < iCnt)
            {
                vPos.push_back({iCnt, DE_TRUSS_R2P(iCurR, obG.m_dAccu)});
                vDebugPos.push_back({iCnt, iCurR});
                iTotalCnt += iCnt;
                iCnt = 0;
            }
            iCurR = iGetR;
        }
        /*printf("SAVE %d (%d, %d) self R: %d\n", iEid,
               obG.findNode(iEid)->paXY.first, obG.findNode(iEid)->paXY.second,
               iGetR);*/
        ++iCnt;
    }
    /* last */
    if (0 < iCnt)
    {
        vPos.push_back({iCnt, DE_TRUSS_R2P(iCurR, obG.m_dAccu)});
        vDebugPos.push_back({iCnt, iCurR});
        iTotalCnt += iCnt;
        iCnt = 0;
    }

    ASSERT(iTotalCnt == vKG.size());

    /* debug */
    //printf("SAVE DEBUG\n");
    verify(obG, vKG, vDebugPos, iCurK, obG.m_dMinP, obG.m_dAccu);

    return 0;
}
/*****************
input:
        myG &mapInitG
        LIST_DECOMP_G &lstDeG
        long double ldMiniP
output:
        LIST_DECOMP_G &lstDeG
description:
        truss decomposition
******************/
int deTruss_update(myG &obG, vector<TPST_INDEX_EDGE> &vQuery, char *pcIndexPath, bool bStore, char *pcSavePath)
{
	struct timeval tv;
	long lStartTime = 0;
	long lCurTime = 0;

    /* eid, p */
	vector<pair<int, double> > vNewIncE;
	vector<pair<int, double> > vNewDecE;
	vector <int> vKTruss;
	vector<pair<int, double> > vPos;

	g_bLock = false;

	for(auto atNode : vQuery)
    {
        TPST_E *pstNode = obG.findNode(atNode.x, atNode.y);
        ASSERT(NULL != pstNode);
        //printf("DE_TRUSS get (%d, %d) %lf old: %lf\n", atNode.x, atNode.y, atNode.p, pstNode->p);
        if (atNode.p > pstNode->p)
        {
            vNewIncE.push_back({pstNode->eid, atNode.p});
        }
        else if (atNode.p < pstNode->p)
        {
            vNewDecE.push_back({pstNode->eid, atNode.p});
        }
    }

    if (!vNewIncE.empty())
    {
        int iMaxK = 0;
        for (auto atN : vNewIncE)
        {
            TPST_E *pstNode = obG.findNode(atN.first);
            int iLcMaxK = pstNode->iTrussness;
            iMaxK = max(iMaxK, iLcMaxK);
            /* printf("DE_TRUSS (%d, %d) old p: %lf new p: %lf max k: %d\n",
                   pstNode->paXY.first, pstNode->paXY.second, pstNode->p, atN.second, iLcMaxK); */
            pstNode->p = atN.second;
        }
        printf("DE_TRUSS max increase affected k: %d\n", iMaxK);
        for (int iCurK = 3; iCurK <= iMaxK; ++iCurK)
        {
            printf("=============INCREMENTAL UPDATE k: %d\n", iCurK);

            /* read index */
            gettimeofday(&tv, NULL);
            lStartTime = tv.tv_sec * 1000 + tv.tv_usec / 1000;

            vKTruss.clear();
            vPos.clear();
            file_readIndex(obG, iCurK, vKTruss, vPos, pcIndexPath);

            gettimeofday(&tv, NULL);
            lCurTime = tv.tv_sec * 1000 + tv.tv_usec / 1000;
            g_lReadTime += lCurTime - lStartTime;

            vector<int> vSeed;
            //printf("DE_TRUSS increase size: %d decrease size: %d\n", vNewIncE.size(), vNewDecE.size());

            /* increase */
            insOneHop(obG, iCurK, vNewIncE, vSeed);
            int iInsFlag = insMaintain(obG, iCurK, vSeed);
            /* decrease */
            //vSeed.clear();
            //decOneHop(obG, iCurK, vNewDecE, vSeed);
            //int iDEcFlag = decMaintain(obG, iCurK, vSeed);

            /* restore */
            for (int iEid : obG.m_vInitE)
            {
                TPST_E *pstNode = obG.findNode(iEid);
                pstNode->bInit = false;
                //printf("DE_TRUSS unset k: %d eid: %d pek: %d\n", iCurK, pstNode->eid, pstNode->vldPek.size());
            }
            obG.m_vInitE.clear();

            /* store */
            if (bStore)
            {
                gettimeofday(&tv, NULL);
                lStartTime = tv.tv_sec * 1000 + tv.tv_usec / 1000;

                vPos.clear();
                save(obG, iCurK, vKTruss, vPos);
                file_saveIndex(obG, iCurK, vKTruss, vPos, pcSavePath);

                gettimeofday(&tv, NULL);
                lCurTime = tv.tv_sec * 1000 + tv.tv_usec / 1000;
                g_lStoreTime += lCurTime - lStartTime;
            }

            printf("DE_TRUSS incremental last R: %d\n", iInsFlag);
            if (1 > iInsFlag)
            {
                /* empty, stop early */
                break;
            }
        }
    }
    if (!vNewDecE.empty())
    {
        int iMaxK = 0;
        for (auto atN : vNewDecE)
        {
            TPST_E *pstNode = obG.findNode(atN.first);
            int iLcMaxK = pstNode->iTrussness;
            iMaxK = max(iMaxK, iLcMaxK);
            //printf("DE_TRUSS old p: %lf max k: %d\n", pstNode->p, iLcMaxK);
            pstNode->p = atN.second;
        }
        printf("DE_TRUSS max increase affected k: %d\n", iMaxK);
        for (int iCurK = 3; iCurK <= iMaxK; ++iCurK)
        {
            printf("=============DECREMENTAL UPDATE k: %d\n", iCurK);

            /* read index */
            gettimeofday(&tv, NULL);
            lStartTime = tv.tv_sec * 1000 + tv.tv_usec / 1000;

            vKTruss.clear();
            vPos.clear();
            file_readIndex(obG, iCurK, vKTruss, vPos, pcIndexPath);

            gettimeofday(&tv, NULL);
            lCurTime = tv.tv_sec * 1000 + tv.tv_usec / 1000;
            g_lReadTime += lCurTime - lStartTime;

            vector<int> vSeed;
            //printf("DE_TRUSS increase size: %d decrease size: %d\n", vNewIncE.size(), vNewDecE.size());

            /* increase */
            //insOneHop(obG, iCurK, vNewIncE, vSeed);
            //int iInsFlag = insMaintain(obG, iCurK, vSeed);
            /* decrease */
            vSeed.clear();
            decOneHop(obG, iCurK, vNewDecE, vSeed);
            int iDEcFlag = decMaintain(obG, iCurK, vSeed);

            /* restore */
            for (int iEid : obG.m_vInitE)
            {
                TPST_E *pstNode = obG.findNode(iEid);
                pstNode->bInit = false;
                //printf("DE_TRUSS unset k: %d eid: %d pek: %d\n", iCurK, pstNode->eid, pstNode->vldPek.size());
            }
            obG.m_vInitE.clear();

            /* store */
            if (bStore)
            {
                gettimeofday(&tv, NULL);
                lStartTime = tv.tv_sec * 1000 + tv.tv_usec / 1000;

                vPos.clear();
                save(obG, iCurK, vKTruss, vPos);
                file_saveIndex(obG, iCurK, vKTruss, vPos, pcSavePath);

                gettimeofday(&tv, NULL);
                lCurTime = tv.tv_sec * 1000 + tv.tv_usec / 1000;
                g_lStoreTime += lCurTime - lStartTime;
            }

            printf("DE_TRUSS decremental last R: %d\n", iDEcFlag);
            if (1 > iDEcFlag)
            {
                /* empty, stop early */
                break;
            }
        }
    }


    return 0;
}

/*****************
input:
        myG &obG
        vector<vector<int> > &vKG   // input
        vector<vector<pair<int, double> > > &vPos   // output
output:
        LIST_DECOMP_G &lstDeG
description:
        truss decomposition
******************/
int deTruss_save(myG &obG, vector<vector<int> > &vKG, vector<vector<pair<int, double> > > &vPos)
{
    vPos.clear();
    vPos.resize(vKG.size());
    for (int iCurK = vKG.size() - 1; iCurK > 2; --iCurK)
    {
        printf("=============SAVE k: %d k-truss: %d\n", iCurK, vKG[iCurK].size());
        sort(vKG[iCurK].begin(), vKG[iCurK].end(),
            [&iCurK, &obG](const int e1, const int e2) {
             int iR1 = obG.findNode(e1)->vSky[iCurK];
             int iR2 = obG.findNode(e2)->vSky[iCurK];
             if (iR1 == iR2) return e1 < e2;
             else return iR1 < iR2;
            });
        save(obG, iCurK, vKG[iCurK], vPos[iCurK]);
    }

    return 0;
}
