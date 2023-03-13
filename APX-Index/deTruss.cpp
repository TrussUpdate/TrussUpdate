/***************
truss decomposition
****************/
#include <sys/time.h>
#include <queue>
#include <numeric>

#include "common.h"
#include "myG.h"
#include "deTruss.h"
#include "calPek.h"
#include "file.h"

#define MAX_R(ldAccu) ((int)(1 / ldAccu))

vector<int> g_vLfE;
vector<int> g_vRtE;
vector<double> g_vLfP;
vector<double> g_vRtP;

vector<pair<int, double> > g_vPos;

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

    if (sum < ldMiniP)
    {
        return 0;
    }


    // DEBUG
    /*for (auto atR : vfDtE)
    {
        printf("%Lf ", atR);
    }
    printf("\nGET_R p: %f real: %Lf number: %d\n", ldP, sum, (int)(sum / ldAccu));
    printf("GET_R ldMiniP: %f size: %d k: %d\n", ldMiniP, vfDtE.size(), iDesK);
    ASSERT(0);*/

    ASSERT(vfDtE.size() - 1 + 2 >= iDesK);

    return (int)(sum / ldAccu);
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
    priority_queue<pair<int, int>, vector<pair<int, int> >, decltype(cmp)> prQ(cmp);
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
    priority_queue<pair<int, int>, vector<pair<int, int> >, decltype(cmp)> prQ(cmp);
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
    priority_queue<pair<int, int>, vector<pair<int, int> >, decltype(cmp)> prQ(cmp);
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
int verify(myG &obG, vector<int> &vKG, vector<pair<int, int> > &vPos, int iCurk, double ldMiniP, double ldAccu)
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
        int iDiff = deTruss_verify(obG, vSubE, iEsR, iCurk, ldMiniP, ldAccu, paXY);
        itBegin += iOffset;

        printf("VERIFY k: %d calculated R: %d difference: %d\n", iCurk, iEsR, iDiff);
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
        myG &obG
        vector<vector<int> > &vKG
        double ldMiniP
        double ldAccu
output:
        LIST_DECOMP_G &lstDeG
description:
        truss decomposition
******************/
int deTruss_PT(myG &obG, const vector<vector<int> > &vKG, double ldMiniP, double ldAccu, bool bStore, char *pcSavePath)
{
	struct timeval tv;
	long lStartTime = 0;
	long lCurTime = 0;

    int iKTrussSize = 0;
    vector<int> vUpG;
    vector<vector<pair<int, int> > > vGPos;
    int iCurK = vKG.size() - 1;
    vGPos.resize(vKG.size());
    //vRetPos.resize(vKG.size());
    printf("=============PEEL top k: %d k-class: %d\n", iCurK, vKG[iCurK].size());
    iKTrussSize = vKG[iCurK].size();
    if (2 >= iCurK)
    {
        //PtoR(obG, vKG[iCurK], vRetPos[iCurK]);
        return 0;
    }
    /* top layer */
    ASSERT(!vKG[iCurK].empty());
    vector <int> vKPlusTruss = vKG[iCurK];
    vector <int> vKTruss;
    deTruss_topL(obG, iCurK, ldMiniP, ldAccu, vKPlusTruss, vGPos[iCurK]);
    --iCurK;
    /* next layers */
    for (; iCurK > 0; --iCurK)
    {
        /* save (k+1)-truss */
        if (bStore)
        {
            gettimeofday(&tv, NULL);
            lStartTime = tv.tv_sec * 1000 + tv.tv_usec / 1000;
            g_vPos.clear();
            g_vPos.reserve(vGPos[iCurK + 1].size());
            for (auto atPos : vGPos[iCurK + 1])
            {
                int iCnt = atPos.first;
                double ldP = recoverR(atPos.second, ldAccu);
                g_vPos.push_back({iCnt, ldP});
            }
            file_saveIndex(obG, iCurK + 1, vKPlusTruss, g_vPos, pcSavePath);
            gettimeofday(&tv, NULL);
            lCurTime = tv.tv_sec * 1000 + tv.tv_usec / 1000;
            g_lStoreTime += lCurTime - lStartTime;
        }

        if (2 >= iCurK)
        {
            //PtoR(obG, vKG[iCurK], vRetPos[iCurK]);
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
            vAffE.insert(vAffE.end(), vKPlusTruss.begin(), vKPlusTruss.begin() + iUpFirstCnt);
            //printf("PEEL vAffEs size: %d\n", vAffE.size());
            //printf("PEEL vKG: %d cnt: %d\n", vKG[iCurK+1].size(), iUpFirstCnt);
            vUpE.insert(vUpE.end(), vKPlusTruss.begin() + iUpFirstCnt, vKPlusTruss.end());
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
            vUpE = vKPlusTruss;
        }

        vector<int> vEdges = vKG[iCurK];
        //vEdges.swap(vKG[iCurK]);
        //vKG[iCurK].clear();
        vKTruss.clear();
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
            vKTruss.insert(vKTruss.end(), vEdges.begin(), vEdges.end() - iCurCnt);
            vEdges.erase(vEdges.begin(), vEdges.end() - iCurCnt);
            vAffE.clear();
            vAffE.insert(vAffE.end(), vKPlusTruss.begin() + iUpCnt, vKPlusTruss.begin() + iUpCnt + iAffCnt);
            //printf("PEEL affE: %d\n", vAffE.size());
            iUpCnt += iAffCnt;
            vUpE.clear();
            if (iUpCnt < vKPlusTruss.size())
            {
                vUpE.insert(vUpE.end(), vKPlusTruss.begin() + iUpCnt, vKPlusTruss.end());
            }
            //printf("PEEL used: %d upE: %d\n", iUpCnt, vUpE.size());
            deTruss_middleL(obG, iCurK, ldMiniP, ldAccu, vEdges, vAffE, vGPos[iCurK], vUpE, iUpR);
            //printf("PEEL second done edges: %d up edges: %d upR: %d\n", vEdges.size(), vUpE.size(), iUpR);
        }

        /* last */
        vKTruss.insert(vKTruss.end(), vEdges.begin(), vEdges.end());
        if (iKTrussSize != vKTruss.size())
        {
            printf("ERROR k-truss: %d get %d\n", iKTrussSize, vKTruss.size());
            ASSERT(0);
        }

        int iLastCnt = vGPos[iCurK].back().first;
        if (0 == iLastCnt)
        {
            vGPos[iCurK].pop_back();
        }

        vKTruss.swap(vKPlusTruss);
        //printf("PEEL k-truss: %d\n", vKG[iCurK].size());
    }

    //printf("PEEL save R size: %d\n", vGPos.size());
    /* save */
    /*vRetPos.resize(vGPos.size());
    for (iCurK = 3; iCurK < vGPos.size(); ++iCurK)
    {
        for (auto atPos : vGPos[iCurK])
        {
            int iCnt = atPos.first;
            double ldP = recoverR(atPos.second, ldAccu);
            vRetPos[iCurK].push_back({atPos.first, ldP});
        }
    }*/
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
