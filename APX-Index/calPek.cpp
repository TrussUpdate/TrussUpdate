/***************
calculate pek
****************/
#include <sys/time.h>
#include <math.h>
#include "common.h"
#include "myG.h"
#include "deTruss.h"
#include "calPek.h"

long g_lInitCnt;
long g_lUptCnt;
vector<int> CAL_g_vLfE;
vector<int> CAL_g_vRtE;
vector<double> CAL_g_vLfP;
vector<double> CAL_g_vRtP;

/*****************
input:
        vector<double> &vLfP
        vector<double> &vRtP
output:
        vector <long double> &vfDtE
description:
        calculate pek
******************/
int calOneEdgePekByNeib(vector<double> &vLfP, vector<double> &vRtP, vector <long double> &vfDtE)
{
    ++g_lInitCnt;

    int iKe = vLfP.size();
	long double *temp;
	long double *pw;
	long double *B;

	temp = new long double [iKe + 1];
	pw = new long double [iKe + 1];
	B = new long double [iKe + 1];

	pw[0] = -1;
	for (int i = 0; i < iKe; ++i)
	{
		pw[i+1] = vLfP[i] * vRtP[i];
	}

	temp[0] = 1;
	for (int i = 1; i <= iKe; i++)
	{
		temp[i] = temp[i - 1] * (1 - pw[i]);
	}

	vfDtE[0] = temp[iKe];
	for (int l=1; l <= iKe; l++)
	{
		B[l-1] = 0;
		for (int i=l; i <= iKe; i++)
		{
			B[i] = temp[i-1]*pw[i] +B[i-1]*(1-pw[i]);
		}
		for (int i=l-1; i <= iKe; i++)
		{
			temp[i] = B[i];
		}
		vfDtE[l] = temp[iKe];
	}

	delete [] temp;
	delete [] pw;
	delete [] B;
    return 0;
}
/*****************
input:
        myG &mapG
        vector<int> &vLfE
        vector<int> &vRtE
output:
        vector <long double> vfDtE
description:
        show MAP_G mapG
******************/
int cal_OneEdgePek(myG &mapG, vector<int> &vLfE, vector<int> &vRtE, vector <long double> &vfDtE)
{
    int iSup = vLfE.size();
    vfDtE.clear();
    vfDtE.resize(iSup + 1);
    CAL_g_vLfP.clear();
    CAL_g_vLfP.resize(iSup);
    CAL_g_vRtP.clear();
    CAL_g_vRtP.resize(iSup);
    for (int i = 0; i < iSup; ++i)
    {
         CAL_g_vLfP[i] = mapG.findNode(vLfE[i])->p;
         CAL_g_vRtP[i] = mapG.findNode(vRtE[i])->p;
    }
    calOneEdgePekByNeib(CAL_g_vLfP, CAL_g_vRtP, vfDtE);
    ASSERT(cal_check(vfDtE));

    return iSup;
}

/*****************
input:
        vector <long double> &vfDtE
output:
        qualified?
description:
        check whether it needs to be recalculated
******************/
bool cal_check(vector <long double> &vfDtE)
{
	long double fSum = 0;
	int iSize = vfDtE.size();
	for (int i=0; i < iSize; ++i)
	{
		fSum+= vfDtE[i];
	}

	if (fabs(fSum-1)>=1e-6)
	{
		return false;
	}

    return true;
}
/*****************
input:
        int x
        int y
        int z
        int ke
        vector <long double> &vldPekXY
        myG &mapG
description:
        update one edge
        x, y is update edge
        x, y, z is vanished triangle
        ke is support number after vanished
******************/
int calPekRm(double ldLf, double ldRt, vector <long double> &vldPekXY)
{
    int iNewSup = vldPekXY.size() - 2;

	long double p_del = 0;

	p_del = ldLf * ldRt;

	if (p_del<0.5)
	{
		vldPekXY[0] = vldPekXY[0]/(1-p_del);

		for (int i=1; i <= iNewSup; ++i)
		{
            vldPekXY[i] = (vldPekXY[i]/(1.0-p_del)-p_del*vldPekXY[i-1]/(1.0-p_del));
		}
	}
	else
    {
        /* new must be 0 */
        /* pet[ke + 1] be 0 */
        long double *pet = new long double[iNewSup+2];
		pet[iNewSup + 1] = 0;
		for (int i = iNewSup + 1; i>=1; --i)
		{
			pet[i-1] = (vldPekXY[i]-pet[i]*(1.0-p_del))/p_del;
		}
		for (int i=0; i<=iNewSup; ++i)
		{
			vldPekXY[i] = pet[i];
		}
        delete [] pet;
	}
    /* clear old max value */
    vldPekXY[iNewSup + 1] = 0;
    vldPekXY.resize(iNewSup + 1);

	return 0;
}
/*****************
input:
        int x
        int y
        int z
        int ke
        vector <long double> &vldPekXY
        myG &mapG
description:
        update one edge
        x, y is update edge
        x, y, z is vanished triangle
        ke is support number after vanished
******************/
int calPekAdd(double ldLf, double ldRt, vector <long double> &vldPekXY)
{
    int iNewSup = vldPekXY.size();
    vldPekXY.resize(iNewSup + 1);

	long double p_add = ldLf * ldRt;

    vldPekXY[iNewSup] = vldPekXY[iNewSup - 1]*p_add;

    for (int i=iNewSup - 1; i > 0; --i)
    {
        vldPekXY[i] = (vldPekXY[i]*(1.0-p_add)+p_add*vldPekXY[i-1]);
    }
    vldPekXY[0] = vldPekXY[0]*(1-p_add);

	return 0;
}
/*****************
input:
        int x
        int y
        int z
        int ke
        vector <long double> &vldPekXY
        myG &mapG
description:
        update one edge
        x, y is update edge
        x, y, z is vanished triangle
        ke is support number after vanished
******************/
int calPek_rm(myG &mapG, int iEid, int iLfEid, int iRtEid, vector <long double> &vldPekXY)
{
    ++g_lUptCnt;
    /* debug */
    //vector <long double> vDebug = vldPekXY;
    calPekRm(mapG.findNode(iLfEid)->p, mapG.findNode(iRtEid)->p, vldPekXY);

    if (cal_check(vldPekXY))
    {
        return 0;
    }
    /* debug */
    /*for (auto atR : vldPekXY)
    {
        printf("%Lf\n", atR);
    }
    printf("ERROR before: \n");
    for (auto atR : vDebug)
    {
        printf("%Lf\n", atR);
    }
    printf("ERROR bRm: %d bUp: %d bNew: %d\n",
           mapG.findNode(iEid)->bRm,
           mapG.findNode(iEid)->bUp,
           mapG.findNode(iEid)->bNew);
    mapG.m_vLfE.clear();
    mapG.m_vRtE.clear();
    mapG.findNotRmNeb(mapG.findNode(iEid)->paXY.first, mapG.findNode(iEid)->paXY.second,
                 mapG.m_vLfE, mapG.m_vRtE);

    mapG.m_vLfP.clear();
    mapG.m_vRtP.clear();
    for (int i = 0; i < mapG.m_vLfE.size(); ++i)
    {
         mapG.m_vLfP.push_back(mapG.findNode(mapG.m_vLfE[i])->p);
         mapG.m_vRtP.push_back(mapG.findNode(mapG.m_vRtE[i])->p);
        if((mapG.m_vLfE[i] == iLfEid) || (mapG.m_vLfE[i] == iRtEid))
        {
            printf("rm %lf, %lf\n", mapG.findNode(mapG.m_vLfE[i])->p,
                   mapG.findNode(mapG.m_vRtE[i])->p);
            continue;
        }
        printf("%lf, %lf\n", mapG.findNode(mapG.m_vLfE[i])->p,
               mapG.findNode(mapG.m_vRtE[i])->p);

    }
    vldPekXY.clear();
    vldPekXY.resize(mapG.m_vLfP.size() + 1);
    calOneEdgePekByNeib(mapG.m_vLfP, mapG.m_vRtP, vldPekXY);
    for (auto atR : vldPekXY)
    {
        printf("%Lf\n", atR);
    }

    calPekRm(mapG.findNode(iLfEid)->p, mapG.findNode(iRtEid)->p, vldPekXY);
    printf("ERROR rm again: \n");
    for (auto atR : vldPekXY)
    {
        printf("%Lf\n", atR);
    }*/

    /* not qualified, recalculate */
    CAL_g_vLfE.clear();
    CAL_g_vRtE.clear();
    mapG.findNotRmNeb(mapG.findNode(iEid)->paXY.first, mapG.findNode(iEid)->paXY.second,
                 CAL_g_vLfE, CAL_g_vRtE);

    int iSup = CAL_g_vLfE.size();
    CAL_g_vLfP.clear();
    CAL_g_vRtP.clear();
    for (int i = 0; i < iSup; ++i)
    {
        if((CAL_g_vLfE[i] == iLfEid) || (CAL_g_vLfE[i] == iRtEid))
        {
            continue;
        }
         CAL_g_vLfP.push_back(mapG.findNode(CAL_g_vLfE[i])->p);
         CAL_g_vRtP.push_back(mapG.findNode(CAL_g_vRtE[i])->p);

    }
    vldPekXY.clear();
    vldPekXY.resize(CAL_g_vLfP.size() + 1);
    calOneEdgePekByNeib(CAL_g_vLfP, CAL_g_vRtP, vldPekXY);

    /* debug */
    /*printf("ERROR recalculate: \n");
    for (auto atR : vldPekXY)
    {
        printf("%Lf\n", atR);
    }
    ASSERT(0);*/

    ASSERT(cal_check(vldPekXY));

	return 0;
}
/*****************
input:
        int x
        int y
        int z
        int ke
        vector <long double> &vldPekXY
        myG &mapG
description:
        update one edge
        x, y is update edge
        x, y, z is vanished triangle
        ke is support number after vanished
******************/
int calPek_add(myG &mapG, int iEid, int iLfEid, int iRtEid, vector <long double> &vldPekXY)
{
    ++g_lUptCnt;
    /* debug */
    //vector <long double> vDebug = vldPekXY;
    calPekAdd(mapG.findNode(iLfEid)->p, mapG.findNode(iRtEid)->p, vldPekXY);

    if (cal_check(vldPekXY))
    {
        return 0;
    }

    /* debug */
    /*for (auto atR : vldPekXY)
    {
        printf("%Lf\n", atR);
    }
    printf("ERROR before: \n");
    for (auto atR : vDebug)
    {
        printf("%Lf\n", atR);
    }
    printf("ERROR bRm: %d bUp: %d bNew: %d\n",
           mapG.findNode(iEid)->bRm,
           mapG.findNode(iEid)->bUp,
           mapG.findNode(iEid)->bNew);
    mapG.m_vLfE.clear();
    mapG.m_vRtE.clear();
    mapG.findNotRmNeb(mapG.findNode(iEid)->paXY.first, mapG.findNode(iEid)->paXY.second,
                 mapG.m_vLfE, mapG.m_vRtE);

    mapG.m_vLfP.clear();
    mapG.m_vRtP.clear();
    for (int i = 0; i < mapG.m_vLfE.size(); ++i)
    {
         mapG.m_vLfP.push_back(mapG.findNode(mapG.m_vLfE[i])->p);
         mapG.m_vRtP.push_back(mapG.findNode(mapG.m_vRtE[i])->p);
        if((mapG.m_vLfE[i] == iLfEid) || (mapG.m_vLfE[i] == iRtEid))
        {
            printf("add %lf, %lf\n", mapG.findNode(mapG.m_vLfE[i])->p,
                   mapG.findNode(mapG.m_vRtE[i])->p);
            continue;
        }
        printf("%lf, %lf\n", mapG.findNode(mapG.m_vLfE[i])->p,
               mapG.findNode(mapG.m_vRtE[i])->p);

    }
    vldPekXY.clear();
    vldPekXY.resize(mapG.m_vLfP.size() + 1);
    calOneEdgePekByNeib(mapG.m_vLfP, mapG.m_vRtP, vldPekXY);
    for (auto atR : vldPekXY)
    {
        printf("%Lf\n", atR);
    }*/
    ASSERT(0);

    /* not qualified, recalculate */

    CAL_g_vLfE.clear();
    CAL_g_vRtE.clear();
    mapG.findNotRmNeb(mapG.findNode(iEid)->paXY.first, mapG.findNode(iEid)->paXY.second,
                 CAL_g_vLfE, CAL_g_vRtE);

    int iSup = CAL_g_vLfE.size();
    CAL_g_vLfP.clear();
    CAL_g_vRtP.clear();
    for (int i = 0; i < iSup; ++i)
    {
         CAL_g_vLfP.push_back(mapG.findNode(CAL_g_vLfE[i])->p);
         CAL_g_vRtP.push_back(mapG.findNode(CAL_g_vRtE[i])->p);

    }
    vldPekXY.clear();
    vldPekXY.resize(CAL_g_vLfP.size() + 1);
    calOneEdgePekByNeib(CAL_g_vLfP, CAL_g_vRtP, vldPekXY);
    ASSERT(cal_check(vldPekXY));

	return 1;
}
/*****************
input:
        MAP_BASIC_G &mapSupG
        myG &OmapG
description:
        calculate myG OmapG part pek, sup, bReadyDel
******************/
int calPek_part(myG &OmapG, vector<int> &vEdges)
{
    for (int eid : vEdges)
    {
        TPST_E *pstNode = OmapG.findNode(eid);
        CAL_g_vLfE.clear();
        CAL_g_vRtE.clear();
        OmapG.findNotRmNeb(pstNode->paXY.first, pstNode->paXY.second, CAL_g_vLfE, CAL_g_vRtE);

        pstNode->iSup = cal_OneEdgePek(OmapG, CAL_g_vLfE, CAL_g_vRtE, pstNode->vldPek);
    }

    return 0;
}
