#pragma once
/***************
my class G
****************/

typedef struct tpstIndexEdge
{
    int x;
    int y;
    float p;
}TPST_INDEX_EDGE;

typedef struct tpstE
{
    int eid;
    double p;
    int iTrussness;  // determined
    int iSup;   // uncertain support
    int iSelfR;  //R calculated by vldPek, varying by k
    int iColR;  //collapse R
    int iPT;  //temp variable
    bool bRm;
    //bool bUp;
    //bool bNew;
    bool bInit;
    bool bInQ;  // maintenance process queue
    bool bBFS;
    bool bLc;   // local peeling process
    pair <int, int> paXY;
    vector <long double> vldPek;    // init
    vector <long double> vldOldPek;    // only for decrement
    //map <int, int> mpDebugPek;   // <eid, eid> only for debug
    //map <int, int> mpDebugOldPek;   // only for debug
    vector <int> vLfE;    // init
    vector <int> vRtE;    // init
    /* k -> R */
    vector <int> vSky;
}TPST_E;

typedef struct tpstAdj
{
    int pid;
    int eid;
}TPST_ADJ;
#define myG_getPair(iX, iY) ((iX < iY) ? (pair<int, int>(iX, iY)) : (pair<int, int>(iY, iX)))

class myG
{
private:
    vector<int> m_vLfE;
    vector<int> m_vRtE;
    vector<double> m_vLfP;
    vector<double> m_vRtP;
    vector<int> m_vTpEdges;
public:
    vector <vector <TPST_ADJ> > m_Adj;
    vector <TPST_E> m_Edges;
    vector <int> m_vInitE;

    int m_iMaxPId;
    int m_iMaxEId;
    int m_iMaxK;
    int m_iMaxPK;
    double m_dMinP;
    double m_dAccu;

    myG();
    ~myG();
    int init();

    bool add(int x, int y, double ldP);
    int findNeb(int x, int y, vector<int> &vLfE, vector<int> &vRtE);
    int findNotRmNeb(int x, int y, vector<int> &vLfE, vector<int> &vRtE);

    TPST_E *findNode(int x, int y);
    TPST_E *findNode(int eid);
};
