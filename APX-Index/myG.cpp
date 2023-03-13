/***************
my class G
****************/

#include "common.h"
#include "myG.h"

/*****************
input:
        none
description:
        init object
        init stage
******************/
myG::myG()
{
    TPST_E stTp = {0};

    m_iMaxPId = 0;
    m_iMaxEId = 0;

    m_Edges.push_back(stTp);
}

/*****************
input:
        myG &oInitG
description:
        copy map
        calculate stage
******************/
myG::~myG()
{
    m_Edges.clear();
    m_Adj.clear();
}

/*****************
input:
        none
description:
        init after add nodes
******************/
int myG::init()
{
    for (int vid = 1; vid <= m_iMaxPId; ++vid)
    {
        sort(m_Adj[vid].begin(), m_Adj[vid].end(),
            [](const TPST_ADJ& e1, const TPST_ADJ& e2) {
             if (e1.pid == e2.pid) return e1.eid < e2.eid;
             else return e1.pid < e2.pid;
            });
        /* debug */
        /*if (!m_Adj[vid].empty())
        {
            for (auto atNode : m_Adj[vid])
            {
                printf("%d(%d) ", atNode.pid, atNode.eid);
            }
            printf("\n");
        }*/
        /* remove useless edges */
        if (!m_Adj[vid].empty())
        {
            for (int iPos = 0; iPos < m_Adj[vid].size() - 1; ++iPos)
            {
                if (m_Adj[vid][iPos].pid == m_Adj[vid][iPos + 1].pid)
                {
                    int iRmEid = m_Adj[vid][iPos + 1].eid;
                    findNode(iRmEid)->eid = 0;
                    findNode(iRmEid)->bRm = true;
                }
            }
        }

        vector <TPST_ADJ>::iterator it = unique(m_Adj[vid].begin(), m_Adj[vid].end(),
            [](const TPST_ADJ& e1, const TPST_ADJ& e2) {
            return (e1.pid == e2.pid);
            });

        /* debug */
        /*if (!m_Adj[vid].empty())
        {
            for (auto atNode : m_Adj[vid])
            {
                printf("%d(%d) ", atNode.pid, atNode.eid);
            }
            printf("\n");
            ASSERT(0);
        }*/

        /*for (auto tpIt = it; tpIt != m_Adj[vid].end(); ++tpIt)
        {
            findNode(tpIt->eid)->eid = 0;
            findNode(tpIt->eid)->bRm = true;
            for (auto itRm = m_Adj[tpIt->pid].begin(); itRm != m_Adj[tpIt->pid].end(); ++itRm)
            {
                if (itRm->eid == tpIt->eid)
                {
                    m_Adj[tpIt->pid].erase(itRm);
                    break;
                }
            }
        }*/
        m_Adj[vid].resize(std::distance(m_Adj[vid].begin(),it));
    }

    return 0;
}
/*****************
input:
        int x
        int y
        long double ldP
description:
        add new edge
        init stage
******************/
bool myG::add(int x, int y, double ldP)
{
    ++m_iMaxEId;
    m_iMaxPId = max(m_iMaxPId, x);
    m_iMaxPId = max(m_iMaxPId, y);
    if (m_iMaxPId + 1 > m_Adj.size())
    {
        m_Adj.resize(m_iMaxPId + 1);
    }
    TPST_ADJ stTp = {0};
    stTp.pid = y;
    stTp.eid = m_iMaxEId;
    m_Adj[x].push_back(stTp);
    stTp.pid = x;
    m_Adj[y].push_back(stTp);


    TPST_E stTpE = {0};
    stTpE.eid = m_iMaxEId;
    stTpE.p = ldP;
    stTpE.bRm = false;
    stTpE.bUp = false;
    stTpE.bNew = false;
    stTpE.iColR = -1;
    stTpE.iPT = -1;
    stTpE.paXY = myG_getPair(x, y);
    m_Edges.push_back(stTpE);
    return true;
}

/*****************
input:
        int x
        int y
        list<int> &lstNeib
description:
        delete edge
        calculate stage
******************/
int myG::findNeb(int x, int y, vector<int> &vLfE, vector<int> &vRtE)
{
    int u = x;
    int v = y;
    if (m_Adj[x].size() > m_Adj[y].size())
    {
        u = y;
        v = x;
    }
    vector<TPST_ADJ>::iterator itNext = m_Adj[v].begin();
    for (int i = 0; i < m_Adj[u].size(); ++i)
    {
        int w = m_Adj[u][i].pid;

        auto prc_info = std::lower_bound(itNext, m_Adj[v].end(), w,
            [](const TPST_ADJ& info, int value)
            {
                return info.pid < value;
            });
        if (prc_info == m_Adj[v].end())
        {
            // end
            break;
        }
        else if (w == prc_info->pid)
        {
            // found
            vLfE.push_back(m_Adj[u][i].eid);
            vRtE.push_back(prc_info->eid);
            itNext = ++prc_info;
        }
        else
        {
            // not found
            itNext = prc_info;
        }
    }

    return vLfE.size();
}
/*****************
input:
        int x
        int y
        list<int> &lstNeib
description:
        delete edge
        calculate stage
******************/
int myG::findNotRmNeb(int x, int y, vector<int> &vLfE, vector<int> &vRtE)
{
    int u = x;
    int v = y;
    if (m_Adj[x].size() > m_Adj[y].size())
    {
        u = y;
        v = x;
    }
    vector<TPST_ADJ>::iterator itNext = m_Adj[v].begin();
    for (int i = 0; i < m_Adj[u].size(); ++i)
    {
        if (findNode(m_Adj[u][i].eid)->bRm)
        {
            continue;
        }
        int w = m_Adj[u][i].pid;

        auto prc_info = std::lower_bound(itNext, m_Adj[v].end(), w,
            [](const TPST_ADJ& info, int value)
            {
                return info.pid < value;
            });
        if (prc_info == m_Adj[v].end())
        {
            // end
            break;
        }
        else if (w == prc_info->pid)
        {
            if (findNode(prc_info->eid)->bRm)
            {
                continue;
            }
            // found
            vLfE.push_back(m_Adj[u][i].eid);
            vRtE.push_back(prc_info->eid);
            itNext = ++prc_info;
        }
        else
        {
            // not found
            itNext = prc_info;
        }
    }

    return vLfE.size();
}
/*****************
input:
        int x
        int y
description:
        find node
******************/
TPST_E *myG::findNode(int x, int y)
{
    int u = x;
    int v = y;
    if (m_Adj[x].size() > m_Adj[y].size())
    {
        u = y;
        v = x;
    }
    auto prc_info = std::lower_bound(m_Adj[u].begin(), m_Adj[u].end(), v,
        [](const TPST_ADJ& info, int value)
        {
            return info.pid < value;
        });
    if ((prc_info != m_Adj[u].end()) && (v == prc_info->pid))
    {
        // found
        return findNode(prc_info->eid);
    }
    return NULL;
}

/*****************
input:
        int x
        int y
description:
        find node
******************/
TPST_E *myG::findNode(int eid)
{
    return &(m_Edges[eid]);
}













