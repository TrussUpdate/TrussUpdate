#pragma once
/***************
truss decomposition
****************/

/* truss decomposition */
#define DE_TRUSS_P2R(gamma, minP, accu) (gamma<minP?0:(int)(gamma/accu))
#define DE_TRUSS_R2P(R, accu) (R>0?(R*accu):0)

#define MAX_R(ldAccu) ((int)(1 / ldAccu))

int deTruss_rePek(myG &obG, int iEid, int iCurK, bool bRmValid);

int deTruss_Detm(myG &obG, vector<vector<int> > &vKG);
int deTruss_ByP(myG &obG, vector<vector<int> > &vKG, double ldMiniP, double ldAccu);
int deTruss_PT(myG &obG, vector<vector<int> > &vKG, vector<vector<pair<int, double> > > &vRetPos, double ldMiniP, double ldAccu);
int deTruss_verify(myG &obG, vector<int> &vSubG, int iResR, int iCurk, double ldMiniP, double ldAccu, pair<int, int> &paXY);

int deTruss_update(myG &obG, vector<TPST_INDEX_EDGE> &vQuery, char *pcIndexPath, bool bStore, char *pcSavePath);
int deTruss_save(myG &obG, vector<vector<int> > &vKG, vector<vector<pair<int, double> > > &vPos);



