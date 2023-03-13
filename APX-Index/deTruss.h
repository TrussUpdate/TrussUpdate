#pragma once
/***************
truss decomposition
****************/

/* truss decomposition */

int deTruss_Detm(myG &obG, vector<vector<int> > &vKG);
int deTruss_ByP(myG &obG, vector<vector<int> > &vKG, double ldMiniP, double ldAccu);
int deTruss_PT(myG &obG, const vector<vector<int> > &vKG, double ldMiniP, double ldAccu, bool bStore, char *pcSavePath);
int deTruss_verify(myG &obG, vector<int> &vSubG, int iResR, int iCurk, double ldMiniP, double ldAccu, pair<int, int> &paXY);




