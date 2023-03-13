#pragma once
/***************
calculate pek
****************/

bool calPek_add(myG &mapG, int iEid, int iLfEid, int iRtEid, vector <long double> &vldPekXY, int iCurK);
bool calPek_rm(myG &mapG, int iEid, int iLfEid, int iRtEid, vector <long double> &vldPekXY, int iCurK);
int cal_OneEdgePek(myG &mapG, vector<int> &vLfE, vector<int> &vRtE, vector <long double> &vfDtE);
bool cal_check(vector <long double> &vfDtE);
int calPek_part(myG &OmapG, vector<int> &vEdges);

