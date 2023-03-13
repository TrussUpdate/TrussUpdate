#pragma once
/***************
calculate pek
****************/

int calPek_add(myG &mapG, int iEid, int iLfEid, int iRtEid, vector <long double> &vldPekXY);
int calPek_rm(myG &mapG, int iEid, int iLfEid, int iRtEid, vector <long double> &vldPekXY);
int cal_OneEdgePek(myG &mapG, vector<int> &vLfE, vector<int> &vRtE, vector <long double> &vfDtE);
bool cal_check(vector <long double> &vfDtE);
int calPek_part(myG &OmapG, vector<int> &vEdges);

