#pragma once
/***************
input file output file function
****************/

#define ONE_LINE_BUFFER 100
#define FIRST_LINE 100
#define FILE_NAME_BUFFER 100

/* read file and fill the list g_lstFindEdgeByP */
int file_fillG(char *szFileName, myG &mpG);
int file_readBlock(myG &mpG, vector<TPST_INDEX_EDGE> &vDataE, vector<pair<int, int> > &vPos, char *pcSavePath);
int file_readQuery(char *szFileName, vector<TPST_INDEX_EDGE> &query);
int file_saveQuery(char *szFileName, vector<TPST_INDEX_EDGE> &query);

//int file_saveIndex(myG &mpG, vector<vector<int> > &vKG, vector<vector<int> > &vPKG, vector<vector<pair<int, double> > > &vRetPos, char *pcSavePath);
//int file_readIndex(myG &mpG, vector<vector<int> > &vKG, vector<vector<int> > &vPKG, vector<vector<pair<int, double> > > &vRetPos, char *pcSavePath);

int file_saveBasic(myG &obG, vector<vector<int> > &vKG, char *pcSavePath);
int file_saveIndex(myG &mpG, int k, vector<int> &vEdges, vector<pair<int, double> > &vPos, char *pcSavePath);

int file_readBasic(myG &mpG, vector<vector<int> > &vKG, char *pcSavePath);
int file_readIndex(myG &mpG, int iCurK, vector<int> &vEdges, vector<pair<int, double> > &vPos, char *pcSavePath);
