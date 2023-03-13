#pragma once
/***************
input file output file function
****************/

#define ONE_LINE_BUFFER 100
#define FIRST_LINE 100
#define FILE_NAME_BUFFER 100

typedef struct tpstIndexEdge
{
    int x;
    int y;
    float p;
}TPST_INDEX_EDGE;
/* read file and fill the list g_lstFindEdgeByP */
int file_fillG(char *szFileName, myG &mpG);

//int file_saveIndex(myG &mpG, vector<vector<int> > &vKG, vector<vector<int> > &vPKG, vector<vector<pair<int, double> > > &vRetPos, char *pcSavePath);
int file_readIndex(myG &mpG, vector<vector<int> > &vKG, vector<vector<pair<int, double> > > &vRetPos, char *pcSavePath);
int file_printIndex(char *pcFile);

int file_saveBasic(myG &obG, vector<vector<int> > &vKG, char *pcSavePath);
int file_saveIndex(myG &mpG, int k, vector<int> &vEdges, vector<pair<int, double> > &vPos, char *pcSavePath);
