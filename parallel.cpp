#include "headers.hpp"

extern int myRank;
extern int nbTasks;

//================================================================================
// Build the list of nodes for MPI communications
//================================================================================

void buildListsNodesMPI(Mesh& m)
{
  if(myRank == 0)
    printf("== buildListsNodesMPI\n");

  //==== Build list with the nodes belonging to current MPI process (i.e. interior + interface)
  
  IntVector maskNodesPart(m.nbOfNodes);
  for(int n=0; n<m.nbOfNodes; n++){
    maskNodesPart(n) = 0;
  }
  for(int iTri=0; iTri<m.nbOfTri; iTri++){
    if(m.triPart(iTri) == myRank){
      int n0 = m.triNodes(iTri, 0);
      int n1 = m.triNodes(iTri, 1);
      int n2 = m.triNodes(iTri, 2);
      maskNodesPart(n0) = 1;
      maskNodesPart(n1) = 1;
      maskNodesPart(n2) = 1;
    }
  }
  
  m.nodesPart.resize(m.nbOfNodes);
  int count = 0;
  for(int n=0; n<m.nbOfNodes; n++){
    if(maskNodesPart(n) == 1){
      m.nodesPart(count) = n;
      count++;
    }
  }
  if(count > 0)
    m.nodesPart.conservativeResize(count);
  m.numNodesPart = count;
  
  //==== Build list with the nodes belonging to both current & neighboring MPI processes (i.e. interface)
  
  IntMatrix maskNodesToExch(nbTasks, m.nbOfNodes);
  for(int nTask=0; nTask<nbTasks; nTask++){
    for(int n=0; n<m.nbOfNodes; n++){
      maskNodesToExch(nTask,n) = 0;
    }
  }
  for(int iTri=0; iTri<m.nbOfTri; iTri++){
    if(m.triPart(iTri) != myRank){
      int n0 = m.triNodes(iTri, 0);
      int n1 = m.triNodes(iTri, 1);
      int n2 = m.triNodes(iTri, 2);
      if(maskNodesPart(n0) == 1)
        maskNodesToExch(m.triPart(iTri),n0) = 1;
      if(maskNodesPart(n1) == 1)
        maskNodesToExch(m.triPart(iTri),n1) = 1;
      if(maskNodesPart(n2) == 1)
        maskNodesToExch(m.triPart(iTri),n2) = 1;
    }
  }
  
  m.numNodesToExch.resize(nbTasks);
  m.nodesToExch.resize(nbTasks,m.nbOfNodes);
  for(int nTask=0; nTask<nbTasks; nTask++){
    m.numNodesToExch(nTask) = 0;
    if(nTask != myRank){
      count = 0;
      for(int n=0; n<m.nbOfNodes; n++){
        if(maskNodesToExch(nTask,n) == 1){
          m.nodesToExch(nTask,count) = n;
          count++;
        }
      }
      m.numNodesToExch(nTask) = count;
      printf("   -> task %i send/recv %i nodes with task %i\n", myRank, m.numNodesToExch(nTask), nTask);
    }
  }
  if(m.numNodesToExch.maxCoeff() > 0)
    m.nodesToExch.conservativeResize(nbTasks, m.numNodesToExch.maxCoeff());
  
}

//================================================================================
// Build the local numbering for the nodes belonging to the current MPI process
//================================================================================

void buildLocalNumbering(Mesh& m)
{
  if(myRank == 0)
    printf("== buildLocalNumbering\n");

  //==== Local numbering for nodes
  
  IntVector localNumNodes(m.nbOfNodes);
  for(int n=0; n<m.nbOfNodes; n++){
    localNumNodes(n) = -1;
  }
  for(int nLoc=0; nLoc<m.numNodesPart; nLoc++){
    int nGlo = m.nodesPart(nLoc);
    localNumNodes(nGlo) = nLoc;
  }
  
  //==== Re-numbering, from nLoc to nGlo
  
  for(int nTask=0; nTask<nbTasks; nTask++){
    for(int n=0; n<m.numNodesToExch(nTask); n++){
      int nGlo = m.nodesToExch(nTask,n);
      m.nodesToExch(nTask,n) = localNumNodes(nGlo);
    }
  }
  
  //==== Build local arrays for nodes/lines/triangles
  
  Matrix    coordsMyRank(m.nbOfNodes,3);
  IntVector linNumMyRank(m.nbOfLin);
  IntVector triNumMyRank(m.nbOfTri);
  IntMatrix linNodesMyRank(m.nbOfLin,2);
  IntMatrix triNodesMyRank(m.nbOfTri,3);
  
  int nNodeLoc = 0;
  int nLinLoc = 0;
  int nTriLoc = 0;
  for(int n=0; n<m.nbOfNodes; n++){
    if(localNumNodes(n) >= 0){
      coordsMyRank(nNodeLoc,0) = m.coords(n,0);
      coordsMyRank(nNodeLoc,1) = m.coords(n,1);
      coordsMyRank(nNodeLoc,2) = m.coords(n,2);
      nNodeLoc++;
    }
  }
  for(int iLin=0; iLin<m.nbOfLin; iLin++){
    if(m.linPart(iLin) == myRank){
      linNumMyRank(nLinLoc) = m.linNum(iLin);
      linNodesMyRank(nLinLoc,0) = localNumNodes(m.linNodes(iLin,0));
      linNodesMyRank(nLinLoc,1) = localNumNodes(m.linNodes(iLin,1));
      nLinLoc++;
    }
  }
  for(int iTri=0; iTri<m.nbOfTri; iTri++){
    if(m.triPart(iTri) == myRank){
      triNumMyRank(nTriLoc) = m.triNum(iTri);
      triNodesMyRank(nTriLoc,0) = localNumNodes(m.triNodes(iTri,0));
      triNodesMyRank(nTriLoc,1) = localNumNodes(m.triNodes(iTri,1));
      triNodesMyRank(nTriLoc,2) = localNumNodes(m.triNodes(iTri,2));
      nTriLoc++;
    }
  }
  
  coordsMyRank.conservativeResize(nNodeLoc,3);
  linNumMyRank.conservativeResize(nLinLoc);
  triNumMyRank.conservativeResize(nTriLoc);
  linNodesMyRank.conservativeResize(nLinLoc,2);
  triNodesMyRank.conservativeResize(nTriLoc,3);
  
  m.nbOfNodes = nNodeLoc;
  m.nbOfLin = nLinLoc;
  m.nbOfTri = nTriLoc;
  
  m.coords = coordsMyRank;
  m.linNum = linNumMyRank;
  m.linNodes = linNodesMyRank;
  m.triNum = triNumMyRank;
  m.triNodes = triNodesMyRank;
}

//================================================================================
// MPI-parallel exchange/add the interface terms
//================================================================================

void exchangeAddInterfMPI(Vector& vec, Mesh& m)
{
  MPI_Request *requestSnd;
  MPI_Request *requestRcv;
  MPI_Status status;
  requestSnd = new MPI_Request[nbTasks];
  requestRcv = new MPI_Request[nbTasks];
  
  double **bufferSnd;
  double **bufferRcv;
  bufferSnd = new double*[nbTasks];
  bufferRcv = new double*[nbTasks];
  
  for(int nTask=0; nTask<nbTasks; nTask++){
    int numToExch = m.numNodesToExch(nTask);
    if(numToExch > 0){
      bufferSnd[nTask] = new double[numToExch];
      bufferRcv[nTask] = new double[numToExch];
      for(int nExch=0; nExch<numToExch; nExch++)
        bufferSnd[nTask][nExch] = vec(m.nodesToExch(nTask,nExch));
      MPI_Isend(bufferSnd[nTask], numToExch, MPI_DOUBLE, nTask, 0, MPI_COMM_WORLD, &requestSnd[nTask]);
      MPI_Irecv(bufferRcv[nTask], numToExch, MPI_DOUBLE, nTask, 0, MPI_COMM_WORLD, &requestRcv[nTask]);
    }
  }
  
  for(int nTask=0; nTask<nbTasks; nTask++){
    int numToExch = m.numNodesToExch(nTask);
    if(numToExch > 0){
      MPI_Wait(&requestRcv[nTask], &status);
      for(int nExch=0; nExch<numToExch; nExch++)
        vec(m.nodesToExch(nTask,nExch)) += bufferRcv[nTask][nExch];
      delete bufferRcv[nTask];
      MPI_Wait(&requestSnd[nTask], &status);
      delete bufferSnd[nTask];
    }
  }
  
  delete[] bufferSnd;
  delete[] bufferRcv;
  delete requestSnd;
  delete requestRcv;
}
