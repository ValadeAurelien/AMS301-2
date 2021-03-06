#include "headers.hpp"

extern int myRank;
extern int nbTasks;

//================================================================================
// Build the list of nodes for MPI communications
//================================================================================

void buildListsNodesMPI(Mesh& m)
{
  if(myRank == 0)
    printf("#== buildListsNodesMPI\n");

  //==== Build list with the nodes belonging to current MPI process (i.e. interior + interface)
  
  IntVector maskNodesPart;
  IntMatrix nodesAssignement;
  
  if(myRank == 0) {
        nodesAssignement.resize(nbTasks - 1, m.nbOfNodes);
        m.nodesToExchFrom0.resize(nbTasks - 1, m.nbOfNodes);
        for (int n = 0; n < m.nbOfNodes; n++) {
            for (int nTask = 0; nTask < nbTasks - 1; ++nTask) {
                nodesAssignement(nTask, n) = 0;
                nodesToExchFrom0(nTask, n) = 0;
            }
        }
    } else {
        maskNodesPart.resize(m.nbOfNodes);
        for (int n = 0; n < m.nbOfNodes; n++) {
            maskNodesPart(n) = 0;
        }
    }

    int n0, n1, n2;
    for (int iTri = 0; iTri < m.nbOfTri; iTri++) {
        n0 = m.triNodes(iTri, 0);
        n1 = m.triNodes(iTri, 1);
        n2 = m.triNodes(iTri, 2);
        if (m.triPart(iTri) == myRank) {
            maskNodesPart(n0) = 1;
            maskNodesPart(n1) = 1;
            maskNodesPart(n2) = 1;
        }else if (myRank == 0) {
            nodesAssignement(m.triPart(iTri)-1, n0) = 1;
            nodesAssignement(m.triPart(iTri)-1, n1) = 1;
            nodesAssignement(m.triPart(iTri)-1, n2) = 1;
        }
    }

  
  //==== Build list with the nodes belonging to both current & neighboring MPI processes (i.e. interface),
  //==== remove these nodes from our part and store it in rank 0
 
  int rankOfTri, count[nbTasks];
  for(int nTask = 0; nTask < nbTasks; ++nTask){
      count[nTask] = 0;
  }
    for (int iTri = 0; iTri < m.nbOfTri; iTri++) {
        rankOfTri = m.triPart(iTri);
        if (rankOfTri != myRank) {
            n0 = m.triNodes(iTri, 0);
            n1 = m.triNodes(iTri, 1);
            n2 = m.triNodes(iTri, 2);
            if (myRank == 0) {
                for (int nTask = 0; nTask < nbTasks; ++nTask) {
                    if (nTask != rankOfTri) {
                        if (nodesAssignement(nTask, n0) == 1) {
                            m.nodesToExchFrom0(nTask, count[nTask]++) = n0;
                        }
                        if (nodesAssignement(nTask, n1) == 1) {
                            m.nodesToExchFrom0(nTask, count[nTask]++) = n1;
                        }
                        if (nodesAssignement(nTask, n2) == 1) {
                            m.nodesToExchFrom0(nTask, count[nTask]++) = n2;
                        }
                    }
                }
            }
            if (maskNodesPart(n0) == 1) maskNodesPart(n0) = 0;
            if (maskNodesPart(n1) == 1) maskNodesPart(n1) = 0;
            if (maskNodesPart(n2) == 1) maskNodesPart(n2) = 0;
        }
    }
  
  
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
      printf("#   -> task %i send/recv %i nodes with task %i\n", myRank, m.numNodesToExch(nTask), nTask);
    }
  }
  
  if(m.numNodesToExch.maxCoeff() > 0)
    m.nodesToExch.conservativeResize(nbTasks, m.numNodesToExch.maxCoeff());
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
}

//================================================================================
// Build the local numbering for the nodes belonging to the current MPI process
//================================================================================

void buildLocalNumbering(Mesh& m)
{
  if(myRank == 0)
    printf("#== buildLocalNumbering\n");

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

void computeL2Err(double& L2_err, Vector& uNum, Vector& uExa, Mesh& m) {
    Vector uErr = (uNum - uExa).cwiseAbs();
    double L2_err_loc = 0;
    
    //==== Compute interface error
    int numToExch, nNode;
    std::vector<int> nNodesToComputeAPriori, nNodesToRemove;
    for (int nTask = 0; nTask < nbTasks; ++nTask) {
        numToExch = m.numNodesToExch(nTask);
        for (int nExch = 0; nExch < numToExch; ++nExch) {
            nNode = m.nodesToExch(nTask, nExch);
            if(nTask < myRank) nNodesToComputeAPriori.push_back(nNode);
            else if(nTask > myRank) nNodesToRemove.push_back(nNode);
        }
    }
    // Computation in O(nlog(n)) (sort + remove all but the 1st elements
    // from every consecutive group of equal elements).
    std::sort(nNodesToComputeAPriori.begin(), nNodesToComputeAPriori.end());
    std::vector<int>::iterator last = std::unique(nNodesToComputeAPriori.begin(),
            nNodesToComputeAPriori.end());
    nNodesToComputeAPriori.erase(last, nNodesToComputeAPriori.end());

    std::sort(nNodesToRemove.begin(), nNodesToRemove.end());
    last = std::unique(nNodesToRemove.begin(), nNodesToRemove.end());
    nNodesToRemove.erase(last, nNodesToRemove.end());

    // Computation in O(n) (remove element present both in nNodesToRemove and 
    // in nNodesToComputeAPriori).
    std::vector<int> nNodesToCompute;
    std::set_difference(nNodesToComputeAPriori.begin(),
            nNodesToComputeAPriori.end(),
            nNodesToRemove.begin(),
            nNodesToRemove.end(),
            std::inserter(nNodesToCompute, nNodesToCompute.begin()));


    for (std::vector<int>::iterator it = nNodesToCompute.begin(); it != nNodesToCompute.end(); ++it) {
        L2_err_loc += pow(uErr(*it), 2);
        //std::cout << m.coords.row(*it) << " " << myRank << std::endl;
    }

    //==== Compute local interior error
    for (int nPart = 0; nPart < m.numNodesPart; nPart++) {
        // Check if the node won't be exchange
        if (std::find(nNodesToComputeAPriori.begin(), 
                      nNodesToComputeAPriori.end(), 
                      m.nodesPart(nPart))           == nNodesToComputeAPriori.end() &&
            std::find(nNodesToRemove.begin(), 
                      nNodesToRemove.end(), 
                      m.nodesPart(nPart))           == nNodesToRemove.end()) {
            
            L2_err_loc += pow(uErr(nPart), 2);
            //std::cout << m.coords.row(m.nodesPart(nPart)) << " " << myRank << std::endl;
        }
    }
      
    MPI_Reduce(&L2_err_loc, &L2_err, 1, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);  
    
    if(myRank == 0){
        L2_err = sqrt(L2_err); 
        printf("#== Affichage Erreur \n#   -> Erreur L2 : %f\n", L2_err); 
    }
}

void sendResTo0(double res2) {
    MPI_Send(&res2, 1, MPI_DOUBLE, 0, 0,
	     MPI_COMM_WORLD);
}

void getAndSumRes(double& res2, Vector& otherRes,
		  std::vector<MPI_Request>& AllRequests) {
    for (int task=1; task<nbTasks; ++task)
	MPI_Irecv(&(otherRes(task-1)), 1, MPI_DOUBLE, task, 0,
		  MPI_COMM_WORLD, &(AllRequests.at(task-1)));
    for (int task=1; task<nbTasks; ++task) {
	MPI_Wait(&(AllRequests.at(task-1)), MPI_STATUS_IGNORE);
	res2 += otherRes(task-1);
    }
}
    

