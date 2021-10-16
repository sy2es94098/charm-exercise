#include <string>
using namespace std;

#if LIVEVIZ_RUN
#include "liveViz.h"
#endif

#include "particleSimulation.decl.h"
#include "custom_rand_gen.h"
/*readonly*/ extern CProxy_Main mainProxy;
/*readonly*/ extern CProxy_Cell cellProxy;
/*readonly*/ extern int particlesPerCell;
/*readonly*/ extern int numCellsPerDim;
/*readonly*/ extern int iterations;
/*readonly*/ extern int lbFreq;
/*readonly*/ extern int reductionFreq;
/*readonly*/ extern double boxMax;
/*readonly*/ extern double boxMin;
/*readonly*/ extern double cellDim;


#include "cell.h"
#include "main.h"

extern CkReduction::reducerType minMaxType;

// Useful function declarations
//void Cell::perturb(Particle* particle);
//void Cell::sendParticles(int xIndex, int yIndex, int iteration,  std::vector<Particle> &outgoing);

struct direction {
  int x;
  int y;
} dir;

//change the position of the particles and send messages to neighbors with their incoming particles
void Cell::updateParticles(int iter) {

  // Variables to use
  // 1. vector<Particle> particles (declared in cell.h)
  // 2. startX, endX (declared in cell.h). Example - The cell (2,3) will have startX = 2.0 and endX= 3.0
  // 3. startY, endY (declared in cell.h). Example - The cell (2,3) will have startY = 3.0 and endY = 4.0
  // 4. thisIndex.x represents my cell's x index (declared in the charm++ runtime system). Example - The cell (2,3) will have thisIndex.x as 2
  // 5. thisIndex.y represents my cell's y index (declared in the charm++ runtime system). Example - The cell (2,3) will have thisIndex.y as 3

  //TODO: Add code for the following
  // 1. Iterate through the particles and call perturb(...) passing each particle. This causes the particle's x and y coordinate to change
  // 2. Identify the new cell that the particle belongs to and construct a vector of particles to be sent to each of the 8 neighbors (topLeft, top, topRight, left, right, bottomLeft, bottom, bottomRight)
  // 3. Make sure that particles that go outside the bounding box are wrapped back i.e for a 2D box consisting of 8 cells in each dimension,
  //    if a particle with the index (7,7) goes to (7, 8), it should be sent back to (7, 0).
  // 4. Call sendParticles(...) to send the 8 different vector of particles to each of the 8 neighbours

  const int dirPerDim = 3;
  const int particlesAllCells = particlesPerCell * numCellsPerDim * numCellsPerDim;
  Particle outgoing[dirPerDim][dirPerDim][particlesAllCells];
  int size[dirPerDim][dirPerDim] = {0};
  vector<Particle>::iterator par = particles.begin();
  while (par != particles.end()) {
    perturb(&(*par));
    
    dir.x = 0;
    dir.y = 0;
    
    if (par->x < startX) {
      --dir.x;
      if (par->y < startY) --dir.y;
      else if (par->y > endY) ++dir.y;
    }
    else if (par->x > endX) {
      ++dir.x;
      if (par->y < startY) --dir.y;
      else if (par->y > endY) ++dir.y;
    }
    else {
      if (par->y < startY) --dir.y;
      else if (par->y > endY) ++dir.y;
      else {
        par++;
        continue;
      }
    }
    
    outgoing[dir.x+1][dir.y+1][size[dir.x+1][dir.y+1]++] = *par;
    par = particles.erase(par);
  }

  int x_out, y_out;

  for (int i = -1; i <= 1; i++) {
    x_out = thisIndex.x + i;
    if (x_out < 0) {
      x_out = numCellsPerDim - 1;
    }
    else if (x_out == numCellsPerDim) {
      x_out = 0;
    }

    for (int j = -1; j <= 1; j++) {
      if (i == 0 && j == 0) continue;

      y_out = thisIndex.y + j;
      if (y_out < 0) {
        y_out = numCellsPerDim - 1;
      }
      else if (y_out == numCellsPerDim) {
        y_out = 0;
      }

      vector<Particle> out(outgoing[i+1][j+1], outgoing[i+1][j+1]+size[i+1][j+1]);
      sendParticles(x_out, y_out, iter, out);
    }
  }
}

#if BONUS_QUESTION
void Main::receiveMinMaxReductionData(CkReductionMsg *data) {
  int *output = (int *) data->getData();

  // TODO: Assign values to maxParticles, maxCellX, maxCellY, minParticles, minCellX, minCellY based
  // on computed reduction values

  maxParticles = output[0];
  maxCellX = output[1];
  maxCellY = output[2];
  minParticles = output[3];
  minCellX = output[4];
  minCellY = output[5];

  CmiPrintf("Max Particles:%d, Cell with Max Particles: (%d, %d)\n", maxParticles, maxCellX, maxCellY);
  CmiPrintf("Min Particles:%d, Cell with Min Particles: (%d, %d)\n", minParticles, minCellX, minCellY);
  readyToOutput();
}

void Cell::contributeToReduction() {
  numParticles = particles.size();

  // TODO: Declare a callback with the function receiveMinMaxReductionData
  // Add code to contribute to reduction for finding the different custom reduction values using the callback
  // declared above

  const int num_data = 6;

  CkCallback cb(CkIndex_Main::receiveMinMaxReductionData(NULL), mainProxy);

  int data[num_data] = {numParticles, (int) startX, (int) startY, 
                        numParticles, (int) startX, (int) startY};

  contribute(num_data*sizeof(int), data, minMaxType, cb);
}

CkReductionMsg *calculateMaxMin(int nMsg, CkReductionMsg **msgs) {
  //TODO: Add code to compute the following:
  // Value of the maximum particles per cell
  // X coordinate of the cell with the maximum particles per cell
  // Y coordinate of the cell with the maximum particles per cell
  // Value of the minimum particles per cell
  // X coordinate of the cell with the minimum particles per cell
  // Y coordinate of the cell with the minimum particles per cell

  const int num_data = 6;
  int *retData = (int *) msgs[0]->getData();

  for (int i = 1; i < nMsg; i++) {
    int *data = (int *) msgs[i]->getData();
    if (data[0] > retData[0]) {
      retData[0] = data[0];
      retData[1] = data[1];
      retData[2] = data[2];
    }
    if (data[3] < retData[3]) {
      retData[3] = data[3];
      retData[4] = data[4];
      retData[5] = data[5];
    }
  }
  return CkReductionMsg::buildNew(num_data*sizeof(int), retData, minMaxType);
}
#endif
