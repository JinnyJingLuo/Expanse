/************************************************************************
 *
 *      Module:      FreeInitArrays
 *      Description: Contains functions needed to release temporary
 *                   most memory allocated for the inData struct.
 *
 *      Includes functions:
 *          FreeInNodeArray()
 *          FreeInitArrays()
 *
 ************************************************************************/
#include "Home.h"
#include "InData.h"
#include "Decomp.h"

void FreeInNodeArray(InData_t *inData, int numNodes) {
  int i;

  if (inData->node == (Node_t *)NULL) {
    return;
  }

  /*
   *      For each node in inData->node, free the node's arm arrays
   */
  for (i = 0; i < numNodes; i++) {
    FreeNodeArms(&inData->node[i]);
  }

  free(inData->node);

  inData->node = (Node_t *)NULL;
}

void FreeInitArrays(Home_t *home, InData_t *inData) {
  /*
   *      Only the domains actively involved in reading nodal data
   *      during initialization will have allocated InData_t
   *      structures, so don't try to free structures that have
   *      not been allocated.
   */
  if (inData == NULL) {
    return;
  }

  /*
   *      Free the inData node array and all arrays associated
   *      with each node's arms.
   */
  FreeInNodeArray(inData, inData->param->nodeCount);

  /*
   *      Memory associated with the domain decomposition is dependent
   *      on the type of decomposition used, so invoke a generic function
   *      that will take the appropriate actions based on the decomposition
   *      type.
   */
  if (inData->decomp != NULL) {
    FreeDecomp(home, inData->decomp);
  }
  inData->decomp = (void *)NULL;
}
