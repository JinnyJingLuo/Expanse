/**************************************************************************
 *
 *      Module:      ForwardEulerIntegrator.c
 *      Description: Contains functions necessary to implement an
 *                   explicit forwardeuler timestep integration method
 *
 *
 ***************************************************************************/
#include "Home.h"
#include "Comm.h"

#define MIN_DELTA_T 1.0e-11

/*------------------------------------------------------------------------
 *
 *      Function:    PreserveNodalData
 *      Description: Copies the current nodal position/velocity
 *                   data to <old*> values for all local and ghost nodes.
 *                   These old values are needed during the timestep
 *                   integration routines and for calculating plastic
 *                   strain.
 *
 *-----------------------------------------------------------------------*/
static void PreserveNodalData(Home_t *home) {
  int i, j;
  Node_t *node;
  int domainIdx, totRemDomCount;
  RemoteDomain_t *remDom;
  Node_t *poNode = NULL;

  for (i = 0; i < home->newNodeKeyPtr; i++) {
    if ((node = home->nodeKeys[i]) == NULL) {
      continue;
    }
    node->oldx = node->x;
    node->oldy = node->y;
    node->oldz = node->z;
  }
  node = home->ghostNodeQ;
  while (node != NULL) {
    node->oldx = node->x;
    node->oldy = node->y;
    node->oldz = node->z;
    node = node->next;
  }
  // junjie: we need to update the position of ghost nodes. Ghost nodes are
  // stored at 2 positions:
  // 1. home->ghostNodeQ 2. home->remoteDomainKeys[remDomID]->nodeKeys
  // the old version only updates ghostNodeQ, However, other subroutine
  // following ForwadEulerIntegrator are all obtaining ghost nodes from 2. So we
  // also need to updates 2.
  totRemDomCount = home->remoteDomainCount + home->secondaryRemoteDomainCount;
  for (i = 0; i < totRemDomCount; i++) {
    domainIdx = home->remoteDomains[i];
    remDom = home->remoteDomainKeys[domainIdx];
    if (remDom == (RemoteDomain_t *)NULL) {
      Fatal("ForwardEulerIntegrator: Missing remote domain struct when "
            "updating ghosts!");
    }
    for (j = 0; j < remDom->maxTagIndex; j++) {
      poNode = remDom->nodeKeys[j];
      if (poNode == NULL) {
        continue;
      }
      poNode->oldx = poNode->x;
      poNode->oldy = poNode->y;
      poNode->oldz = poNode->z;
    }
  }
  // junjie
}

/*---------------------------------------------------------------------------
 *
 *      Function:       AdvanceAllNodes
 *      Description:    Advance all nodes (local and ghost) in time
 *                      and space by a specified amount.  Note: the
 *                      function assumes current nodal positions and
 *                      velocities have already been preserved in
 *                      the old* variables by a call to PreserveNodalData().
 *
 *-------------------------------------------------------------------------*/
static void AdvanceAllNodes(Home_t *home) {
  int i, j;
  real8 dt;
  Node_t *node;
  int domainIdx, totRemDomCount;
  RemoteDomain_t *remDom;
  Node_t *poNode = NULL;

  dt = home->param->deltaTT;

  for (i = 0; i < home->newNodeKeyPtr; i++) {
    if ((node = home->nodeKeys[i]) == NULL) {
      continue;
    }
    // for local nodes, compute the plastic strain due to this displacement
    // we do this by calculating the swept area, which is the cross product of
    // the velocity vector and the CURRENT line direction of the segment (before
    // the displacement) do this only once per segment

    // if the segment is remote-remote, this calculation won't be carried out
    // anyway if it is native native, it will be carried out twice, once per
    // node, the first will compute the area swept by moving the first node, the
    // second will compute the area swept by moving the second node. if it is
    // native-remote, it will be carried out once to compute the area swept due
    // to the motion of this node, if the remote node position hasn't been
    // updated yet, that's fine, otherwise, we are calculating only one half of
    // the area swept, the other half will be calculated by the processor owning
    // the remote node and having the native node as a ghost node overthere.

    node->x = node->oldx + node->vX * dt;
    node->y = node->oldy + node->vY * dt;
    node->z = node->oldz + node->vZ * dt;
    GetPrimaryImage(home->param, &node->x, &node->y, &node->z);
  }

  /*
   *      Also need to move Ghost nodes
   */
  node = home->ghostNodeQ;

  while (node) {
    node->x = node->oldx + node->vX * dt;
    node->y = node->oldy + node->vY * dt;
    node->z = node->oldz + node->vZ * dt;
    GetPrimaryImage(home->param, &node->x, &node->y, &node->z);
    node = node->next;
  }

  // junjie: we need to update the position of ghost nodes. Ghost nodes are
  // stored at 2 positions:
  // 1. home->ghostNodeQ 2. home->remoteDomainKeys[remDomID]->nodeKeys
  // the old version only updates ghostNodeQ, However, other subroutine
  // following ForwadEulerIntegrator are all obtaining ghost nodes from 2. So we
  // also need to updates 2.
  totRemDomCount = home->remoteDomainCount + home->secondaryRemoteDomainCount;
  for (i = 0; i < totRemDomCount; i++) {
    domainIdx = home->remoteDomains[i];
    remDom = home->remoteDomainKeys[domainIdx];
    if (remDom == (RemoteDomain_t *)NULL) {
      Fatal("ForwardEulerIntegrator: Missing remote domain struct when "
            "updating ghosts!");
    }
    for (j = 0; j < remDom->maxTagIndex; j++) {
      poNode = remDom->nodeKeys[j];
      if (poNode == NULL) {
        continue;
      }
      poNode->x = poNode->oldx + poNode->vX * dt;
      poNode->y = poNode->oldy + poNode->vY * dt;
      poNode->z = poNode->oldz + poNode->vZ * dt;
      GetPrimaryImage(home->param, &poNode->x, &poNode->y, &poNode->z);
    }
  }
  // junjie
}

/*------------------------------------------------------------------------
 *
 *      Function:    ForwardEulerIntegrator()
 *      Description: Use the current nodal velocities, the maximum flight
 *                   distance (param->rmax), and the change in nodal
 *                   velocities since the previous step to determine
 *                   the duration of the next timestep.
 *
 *                   Before returning to the caller, the subroutine
 *                   will advance all nodes to their new positions
 *                   and recalculate the force/velocity of all nodes
 *                   at their new positions.
 *
 *------------------------------------------------------------------------*/
void ForwardEulerIntegrator(Home_t *home) {
  int i, doAll = 1;
  real8 vx, vy, vz, v2, vmax, vmax2;
  real8 deltaTT;
  Param_t *param;
  Node_t *node;

  param = home->param;

  vmax = 0.0;
  vmax2 = 0.0;

  /*
   *      Loop over all the nodes, obtain the maximum velocity from among
   *      all nodes, and calculate what the new timestep would be based
   *      solely on the velocty changes since the previous timestep.
   */
  for (i = 0; i < home->newNodeKeyPtr; i++) {
    node = home->nodeKeys[i];
    if (node == NULL) {
      continue;
    }
    if (node->constraint == SURFACE_NODE) {
      continue;
    }
    vx = node->vX;
    vy = node->vY;
    vz = node->vZ;
    v2 = vx * vx + vy * vy + vz * vz;
    if (vmax2 < v2) {
      vmax2 = v2;
    }
  }
  vmax = sqrt(vmax2);

  /*
   *      No node is permitted to move more than a distance of
   *      param->rmax in a single timestep.  Use the highest nodal
   *      velocity to calculate the maximum time step delta that
   *      can be used while still limiting the flight distance of
   *      the fastest node to rmax.
   */
  double dMaxTimeStep = 1.0e-7;
  if (vmax > 0.0) {
    deltaTT = param->rmax / vmax;
  } else {
    deltaTT = dMaxTimeStep;
  }

  if (deltaTT > dMaxTimeStep) {
    deltaTT = dMaxTimeStep;
  }

  // the time step cannot suddenly increase, the new time step is at most
  // 1.2 times the old time step
  double dMaxTimeStepIncrease = 1.2;
  if (deltaTT / param->deltaTT > dMaxTimeStepIncrease) {
    deltaTT = dMaxTimeStepIncrease * param->deltaTT;
  }

#ifdef PARALLEL
  MPI_Allreduce(&deltaTT, &param->deltaTT, 1, MPI_DOUBLE, MPI_MIN,
                MPI_COMM_WORLD);
#else
  param->deltaTT = deltaTT;
#endif
  /*
   *      Update the force and velocity for ghost nodes, and then reposition
   *      all nodes based on their velocities and the new timestep.
   */
  PreserveNodalData(home);
  AdvanceAllNodes(home);
}
