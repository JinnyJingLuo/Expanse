/***************************************************************************
 *
 *      Mobility.h      Declare functioin prototypes and other relevant
 *			mobility related data needed by modules invoking
 *			the mobility laws.
 *
 **************************************************************************/

#ifndef _MOBILITY_H
#define _MOBILITY_H

#include "Home.h"

/*
 *      Define some values for indicating the type of crystal
 *      structure the material has (determined by the mobility
 *      function selected).  Used mainly in sections of the
 *      code that behave differently for BCC and FCC materials.
 *      By grouping things this way we don't have to check for every
 *      single BCC or FCC mobility in those areas, and in
 *      particular don't have to modify a dozen pieces of
 *      code to add a new mobility function.
 */
typedef enum {
  MAT_TYPE_BCC = 0,
  MAT_TYPE_FCC,
} MatType_t;

/*
 *      Define a set of integer values corresponding to the
 *      various types of available mobility laws.  During
 *      initialization the mobilityType parameter will be
 *      set to one of these values which can be used throughout
 *      the code to determine the current mobility in use.
 *      this is easier than using strcmp() everywhere we want
 *      to know make decisions based on the mobility law.
 */
typedef enum {
  MOB_FCC_0,
} MobType_t;

int Mobility_FCC_0(Home_t *home, Node_t *node);

#endif /* _MOBILITY_H */
