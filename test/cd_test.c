#include <roki/roki.h>

/* vertex-based collision detection: the contact frame of a vertex that
 * persists between two checks is restored from a cache expressed in
 * the link frame of the penetrated body.  When that body moves between
 * the checks, the restored normal and axes must still be the unit,
 * mutually orthogonal directions of the contact, and the reference
 * point must move rigidly with the body. */

#define TOL 1.0e-9
#define NUM_VERT 4

typedef struct{
  zVec3D *vert;
  zVec3D ref;
} record_t;

static void set_pose(rkChain *chain, double x, double y, double z, double yaw)
{
  zVec dis;

  dis = zVecAlloc( 6 );
  zVecSetElem( dis, 0, x );
  zVecSetElem( dis, 1, y );
  zVecSetElem( dis, 2, z );
  zVecSetElem( dis, 5, yaw );
  rkChainFK( chain, dis );
  zVecFree( dis );
}

/* the contact frames of all vertices: count them, check the invariants,
 * and record the vertices and their reference points */
static int check_frames(rkCD *cd, rkChain *box, bool *frame_ok, bool *normal_ok, bool *on_box, record_t rec[])
{
  rkCDPair *cp;
  rkCDVert *v;
  zVec3D z = { { 0, 0, 1 } };
  int n = 0, i;

  *frame_ok = *normal_ok = *on_box = true;
  zListForEach( &cd->plist, cp ){
    if( !cp->data.is_col ) continue;
    zListForEach( &cp->data.vlist, v ){
      if( v->data.cell->data.chain != box ) *on_box = false;
      for( i=0; i<3; i++ )
        if( !zIsTol( zVec3DNorm( &v->data.axis[i] ) - 1, TOL ) ) *frame_ok = false;
      if( !zIsTol( zVec3DInnerProd( &v->data.axis[0], &v->data.axis[1] ), TOL ) ||
          !zIsTol( zVec3DInnerProd( &v->data.axis[1], &v->data.axis[2] ), TOL ) ||
          !zIsTol( zVec3DInnerProd( &v->data.axis[2], &v->data.axis[0] ), TOL ) ) *frame_ok = false;
      if( !zVec3DEqual( &v->data.axis[0], &v->data.norm ) ) *frame_ok = false;
      if( !zIsTol( zVec3DNorm( &v->data.norm ) - 1, TOL ) || !zVec3DEqual( &v->data.norm, &z ) ) *normal_ok = false;
      if( n < NUM_VERT ){
        rec[n].vert = v->data.vert;
        zVec3DCopy( &v->data.ref, &rec[n].ref );
      }
      n++;
    }
  }
  return n;
}

int main(int argc, char *argv[])
{
  rkChain floor, box;
  rkCD cd;
  record_t before[NUM_VERT], after[NUM_VERT];
  zVec3D d = { { 0.3, -0.2, 0.05 } }, expected, moved;
  zMat3D rot;
  double yaw = zDeg2Rad( 30.0 );
  bool frame_ok, normal_ok, on_box, ref_ok = true;
  int n1, n2, i, j;

  rkChainReadZTK( &floor, "cd_floor.ztk" );
  rkChainReadZTK( &box, "cd_box.ztk" );
  rkCDCreate( &cd );
  rkCDChainReg( &cd, &floor, RK_CD_CELL_MOVE );
  rkCDChainReg( &cd, &box, RK_CD_CELL_MOVE );

  /* the cube rests on the slab, its bottom face 1 mm into it */
  set_pose( &floor, 0, 0, 0, 0 );
  set_pose( &box, 0, 0, 0.049, 0 );
  rkCDColChkVert( &cd );
  n1 = check_frames( &cd, &box, &frame_ok, &normal_ok, &on_box, before );
  zAssert( rkCDColChkVert (four vertices of the cube penetrate the slab), n1 == NUM_VERT && on_box );
  zAssert( rkCDColChkVert (contact frames at registration), frame_ok && normal_ok );

  /* both bodies move rigidly together: the vertices persist, and the
     contact frames are restored from the cache in the slab's frame */
  zMat3DFromZYX( &rot, yaw, 0, 0 );
  set_pose( &floor, d.c.x, d.c.y, d.c.z, yaw );
  set_pose( &box, d.c.x, d.c.y, d.c.z + 0.049, yaw );
  rkCDColChkVert( &cd );
  n2 = check_frames( &cd, &box, &frame_ok, &normal_ok, &on_box, after );
  zAssert( rkCDColChkVert (the vertices persist), n2 == NUM_VERT && on_box );
  zAssert( rkCDColChkVert (restored contact axes are unit and orthogonal), frame_ok );
  zAssert( rkCDColChkVert (restored contact normal), normal_ok );
  /* the reference points move with the slab */
  for( i=0; i<NUM_VERT; i++ ){
    for( j=0; j<NUM_VERT; j++ )
      if( after[j].vert == before[i].vert ) break;
    if( j == NUM_VERT ){ ref_ok = false; continue; }
    zMulMat3DVec3D( &rot, &before[i].ref, &moved );
    zVec3DAdd( &moved, &d, &expected );
    if( !zVec3DEqual( &after[j].ref, &expected ) ) ref_ok = false;
  }
  zAssert( rkCDColChkVert (restored reference points move with the body), ref_ok );

  rkCDDestroy( &cd );
  rkChainDestroy( &floor );
  rkChainDestroy( &box );
  return 0;
}
