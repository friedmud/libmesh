#include <libmesh/lumped_matrix.h>

// test includes
#include "test_comm.h"

#include <libmesh/ignore_warnings.h>
#include <cppunit/extensions/HelperMacros.h>
#include <cppunit/TestCase.h>
#include <libmesh/restore_warnings.h>
#include <libmesh/dense_matrix.h>

using namespace libMesh;

class LumpedMatrixTest : public CppUnit::TestCase
{
public:
  void setUp() {}

  void tearDown() {}

  CPPUNIT_TEST_SUITE( LumpedMatrixTest );

  CPPUNIT_TEST( testAdd );
  CPPUNIT_TEST( testAddMatrix );

  CPPUNIT_TEST_SUITE_END();

private:
  void testAdd()
  {
    LumpedMatrix<double> lm(*TestCommWorld);

    lm.init(4,4,4,4);

    lm.add(1,2, 4.3);

    lm.add(1,1, 7.2);

    lm.close();

    CPPUNIT_ASSERT_DOUBLES_EQUAL(11.5,
                                 lm(1,1),
                                 libMesh::TOLERANCE*libMesh::TOLERANCE);
  }

  void testAddMatrix()
  {
    LumpedMatrix<double> lm(*TestCommWorld);

    lm.init(8,8,8,8);

    DenseMatrix<double> dm(2, 2);

    dm(0,0) = 1.2;
    dm(0,1) = 2.3;
    dm(1,0) = 3.4;
    dm(1,1) = 4.5;

    lm.add_matrix(dm, {2,4}, {2,5});

    lm.close();

    CPPUNIT_ASSERT_DOUBLES_EQUAL(3.5,
                                 lm(2,2),
                                 libMesh::TOLERANCE*libMesh::TOLERANCE);

    CPPUNIT_ASSERT_DOUBLES_EQUAL(7.9,
                                 lm(4,4),
                                 libMesh::TOLERANCE*libMesh::TOLERANCE);

  }
};

CPPUNIT_TEST_SUITE_REGISTRATION( LumpedMatrixTest );
