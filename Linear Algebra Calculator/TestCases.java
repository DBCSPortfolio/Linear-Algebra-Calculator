public class TestCases
{
   public static void main(String args[]) throws Exception
   {
      // test for real eigenvalues of a 2x2 matrix
      ///*
      System.out.println("\n----------------TEST CASE 1----------------\n");
      double[][] exampleArray = new double[2][2];
      exampleArray[0][0] = 4;
      exampleArray[0][1] = 7;
      exampleArray[1][0] = 1.5;
      exampleArray[1][1] = 2;
      Matrix exampleMatrix = new Matrix(exampleArray);
      
      exampleMatrix.findEigenvalues();
      System.out.println(exampleMatrix.eigenvaluesToString());
      System.out.println("\n--------------------------------\n");
      //*/
      
      ///*
      System.out.println("\n----------------TEST CASE 2----------------\n");
      double[][] exampleRef = new double[3][3];
      exampleRef[0][0] = 4.0;
      exampleRef[0][1] = 0.0;
      exampleRef[0][2] = 1.0;
      exampleRef[1][0] = 2.0;
      exampleRef[1][1] = 1.0;
      exampleRef[1][2] = -2.0;
      exampleRef[2][0] = 3.0;
      exampleRef[2][1] = -3.0;
      exampleRef[2][2] = -5.0;
      Matrix exampleRefMatrix = new Matrix(exampleRef);
      Matrix echelonForm = exampleRefMatrix.ref();
      System.out.println(echelonForm.refToString());
      System.out.println("\n--------------------------------\n");
      //*/
      
      ///*
      System.out.println("\n----------------TEST CASE 3----------------\n");
      double[][] exampleRef2 = new double[4][4];
      exampleRef2[0][0] = 0;
      exampleRef2[0][1] = 1;
      exampleRef2[0][2] = 4;
      exampleRef2[0][3] = 5;
      exampleRef2[1][0] = -2;
      exampleRef2[1][1] = 6;
      exampleRef2[1][2] = 7;
      exampleRef2[1][3] = 3;
      exampleRef2[2][0] = -3;
      exampleRef2[2][1] = 1;
      exampleRef2[2][2] = -1;
      exampleRef2[2][3] = 2;
      exampleRef2[3][0] = 1;
      exampleRef2[3][1] = -1;
      exampleRef2[3][2] = 4;
      exampleRef2[3][3] = 0;
      Matrix exampleRefMatrix2 = new Matrix(exampleRef2);
      Matrix echelonForm2 = exampleRefMatrix2.ref();
      System.out.println(echelonForm2.refToString());
      System.out.println("\n--------------------------------\n");
      //*/
      
      ///*
      System.out.println("\n----------------TEST CASE 4----------------\n");
      double[][] exampleRef3 = new double[3][6];
      exampleRef3[0][0] = 0;
      exampleRef3[0][1] = 3;
      exampleRef3[0][2] = -6;
      exampleRef3[0][3] = 6;
      exampleRef3[0][4] = 4;
      exampleRef3[0][5] = -5;
      exampleRef3[1][0] = 3;
      exampleRef3[1][1] = -7;
      exampleRef3[1][2] = 8;
      exampleRef3[1][3] = -5;
      exampleRef3[1][4] = 8;
      exampleRef3[1][5] = 9;
      exampleRef3[2][0] = 3;
      exampleRef3[2][1] = -9;
      exampleRef3[2][2] = 12;
      exampleRef3[2][3] = -9;
      exampleRef3[2][4] = 6;
      exampleRef3[2][5] = 15;
      Matrix exampleRefMatrix3 = new Matrix(exampleRef);
      Matrix echelonForm3 = exampleRefMatrix3.ref();
      System.out.println(echelonForm3.refToString());
      System.out.println("\n--------------------------------\n");
      //*/
      
      ///*
      System.out.println("\n----------------TEST CASE 5----------------\n");
      double[][] exampleRRef2 = new double[4][4];
      exampleRRef2[0][0] = 0;
      exampleRRef2[0][1] = 1;
      exampleRRef2[0][2] = 4;
      exampleRRef2[0][3] = 5;
      exampleRRef2[1][0] = -2;
      exampleRRef2[1][1] = 6;
      exampleRRef2[1][2] = 7;
      exampleRRef2[1][3] = 3;
      exampleRRef2[2][0] = -3;
      exampleRRef2[2][1] = 1;
      exampleRRef2[2][2] = -1;
      exampleRRef2[2][3] = 2;
      exampleRRef2[3][0] = 1;
      exampleRRef2[3][1] = -1;
      exampleRRef2[3][2] = 4;
      exampleRRef2[3][3] = 0;
      Matrix exampleRRefMatrix2 = new Matrix(exampleRRef2);
      Matrix reducedEchelonForm2 = exampleRRefMatrix2.rref();
      
      System.out.println(reducedEchelonForm2.rRefToString());
      System.out.println(reducedEchelonForm2.generalSolutionZero());
      
      //Matrix transposeExample = ((new Matrix(exampleRRef2)).transpose());
      System.out.println(exampleRRefMatrix2.transposeToString());
      
      // this shows 11 matrices, but there are 14 steps in the output.
      System.out.println("\n--------------------------------\n");
      //*/
      
      ///*
      System.out.println("\n----------------TEST CASE 6----------------\n");
      double[][] exampleRRef3 = new double[3][4];
      exampleRRef3[0][0] = 1;
      exampleRRef3[0][1] = -2;
      exampleRRef3[0][2] = -1;
      exampleRRef3[0][3] = 3;
      exampleRRef3[1][0] = -2;
      exampleRRef3[1][1] = 4;
      exampleRRef3[1][2] = 5;
      exampleRRef3[1][3] = -5;
      exampleRRef3[2][0] = 3;
      exampleRRef3[2][1] = -6;
      exampleRRef3[2][2] = -6;
      exampleRRef3[2][3] = 8;
      Matrix exampleRRefMatrix3 = new Matrix(exampleRRef3);
      Matrix reducedEchelonForm3 = exampleRRefMatrix3.rref();
      System.out.println(reducedEchelonForm3.rRefToString());
      System.out.println("\n--------------------------------\n");
      //*/
      
      ///*
      // testing inverse for 2x2 matrix
      System.out.println("\n----------------TEST CASE 7----------------\n");
      double[][] exampleInverseArray = new double[2][2];
      exampleInverseArray[0][0] = 4;
      exampleInverseArray[0][1] = 7;
      exampleInverseArray[1][0] = 1.5;
      exampleInverseArray[1][1] = 2;
      Matrix exampleInverseMatrix = new Matrix(exampleInverseArray);
      
      Matrix inverseMatrix = exampleInverseMatrix.inverse2();
      System.out.println(inverseMatrix.inverseToString());
      System.out.println("\n--------------------------------\n");
      //*/
      
      ///*
      // testing transpose of a non-square matrix
      System.out.println("\n----------------TEST CASE 8----------------\n");
      double[][] exampleTranspose = new double[2][4];
      exampleTranspose[0][0] = 4;
      exampleTranspose[0][1] = -2;
      exampleTranspose[0][2] = 3;
      exampleTranspose[0][3] = 6;
      exampleTranspose[1][0] = 1;
      exampleTranspose[1][1] = -1;
      exampleTranspose[1][2] = 5;
      exampleTranspose[1][3] = 7;
      
      Matrix newTransposeMatrix = new Matrix(exampleTranspose);
      //Matrix transposeMatrix = newTransposeMatrix.transpose();
      System.out.println(newTransposeMatrix.transposeToString());
      System.out.println("\n--------------------------------\n");
      //*/
      
      ///*
      // testing Cramer's Rule for 2x2 matrix
      System.out.println("\n----------------TEST CASE 9----------------\n");
      double[][] exampleCramer = new double[2][2];
      exampleCramer[0][0] = 3;
      exampleCramer[0][1] = -2;
      exampleCramer[1][0] = -5;
      exampleCramer[1][1] = 4;
      double[] exampleEqualsVector = new double[2];
      exampleEqualsVector[0] = 6;
      exampleEqualsVector[1] = 8;
      Matrix cramerMatrix = new Matrix(exampleCramer);
      Vector cramerVector = new Vector(exampleEqualsVector);
      
      //List<Double> exampleValues = cramerMatrix.cramersRule(cramerVector);
      System.out.println(cramerMatrix.cramersRuleToString(cramerVector));
      System.out.println("\n--------------------------------\n");
      //*/
      
      ///*
      // testing Cramer's Rule for 3x3 matrix
      System.out.println("\n----------------TEST CASE 10----------------\n");
      double[][] exampleCramer2 = new double[3][3];
      exampleCramer2[0][0] = 1;
      exampleCramer2[0][1] = 3;
      exampleCramer2[0][2] = 1;
      exampleCramer2[1][0] = -1;
      exampleCramer2[1][1] = 0;
      exampleCramer2[1][2] = 2;
      exampleCramer2[2][0] = 3;
      exampleCramer2[2][1] = 1;
      exampleCramer2[2][2] = 0;
      double[] exampleEqualsVector2 = new double[3];
      exampleEqualsVector2[0] = 4;
      exampleEqualsVector2[1] = 2;
      exampleEqualsVector2[2] = 2;
      Matrix cramerMatrix2 = new Matrix(exampleCramer2);
      Vector cramerVector2 = new Vector(exampleEqualsVector2);
      
      //List<Double> exampleValues2 = cramerMatrix2.cramersRule(cramerVector2);
      System.out.println(cramerMatrix2.cramersRuleToString(cramerVector2));
      System.out.println("\n--------------------------------\n");
      //*/
      
      
      ///*
      // testing for determinant of matrix of dimension >= 4x4
      System.out.println("\n----------------TEST CASE 11----------------\n");
      double[][] exampleDetFourPlus = new double[7][7];
      exampleDetFourPlus[0][0] = -1;
      exampleDetFourPlus[0][1] = 2;
      exampleDetFourPlus[0][2] = 4;
      exampleDetFourPlus[0][3] = -3;
      exampleDetFourPlus[0][4] = 5;
      exampleDetFourPlus[0][5] = 0;
      exampleDetFourPlus[0][6] = 6;
      exampleDetFourPlus[1][0] = 6;
      exampleDetFourPlus[1][1] = 7;
      exampleDetFourPlus[1][2] = -8;
      exampleDetFourPlus[1][3] = 1;
      exampleDetFourPlus[1][4] = -1;
      exampleDetFourPlus[1][5] = -1;
      exampleDetFourPlus[1][6] = 2;
      exampleDetFourPlus[2][0] = -3;
      exampleDetFourPlus[2][1] = 0;
      exampleDetFourPlus[2][2] = 1;
      exampleDetFourPlus[2][3] = 1;
      exampleDetFourPlus[2][4] = -4;
      exampleDetFourPlus[2][5] = -2;
      exampleDetFourPlus[2][6] = 3;
      exampleDetFourPlus[3][0] = 0;
      exampleDetFourPlus[3][1] = 0;
      exampleDetFourPlus[3][2] = 4;
      exampleDetFourPlus[3][3] = 6;
      exampleDetFourPlus[3][4] = -7;
      exampleDetFourPlus[3][5] = 1;
      exampleDetFourPlus[3][6] = 4;
      exampleDetFourPlus[4][0] = 2;
      exampleDetFourPlus[4][1] = -2;
      exampleDetFourPlus[4][2] = -1;
      exampleDetFourPlus[4][3] = -3;
      exampleDetFourPlus[4][4] = 3;
      exampleDetFourPlus[4][5] = 8;
      exampleDetFourPlus[4][6] = -4;
      exampleDetFourPlus[5][0] = 6;
      exampleDetFourPlus[5][1] = 0;
      exampleDetFourPlus[5][2] = 3;
      exampleDetFourPlus[5][3] = 1;
      exampleDetFourPlus[5][4] = 0;
      exampleDetFourPlus[5][5] = -2;
      exampleDetFourPlus[5][6] = -5;
      exampleDetFourPlus[6][0] = 4;
      exampleDetFourPlus[6][1] = -5;
      exampleDetFourPlus[6][2] = -1;
      exampleDetFourPlus[6][3] = 0;
      exampleDetFourPlus[6][4] = -1;
      exampleDetFourPlus[6][5] = 0;
      exampleDetFourPlus[6][6] = 3;
      
      Matrix detFourPlus = new Matrix(exampleDetFourPlus);
      System.out.println(detFourPlus.determinantToString());
      System.out.println("\n--------------------------------\n");
      //*/
      
      ///*
      // test for real eigenvalues of a 3x3 matrix
      // should be 3, -5, and 6
      System.out.println("\n----------------TEST CASE 12----------------\n");
      double[][] exampleEigen3Array = new double[3][3];
      exampleEigen3Array[0][0] = -2;
      exampleEigen3Array[0][1] = -4;
      exampleEigen3Array[0][2] = 2;
      exampleEigen3Array[1][0] = -2;
      exampleEigen3Array[1][1] = 1;
      exampleEigen3Array[1][2] = 2;
      exampleEigen3Array[2][0] = 4;
      exampleEigen3Array[2][1] = 2;
      exampleEigen3Array[2][2] = 5;
      Matrix exampleEigen3Matrix = new Matrix(exampleEigen3Array);
      
      exampleEigen3Matrix.findEigenvalues();
      System.out.println(exampleEigen3Matrix.eigenvaluesToString());
      System.out.println("\n--------------------------------\n");
      //*/
      
      /*
      // test for real/complex eigenvalues of a 4x4 matrix
      System.out.println("\n----------------TEST CASE 13----------------\n");
      double[][] exampleEigen4Array = new double[4][4];
      exampleEigen4Array[0][0] = -2;
      exampleEigen4Array[0][1] = 3;
      exampleEigen4Array[0][2] = -1;
      exampleEigen4Array[0][3] = 6;
      exampleEigen4Array[1][0] = 4;
      exampleEigen4Array[1][1] = 1;
      exampleEigen4Array[1][2] = -1;
      exampleEigen4Array[1][3] = 0;
      exampleEigen4Array[2][0] = 1;
      exampleEigen4Array[2][1] = 6;
      exampleEigen4Array[2][2] = 7;
      exampleEigen4Array[2][3] = 8;
      exampleEigen4Array[3][0] = -3;
      exampleEigen4Array[3][1] = -3;
      exampleEigen4Array[3][2] = -4;
      exampleEigen4Array[3][3] = -1;
      Matrix exampleEigen4Matrix = new Matrix(exampleEigen4Array);
      
      exampleEigen4Matrix.findEigenvalues();
      System.out.println(exampleEigen4Matrix.eigenvaluesToString());
      System.out.println("\n--------------------------------\n");
      //*/
      
      ///*
      // test for finding 4 real eigenvalues of a 4x4 matrix
      System.out.println("\n----------------TEST CASE 14----------------\n");
      double[][] exampleEigen4Array2 = new double[4][4];
      exampleEigen4Array2[0][0] = 2;
      exampleEigen4Array2[0][1] = -1;
      exampleEigen4Array2[0][2] = -1;
      exampleEigen4Array2[0][3] = 0;
      exampleEigen4Array2[1][0] = -1;
      exampleEigen4Array2[1][1] = 3;
      exampleEigen4Array2[1][2] = -1;
      exampleEigen4Array2[1][3] = -1;
      exampleEigen4Array2[2][0] = -1;
      exampleEigen4Array2[2][1] = -1;
      exampleEigen4Array2[2][2] = 3;
      exampleEigen4Array2[2][3] = -1;
      exampleEigen4Array2[3][0] = 0;
      exampleEigen4Array2[3][1] = -1;
      exampleEigen4Array2[3][2] = -1;
      exampleEigen4Array2[3][3] = 2;
      Matrix exampleEigen4Matrix2 = new Matrix(exampleEigen4Array2);
      
      exampleEigen4Matrix2.findEigenvalues();
      System.out.println(exampleEigen4Matrix2.eigenvaluesToString());
      System.out.println("\n--------------------------------\n");
      //*/
      
      ///*
      System.out.println("\n----------------TEST CASE 15----------------\n");
      double[][] exampleSecondEig3Array = new double[4][4];
      exampleSecondEig3Array[0][0] = 1;
      exampleSecondEig3Array[0][1] = 2;
      exampleSecondEig3Array[0][2] = 3;
      exampleSecondEig3Array[0][3] = 4;
      exampleSecondEig3Array[1][0] = 2;
      exampleSecondEig3Array[1][1] = 3;
      exampleSecondEig3Array[1][2] = 4;
      exampleSecondEig3Array[1][3] = 5;
      exampleSecondEig3Array[2][0] = 3;
      exampleSecondEig3Array[2][1] = 4;
      exampleSecondEig3Array[2][2] = 5;
      exampleSecondEig3Array[2][3] = 6;
      exampleSecondEig3Array[3][0] = 4;
      exampleSecondEig3Array[3][1] = 5;
      exampleSecondEig3Array[3][2] = 6;
      exampleSecondEig3Array[3][3] = 7;
      Matrix exampleSecondEig3Matrix = new Matrix(exampleSecondEig3Array);
      
      exampleSecondEig3Matrix.findEigenvalues();
      System.out.println(exampleSecondEig3Matrix.eigenvaluesToString());
      System.out.println("\n--------------------------------\n");
      //*/
      
      ///*
      System.out.println("\n----------------TEST CASE 16----------------\n");
      double[][] exampleAnotherEig3Array = new double[4][4];
      exampleAnotherEig3Array[0][0] = 1;
      exampleAnotherEig3Array[0][1] = 2;
      exampleAnotherEig3Array[0][2] = 3;
      exampleAnotherEig3Array[0][3] = 4;
      exampleAnotherEig3Array[1][0] = 2;
      exampleAnotherEig3Array[1][1] = 4;
      exampleAnotherEig3Array[1][2] = 6;
      exampleAnotherEig3Array[1][3] = 8;
      exampleAnotherEig3Array[2][0] = 3;
      exampleAnotherEig3Array[2][1] = 6;
      exampleAnotherEig3Array[2][2] = 9;
      exampleAnotherEig3Array[2][3] = 12;
      exampleAnotherEig3Array[3][0] = 4;
      exampleAnotherEig3Array[3][1] = 8;
      exampleAnotherEig3Array[3][2] = 12;
      exampleAnotherEig3Array[3][3] = 16;
      Matrix exampleAnotherEig3Matrix = new Matrix(exampleAnotherEig3Array);
      
      exampleAnotherEig3Matrix.findEigenvalues();
      System.out.println(exampleAnotherEig3Matrix.eigenvaluesToString());
      System.out.println("\n--------------------------------\n");
      //*/
      
      ///*
      System.out.println("\n----------------TEST CASE 17----------------\n");
      double[][] exampleAnotherEig4Array = new double[4][4];
      exampleAnotherEig4Array[0][0] = 0;
      exampleAnotherEig4Array[0][1] = 1;
      exampleAnotherEig4Array[0][2] = 1;
      exampleAnotherEig4Array[0][3] = 0;
      exampleAnotherEig4Array[1][0] = -1;
      exampleAnotherEig4Array[1][1] = 0;
      exampleAnotherEig4Array[1][2] = 0;
      exampleAnotherEig4Array[1][3] = 1;
      exampleAnotherEig4Array[2][0] = 0;
      exampleAnotherEig4Array[2][1] = 0;
      exampleAnotherEig4Array[2][2] = 0;
      exampleAnotherEig4Array[2][3] = 1;
      exampleAnotherEig4Array[3][0] = 0;
      exampleAnotherEig4Array[3][1] = 0;
      exampleAnotherEig4Array[3][2] = -1;
      exampleAnotherEig4Array[3][3] = 0;
      Matrix exampleAnotherEig4Matrix = new Matrix(exampleAnotherEig4Array);
      
      exampleAnotherEig4Matrix.findEigenvalues();
      System.out.println(exampleAnotherEig4Matrix.eigenvaluesToString());
      System.out.println("\n--------------------------------\n");
      //*/
      
      ///*
      System.out.println("\n----------------TEST CASE 18----------------\n");
      double[][] exampleYETAnotherEig4Array = new double[4][4];
      exampleYETAnotherEig4Array[0][0] = 1;
      exampleYETAnotherEig4Array[0][1] = 3;
      exampleYETAnotherEig4Array[0][2] = 0;
      exampleYETAnotherEig4Array[0][3] = 3;
      exampleYETAnotherEig4Array[1][0] = 1;
      exampleYETAnotherEig4Array[1][1] = 1;
      exampleYETAnotherEig4Array[1][2] = 1;
      exampleYETAnotherEig4Array[1][3] = 1;
      exampleYETAnotherEig4Array[2][0] = 0;
      exampleYETAnotherEig4Array[2][1] = 4;
      exampleYETAnotherEig4Array[2][2] = 2;
      exampleYETAnotherEig4Array[2][3] = 8;
      exampleYETAnotherEig4Array[3][0] = 2;
      exampleYETAnotherEig4Array[3][1] = 0;
      exampleYETAnotherEig4Array[3][2] = 3;
      exampleYETAnotherEig4Array[3][3] = 1;
      Matrix exampleYETAnotherEig4Matrix = new Matrix(exampleYETAnotherEig4Array);
      
      exampleYETAnotherEig4Matrix.findEigenvalues();
      System.out.println(exampleYETAnotherEig4Matrix.eigenvaluesToString());
      System.out.println("\n--------------------------------\n");
      //*/
      
      ///*
      System.out.println("\n----------------TEST CASE 19----------------\n");
      double[][] exampleDet2Array = new double[2][2];
      exampleDet2Array[0][0] = 4;
      exampleDet2Array[0][1] = 7;
      exampleDet2Array[1][0] = 1.5;
      exampleDet2Array[1][1] = 2;
      Matrix exampleDet2Matrix = new Matrix(exampleDet2Array);
      
      System.out.println(exampleDet2Matrix.determinantToString());
      System.out.println("\n--------------------------------\n");
      //*/
      
      ///*
      System.out.println("\n----------------TEST CASE 20----------------\n");
      double[][] example3Ref = new double[3][3];
      example3Ref[0][0] = 4.0;
      example3Ref[0][1] = 0.0;
      example3Ref[0][2] = 1.0;
      example3Ref[1][0] = 2.0;
      example3Ref[1][1] = 1.0;
      example3Ref[1][2] = -2.0;
      example3Ref[2][0] = 3.0;
      example3Ref[2][1] = -3.0;
      example3Ref[2][2] = -5.0;
      Matrix example3RefMatrix = new Matrix(example3Ref);
      
      System.out.println(example3RefMatrix.determinantToString());
      System.out.println("\n--------------------------------\n");
      //*/
      
      ///*
      System.out.println("\n----------------TEST CASE 21----------------\n");
      double[][] exampleDet4Array = new double[4][4];
      exampleDet4Array[0][0] = 1;
      exampleDet4Array[0][1] = 3;
      exampleDet4Array[0][2] = 0;
      exampleDet4Array[0][3] = 3;
      exampleDet4Array[1][0] = 1;
      exampleDet4Array[1][1] = 1;
      exampleDet4Array[1][2] = 1;
      exampleDet4Array[1][3] = 1;
      exampleDet4Array[2][0] = 0;
      exampleDet4Array[2][1] = 4;
      exampleDet4Array[2][2] = 2;
      exampleDet4Array[2][3] = 8;
      exampleDet4Array[3][0] = 2;
      exampleDet4Array[3][1] = 0;
      exampleDet4Array[3][2] = 3;
      exampleDet4Array[3][3] = 1;
      Matrix exampleDet4Matrix = new Matrix(exampleDet4Array);
      
      System.out.println(exampleDet4Matrix.determinantToString());
      System.out.println("\n--------------------------------\n");
      //*/
      
      ///* IT WORKS!!! No zeros - it works!!!
      System.out.println("\n----------------TEST CASE 22----------------\n");
      double[][] exampleDet5Array = new double[5][5];
      exampleDet5Array[0][0] = 5;
      exampleDet5Array[0][1] = -1;
      exampleDet5Array[0][2] = 2;
      exampleDet5Array[0][3] = -3;
      exampleDet5Array[0][4] = -2;
      exampleDet5Array[1][0] = -7;
      exampleDet5Array[1][1] = 1;
      exampleDet5Array[1][2] = 1;
      exampleDet5Array[1][3] = 4;
      exampleDet5Array[1][4] = -3;
      exampleDet5Array[2][0] = -1;
      exampleDet5Array[2][1] = -2;
      exampleDet5Array[2][2] = -5;
      exampleDet5Array[2][3] = 3;
      exampleDet5Array[2][4] = 3;
      exampleDet5Array[3][0] = 3;
      exampleDet5Array[3][1] = 1;
      exampleDet5Array[3][2] = -1;
      exampleDet5Array[3][3] = 1;
      exampleDet5Array[3][4] = -4;
      exampleDet5Array[4][0] = 6;
      exampleDet5Array[4][1] = 2;
      exampleDet5Array[4][2] = -2;
      exampleDet5Array[4][3] = -2;
      exampleDet5Array[4][4] = 4;
      Matrix exampleDet5Matrix = new Matrix(exampleDet5Array);
      
      System.out.println(exampleDet5Matrix.determinantToString());
      System.out.println("\n--------------------------------\n");
      //*/
      
      ///*
      System.out.println("\n----------------TEST CASE 23----------------\n");
      double[][] exampleDet5Array2 = new double[5][5];
      exampleDet5Array2[0][0] = 2;
      exampleDet5Array2[0][1] = 5;
      exampleDet5Array2[0][2] = -7;
      exampleDet5Array2[0][3] = -2;
      exampleDet5Array2[0][4] = -1;
      exampleDet5Array2[1][0] = -1;
      exampleDet5Array2[1][1] = 0;
      exampleDet5Array2[1][2] = 4;
      exampleDet5Array2[1][3] = 0;
      exampleDet5Array2[1][4] = 3;
      exampleDet5Array2[2][0] = 1;
      exampleDet5Array2[2][1] = 6;
      exampleDet5Array2[2][2] = 0;
      exampleDet5Array2[2][3] = -3;
      exampleDet5Array2[2][4] = -2;
      exampleDet5Array2[3][0] = -2;
      exampleDet5Array2[3][1] = 4;
      exampleDet5Array2[3][2] = 7;
      exampleDet5Array2[3][3] = 8;
      exampleDet5Array2[3][4] = -1;
      exampleDet5Array2[4][0] = -3;
      exampleDet5Array2[4][1] = -1;
      exampleDet5Array2[4][2] = -1;
      exampleDet5Array2[4][3] = 1;
      exampleDet5Array2[4][4] = 4;
      Matrix exampleDet5Matrix2 = new Matrix(exampleDet5Array2);
      
      System.out.println(exampleDet5Matrix2.determinantToString());
      System.out.println("\n--------------------------------\n");
      //*/
      
      ///* even this works! 10x10 with no zeros in about ~3 seconds
      System.out.println("\n----------------TEST CASE 24----------------\n");
      double[][] exampleDet10Array = new double[10][10];
      exampleDet10Array[0][0] = -4;
      exampleDet10Array[0][1] = 3;
      exampleDet10Array[0][2] = -1;
      exampleDet10Array[0][3] = 9;
      exampleDet10Array[0][4] = -8;
      exampleDet10Array[0][5] = 12;
      exampleDet10Array[0][6] = 6;
      exampleDet10Array[0][7] = -2;
      exampleDet10Array[0][8] = -7;
      exampleDet10Array[0][9] = 1;
      exampleDet10Array[1][0] = 5;
      exampleDet10Array[1][1] = 2;
      exampleDet10Array[1][2] = -4;
      exampleDet10Array[1][3] = 3;
      exampleDet10Array[1][4] = 7;
      exampleDet10Array[1][5] = 1;
      exampleDet10Array[1][6] = 1;
      exampleDet10Array[1][7] = -1;
      exampleDet10Array[1][8] = 4;
      exampleDet10Array[1][9] = 9;
      exampleDet10Array[2][0] = 1;
      exampleDet10Array[2][1] = 3;
      exampleDet10Array[2][2] = -2;
      exampleDet10Array[2][3] = -7;
      exampleDet10Array[2][4] = 7;
      exampleDet10Array[2][5] = 6;
      exampleDet10Array[2][6] = 8;
      exampleDet10Array[2][7] = 3;
      exampleDet10Array[2][8] = -3;
      exampleDet10Array[2][9] = -3;
      exampleDet10Array[3][0] = 4;
      exampleDet10Array[3][1] = 1;
      exampleDet10Array[3][2] = -2;
      exampleDet10Array[3][3] = -1;
      exampleDet10Array[3][4] = 2;
      exampleDet10Array[3][5] = 3;
      exampleDet10Array[3][6] = -3;
      exampleDet10Array[3][7] = 5;
      exampleDet10Array[3][8] = 5;
      exampleDet10Array[3][9] = -5;
      exampleDet10Array[4][0] = -7;
      exampleDet10Array[4][1] = 6;
      exampleDet10Array[4][2] = 3;
      exampleDet10Array[4][3] = -3;
      exampleDet10Array[4][4] = 4;
      exampleDet10Array[4][5] = 1;
      exampleDet10Array[4][6] = -1;
      exampleDet10Array[4][7] = -1;
      exampleDet10Array[4][8] = -2;
      exampleDet10Array[4][9] = -8;
      exampleDet10Array[5][0] = -8;
      exampleDet10Array[5][1] = 1;
      exampleDet10Array[5][2] = 4;
      exampleDet10Array[5][3] = -3;
      exampleDet10Array[5][4] = 4;
      exampleDet10Array[5][5] = 5;
      exampleDet10Array[5][6] = -7;
      exampleDet10Array[5][7] = 1;
      exampleDet10Array[5][8] = -1;
      exampleDet10Array[5][9] = -2;
      exampleDet10Array[6][0] = 3;
      exampleDet10Array[6][1] = -2;
      exampleDet10Array[6][2] = 1;
      exampleDet10Array[6][3] = 4;
      exampleDet10Array[6][4] = 6;
      exampleDet10Array[6][5] = -5;
      exampleDet10Array[6][6] = -1;
      exampleDet10Array[6][7] = -3;
      exampleDet10Array[6][8] = 2;
      exampleDet10Array[6][9] = 8;
      exampleDet10Array[7][0] = 2;
      exampleDet10Array[7][1] = -1;
      exampleDet10Array[7][2] = -3;
      exampleDet10Array[7][3] = -4;
      exampleDet10Array[7][4] = -2;
      exampleDet10Array[7][5] = -7;
      exampleDet10Array[7][6] = 8;
      exampleDet10Array[7][7] = 6;
      exampleDet10Array[7][8] = -3;
      exampleDet10Array[7][9] = 4;
      exampleDet10Array[8][0] = 2;
      exampleDet10Array[8][1] = 7;
      exampleDet10Array[8][2] = 3;
      exampleDet10Array[8][3] = 8;
      exampleDet10Array[8][4] = -5;
      exampleDet10Array[8][5] = -1;
      exampleDet10Array[8][6] = -4;
      exampleDet10Array[8][7] = -1;
      exampleDet10Array[8][8] = 1;
      exampleDet10Array[8][9] = -3;
      exampleDet10Array[9][0] = 1;
      exampleDet10Array[9][1] = 6;
      exampleDet10Array[9][2] = -6;
      exampleDet10Array[9][3] = -5;
      exampleDet10Array[9][4] = -4;
      exampleDet10Array[9][5] = 1;
      exampleDet10Array[9][6] = -1;
      exampleDet10Array[9][7] = 3;
      exampleDet10Array[9][8] = -2;
      exampleDet10Array[9][9] = 8;
      Matrix exampleDet10Matrix = new Matrix(exampleDet10Array);
      
      System.out.println(exampleDet10Matrix.determinantToString());
      System.out.println("\n--------------------------------\n");
      //*/
      
      ///*
      System.out.println("\n----------------TEST CASE 25----------------\n");
      double[][] anotherExampleRef = new double[3][3];
      anotherExampleRef[0][0] = 1;
      anotherExampleRef[0][1] = 2;
      anotherExampleRef[0][2] = 3;
      anotherExampleRef[1][0] = 4;
      anotherExampleRef[1][1] = 8;
      anotherExampleRef[1][2] = 9;
      anotherExampleRef[2][0] = 3;
      anotherExampleRef[2][1] = 6;
      anotherExampleRef[2][2] = 7;
      Matrix anotherExampleRefMatrix = new Matrix(anotherExampleRef);
      Matrix anotherEchelonForm = anotherExampleRefMatrix.ref();
      
      //List<Matrix> publicRefSteps = anotherExampleRefMatrix.getRefSteps();
      System.out.println(anotherEchelonForm.refToString());
      System.out.println("\n--------------------------------\n");
      //*/
      
      ///*
      System.out.println("\n----------------TEST CASE 26----------------\n");
      double[][] extraExampleRREF = new double[3][3];
      extraExampleRREF[0][0] = 1;
      extraExampleRREF[0][1] = 2;
      extraExampleRREF[0][2] = 3;
      extraExampleRREF[1][0] = 4;
      extraExampleRREF[1][1] = 8;
      extraExampleRREF[1][2] = 12;
      extraExampleRREF[2][0] = 3;
      extraExampleRREF[2][1] = 5;
      extraExampleRREF[2][2] = 8;
      Matrix extraExampleRREFMatrix = new Matrix(extraExampleRREF);
      //Matrix echelonForm = exampleRREFMatrix.ref();
      //System.out.println(echelonForm.getCurrentRow());
      Matrix extraReducedEchelonForm = extraExampleRREFMatrix.rref();
      
      //List<Matrix> publicRREFSteps = exampleRREFMatrix.getRREFSteps();
      System.out.println(extraReducedEchelonForm.refToString());
      
      System.out.println(extraReducedEchelonForm.generalSolutionZero());
      System.out.println("\n--------------------------------\n");
      //*/
      
      ///*
      System.out.println("\n----------------TEST CASE 27----------------\n");
      double[][] exampleUnevenRREF = new double[4][5];
      exampleUnevenRREF[0][0] = 8;
      exampleUnevenRREF[0][1] = 11;
      exampleUnevenRREF[0][2] = -6;
      exampleUnevenRREF[0][3] = -7;
      exampleUnevenRREF[0][4] = 13;
      exampleUnevenRREF[1][0] = -7;
      exampleUnevenRREF[1][1] = -8;
      exampleUnevenRREF[1][2] = 5;
      exampleUnevenRREF[1][3] = 6;
      exampleUnevenRREF[1][4] = -9;
      exampleUnevenRREF[2][0] = 11;
      exampleUnevenRREF[2][1] = 7;
      exampleUnevenRREF[2][2] = -7;
      exampleUnevenRREF[2][3] = -9;
      exampleUnevenRREF[2][4] = -6;
      exampleUnevenRREF[3][0] = -3;
      exampleUnevenRREF[3][1] = 4;
      exampleUnevenRREF[3][2] = 1;
      exampleUnevenRREF[3][3] = 8;
      exampleUnevenRREF[3][4] = 7;
      Matrix exampleUnevenRREFMatrix = new Matrix(exampleUnevenRREF);
      //Matrix echelonForm = exampleRREFMatrix.ref();
      //System.out.println(echelonForm.getCurrentRow());
      //System.out.println(echelonForm.getNumPivots());
      Matrix reducedEchelonFormUneven = exampleUnevenRREFMatrix.rref();
      
      //List<Matrix> publicRREFSteps = exampleRREFMatrix.getRREFSteps();
      System.out.println(reducedEchelonFormUneven.refToString());
      System.out.println("\n--------------------------------\n");
      //*/
      
      ///*
      System.out.println("\n----------------TEST CASE 28----------------\n");
      double[][] exampleRREF = new double[7][7];
      exampleRREF[0][0] = -1;
      exampleRREF[0][1] = 2;
      exampleRREF[0][2] = 4;
      exampleRREF[0][3] = -3;
      exampleRREF[0][4] = 5;
      exampleRREF[0][5] = 0;
      exampleRREF[0][6] = 6;
      exampleRREF[1][0] = 6;
      exampleRREF[1][1] = 7;
      exampleRREF[1][2] = -8;
      exampleRREF[1][3] = 1;
      exampleRREF[1][4] = -1;
      exampleRREF[1][5] = -1;
      exampleRREF[1][6] = 2;
      exampleRREF[2][0] = -3;
      exampleRREF[2][1] = 0;
      exampleRREF[2][2] = 1;
      exampleRREF[2][3] = 1;
      exampleRREF[2][4] = -4;
      exampleRREF[2][5] = -2;
      exampleRREF[2][6] = 3;
      exampleRREF[3][0] = 0;
      exampleRREF[3][1] = 0;
      exampleRREF[3][2] = 4;
      exampleRREF[3][3] = 6;
      exampleRREF[3][4] = -7;
      exampleRREF[3][5] = 1;
      exampleRREF[3][6] = 4;
      exampleRREF[4][0] = 2;
      exampleRREF[4][1] = -2;
      exampleRREF[4][2] = -1;
      exampleRREF[4][3] = -3;
      exampleRREF[4][4] = 3;
      exampleRREF[4][5] = 8;
      exampleRREF[4][6] = -4;
      exampleRREF[5][0] = 6;
      exampleRREF[5][1] = 0;
      exampleRREF[5][2] = 3;
      exampleRREF[5][3] = 1;
      exampleRREF[5][4] = 0;
      exampleRREF[5][5] = -2;
      exampleRREF[5][6] = -5;
      exampleRREF[6][0] = 4;
      exampleRREF[6][1] = -5;
      exampleRREF[6][2] = -1;
      exampleRREF[6][3] = 0;
      exampleRREF[6][4] = -1;
      exampleRREF[6][5] = 0;
      exampleRREF[6][6] = 3;
      Matrix exampleRREFMatrix7 = new Matrix(exampleRREF);
      //Matrix echelonForm = exampleRREFMatrix.ref();
      //System.out.println(echelonForm.getCurrentRow());
      //System.out.println(echelonForm.getNumPivots());
      Matrix reducedEchelonForm7 = exampleRREFMatrix7.rref();
      
      //List<Matrix> publicRREFSteps = exampleRREFMatrix.getRREFSteps();
      System.out.println(reducedEchelonForm7.refToString());
      System.out.println("\n--------------------------------\n");
      //*/
      
      ///*
      System.out.println("\n----------------TEST CASE 29----------------\n");
      double[][] exampleRREF10 = new double[10][10];
      exampleRREF10[0][0] = -4;
      exampleRREF10[0][1] = 3;
      exampleRREF10[0][2] = -1;
      exampleRREF10[0][3] = 9;
      exampleRREF10[0][4] = -8;
      exampleRREF10[0][5] = 12;
      exampleRREF10[0][6] = 6;
      exampleRREF10[0][7] = -2;
      exampleRREF10[0][8] = -7;
      exampleRREF10[0][9] = 1;
      exampleRREF10[1][0] = 5;
      exampleRREF10[1][1] = 2;
      exampleRREF10[1][2] = -4;
      exampleRREF10[1][3] = 3;
      exampleRREF10[1][4] = 7;
      exampleRREF10[1][5] = 1;
      exampleRREF10[1][6] = 1;
      exampleRREF10[1][7] = -1;
      exampleRREF10[1][8] = 4;
      exampleRREF10[1][9] = 9;
      exampleRREF10[2][0] = 1;
      exampleRREF10[2][1] = 3;
      exampleRREF10[2][2] = -2;
      exampleRREF10[2][3] = -7;
      exampleRREF10[2][4] = 7;
      exampleRREF10[2][5] = 6;
      exampleRREF10[2][6] = 8;
      exampleRREF10[2][7] = 3;
      exampleRREF10[2][8] = -3;
      exampleRREF10[2][9] = -3;
      exampleRREF10[3][0] = 4;
      exampleRREF10[3][1] = 1;
      exampleRREF10[3][2] = -2;
      exampleRREF10[3][3] = -1;
      exampleRREF10[3][4] = 2;
      exampleRREF10[3][5] = 3;
      exampleRREF10[3][6] = -3;
      exampleRREF10[3][7] = 5;
      exampleRREF10[3][8] = 5;
      exampleRREF10[3][9] = -5;
      exampleRREF10[4][0] = -7;
      exampleRREF10[4][1] = 6;
      exampleRREF10[4][2] = 3;
      exampleRREF10[4][3] = -3;
      exampleRREF10[4][4] = 4;
      exampleRREF10[4][5] = 1;
      exampleRREF10[4][6] = -1;
      exampleRREF10[4][7] = -1;
      exampleRREF10[4][8] = -2;
      exampleRREF10[4][9] = -8;
      exampleRREF10[5][0] = -8;
      exampleRREF10[5][1] = 1;
      exampleRREF10[5][2] = 4;
      exampleRREF10[5][3] = -3;
      exampleRREF10[5][4] = 4;
      exampleRREF10[5][5] = 5;
      exampleRREF10[5][6] = -7;
      exampleRREF10[5][7] = 1;
      exampleRREF10[5][8] = -1;
      exampleRREF10[5][9] = -2;
      exampleRREF10[6][0] = 3;
      exampleRREF10[6][1] = -2;
      exampleRREF10[6][2] = 1;
      exampleRREF10[6][3] = 4;
      exampleRREF10[6][4] = 6;
      exampleRREF10[6][5] = -5;
      exampleRREF10[6][6] = -1;
      exampleRREF10[6][7] = -3;
      exampleRREF10[6][8] = 2;
      exampleRREF10[6][9] = 8;
      exampleRREF10[7][0] = 2;
      exampleRREF10[7][1] = -1;
      exampleRREF10[7][2] = -3;
      exampleRREF10[7][3] = -4;
      exampleRREF10[7][4] = -2;
      exampleRREF10[7][5] = -7;
      exampleRREF10[7][6] = 8;
      exampleRREF10[7][7] = 6;
      exampleRREF10[7][8] = -3;
      exampleRREF10[7][9] = 4;
      exampleRREF10[8][0] = 2;
      exampleRREF10[8][1] = 7;
      exampleRREF10[8][2] = 3;
      exampleRREF10[8][3] = 8;
      exampleRREF10[8][4] = -5;
      exampleRREF10[8][5] = -1;
      exampleRREF10[8][6] = -4;
      exampleRREF10[8][7] = -1;
      exampleRREF10[8][8] = 1;
      exampleRREF10[8][9] = -3;
      exampleRREF10[9][0] = 1;
      exampleRREF10[9][1] = 6;
      exampleRREF10[9][2] = -6;
      exampleRREF10[9][3] = -5;
      exampleRREF10[9][4] = -4;
      exampleRREF10[9][5] = 1;
      exampleRREF10[9][6] = -1;
      exampleRREF10[9][7] = 3;
      exampleRREF10[9][8] = -2;
      exampleRREF10[9][9] = 8;
      Matrix exampleRREFMatrix10 = new Matrix(exampleRREF10);
      Matrix echelonForm10 = exampleRREFMatrix10.ref();
      //System.out.println(echelonForm.getCurrentRow());
      //Matrix reducedEchelonForm = exampleRREFMatrix.rref();
      
      //List<Matrix> publicRefSteps = exampleRREFMatrix.getRefSteps();
      System.out.println(echelonForm10.refToString());
      //List<Matrix> publicRREFSteps = exampleRREFMatrix.getRREFSteps();
      //System.out.println(rRefToString(publicRREFSteps));
      System.out.println("\n--------------------------------\n");
      //*/
      
      ///*
      System.out.println("\n----------------TEST CASE 30----------------\n");
      double[][] originalGJArray = new double[2][2];
      originalGJArray[0][0] = 4;
      originalGJArray[0][1] = 3;
      originalGJArray[1][0] = 7;
      originalGJArray[1][1] = 6;
      Matrix originalGJMatrix = new Matrix(originalGJArray);
      
      Matrix inverseGJMatrix = originalGJMatrix.inverseGaussJordan();
      System.out.println(inverseGJMatrix.inverseToString());
      System.out.println("\n--------------------------------\n");
      //*/
      
      ///*
      System.out.println("\n----------------TEST CASE 31----------------\n");
      double[][] originalGJArray2 = new double[4][4];
      originalGJArray2[0][0] = -1;
      originalGJArray2[0][1] = 4;
      originalGJArray2[0][2] = 6;
      originalGJArray2[0][3] = 8;
      originalGJArray2[1][0] = 2;
      originalGJArray2[1][1] = 7;
      originalGJArray2[1][2] = 3;
      originalGJArray2[1][3] = -1;
      originalGJArray2[2][0] = 9;
      originalGJArray2[2][1] = 5;
      originalGJArray2[2][2] = -3;
      originalGJArray2[2][3] = 3;
      originalGJArray2[3][0] = 1;
      originalGJArray2[3][1] = 8;
      originalGJArray2[3][2] = 17;
      originalGJArray2[3][3] = 2;
      Matrix originalGJMatrix2 = new Matrix(originalGJArray2);
      
      Matrix inverseGJMatrix2 = originalGJMatrix2.inverseGaussJordan();
      System.out.println(inverseGJMatrix2.inverseToString());
      System.out.println("\n--------------------------------\n");
      //*/
      
      ///*
      System.out.println("\n----------------TEST CASE 32----------------\n");
      double[][] exampleRREFMore = new double[3][3];
      exampleRREFMore[0][0] = 2;
      exampleRREFMore[0][1] = -1;
      exampleRREFMore[0][2] = 6;
      exampleRREFMore[1][0] = 2;
      exampleRREFMore[1][1] = -1;
      exampleRREFMore[1][2] = 6;
      exampleRREFMore[2][0] = 2;
      exampleRREFMore[2][1] = -1;
      exampleRREFMore[2][2] = 6;
      Matrix exampleRREFMatrixMore = new Matrix(exampleRREFMore);
      //Matrix echelonForm = exampleRREFMatrix.ref();
      //System.out.println(echelonForm.getCurrentRow());
      Matrix reducedEchelonFormMore = exampleRREFMatrixMore.rref();
      
      //List<Matrix> publicRREFSteps = exampleRREFMatrix.getRREFSteps();
      System.out.println(reducedEchelonFormMore.refToString());
      
      System.out.println(reducedEchelonFormMore.generalSolutionZero());
      System.out.println("\n--------------------------------\n");
      //*/
      
      ///*
      System.out.println("\n----------------TEST CASE 33----------------\n");
      double[][] exampleRREFUnevenAgain = new double[3][5];
      exampleRREFUnevenAgain[0][0] = 1;
      exampleRREFUnevenAgain[0][1] = 6;
      exampleRREFUnevenAgain[0][2] = 2;
      exampleRREFUnevenAgain[0][3] = -5;
      exampleRREFUnevenAgain[0][4] = -2;
      exampleRREFUnevenAgain[1][0] = 0;
      exampleRREFUnevenAgain[1][1] = 0;
      exampleRREFUnevenAgain[1][2] = 2;
      exampleRREFUnevenAgain[1][3] = -8;
      exampleRREFUnevenAgain[1][4] = -1;
      exampleRREFUnevenAgain[2][0] = 0;
      exampleRREFUnevenAgain[2][1] = 0;
      exampleRREFUnevenAgain[2][2] = 0;
      exampleRREFUnevenAgain[2][3] = 0;
      exampleRREFUnevenAgain[2][4] = 1;
      Matrix exampleRREFMatrixUnevenAgain = new Matrix(exampleRREFUnevenAgain);
      //Matrix echelonForm = exampleRREFMatrix.ref();
      //System.out.println(echelonForm.getCurrentRow());
      Matrix reducedEchelonFormUnevenAgain = exampleRREFMatrixUnevenAgain.rref();
      
      //List<Matrix> publicRREFSteps = exampleRREFMatrix.getRREFSteps();
      System.out.println(reducedEchelonFormUnevenAgain.refToString());
      
      System.out.println(reducedEchelonFormUnevenAgain.generalSolutionZero());
      System.out.println("\n--------------------------------\n");
      //*/
      
      ///*
      System.out.println("\n----------------TEST CASE 34----------------\n");
      double[][] eigenFourArrayExtra = new double[4][4];
      eigenFourArrayExtra[0][0] = 5;
      eigenFourArrayExtra[0][1] = 0;
      eigenFourArrayExtra[0][2] = 0;
      eigenFourArrayExtra[0][3] = 0;
      eigenFourArrayExtra[1][0] = 0;
      eigenFourArrayExtra[1][1] = 5;
      eigenFourArrayExtra[1][2] = 0;
      eigenFourArrayExtra[1][3] = 0;
      eigenFourArrayExtra[2][0] = 1;
      eigenFourArrayExtra[2][1] = 4;
      eigenFourArrayExtra[2][2] = -3;
      eigenFourArrayExtra[2][3] = 0;
      eigenFourArrayExtra[3][0] = -1;
      eigenFourArrayExtra[3][1] = -2;
      eigenFourArrayExtra[3][2] = 0;
      eigenFourArrayExtra[3][3] = -3;
      Matrix eigenFourMatrixExtra = new Matrix(eigenFourArrayExtra);
      
      eigenFourMatrixExtra.findEigenvalues();
      System.out.println(eigenFourMatrixExtra.eigenvaluesToString());
      System.out.println(eigenFourMatrixExtra.eigenvectorsToString());
      System.out.println("\n--------------------------------\n");
      //*/
      
      ///*
      System.out.println("\n----------------TEST CASE 35----------------\n");
      double[][] eigenFourArray = new double[4][4];
      eigenFourArray[0][0] = 0;
      eigenFourArray[0][1] = 0;
      eigenFourArray[0][2] = 0;
      eigenFourArray[0][3] = 0;
      eigenFourArray[1][0] = 0;
      eigenFourArray[1][1] = 0;
      eigenFourArray[1][2] = 0;
      eigenFourArray[1][3] = 0;
      eigenFourArray[2][0] = 1;
      eigenFourArray[2][1] = 4;
      eigenFourArray[2][2] = -8;
      eigenFourArray[2][3] = 0;
      eigenFourArray[3][0] = -1;
      eigenFourArray[3][1] = -2;
      eigenFourArray[3][2] = 0;
      eigenFourArray[3][3] = -8;
      Matrix eigenFourMatrix = new Matrix(eigenFourArray);
      Matrix testRREFEigen4 = eigenFourMatrix.rref();
      
      //List<Matrix> testList = testRREF.getRefSteps();
      System.out.println(testRREFEigen4.rRefToString());
      System.out.println("\n--------------------------------\n");
      //*/
      
      ///*
      System.out.println("\n----------------TEST CASE 36----------------\n");
      double[][] exampleArrayAgain = new double[4][4];
      exampleArrayAgain[0][0] = -2;
      exampleArrayAgain[0][1] = 3;
      exampleArrayAgain[0][2] = -1;
      exampleArrayAgain[0][3] = 6;
      exampleArrayAgain[1][0] = 4;
      exampleArrayAgain[1][1] = 1;
      exampleArrayAgain[1][2] = -1;
      exampleArrayAgain[1][3] = 0;
      exampleArrayAgain[2][0] = 1;
      exampleArrayAgain[2][1] = 6;
      exampleArrayAgain[2][2] = 7;
      exampleArrayAgain[2][3] = 8;
      exampleArrayAgain[3][0] = -3;
      exampleArrayAgain[3][1] = -3;
      exampleArrayAgain[3][2] = -4;
      exampleArrayAgain[3][3] = -1;
      Matrix exampleMatrixAgain = new Matrix(exampleArrayAgain);
      //Matrix testRef = exampleMatrix.ref();
      //System.out.println(testRef.refToString());
      
      Matrix testRREFAgain = exampleMatrixAgain.rref();
      System.out.println(testRREFAgain.rRefToString());
      System.out.println("\n--------------------------------\n");
      //*/
      
      ///*
      System.out.println("\n----------------TEST CASE 37----------------\n");
      double[][] eigenThreeArray = new double[3][3];
      eigenThreeArray[0][0] = 4;
      eigenThreeArray[0][1] = -1;
      eigenThreeArray[0][2] = 6;
      eigenThreeArray[1][0] = 2;
      eigenThreeArray[1][1] = 1;
      eigenThreeArray[1][2] = 6;
      eigenThreeArray[2][0] = 2;
      eigenThreeArray[2][1] = -1;
      eigenThreeArray[2][2] = 8;
      Matrix eigenThreeMatrix = new Matrix(eigenThreeArray);
      
      eigenThreeMatrix.findEigenvalues();
      System.out.println(eigenThreeMatrix.eigenvaluesToString());
      System.out.println(eigenThreeMatrix.eigenvectorsToString());
      System.out.println("\n--------------------------------\n");
      //*/
      
      ///*
      System.out.println("\n----------------TEST CASE 38----------------\n");
      double[][] exampleArrayIGJ4 = new double[4][4];
      exampleArrayIGJ4[0][0] = -2;
      exampleArrayIGJ4[0][1] = 3;
      exampleArrayIGJ4[0][2] = -1;
      exampleArrayIGJ4[0][3] = 6;
      exampleArrayIGJ4[1][0] = 4;
      exampleArrayIGJ4[1][1] = 1;
      exampleArrayIGJ4[1][2] = -1;
      exampleArrayIGJ4[1][3] = 0;
      exampleArrayIGJ4[2][0] = 1;
      exampleArrayIGJ4[2][1] = 6;
      exampleArrayIGJ4[2][2] = 7;
      exampleArrayIGJ4[2][3] = 8;
      exampleArrayIGJ4[3][0] = -3;
      exampleArrayIGJ4[3][1] = -3;
      exampleArrayIGJ4[3][2] = -4;
      exampleArrayIGJ4[3][3] = -1;
      Matrix exampleMatrixIGJ4 = new Matrix(exampleArrayIGJ4);
      
      Matrix inverseGJMatrixIGJ4 = exampleMatrixIGJ4.inverseGaussJordan();
      System.out.println(inverseGJMatrixIGJ4.inverseToString());
      System.out.println("\n--------------------------------\n");
      //*/
      
      ///*
      System.out.println("\n----------------TEST CASE 39----------------\n");
      double[][] exampleRREFIGJ7 = new double[7][7];
      exampleRREFIGJ7[0][0] = -1;
      exampleRREFIGJ7[0][1] = 2;
      exampleRREFIGJ7[0][2] = 4;
      exampleRREFIGJ7[0][3] = -3;
      exampleRREFIGJ7[0][4] = 5;
      exampleRREFIGJ7[0][5] = 0;
      exampleRREFIGJ7[0][6] = 6;
      exampleRREFIGJ7[1][0] = 6;
      exampleRREFIGJ7[1][1] = 7;
      exampleRREFIGJ7[1][2] = -8;
      exampleRREFIGJ7[1][3] = 1;
      exampleRREFIGJ7[1][4] = -1;
      exampleRREFIGJ7[1][5] = -1;
      exampleRREFIGJ7[1][6] = 2;
      exampleRREFIGJ7[2][0] = -3;
      exampleRREFIGJ7[2][1] = 0;
      exampleRREFIGJ7[2][2] = 1;
      exampleRREFIGJ7[2][3] = 1;
      exampleRREFIGJ7[2][4] = -4;
      exampleRREFIGJ7[2][5] = -2;
      exampleRREFIGJ7[2][6] = 3;
      exampleRREFIGJ7[3][0] = 0;
      exampleRREFIGJ7[3][1] = 0;
      exampleRREFIGJ7[3][2] = 4;
      exampleRREFIGJ7[3][3] = 6;
      exampleRREFIGJ7[3][4] = -7;
      exampleRREFIGJ7[3][5] = 1;
      exampleRREFIGJ7[3][6] = 4;
      exampleRREFIGJ7[4][0] = 2;
      exampleRREFIGJ7[4][1] = -2;
      exampleRREFIGJ7[4][2] = -1;
      exampleRREFIGJ7[4][3] = -3;
      exampleRREFIGJ7[4][4] = 3;
      exampleRREFIGJ7[4][5] = 8;
      exampleRREFIGJ7[4][6] = -4;
      exampleRREFIGJ7[5][0] = 6;
      exampleRREFIGJ7[5][1] = 0;
      exampleRREFIGJ7[5][2] = 3;
      exampleRREFIGJ7[5][3] = 1;
      exampleRREFIGJ7[5][4] = 0;
      exampleRREFIGJ7[5][5] = -2;
      exampleRREFIGJ7[5][6] = -5;
      exampleRREFIGJ7[6][0] = 4;
      exampleRREFIGJ7[6][1] = -5;
      exampleRREFIGJ7[6][2] = -1;
      exampleRREFIGJ7[6][3] = 0;
      exampleRREFIGJ7[6][4] = -1;
      exampleRREFIGJ7[6][5] = 0;
      exampleRREFIGJ7[6][6] = 3;
      Matrix exampleMatrixIGJ7 = new Matrix(exampleRREFIGJ7);
      
      Matrix inverseGJMatrixIGJ7 = exampleMatrixIGJ7.inverseGaussJordan();
      System.out.println(inverseGJMatrixIGJ7.inverseToString());
      System.out.println("\n--------------------------------\n");
      //*/
      
      ///*
      System.out.println("\n----------------TEST CASE 40----------------\n");
      double[][] exampleRREFIGJ10 = new double[10][10];
      exampleRREFIGJ10[0][0] = -4;
      exampleRREFIGJ10[0][1] = 3;
      exampleRREFIGJ10[0][2] = -1;
      exampleRREFIGJ10[0][3] = 9;
      exampleRREFIGJ10[0][4] = -8;
      exampleRREFIGJ10[0][5] = 12;
      exampleRREFIGJ10[0][6] = 6;
      exampleRREFIGJ10[0][7] = -2;
      exampleRREFIGJ10[0][8] = -7;
      exampleRREFIGJ10[0][9] = 1;
      exampleRREFIGJ10[1][0] = 5;
      exampleRREFIGJ10[1][1] = 2;
      exampleRREFIGJ10[1][2] = -4;
      exampleRREFIGJ10[1][3] = 3;
      exampleRREFIGJ10[1][4] = 7;
      exampleRREFIGJ10[1][5] = 1;
      exampleRREFIGJ10[1][6] = 1;
      exampleRREFIGJ10[1][7] = -1;
      exampleRREFIGJ10[1][8] = 4;
      exampleRREFIGJ10[1][9] = 9;
      exampleRREFIGJ10[2][0] = 1;
      exampleRREFIGJ10[2][1] = 3;
      exampleRREFIGJ10[2][2] = -2;
      exampleRREFIGJ10[2][3] = -7;
      exampleRREFIGJ10[2][4] = 7;
      exampleRREFIGJ10[2][5] = 6;
      exampleRREFIGJ10[2][6] = 8;
      exampleRREFIGJ10[2][7] = 3;
      exampleRREFIGJ10[2][8] = -3;
      exampleRREFIGJ10[2][9] = -3;
      exampleRREFIGJ10[3][0] = 4;
      exampleRREFIGJ10[3][1] = 1;
      exampleRREFIGJ10[3][2] = -2;
      exampleRREFIGJ10[3][3] = -1;
      exampleRREFIGJ10[3][4] = 2;
      exampleRREFIGJ10[3][5] = 3;
      exampleRREFIGJ10[3][6] = -3;
      exampleRREFIGJ10[3][7] = 5;
      exampleRREFIGJ10[3][8] = 5;
      exampleRREFIGJ10[3][9] = -5;
      exampleRREFIGJ10[4][0] = -7;
      exampleRREFIGJ10[4][1] = 6;
      exampleRREFIGJ10[4][2] = 3;
      exampleRREFIGJ10[4][3] = -3;
      exampleRREFIGJ10[4][4] = 4;
      exampleRREFIGJ10[4][5] = 1;
      exampleRREFIGJ10[4][6] = -1;
      exampleRREFIGJ10[4][7] = -1;
      exampleRREFIGJ10[4][8] = -2;
      exampleRREFIGJ10[4][9] = -8;
      exampleRREFIGJ10[5][0] = -8;
      exampleRREFIGJ10[5][1] = 1;
      exampleRREFIGJ10[5][2] = 4;
      exampleRREFIGJ10[5][3] = -3;
      exampleRREFIGJ10[5][4] = 4;
      exampleRREFIGJ10[5][5] = 5;
      exampleRREFIGJ10[5][6] = -7;
      exampleRREFIGJ10[5][7] = 1;
      exampleRREFIGJ10[5][8] = -1;
      exampleRREFIGJ10[5][9] = -2;
      exampleRREFIGJ10[6][0] = 3;
      exampleRREFIGJ10[6][1] = -2;
      exampleRREFIGJ10[6][2] = 1;
      exampleRREFIGJ10[6][3] = 4;
      exampleRREFIGJ10[6][4] = 6;
      exampleRREFIGJ10[6][5] = -5;
      exampleRREFIGJ10[6][6] = -1;
      exampleRREFIGJ10[6][7] = -3;
      exampleRREFIGJ10[6][8] = 2;
      exampleRREFIGJ10[6][9] = 8;
      exampleRREFIGJ10[7][0] = 2;
      exampleRREFIGJ10[7][1] = -1;
      exampleRREFIGJ10[7][2] = -3;
      exampleRREFIGJ10[7][3] = -4;
      exampleRREFIGJ10[7][4] = -2;
      exampleRREFIGJ10[7][5] = -7;
      exampleRREFIGJ10[7][6] = 8;
      exampleRREFIGJ10[7][7] = 6;
      exampleRREFIGJ10[7][8] = -3;
      exampleRREFIGJ10[7][9] = 4;
      exampleRREFIGJ10[8][0] = 2;
      exampleRREFIGJ10[8][1] = 7;
      exampleRREFIGJ10[8][2] = 3;
      exampleRREFIGJ10[8][3] = 8;
      exampleRREFIGJ10[8][4] = -5;
      exampleRREFIGJ10[8][5] = -1;
      exampleRREFIGJ10[8][6] = -4;
      exampleRREFIGJ10[8][7] = -1;
      exampleRREFIGJ10[8][8] = 1;
      exampleRREFIGJ10[8][9] = -3;
      exampleRREFIGJ10[9][0] = 1;
      exampleRREFIGJ10[9][1] = 6;
      exampleRREFIGJ10[9][2] = -6;
      exampleRREFIGJ10[9][3] = -5;
      exampleRREFIGJ10[9][4] = -4;
      exampleRREFIGJ10[9][5] = 1;
      exampleRREFIGJ10[9][6] = -1;
      exampleRREFIGJ10[9][7] = 3;
      exampleRREFIGJ10[9][8] = -2;
      exampleRREFIGJ10[9][9] = 8;
      
      Matrix exampleMatrixIGJ10 = new Matrix(exampleRREFIGJ10);
      
      Matrix inverseGJMatrixIGJ10 = exampleMatrixIGJ10.inverseGaussJordan();
      System.out.println(inverseGJMatrixIGJ10.inverseToString());
      System.out.println("\n--------------------------------\n");
      //*/
      
      
      /*
      // FIX THIS. An exception should pop up instead of an NPE.
      System.out.println("\n----------------TEST CASE 41----------------\n");
      double[][] exampleArrayNoInverse = new double[3][3];
      exampleArrayNoInverse[0][0] = 1;
      exampleArrayNoInverse[0][1] = -2;
      exampleArrayNoInverse[0][2] = 1;
      exampleArrayNoInverse[1][0] = 4;
      exampleArrayNoInverse[1][1] = -7;
      exampleArrayNoInverse[1][2] = 3;
      exampleArrayNoInverse[2][0] = -2;
      exampleArrayNoInverse[2][1] = 6;
      exampleArrayNoInverse[2][2] = -4;
      Matrix exampleMatrixNoInverse = new Matrix(exampleArrayNoInverse);
      
      Matrix inverseGJMatrixNoInverse = exampleMatrixNoInverse.inverseGaussJordan();
      System.out.println(inverseGJMatrixNoInverse.inverseToString());
      System.out.println("\n--------------------------------\n");
      //*/
      
      ///*
      System.out.println("\n----------------TEST CASE 42----------------\n");
      double[][] exampleArrayCramer = new double[4][4];
      exampleArrayCramer[0][0] = -2;
      exampleArrayCramer[0][1] = 3;
      exampleArrayCramer[0][2] = -1;
      exampleArrayCramer[0][3] = 6;
      exampleArrayCramer[1][0] = 4;
      exampleArrayCramer[1][1] = 1;
      exampleArrayCramer[1][2] = -1;
      exampleArrayCramer[1][3] = 0;
      exampleArrayCramer[2][0] = 1;
      exampleArrayCramer[2][1] = 6;
      exampleArrayCramer[2][2] = 7;
      exampleArrayCramer[2][3] = 8;
      exampleArrayCramer[3][0] = -3;
      exampleArrayCramer[3][1] = -3;
      exampleArrayCramer[3][2] = -4;
      exampleArrayCramer[3][3] = -1;
      Matrix exampleMatrixCramer = new Matrix(exampleArrayCramer);
      
      Matrix inverseMatrixCramer = exampleMatrixCramer.inverseCramer();
      System.out.println(inverseMatrixCramer.inverseToString());
      System.out.println("\n--------------------------------\n");
      //*/
      
      ///*
      // testing for determinant of matrix of dimension >= 4x4
      System.out.println("\n----------------TEST CASE 43----------------\n");
      double[][] exampleDetFive = new double[5][5];
      exampleDetFive[0][0] = 4;
      exampleDetFive[0][1] = 0;
      exampleDetFive[0][2] = -3;
      exampleDetFive[0][3] = 6;
      exampleDetFive[0][4] = 1;
      exampleDetFive[1][0] = -1;
      exampleDetFive[1][1] = 4;
      exampleDetFive[1][2] = -2;
      exampleDetFive[1][3] = 8;
      exampleDetFive[1][4] = 0;
      exampleDetFive[2][0] = 3;
      exampleDetFive[2][1] = 1;
      exampleDetFive[2][2] = 4;
      exampleDetFive[2][3] = -7;
      exampleDetFive[2][4] = 9;
      exampleDetFive[3][0] = 0;
      exampleDetFive[3][1] = 8;
      exampleDetFive[3][2] = 0;
      exampleDetFive[3][3] = -1;
      exampleDetFive[3][4] = -3;
      exampleDetFive[4][0] = -5;
      exampleDetFive[4][1] = -2;
      exampleDetFive[4][2] = -1;
      exampleDetFive[4][3] = 3;
      exampleDetFive[4][4] = 4;
           
      Matrix detFive = new Matrix(exampleDetFive);
      System.out.println(detFive.determinantToString());
      System.out.println("\n--------------------------------\n");
      //*/
      
      /*
      System.out.println("\n----------------TEST CASE 44----------------\n");
      double[][] exampleDetTenPlus = new double[10][10];
      //double[][] exampleDetFourPlus = new double[15][15];
      exampleDetTenPlus[0][0] = 4;
      exampleDetTenPlus[0][1] = 3;
      exampleDetTenPlus[0][2] = -6;
      exampleDetTenPlus[0][3] = 7;
      exampleDetTenPlus[0][4] = 1;
      exampleDetTenPlus[0][5] = 0;
      exampleDetTenPlus[0][6] = 8;
      exampleDetTenPlus[0][7] = 9;
      exampleDetTenPlus[0][8] = 1;
      exampleDetTenPlus[0][9] = -4;
      //exampleDetTenPlus[0][10] = -3;
      //exampleDetTenPlus[0][11] = 0;
      //exampleDetTenPlus[0][12] = 2;
      //exampleDetTenPlus[0][13] = -8;
      //exampleDetTenPlus[0][14] = 1;
      exampleDetTenPlus[1][0] = 7;
      exampleDetTenPlus[1][1] = 6;
      exampleDetTenPlus[1][2] = -5;
      exampleDetTenPlus[1][3] = -1;
      exampleDetTenPlus[1][4] = 0;
      exampleDetTenPlus[1][5] = -3;
      exampleDetTenPlus[1][6] = -2;
      exampleDetTenPlus[1][7] = 0;
      exampleDetTenPlus[1][8] = 8;
      exampleDetTenPlus[1][9] = 1;
      //exampleDetTenPlus[1][10] = 4;
      //exampleDetTenPlus[1][11] = -6;
      //exampleDetTenPlus[1][12] = 0;
      //exampleDetTenPlus[1][13] = 3;
      //exampleDetTenPlus[1][14] = 8;
      exampleDetTenPlus[2][0] = 9;
      exampleDetTenPlus[2][1] = 1;
      exampleDetTenPlus[2][2] = -1;
      exampleDetTenPlus[2][3] = -3;
      exampleDetTenPlus[2][4] = -8;
      exampleDetTenPlus[2][5] = -4;
      exampleDetTenPlus[2][6] = 1;
      exampleDetTenPlus[2][7] = 4;
      exampleDetTenPlus[2][8] = -5;
      exampleDetTenPlus[2][9] = 7;
      //exampleDetTenPlus[2][10] = 0;
      //exampleDetTenPlus[2][11] = 2;
      //exampleDetTenPlus[2][12] = -2;
      //exampleDetTenPlus[2][13] = -3;
      //exampleDetTenPlus[2][14] = -9;
      exampleDetTenPlus[3][0] = -2;
      exampleDetTenPlus[3][1] = -8;
      exampleDetTenPlus[3][2] = 3;
      exampleDetTenPlus[3][3] = 4;
      exampleDetTenPlus[3][4] = -1;
      exampleDetTenPlus[3][5] = 0;
      exampleDetTenPlus[3][6] = 6;
      exampleDetTenPlus[3][7] = -7;
      exampleDetTenPlus[3][8] = -3;
      exampleDetTenPlus[3][9] = 4;
      //exampleDetTenPlus[3][10] = -1;
      //exampleDetTenPlus[3][11] = -2;
      //exampleDetTenPlus[3][12] = 3;
      //exampleDetTenPlus[3][13] = 6;
      //exampleDetTenPlus[3][14] = -4;
      exampleDetTenPlus[4][0] = 1;
      exampleDetTenPlus[4][1] = 8;
      exampleDetTenPlus[4][2] = 8;
      exampleDetTenPlus[4][3] = 0;
      exampleDetTenPlus[4][4] = 8;
      exampleDetTenPlus[4][5] = -3;
      exampleDetTenPlus[4][6] = 4;
      exampleDetTenPlus[4][7] = 0;
      exampleDetTenPlus[4][8] = -6;
      exampleDetTenPlus[4][9] = 7;
      //exampleDetTenPlus[4][10] = 9;
      //exampleDetTenPlus[4][11] = -5;
      //exampleDetTenPlus[4][12] = -3;
      //exampleDetTenPlus[4][13] = 4;
      //exampleDetTenPlus[4][14] = -7;
      exampleDetTenPlus[5][0] = 2;
      exampleDetTenPlus[5][1] = 2;
      exampleDetTenPlus[5][2] = 3;
      exampleDetTenPlus[5][3] = -1;
      exampleDetTenPlus[5][4] = -7;
      exampleDetTenPlus[5][5] = -9;
      exampleDetTenPlus[5][6] = -3;
      exampleDetTenPlus[5][7] = -4;
      exampleDetTenPlus[5][8] = 5;
      exampleDetTenPlus[5][9] = -1;
      //exampleDetTenPlus[5][10] = 1;
      //exampleDetTenPlus[5][11] = -5;
      //exampleDetTenPlus[5][12] = -6;
      //exampleDetTenPlus[5][13] = 6;
      //exampleDetTenPlus[5][14] = 10;
      exampleDetTenPlus[6][0] = 0;
      exampleDetTenPlus[6][1] = 0;
      exampleDetTenPlus[6][2] = -3;
      exampleDetTenPlus[6][3] = 4;
      exampleDetTenPlus[6][4] = 1;
      exampleDetTenPlus[6][5] = 5;
      exampleDetTenPlus[6][6] = 6;
      exampleDetTenPlus[6][7] = -8;
      exampleDetTenPlus[6][8] = -3;
      exampleDetTenPlus[6][9] = -7;
      //exampleDetTenPlus[6][10] = 9;
      //exampleDetTenPlus[6][11] = 7;
      //exampleDetTenPlus[6][12] = 1;
      //exampleDetTenPlus[6][13] = 4;
      //exampleDetTenPlus[6][14] = -8;
      exampleDetTenPlus[7][0] = 2;
      exampleDetTenPlus[7][1] = 1;
      exampleDetTenPlus[7][2] = -2;
      exampleDetTenPlus[7][3] = 3;
      exampleDetTenPlus[7][4] = -6;
      exampleDetTenPlus[7][5] = 5;
      exampleDetTenPlus[7][6] = -4;
      exampleDetTenPlus[7][7] = -2;
      exampleDetTenPlus[7][8] = -1;
      exampleDetTenPlus[7][9] = -1;
      //exampleDetTenPlus[7][10] = 2;
      //exampleDetTenPlus[7][11] = 3;
      //exampleDetTenPlus[7][12] = -9;
      //exampleDetTenPlus[7][13] = 7;
      //exampleDetTenPlus[7][14] = -6;
      exampleDetTenPlus[8][0] = 3;
      exampleDetTenPlus[8][1] = 0;
      exampleDetTenPlus[8][2] = 3;
      exampleDetTenPlus[8][3] = -3;
      exampleDetTenPlus[8][4] = 8;
      exampleDetTenPlus[8][5] = 1;
      exampleDetTenPlus[8][6] = -1;
      exampleDetTenPlus[8][7] = 7;
      exampleDetTenPlus[8][8] = 8;
      exampleDetTenPlus[8][9] = 3;
      //exampleDetTenPlus[8][10] = -2;
      //exampleDetTenPlus[8][11] = 4;
      //exampleDetTenPlus[8][12] = 2;
      //exampleDetTenPlus[8][13] = -6;
      //exampleDetTenPlus[8][14] = 5;
      exampleDetTenPlus[9][0] = -5;
      exampleDetTenPlus[9][1] = 4;
      exampleDetTenPlus[9][2] = -2;
      exampleDetTenPlus[9][3] = -1;
      exampleDetTenPlus[9][4] = 6;
      exampleDetTenPlus[9][5] = -7;
      exampleDetTenPlus[9][6] = 3;
      exampleDetTenPlus[9][7] = 6;
      exampleDetTenPlus[9][8] = -1;
      exampleDetTenPlus[9][9] = 0;
      //exampleDetTenPlus[9][10] = 3;
      //exampleDetTenPlus[9][11] = 4;
      //exampleDetTenPlus[9][12] = -4;
      //exampleDetTenPlus[9][13] = -1;
      //exampleDetTenPlus[9][14] = -2;
      /*
      exampleDetTenPlus[10][0] = -1;
      exampleDetTenPlus[10][1] = 3;
      exampleDetTenPlus[10][2] = 7;
      exampleDetTenPlus[10][3] = 6;
      exampleDetTenPlus[10][4] = 4;
      exampleDetTenPlus[10][5] = -2;
      exampleDetTenPlus[10][6] = -5;
      exampleDetTenPlus[10][7] = -1;
      exampleDetTenPlus[10][8] = 0;
      exampleDetTenPlus[10][9] = -8;
      exampleDetTenPlus[10][10] = -9;
      exampleDetTenPlus[10][11] = -3;
      exampleDetTenPlus[10][12] = 3;
      exampleDetTenPlus[10][13] = -2;
      exampleDetTenPlus[10][14] = 4;
      exampleDetTenPlus[11][0] = -8;
      exampleDetTenPlus[11][1] = 3;
      exampleDetTenPlus[11][2] = -4;
      exampleDetTenPlus[11][3] = 2;
      exampleDetTenPlus[11][4] = -1;
      exampleDetTenPlus[11][5] = 0;
      exampleDetTenPlus[11][6] = -2;
      exampleDetTenPlus[11][7] = -3;
      exampleDetTenPlus[11][8] = 1;
      exampleDetTenPlus[11][9] = 7;
      exampleDetTenPlus[11][10] = 4;
      exampleDetTenPlus[11][11] = -4;
      exampleDetTenPlus[11][12] = 6;
      exampleDetTenPlus[11][13] = 0;
      exampleDetTenPlus[11][14] = -3;
      exampleDetTenPlus[12][0] = -3;
      exampleDetTenPlus[12][1] = -7;
      exampleDetTenPlus[12][2] = 8;
      exampleDetTenPlus[12][3] = 1;
      exampleDetTenPlus[12][4] = 0;
      exampleDetTenPlus[12][5] = -1;
      exampleDetTenPlus[12][6] = -9;
      exampleDetTenPlus[12][7] = -1;
      exampleDetTenPlus[12][8] = 5;
      exampleDetTenPlus[12][9] = 3;
      exampleDetTenPlus[12][10] = 2;
      exampleDetTenPlus[12][11] = -4;
      exampleDetTenPlus[12][12] = -1;
      exampleDetTenPlus[12][13] = 5;
      exampleDetTenPlus[12][14] = 0;
      exampleDetTenPlus[13][0] = -6;
      exampleDetTenPlus[13][1] = 1;
      exampleDetTenPlus[13][2] = 3;
      exampleDetTenPlus[13][3] = 9;
      exampleDetTenPlus[13][4] = -2;
      exampleDetTenPlus[13][5] = -2;
      exampleDetTenPlus[13][6] = 2;
      exampleDetTenPlus[13][7] = -2;
      exampleDetTenPlus[13][8] = 4;
      exampleDetTenPlus[13][9] = -7;
      exampleDetTenPlus[13][10] = 0;
      exampleDetTenPlus[13][11] = 3;
      exampleDetTenPlus[13][12] = -5;
      exampleDetTenPlus[13][13] = 4;
      exampleDetTenPlus[13][14] = 1;
      exampleDetTenPlus[14][0] = 1;
      exampleDetTenPlus[14][1] = 0;
      exampleDetTenPlus[14][2] = 8;
      exampleDetTenPlus[14][3] = -2;
      exampleDetTenPlus[14][4] = 7;
      exampleDetTenPlus[14][5] = 9;
      exampleDetTenPlus[14][6] = 0;
      exampleDetTenPlus[14][7] = -4;
      exampleDetTenPlus[14][8] = 5;
      exampleDetTenPlus[14][9] = -3;
      exampleDetTenPlus[14][10] = -3;
      exampleDetTenPlus[14][11] = 4;
      exampleDetTenPlus[14][12] = 3;
      exampleDetTenPlus[14][13] = -1;
      exampleDetTenPlus[14][14] = 8;
      */
      
      /*
      Matrix detTenPlus = new Matrix(exampleDetTenPlus);
      //Matrix testRef = detFourPlus.ref();
      //System.out.println(testRef.refToString());
      //System.out.println(testRef.isTriangular());
      //System.out.println(testRef.determinantToString());
      //Matrix testRREF = detFourPlus.rref();
      //System.out.println(testRREF.rRefToString());
      System.out.println(detTenPlus.determinantToString());
      System.out.println("\n--------------------------------\n");
      //*/
      
      ///*
      // M == N
      System.out.println("\n----------------TEST CASE 45----------------\n");
      double[][] firstMatrixArray = new double[3][3];
      double[][] secondMatrixArray = new double[3][3];
      firstMatrixArray[0][0] = 3;
      firstMatrixArray[0][1] = 12;
      firstMatrixArray[0][2] = 7;
      firstMatrixArray[1][0] = 6;
      firstMatrixArray[1][1] = 18;
      firstMatrixArray[1][2] = 4;
      firstMatrixArray[2][0] = 9;
      firstMatrixArray[2][1] = 3;
      firstMatrixArray[2][2] = 1;
      secondMatrixArray[0][0] = 2;
      secondMatrixArray[0][1] = 4;
      secondMatrixArray[0][2] = 6;
      secondMatrixArray[1][0] = 8;
      secondMatrixArray[1][1] = 10;
      secondMatrixArray[1][2] = 12;
      secondMatrixArray[2][0] = 1;
      secondMatrixArray[2][1] = 2;
      secondMatrixArray[2][2] = 3;
      Matrix firstMatrix = new Matrix(firstMatrixArray);
      Matrix secondMatrix = new Matrix(secondMatrixArray);
      
      Matrix newMatrix = firstMatrix.multiplyMatrices(secondMatrix);
      System.out.println(newMatrix.multiplyToString());
      System.out.println("\n--------------------------------\n");
      //*/
      
      ///*
      // M < N
      System.out.println("\n----------------TEST CASE 46----------------\n");
      double[][] firstMatrixArray2 = new double[1][2];
      double[][] secondMatrixArray2 = new double[2][3];
      firstMatrixArray2[0][0] = 3;
      firstMatrixArray2[0][1] = 12;
      secondMatrixArray2[0][0] = 2;
      secondMatrixArray2[0][1] = 4;
      secondMatrixArray2[0][2] = 6;
      secondMatrixArray2[1][0] = 8;
      secondMatrixArray2[1][1] = 10;
      secondMatrixArray2[1][2] = 12;
      Matrix firstMatrix2 = new Matrix(firstMatrixArray2);
      Matrix secondMatrix2 = new Matrix(secondMatrixArray2);
      
      Matrix newMatrix2 = firstMatrix2.multiplyMatrices(secondMatrix2);
      System.out.println(newMatrix2.multiplyToString());
      System.out.println("\n--------------------------------\n");
      //*/
      
      ///*
      // M > N
      System.out.println("\n----------------TEST CASE 47----------------\n");
      double[][] firstMatrixArray3 = new double[3][3];
      double[][] secondMatrixArray3 = new double[3][2];
      firstMatrixArray3[0][0] = 3;
      firstMatrixArray3[0][1] = 12;
      firstMatrixArray3[0][2] = 7;
      firstMatrixArray3[1][0] = 6;
      firstMatrixArray3[1][1] = 18;
      firstMatrixArray3[1][2] = 4;
      firstMatrixArray3[2][0] = 9;
      firstMatrixArray3[2][1] = 3;
      firstMatrixArray3[2][2] = -1;
      secondMatrixArray3[0][0] = 2;
      secondMatrixArray3[0][1] = 4;
      secondMatrixArray3[1][0] = 8;
      secondMatrixArray3[1][1] = 10;
      secondMatrixArray3[2][0] = 1;
      secondMatrixArray3[2][1] = 2;
      
      Matrix firstMatrix3 = new Matrix(firstMatrixArray3);
      Matrix secondMatrix3 = new Matrix(secondMatrixArray3);
      
      Matrix newMatrix3 = firstMatrix3.multiplyMatrices(secondMatrix3);
      System.out.println(newMatrix3.multiplyToString());
      System.out.println("\n--------------------------------\n");
      //*/
      
      /*
      double[][] matrixArray = new double[3][3];
      matrixArray[0][0] = 1;
      matrixArray[0][1] = 2;
      matrixArray[0][2] = 3;
      matrixArray[1][0] = 2;
      matrixArray[1][1] = 4;
      matrixArray[1][2] = 6;
      matrixArray[2][0] = -3;
      matrixArray[2][1] = 8;
      matrixArray[2][2] = 4;
      Matrix AXMatrix = new Matrix(matrixArray);
      
      double[] vectorArray = new double[3];
      vectorArray[0] = 12;
      vectorArray[1] = 12;
      vectorArray[2] = 1;
      Vector BVector = new Vector(vectorArray);
      */
   }
}
