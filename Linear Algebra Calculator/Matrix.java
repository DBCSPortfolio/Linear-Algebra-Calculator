/****************************************
 * Linear Algebra Calculator
 * 
 * Supports representation of complex numbers, matrix multiplication
 *    and addition, row reduction, transposing, methods for finding
 *    determinants, several methods for inverses, eigenvalues, eigenvectors,
 *    diagonalization, and operations related to orthogonality.
 * 
 * 
 * Matrix.java
 * 
 * API:
 *    public Matrix(int m, int n)
 *    public Matrix(double[][] matrixEntries)
 *    private Matrix copy()
 *    private void setNumberOfZeros(int number)
 *    private Matrix multiplyByScalar(Matrix matrix, double scalar)
 *    private Matrix add(Matrix matrixTwo)
 *    public Matrix multiplyMatrices(Matrix matrixTwo)
 *    public String multiplyToString()
 *    private boolean validateSameSize(Matrix matrixOne, Matrix matrixTwo)
 *    private boolean squareMatrix()
 *    private boolean equalMatrices(Matrix matrixOne, Matrix matrixTwo)
 *    public boolean isTriangular()
 *    private boolean isUpperTriangular()
 *    private boolean isLowerTriangular()
 *    public Matrix rref()
 *    public String rRefToString()
 *    private void findPivotColumn()
 *    private Matrix reduceUpward()
 *    private Matrix scale()
 *    public String generalSolutionZero()
 *    public Matrix ref()
 *    public String refToString()
 *    private void findPivotRow()
 *    private Matrix reduceDownward()
 *    private void interchange(int pivotRow, int otherRow)
 *    public Matrix refAXB(Vector bVector)
 *    private void findPivotRowAXB()
 *    private Matrix reduceDownwardAXB(Vector bVector)
 *    private Matrix interchangeAXB(int pivotRow, int otherRow)
 *    private Matrix transpose()
 *    public String transposeToString()
 *    public Matrix inverse2()
 *    public String inverseToString()
 *    public Matrix inverseGaussJordan()
 *    private Matrix identityMatrix(int dimension)
 *    private List<Double> cramersRule(Vector equalsVector)
 *    public String cramersRuleToString(Vector equalsVector)
 *    private double[][] getCofactorMatrix(int mostZerosLine, boolean isRow,
 *                            int currentOther)
 *    private double findCofactor(int row, int col)
 *    private Matrix findAdjugate()
 *    public Matrix inverseCramer()
 *    public double determinant2()
 *    public double determinant3()
 *    public double determinantFourPlus()
 *    private void setZeroPositions(Matrix cofactorMatrix, int mostZerosLine,
 *                      boolean isRow)
 *    public String determinantToString()
 *    public void findEigenvalues()
 *    public String eigenvaluesToString()
 *    private double calculateDiscriminant()
 *    private boolean findTypeEigenvaluesQuadratic(double discriminant)
 *    private void findRealEigenvaluesQuadratic()
 *    private void findComplexEigenvaluesQuadratic()
 *    private double calculateHCubic()
 *    private static double calculateHCubic(double a, double b, double c,
 *                               double d)
 *    private static boolean findTypeEigenvaluesCubic(double h)
 *    private void findRealEigenvaluesCubic()
 *    private static List<Double> findRealEigenvaluesCubic(double a, double b,
 *                                     double c, double d, double h)
 *    private void findComplexEigenvaluesCubic()
 *    private static List<ComplexNumber> findComplexEigenvaluesCubic(double a,
 *                                           double b, double c, double d,
 *                                           double h)
 *    private boolean findTypeEigenvaluesQuartic()
 *    private void findRealEigenvaluesQuartic()
 *    private void findComplexEigenvaluesQuartic()
 *    private String eigenDeterminant2()
 *    private String eigenDeterminant3()
 *    private String eigenDeterminant4Plus()
 *    public List<Vector> findEigenvectors()
 *    public String eigenvectorsToString()
 *    public List<Matrix> diagonalize(List<Double> eigenvalues,
 *                            List<Vector> eigenvectors)
 *    private void diagonalizeBothReal()
 *    private void diagonalizeValuesComplex()
 *    private void diagonalizeVectorsComplex()
 *    private void diagonalizeBothComplex()
 *    public String diagonalizeToString()
 *    
 *    
 * Dependencies:
 *    -ComplexNumber.java
 *    -Vector.java
 *    -MatrixZeros.java
 * 
 ****************************************/

import java.lang.Math;
import java.util.List;
import java.util.ArrayList;
import java.util.Collections;
import java.lang.Exception;
import java.text.DecimalFormat;
import java.math.RoundingMode;

public class Matrix extends Exception
{
   // final constants
   public static final String LAMBDA_STRING = "\u03BB"; //"lambda";
   public static final String MATRIX_SPACING = "             ";
      // made it 12 instead of 9 for large numbers (hundreds, thousands, etc.)
   
   // private data tied to Matrix objects
   private int                   m;
   private int                   n;
   private double[][]            matrix;
   private int                   numberOfZeros;
   
   private double[]              coefficientArray;
   private List<MatrixZeros>     zeroPositionMap;
   private boolean               zerosSet;
   private List<Double>          realEigenvalueList;
   private List<ComplexNumber>   complexEigenvalueList;
   private List<Vector>          realEigenvectors;
   private List<ComplexVector>   complexEigenvectors;
   private List<Double>          internalCubicCoeffs;
   private List<Double>          internalCubicRootsReal;
   private List<ComplexNumber>   internalCubicRootsComplex;
   
   // private data for use in multiple methods
   private int                   pivotColumn;
   private int                   currentRow;
   private int                   numPivots;
   private int                   lastRowWithPivot;
   private boolean               refDoneFlag;
   private List<Matrix>          refStepsList;
   private List<Matrix>          rRefStepsList;
   private List<String>          refExplanationList;
   private List<String>          rRefExplanationList;
   private Matrix                matrixP;
   private Matrix                matrixD;
   private Matrix                matrixPInverse;
   private static DecimalFormat  df = new DecimalFormat("#.#####");
   
   
   // ************************** NOTE!!! ***************************
   // EVENTUALLY SEPARATE THIS WHOLE THING INTO AN MVC STRUCTURE.
   // **************************************************************
   
   // BUILD IN EXCEPTIONS FOR EVERY METHOD!
   
   // include method for testing roots to make sure that they actually make
   //       an equation equal to zero!
   
   
   /*
   // constructor with dimensions
   public Matrix(int m, int n)
   {
      this.m = m;
      this.n = n;
      matrix = new double[m][n];
      numberOfZeros = -1;
      if (m == n)
      {
         realEigenvalueList = new ArrayList<Double>();
         complexEigenvalueList = new ArrayList<ComplexNumber>();
         coefficientArray = new double[m + 1];
         
         if (m == 2)
         {
            //
         }
         else if (m == 3)
         {
            
         }
      }
   }
   */
   
   // constructor with actual matrix entries
   public Matrix(double[][] matrixEntries)
   {
      m = matrixEntries.length;
      n = matrixEntries[0].length;
      this.matrix = new double[m][n];
      for (int i = 0; i < m; i++)
      {
         for (int j = 0; j < n; j++)
         {
            this.matrix[i][j] = matrixEntries[i][j];
         }
      }
      numberOfZeros = -1; // this is terrible but it should work
      
      pivotColumn = 0;
      currentRow = 0;
      numPivots = 0;
      lastRowWithPivot = -1;
      refDoneFlag = false;
      
      refStepsList = new ArrayList<Matrix>();
      rRefStepsList = new ArrayList<Matrix>();
      
      refExplanationList = new ArrayList<String>();
      rRefExplanationList = new ArrayList<String>();
      
      if (m == n)
      {
         zeroPositionMap = new ArrayList<MatrixZeros>();
         zerosSet = false;
         realEigenvalueList = new ArrayList<Double>();
         complexEigenvalueList = new ArrayList<ComplexNumber>();
         realEigenvectors = new ArrayList<Vector>();
         complexEigenvectors = new ArrayList<ComplexVector>();
            // this doesn't necessarily need to be initiated?
         
         coefficientArray = new double[m + 1];
         
         if (m == 2)
         {
            coefficientArray[0] = 1;
            coefficientArray[1] = // -(a + d)
                  -(this.matrix[0][0] + this.matrix[1][1]);
            coefficientArray[2] = // (ad - bc)
                  ((this.matrix[0][0] * this.matrix[1][1])
                  - (this.matrix[0][1] * this.matrix[1][0]));
         }
         else if (m == 3)
         {
            coefficientArray[0] = -1;
            coefficientArray[1] = // a + e + i
                  (this.matrix[0][0] + this.matrix[1][1] +
                        this.matrix[2][2]);
            coefficientArray[2] = // cg + fh + bd - ae - ai - ei
                  (this.matrix[0][2] * this.matrix[2][0]) +
                  (this.matrix[1][2] * this.matrix[2][1]) +
                  (this.matrix[0][1] * this.matrix[1][0]) -
                  (this.matrix[0][0] * this.matrix[1][1]) -
                  (this.matrix[0][0] * this.matrix[2][2]) -
                  (this.matrix[1][1] * this.matrix[2][2]);
            coefficientArray[3] = // aei + bfg + cdh - afh - bdi - ceg
                  (this.matrix[0][0] * this.matrix[1][1] *
                        this.matrix[2][2]) +
                  (this.matrix[0][1] * this.matrix[1][2] *
                        this.matrix[2][0]) +
                  (this.matrix[0][2] * this.matrix[1][0] *
                        this.matrix[2][1]) -
                  (this.matrix[0][0] * this.matrix[1][2] *
                        this.matrix[2][1]) -
                  (this.matrix[0][1] * this.matrix[1][0] *
                        this.matrix[2][2]) -
                  (this.matrix[0][2] * this.matrix[1][1] *
                        this.matrix[2][0]);
         }
         else if (m == 4)
         {
            internalCubicCoeffs = new ArrayList<Double>();
            internalCubicRootsReal = new ArrayList<Double>();
            internalCubicRootsComplex = new ArrayList<ComplexNumber>();
            
            double a = this.matrix[0][0];
            double b = this.matrix[0][1];
            double c = this.matrix[0][2];
            double d = this.matrix[0][3];
            double e = this.matrix[1][0];
            double f = this.matrix[1][1];
            double g = this.matrix[1][2];
            double h = this.matrix[1][3];
            double i = this.matrix[2][0];
            double j = this.matrix[2][1];
            double k = this.matrix[2][2];
            double l = this.matrix[2][3];
            double m = this.matrix[3][0];
            double n = this.matrix[3][1];
            double o = this.matrix[3][2];
            double p = this.matrix[3][3];
            
            // a = 1
            // b = -a + -f + -k + -p
            // c = af + ak + ap + fk + fp + kp
            //          - be - ci - dm - gj - hn - lo
            // d = agj + ahn + alo + bek + bep + cfi + cip + dfm + dkm
            //          + flo + gjp + hkn - afk - afp - akp - bgi - bhm
            //          - cej - clm - den - dio - fkp - gln - hjo
            // e = afkp + agln + ahjo + belo + bgip + bhkm + cejp + cflm
            //          + chin + dekn + dfio + dgjm - aflo - agjp - ahkn
            //          - bekp - bglm - bhio - celn - cfip - chjm - dejo
            //          - dfkm - dgin
            
            coefficientArray[0] = 1;
            coefficientArray[1] = - (a + f + k + p);
            coefficientArray[2] =
                  (a * f) + (a * k) + (a * p) + (f * k) + (f * p) + (k * p)
                  - (b * e) - (c * i) - (d * m) - (g * j) - (h * n) - (l * o);
            coefficientArray[3] =
                  (a * g * j) + (a * h * n) + (a * l * o) + (b * e * k)
                  + (b * e * p) + (c * f * i) + (c * i * p) + (d * f * m)
                  + (d * k * m) + (f * l * o) + (g * j * p) + (h * k * n)
                  - (a * f * k) - (a * f * p) - (a * k * p) - (b * g * i)
                  - (b * h * m) - (c * e * j) - (c * l * m) - (d * e * n)
                  - (d * i * o) - (f * k * p) - (g * l * n) - (h * j * o);
            coefficientArray[4] =
                  (a * f * k * p) + (a * g * l * n) + (a * h * j * o)
                  + (b * e * l * o) + (b * g * i * p) + (b * h * k * m)
                  + (c * e * j * p) + (c * f * l * m) + (c * h * i * n)
                  + (d * e * k * n) + (d * f * i * o) + (d * g * j * m)
                  - (a * f * l * o) - (a * g * j * p) - (a * h * k * n)
                  - (b * e * k * p) - (b * g * l * m) - (b * h * i * o)
                  - (c * e * l * n) - (c * f * i * p) - (c * h * j * m)
                  - (d * e * j * o) - (d * f * k * m) - (d * g * i * n);
         }
      }
   }
   
   // returns a copy of the instance Matrix that calls this method
   private Matrix copy()
   {
      double[][] copyArray = new double[this.m][this.n];
      for (int i = 0; i < this.m; i++)
      {
         for (int j = 0; j < this.n; j++)
         {
            copyArray[i][j] = this.matrix[i][j];
         }
      }
      Matrix matrixCopy = new Matrix(copyArray);
      return matrixCopy;
   }
   
   private void setNumberOfZeros(int number)
   {
      this.numberOfZeros = number;
   }
   
   // multiplies a matrix by a scalar quantity, returning the result
   private Matrix multiplyByScalar(double scalar)
   {
      for (int i = 0; i < this.m; i++)
      {
         for (int j = 0; j < this.n; j++)
         {
            this.matrix[i][j] = this.matrix[i][j] * scalar;
         }
      }
      return this;
   }
   
   // adds two matrices together, returning the result
   private Matrix add(Matrix matrixTwo)
   {
      if (validateSameSize(this, matrixTwo))
      {
         // modify max for particular matrices
         for (int i = 0; i < this.m; i++)
         {
            for (int j = 0; j < this.n; j++)
            {
               this.matrix[i][j] = 
                     this.matrix[i][j] + matrixTwo.matrix[i][j];
            }
         }
      }
      return this;
   }
   
   public Matrix multiplyMatrices(Matrix matrixTwo) throws Exception
   {
      if (this.n != matrixTwo.m)
      {
         throw new Exception("Matrices cannot be multiplied.");
      }
      
      double[][] multipliedArray = new double[this.m][matrixTwo.n];
      
      for (int i = 0; i < this.m; i++)
      {
         for (int j = 0; j < matrixTwo.n; j++)
         {
            for (int k = 0; k < this.n; k++) // need to change upper bounds!!!
            {
               // check to see if this is even right.
               multipliedArray[i][j] = multipliedArray[i][j] +
                     (this.matrix[i][k] * matrixTwo.matrix[k][j]);
            }
         }
      }
      Matrix multResult = new Matrix(multipliedArray);
      return multResult;
   }
   
   public String multiplyToString()
   {  
      String multipliedString =
            "The result of the matrix multiplication is as follows: \n\n";
      for (int i = 0; i < this.m; i++)
      {
         for (int j = 0; j < this.n; j++)
         {
            multipliedString = multipliedString + df.format(this.matrix[i][j])
            + MATRIX_SPACING.substring(0, (MATRIX_SPACING.length() - 
                  (df.format(this.matrix[i][j]).length()))) + "   ";
         }
         multipliedString = multipliedString + "\n";
      }
      return multipliedString;
   }
   
   // makes sure that two matrices/vectors are same size and can be added
   private boolean validateSameSize(Matrix matrixOne, Matrix matrixTwo)
   {
      if (matrixOne.m == matrixTwo.m && matrixOne.n == matrixTwo.n)
      {
         return true;
      }
      return false;
   }
   
   // determines whether or not a matrix is square
   private boolean squareMatrix()
   {
      if (this.m == this.n)
      {
         return true;
      }
      return false;
   }
   
   // checks if two matrices/vectors are perfectly equal
   private boolean equalMatrices(Matrix matrixOne, Matrix matrixTwo)
   {
      if (!validateSameSize(matrixOne, matrixTwo))
      {
         return false;
      }
      else
      {
         for (int i = 0; i < matrixOne.n; i++)
         {
            for (int j = 0; j < matrixOne.m; j++)
            {
               if (!(matrixOne.matrix[i][j] == matrixTwo.matrix[i][j]))
               {
                  return false;
               }
            }
         }
      }
      return true;
   }
   
   // determines whether or not a matrix is triangular (upper or lower)
   public boolean isTriangular() throws Exception // PRIVATE!
   {
      if (!this.squareMatrix())
      {
         throw new Exception("Matrix must be square to be triangular.");
      }
      
      if (isUpperTriangular() || isLowerTriangular())
      {
         return true;
      }
      return false;
   }
   
   // determines whether or not a matrix is upper triangular
   private boolean isUpperTriangular() throws Exception
   {
      if (!this.squareMatrix())
      {
         throw new Exception("Matrix must be square to be triangular.");
      }
      
      boolean upperTriangular = true;
      
      for (int i = 0; i < this.m; i++)
      {
         for (int j = i + 1; j < this.n; j++)
         {
            if (this.matrix[j][i] != 0)
            {
               upperTriangular = false;
            }
         }
      }
      if (!upperTriangular)
      {
         return false;
      }
      return true;
   }
   
   // determines whether or not a matrix is lower triangular
   private boolean isLowerTriangular() throws Exception
   {
      if (!this.squareMatrix())
      {
         throw new Exception("Matrix must be square to be triangular.");
      }
      
      boolean lowerTriangular = true;
      
      for (int k = 0; k < this.m - 1; k++)
      {
         for (int l = k + 1; l < this.n; l++)
         {
            if (this.matrix[k][l] != 0)
            {
               lowerTriangular = false;
            }
         }
      }
      if (!lowerTriangular)
      {
         return false;
      }
      return true;
   }
   
   // determines the reduced row echelon form of the matrix in question
   public Matrix rref()
   {    
      Matrix refPart = this.ref();
      refPart.rRefStepsList = refPart.refStepsList;
      refPart.rRefExplanationList = refPart.refExplanationList;
      
      refPart.currentRow = refPart.numPivots - 1; // what if there's one pivot?
      
      while (refPart.currentRow > 0)
      {
         refPart.findPivotColumn();
         if (refPart.matrix[refPart.currentRow][refPart.pivotColumn] != 1)
         {
            refPart = refPart.scale();
         }
         refPart = refPart.reduceUpward();
      }
      
      refPart.findPivotColumn();
      if (refPart.matrix[refPart.currentRow][refPart.pivotColumn] != 1)
      {
         refPart = refPart.scale();
      }
      return refPart;
   }
   
   // returns the steps required to find a matrix's reduced row echelon form
   public String rRefToString()
   {
      String rRefOutput = "The steps leading to a *reduced* row echelon form "
            + "of the given matrix are: \n";
      
      // this runs in cubic time!
      for (Matrix rS : this.rRefStepsList)
      {
         rRefOutput = rRefOutput + "\n" +
               this.refExplanationList.get(this.refStepsList.indexOf(rS)) +
               "\n";
         for (int i = 0; i < rS.m; i++)
         {
            for (int j = 0; j < rS.n; j++)
            {
               rRefOutput = rRefOutput + df.format(rS.matrix[i][j])
                     + MATRIX_SPACING.substring(0, (MATRIX_SPACING.length() - 
                           (df.format(rS.matrix[i][j]).length()))) + "   ";
               
               /*
               if (rS.matrix[i][j] >= 0)
               {
                  rRefOutput = rRefOutput + " " +
                        df.format(rS.matrix[i][j]) + "   ";
               }
               else
               {
                  rRefOutput = rRefOutput +
                     df.format(rS.matrix[i][j]) + "   ";
               }
               */
            }
            rRefOutput = rRefOutput + "\n";
         }
      }
      rRefOutput = rRefOutput +
            "This is the reduced row echelon form of the original matrix.";
      return rRefOutput;
   }
   
   // finds the next column, in backwards order, containing a pivot (for RREF)
   private void findPivotColumn()
   {
      // finds pivot position for use when reducing upward
      for (int j = 0; j < this.n; j++)
      {
         if (this.matrix[this.currentRow][j] != 0)
         {
            this.pivotColumn = j;
            return;
         }
      }
   }
   
   // makes all rows above a pivot zero (for RREF)
   private Matrix reduceUpward()
   {
      for (int i = 0; i < this.currentRow; i++)
      {
         double currAddFactor = -(this.matrix[i][this.pivotColumn]);
         
         for (int j = 0; j < this.n; j++)
         {
            this.matrix[i][j] = this.matrix[i][j] +
                  (currAddFactor * this.matrix[this.currentRow][j]);
            if (this.matrix[i][j] == -0.0)
            {
               this.matrix[i][j] = 0.0;
            }
            else if (this.matrix[i][j] < 0.00001 &&
                  this.matrix[i][j] > -0.00001)
            {
               this.matrix[i][j] = 0.0;
            }
         }
      }
      // so pivot finding and reduction continues upward in matrix
      this.currentRow--;
      
      rRefExplanationList.add("Reducing the matrix upward from row " +
            Integer.toString(currentRow + 1) + " yields: ");
      rRefStepsList.add(this.copy());
      return this;
   }
   
   // multiplies a row of a matrix by a scalar quantity to make its pivot 1
   private Matrix scale()
   {
      double scaleBase = this.matrix[this.currentRow][this.pivotColumn];
      if (scaleBase != 1)
      {
         for (int j = this.pivotColumn; j < this.n; j++)
         {
            
            this.matrix[this.currentRow][j] =
                  this.matrix[this.currentRow][j] / scaleBase;

            if (this.matrix[this.currentRow][j] == -0.0)
            {
               this.matrix[this.currentRow][j] = 0.0;
            }
            else if (this.matrix[this.currentRow][j] < 0.00001 &&
                  this.matrix[this.currentRow][j] > -0.00001)
            {
               this.matrix[this.currentRow][j] = 0.0;
            }
         }
      }
      rRefExplanationList.add("Scaling row " + 
            Integer.toString(currentRow + 1) + " by a factor of " +
            df.format(scaleBase) + " yields: ");
      rRefStepsList.add(this.copy());
      return this;
   }
   
   // returns the general solution for a matrix
   public String generalSolutionZero()
   {
      String generalSolution = "The general solution for the original matrix "
            + "is as follows:\n\n";
      String equalsSomething;
      
      for (int i = 0; i < this.m; i++)
      {
         boolean pivotFound = false;
         equalsSomething = "";
         
         for (int j = 0; j < this.n; j++)
         {
            if (pivotFound)
            {
               if (this.matrix[i][j] == 0)
               {
                  continue;
               }
               else if (this.matrix[i][j] != 0 &&
                     equalsSomething.equals(""))
               {
                  if (this.matrix[i][j] == -1)
                  {
                     equalsSomething = "x" + Integer.toString(j + 1);
                  }
                  else if (this.matrix[i][j] == 1)
                  {
                     equalsSomething = "-x" + Integer.toString(j + 1);
                  }
                  else
                  {
                     if (this.matrix[i][j] % 1 == 0)
                     {
                        Double tempDouble = new Double(this.matrix[i][j]);
                        int tempInt = tempDouble.intValue();
                        equalsSomething = Integer.toString(-(tempInt)) +
                              "x" + Integer.toString(j + 1);
                     }
                     else
                     {
                        equalsSomething =
                              Double.toString(-(this.matrix[i][j])) +
                              "x" + Integer.toString(j + 1);
                     }
                  }
               }
               else if (this.matrix[i][j] != 0 &&
                     !equalsSomething.equals(""))
               {
                  if (this.matrix[i][j] == -1)
                  {
                     equalsSomething = equalsSomething + " + " + "x" +
                           Integer.toString(j + 1);
                  }
                  else if (this.matrix[i][j] == 1)
                  {
                     equalsSomething = equalsSomething + " + " + "-x" +
                           Integer.toString(j + 1);
                  }
                  else
                  {
                     if (this.matrix[i][j] % 1 == 0)
                     {
                        Double tempDouble = new Double(this.matrix[i][j]);
                        int tempInt = tempDouble.intValue();
                        equalsSomething = equalsSomething + " + " +
                              Integer.toString(-(tempInt)) + "x" +
                              Integer.toString(j + 1);
                     }
                     else
                     {
                        equalsSomething = equalsSomething + " + " + 
                              Double.toString(-(this.matrix[i][j])) +
                              "x" + Integer.toString(j + 1);
                     }
                  }
               }
            }
            else if (this.matrix[i][j] != 0)
            {
               generalSolution =
                     generalSolution + "x" + Integer.toString(j + 1) + " = ";
               pivotFound = true;
            }
         }
         if (equalsSomething.equals("") && pivotFound)
         {
            generalSolution = generalSolution + "0\n";
         }
         if (!equalsSomething.equals(""))
         {
            generalSolution = generalSolution + equalsSomething + "\n";
         }
      }
      
      return generalSolution;
   }
   
   // determines one reduced echelon form (of many possible) of a matrix
   public Matrix ref()
   {
      Matrix reductionStep = this.copy();
      reductionStep.refStepsList.add(reductionStep.copy());
      reductionStep.refExplanationList.add("Original matrix: ");
      
      while (!reductionStep.refDoneFlag)
      {
         reductionStep.findPivotRow();
         reductionStep = reductionStep.reduceDownward();
         
         if (reductionStep.currentRow == reductionStep.m)
         {
            reductionStep.refDoneFlag = true;
         }
      }
      return reductionStep;
   }
   
   // returns the steps required to find a row echelon form of a matrix
   public String refToString()
   {
      String refOutput = "The steps leading to a row echelon form of the "
            + "given matrix are: \n";
      df.setRoundingMode(RoundingMode.HALF_UP);
      
      // this runs in cubic time!
      for (Matrix rS : this.refStepsList)
      {
         refOutput = refOutput + "\n" +
               this.refExplanationList.get(this.refStepsList.indexOf(rS)) +
               "\n";
         for (int i = 0; i < rS.m; i++)
         {
            for (int j = 0; j < rS.n; j++)
            {
               /*
               refOutput = refOutput + df.format(rS.matrix[i][j])
               + MATRIX_SPACING.substring(0, (MATRIX_SPACING.length() - 
                     (df.format(rS.matrix[i][j]).length()))); //+ "   ";
               */
               
               
               if (rS.matrix[i][j] >= 0)
               {
                  refOutput = refOutput + " " +
                        df.format(rS.matrix[i][j]) + "   ";
               }
               else
               {
                  refOutput = refOutput +
                        df.format(rS.matrix[i][j]) + "   ";
               }
               
            }
            refOutput = refOutput + "\n";
         }
      }
      refOutput = refOutput +
            "This is a row echelon form of the original matrix.";
      return refOutput;
   }
   
   // finds the next row in the current matrix that contains a pivot point
   private void findPivotRow()
   {  
      // goes downward by row, not the typical across-by-column way
      for (int j = 0; j < this.n; j++) // shouldn't start at zero!!!!!!!!
      {
         for (int i = this.currentRow; i < this.m; i++)
         {
            if (this.matrix[i][j] != 0)
            {
               this.pivotColumn = j;
               this.numPivots++;
               this.lastRowWithPivot++;
               if (i > this.currentRow)
               {
                  this.interchange(i, this.currentRow);
               }
               return;
            }
         }
      }
   }
   
   // performs a row reduction operation on one row, based on another row
   private Matrix reduceDownward()
   {      
      // SHOULD ONLY ACTUALLY BE USED IF A ROW BELOW ALSO CONTAINS
      // NONZERO VALUE IN SAME COLUMN. WASTE OF TIME OTHERWISE.
      
      boolean reductionOccurred = false;
      
      for (int i = this.currentRow + 1; i < this.m; i++)
      {
         // only if matrixArray[i][pivotColumn] != 0
         // only reduces a row if the value in the pivot column is nonzero
         if (this.matrix[i][this.pivotColumn] == 0)
         {
            continue;
         }
         else
         {
            reductionOccurred = true;
            // changes with each following row
            double currAddFactor =
                  (-1) * ((this.matrix[i][this.pivotColumn]) / 
                        (this.matrix[this.currentRow][this.pivotColumn]));
            
            for (int j = this.pivotColumn; j < this.n; j++)
            {
               this.matrix[i][j] = this.matrix[i][j] +
                     (currAddFactor * this.matrix[this.currentRow][j]);
               if (this.matrix[i][j] == -0.0)
               {
                  this.matrix[i][j] = 0.0;
               }
               else if (this.matrix[i][j] < 0.00001 &&
                     this.matrix[i][j] > -0.00001)
               {
                  this.matrix[i][j] = 0.0;
               }
            }
         }
      }
      
      if (reductionOccurred)
      {
         this.refStepsList.add(this.copy());
      }
      
      refExplanationList.add("Performing row replacement based on row "
            + Integer.toString(this.currentRow + 1) + " yields: ");
      
      // so pivot finding and reduction continues downward in matrix
      this.currentRow++;
      
      return this;
   }
   
   // moves row with earlier pivot above other rows
   // pivot row first, less significant row second
   private Matrix interchange(int pivotRow, int otherRow)
   {
      double[] tempArray = new double[this.n];
      for (int i = 0; i < this.n; i++)
      {
         tempArray[i] = this.matrix[otherRow][i];
         this.matrix[otherRow][i] = this.matrix[pivotRow][i];
         this.matrix[pivotRow][i] = tempArray[i];
      }
      String explanationString =
            "Swapping rows " + Integer.toString(otherRow + 1)
            + " and " + Integer.toString(pivotRow + 1) + " yields: ";
      this.refExplanationList.add(explanationString);
      this.refStepsList.add(this.copy());
      return this;
   }
   
   // determines one reduced echelon form (of many possible) of a matrix
   public Matrix refAXB(Vector bVector)
   {
      Matrix reductionStep = this.copy();
      reductionStep.refStepsList.add(reductionStep.copy());
      reductionStep.refExplanationList.add("Original matrix: ");
      
      while (!reductionStep.refDoneFlag)
      {
         reductionStep.findPivotRow();
         reductionStep = reductionStep.reduceDownward();
         
         if (reductionStep.currentRow == reductionStep.m)
         {
            reductionStep.refDoneFlag = true;
         }
      }
      return reductionStep;
   }
   
   // finds the next row in the current matrix that contains a pivot point
   private void findPivotRowAXB()
   {  
      // goes downward by row, not the typical across-by-column way
      for (int j = 0; j < this.n; j++) // shouldn't start at zero!!!!!!!!
      {
         for (int i = this.currentRow; i < this.m; i++)
         {
            if (this.matrix[i][j] != 0)
            {
               this.pivotColumn = j;
               this.numPivots++;
               this.lastRowWithPivot++;
               if (i > this.currentRow)
               {
                  this.interchange(i, this.currentRow);
               }
               return;
            }
         }
      }
   }
   
   // performs a row reduction operation on one row, based on another row
   private Matrix reduceDownwardAXB(Vector bVector) throws Exception
   {      
      // SHOULD ONLY ACTUALLY BE USED IF A ROW BELOW ALSO CONTAINS
      // NONZERO VALUE IN SAME COLUMN. WASTE OF TIME OTHERWISE.
      
      boolean reductionOccurred = false;
      
      for (int i = this.currentRow + 1; i < this.m; i++)
      {
         boolean inconsistentRow = true;
         
         // only if matrixArray[i][pivotColumn] != 0
         // only reduces a row if the value in the pivot column is nonzero
         if (this.matrix[i][this.pivotColumn] == 0)
         {
            continue;
         }
         else
         {
            //inconsistentRow = false;
            reductionOccurred = true;
            // changes with each following row
            double currAddFactor =
                  (-1) * ((this.matrix[i][this.pivotColumn]) / 
                        (this.matrix[this.currentRow][this.pivotColumn]));
            
            for (int j = this.pivotColumn; j < this.n; j++)
            {
               this.matrix[i][j] = this.matrix[i][j] +
                     (currAddFactor * this.matrix[this.currentRow][j]);
               if (this.matrix[i][j] == -0.0)
               {
                  this.matrix[i][j] = 0.0;
               }
               else if (this.matrix[i][j] < 0.00001 &&
                     this.matrix[i][j] > -0.00001)
               {
                  this.matrix[i][j] = 0.0;
               }
               if (this.matrix[i][j] != 0 && inconsistentRow)
               {
                  inconsistentRow = false;
               }
            }
            
            // can't easily access vector entries... figure out another way
            bVector.setValue((bVector.getValue(i) +
                  (currAddFactor * bVector.getValue(this.currentRow))), i);
            // where is inconsistentRow ever set though?
            if (inconsistentRow && bVector.getValue(i) != 0)
            {
               throw new Exception("The system is inconsistent.");
            }
         }
      }
      
      if (reductionOccurred)
      {
         this.refStepsList.add(this.copy());
      }
      
      /*
      refExplanationList.add("Performing row replacement based on row "
            + Integer.toString(this.currentRow + 1) + " yields: ");
      */
      
      // so pivot finding and reduction continues downward in matrix
      this.currentRow++;
      
      return this;
   }
   
   // moves row with earlier pivot above other rows
   // pivot row first, less significant row second
   private Matrix interchangeAXB(int pivotRow, int otherRow)
   {
      double[] tempArray = new double[this.n];
      for (int i = 0; i < this.n; i++)
      {
         tempArray[i] = this.matrix[otherRow][i];
         this.matrix[otherRow][i] = this.matrix[pivotRow][i];
         this.matrix[pivotRow][i] = tempArray[i];
      }
      String explanationString =
            "Swapping rows " + Integer.toString(otherRow + 1)
            + " and " + Integer.toString(pivotRow + 1) + " yields: ";
      this.refExplanationList.add(explanationString);
      this.refStepsList.add(this.copy());
      return this;
   }
   
   
   
   // *** FULLY TESTED ***
   // determines the transpose of a passed matrix
   // should run in n^2 - n time for square matrices
   // n^2 - n array accesses, BUT n^2 comparisons
   private Matrix transpose()
   {  
      double[][] auxArray = new double[this.n][this.m];
      //Matrix auxMatrix = new Matrix(this.n, this.m);
      double tempValue;
      
      if (this.m != this.n)
      {
         for (int i = 0; i < this.n; i++)
         {
            for (int j = 0; j < this.m; j++)
            {
               if (i != j)
               {
                  auxArray[i][j] = this.matrix[j][i];
               }
               else
               {
                  auxArray[i][j] = this.matrix[i][j];
               }
            }
         }
         Matrix auxMatrix = new Matrix(auxArray);
         return auxMatrix;
      }
      else
      {
         // see isTriangular method for a potential improvement on running time
         for (int i = 0; i < this.n; i++)
         {
            for (int j = 0; j < this.m; j++)
            {
               if (i > j) // using the upper triangle above main diagonal
               {
                   // sets two places at once!
                   tempValue = this.matrix[i][j];
                   this.matrix[i][j] = this.matrix[j][i];
                   this.matrix[j][i] = tempValue;
               }
            }
         }
         return this;
      }
   }
   
   // produces String output for the transpose of a matrix
   public String transposeToString()
   {
      Matrix transposeMatrix = this.transpose();
      
      String transposeString = "The transpose of the given matrix is: ";
      df.setRoundingMode(RoundingMode.HALF_UP);
      
      
      transposeString = transposeString + "\n\n";
      for (int i = 0; i < transposeMatrix.m; i++)
      {
         for (int j = 0; j < transposeMatrix.n; j++)
         {
            if (transposeMatrix.matrix[i][j] >= 0)
            {
               transposeString = transposeString + " " +
                     df.format(transposeMatrix.matrix[i][j]) + "   ";
            }
            else
            {
               transposeString = transposeString +
                     df.format(transposeMatrix.matrix[i][j]) + "   ";
            }
         }
         transposeString = transposeString + "\n";
      }
      return transposeString;
   }
   
   // returns the inverse of a 2x2 matrix, if it is invertible
   public Matrix inverse2() throws Exception
   {
      if (!this.squareMatrix())
      {
         throw new Exception("A non-square matrix is not invertible.");
      }
      
      double determinant = this.determinant2();
      if (determinant == 0)
      {
         // fix this later to make easier for user, but works as placeholder
         throw new Exception("This matrix is not invertible.");
      }
      
      Matrix inverseMatrix = this;
      inverseMatrix.matrix[0][1] = -(inverseMatrix.matrix[0][1]);
      inverseMatrix.matrix[1][0] = -(inverseMatrix.matrix[1][0]);
      double temp = inverseMatrix.matrix[0][0];
      inverseMatrix.matrix[0][0] = inverseMatrix.matrix[1][1];
      inverseMatrix.matrix[1][1] = temp;
            
      inverseMatrix = inverseMatrix.multiplyByScalar(1 / determinant);
      return inverseMatrix;
   }
   
   // returns the inverse of the matrix in question, formatted as a string
   public String inverseToString() throws Exception
   {
      if (this.equals(null))
      {
         throw new Exception("Inverse of the given matrix does not exist.");
      }
      
      String inverseString = "The inverse of the given matrix is: ";
      df.setRoundingMode(RoundingMode.HALF_UP);
      
      inverseString = inverseString + "\n\n";
      for (int i = 0; i < this.m; i++)
      {
         for (int j = 0; j < this.n; j++)
         {
            if (this.matrix[i][j] >= 0)
            {
               inverseString = inverseString + " " +
                     df.format(this.matrix[i][j]) + "   ";
            }
            else
            {
               inverseString = inverseString +
                     df.format(this.matrix[i][j]) + "   ";
            }
         }
         inverseString = inverseString + "\n";
      }
      return inverseString;
   }
   
   // MAKE SURE TO INCLUDE CHECKS FOR WHETHER THE MATRIX CAN EVEN BE INVERTED.
   // RELY ON EVERY PIECE OF THE INVERTIBLE MATRIX THEOREM.
   
   /* Invertible Matrix Theorem:
    *    Let A be a square n x n matrix. Then the following statements are
    *    equivalent. That is, for a given A, the statements are either all true
    *    or all false.
    *    
    *    (a) A is an invertible matrix.
    *    (c) A has n pivot positions.
    *    (d) The equation Ax = 0 has only the trivial solution.
    *    (g) The equation Ax = b has at least one solution for each b in R^n.
    * 
    */
   public Matrix inverseGaussJordan() throws Exception
   {
      if (!this.squareMatrix())
      {
         throw new Exception("A non-square matrix is not invertible.");
      }
      
      int length = matrix.length;
      double[][] newMatrixArray = new double[length][length * 2];
      Matrix sizedIdentity = identityMatrix(length);
      
      for (int h = 0; h < length; h++)
      {
         for (int i = 0; i < length; i++)
         {
            // better way to do this?
            newMatrixArray[h][i] = this.matrix[h][i];
            newMatrixArray[h][i + length] =
                  sizedIdentity.matrix[h][i];
         }
      }
      
      Matrix newMatrix = new Matrix(newMatrixArray);
      
      Matrix newMatrixInverse = newMatrix.rref();
      
      for (int j = 0; j < length; j++)
      {
         if (newMatrixInverse.matrix[j][j] == 0)
         {
            return null;
         }
      }
      
      // make sure to check here that the rref is actually invertible!
      // this is where it could get complicated if a pivot is lost on the
      //    left-hand side of the matrix!
      
      double[][] inverseGJArray = new double[length][length];
      
      for (int p = 0; p < length; p++)
      {
         for (int q = 0; q < length; q++)
         {
            inverseGJArray[p][q] =
                  newMatrixInverse.matrix[p][q + length];
         }
      }
      Matrix inverseGJ = new Matrix(inverseGJArray);
      
      return inverseGJ;
   }
   
   
   private Matrix identityMatrix(int dimension)
   {
      double[][] identityArray = new double[dimension][dimension];
      for (int i = 0; i < dimension; i++)
      {
         for (int j = 0; j < dimension; j++)
         {
            if (i == j)
            {
               identityArray[i][j] = 1;
            }
            else
            {
               identityArray[i][j] = 0;
            }
         }
      }
      Matrix identityMatrix = new Matrix(identityArray);
      return identityMatrix;
   }
   
   // Cramer's is really inefficient - rarely used and only on small matrices
   // solves systems of equations (passedMatrix = equalsVector) w/Cramer's rule
   // must be a square matrix - build in exceptions!
   private List<Double> cramersRule(Vector equalsVector) throws Exception
   {
      if (!this.squareMatrix())
      {
         throw new Exception("A non-square matrix is not invertible.");
      }
      
      int numberOfVariables = this.n;
      double determinantPassed;
      double determinantCramer;
      double currentVariable;
      List<Double> variableValues = new ArrayList<Double>();
      
      if (numberOfVariables == 2)
      {
         determinantPassed = this.determinant2();
      }
      else if (numberOfVariables == 3)
      {
         determinantPassed = this.determinant3();
      }
      else
      {
         determinantPassed = this.determinantFourPlus();
      }
      
      for (int j = 0; j < numberOfVariables; j++)
      {
         Matrix newMatrix = this.copy();
         
         for (int i = 0; i < this.m; i++)
         {
            newMatrix.matrix[i][j] = equalsVector.getValue(i);
         }
         
         
         if (numberOfVariables == 2)
         {
            determinantCramer = newMatrix.determinant2();
         }
         else if (numberOfVariables == 3)
         {
            determinantCramer = newMatrix.determinant3();
         }
         else
         {
            determinantCramer = newMatrix.determinantFourPlus();
         }
         
         currentVariable = (determinantCramer / determinantPassed);
         variableValues.add(currentVariable);
      }
      return variableValues;
   }
   
   public String cramersRuleToString(Vector equalsVector) throws Exception
   {
      List<Double> variables = this.cramersRule(equalsVector);
      
      String cramerOutput = "By Cramer's Rule, the values of the variables "
            + "in the original system are: \n\n";
      
      for (int i = 0; i < variables.size(); i++)
      {
         cramerOutput = cramerOutput + "x" + Integer.toString(i + 1) + " = "
               + Double.toString(variables.get(i)) + "\n";
      }
      return cramerOutput;
   }
   
   /* From Wikipedia: "The Laplace [cofactor] expansion is computationally
    * inefficient for high dimension matrices, with a time complexity in big O
    * notation of O(n!). Alternatively, using a decomposition into triangular
    * matrices as in the LU decomposition can yield determinants with a time
    * complexity of O(n^3)."
    */
   private double[][] getCofactorMatrix (int mostZerosLine, boolean isRow,
         int currentOther) throws Exception
   {
      // technically methods using this should already check for squareness
      if (!this.squareMatrix())
      {
         throw new Exception("Matrix must be square for this operation.");
      }
      
      double[][] cofactorMatrixArray;
      //Matrix cofactorMatrix;
      
      /*
      if (passedMatrix.m == 3)
      {
         cofactorMatrixArray = new double[2][2];
         
         if (isRow)
         {
            // for loop to cycle through each column
         }
         else
         {
            // for loop to cycle through each row
         }
         
         for (int i = 0; i < 2; i++)
         {
            if (i == row)
            {
               continue;
            }
            else if (i < row)
            {
               for (int j = 0; j < 2; j++)
               {
                  if (j == col)
                  {
                     continue;
                  }
                  else if (j < col)
                  {
                     cofactorMatrixArray[i][j] = passedMatrix.matrix[i][j];
                  }
                  else // if (j > col)
                  {
                     cofactorMatrixArray[i][j] = passedMatrix.matrix[i][j - 1];
                  }
               }
            }
            else // if (i > row)
            {
               for (int j = 0; j < 2; j++)
               {
                  if (j == col)
                  {
                     continue;
                  }
                  else if (j < col)
                  {
                     cofactorMatrixArray[i][j] = passedMatrix.matrix[i - 1][j];
                  }
                  else // if (j > col)
                  {
                     cofactorMatrixArray[i][j] =
                           passedMatrix.matrix[i - 1][j - 1];
                  }
               }
            }
         }
         
         //cofactorMatrix = new Matrix(cofactorMatrixArray);
      }
      else if (passedMatrix.m == 4)
      {
         cofactorMatrixArray = new double[3][3];
         
         for (int i = 0; i < 4; i++)
         {
            if (i == row)
            {
               continue;
            }
            else if (i < row)
            {
               for (int j = 0; j < 4; j++)
               {
                  if (j == col)
                  {
                     continue;
                  }
                  else if (j < col)
                  {
                     cofactorMatrixArray[i][j] = passedMatrix.matrix[i][j];
                  }
                  else // if (j > col)
                  {
                     cofactorMatrixArray[i][j - 1] = passedMatrix.matrix[i][j];
                  }
               }
            }
            else // if (i > row)
            {
               for (int j = 0; j < 4; j++)
               {
                  if (j == col)
                  {
                     continue;
                  }
                  else if (j < col)
                  {
                     cofactorMatrixArray[i - 1][j] = passedMatrix.matrix[i][j];
                  }
                  else // if (j > col)
                  {
                     cofactorMatrixArray[i - 1][j - 1] =
                           passedMatrix.matrix[i][j];
                  }
               }
            }
         }
         
         //cofactorMatrix = new Matrix(cofactorMatrixArray);
      }
      */
      
      //else if (passedMatrix.m > 4)
      if (this.m >= 3)
      {
         int dimension = this.m;
         
         cofactorMatrixArray = 
               new double[dimension - 1][dimension - 1];
         
         // runs in cubic time!
         
         if (isRow)
         {
            // for loop to cycle through each column
            //for (int currentCol = 0; currentCol < dimension; currentCol++)
            //{
               
               
            for (int i = 0; i < dimension; i++)
            {
               if (i == mostZerosLine)
               {
                  continue;
               }
               else if (i < mostZerosLine)
               {
                  for (int j = 0; j < dimension; j++)
                  {
                     if (j == currentOther)
                     {
                        continue;
                     }
                     else if (j < currentOther)
                     {
                        cofactorMatrixArray[i][j] =
                              this.matrix[i][j];
                     }
                     else // if (j > currentCol)
                     {
                        cofactorMatrixArray[i][j - 1] =
                              this.matrix[i][j];
                     }
                  }
               }
               else // if (i > mostZerosLine)
               {
                  for (int j = 0; j < dimension; j++)
                  {
                     if (j == currentOther)
                     {
                        continue;
                     }
                     else if (j < currentOther)
                     {
                        cofactorMatrixArray[i - 1][j] =
                              this.matrix[i][j];
                     }
                     else // if (j > currentCol)
                     {
                        cofactorMatrixArray[i - 1][j - 1] =
                              this.matrix[i][j];
                     }
                  }
               }
            }
         }
         else
         {
            // for loop to cycle through each row
            // ADJUST THIS TO BE TRUE FOR COLUMN COFACTOR EXPANSION
            
            //for (int currentCol = 0; currentCol < dimension; currentCol++)
            //{
            
            for (int j = 0; j < dimension; j++)
            {
               if (j == mostZerosLine)
               {
                  continue;
               }
               else if (j < mostZerosLine)
               {
                  for (int i = 0; i < dimension; i++)
                  {
                     if (i == currentOther)
                     {
                        continue;
                     }
                     else if (i < currentOther)
                     {
                        cofactorMatrixArray[i][j] =
                              this.matrix[i][j];
                     }
                     else // if (i > currentCol)
                     {
                        cofactorMatrixArray[i - 1][j] =
                              this.matrix[i][j];
                     }
                  }
               }
               else // if (j > mostZerosLine)
               {
                  for (int i = 0; i < dimension; i++)
                  {
                     if (i == currentOther)
                     {
                        continue;
                     }
                     else if (i < currentOther)
                     {
                        cofactorMatrixArray[i][j - 1] =
                              this.matrix[i][j];
                     }
                     else // if (i > currentCol)
                     {
                        cofactorMatrixArray[i - 1][j - 1] =
                              this.matrix[i][j];
                     }
                  }
               }
            }
         }
      }
      else
      {
         throw new Exception("Passed matrix is not large enough.");
      }
      return cofactorMatrixArray;
   }
   
   private double findCofactor(int row, int col) throws Exception
   {  
      // technically methods using this should already check for squareness
      if (!this.squareMatrix())
      {
         throw new Exception("Matrix must be square for this operation.");
      }
      
      double cofactor;
      int negativeOneResult;
      double[][] cofactorMatrixArray;
      Matrix cofactorMatrix;
      
      if (((row + 1) + (col + 1)) % 2 == 0)
      {
         negativeOneResult = 1;
      }
      else
      {
         negativeOneResult = -1;
      }
      
      if (this.m == 3)
      {
         cofactorMatrixArray = new double[2][2];
         
         for (int i = 0; i < 3; i++)
         {
            if (i == row)
            {
               continue;
            }
            else if (i < row)
            {
               for (int j = 0; j < 3; j++)
               {
                  if (j == col)
                  {
                     continue;
                  }
                  else if (j < col)
                  {
                     cofactorMatrixArray[i][j] = this.matrix[i][j];
                  }
                  else // if (j > col)
                  {
                     cofactorMatrixArray[i][j - 1] = this.matrix[i][j];
                  }
               }
            }
            else // if (i > row)
            {
               for (int j = 0; j < 3; j++)
               {
                  if (j == col)
                  {
                     continue;
                  }
                  else if (j < col)
                  {
                     cofactorMatrixArray[i - 1][j] = this.matrix[i][j];
                  }
                  else // if (j > col)
                  {
                     cofactorMatrixArray[i - 1][j - 1] =
                           this.matrix[i][j];
                  }
               }
            }
         }
         
         cofactorMatrix = new Matrix(cofactorMatrixArray);
         
         cofactor = /*this.matrix[row][col] **/ negativeOneResult *
               cofactorMatrix.determinant2();
      }
      else if (this.m == 4)
      {
         cofactorMatrixArray = new double[3][3];
         
         for (int i = 0; i < 4; i++)
         {
            if (i == row)
            {
               continue;
            }
            else if (i < row)
            {
               for (int j = 0; j < 4; j++)
               {
                  if (j == col)
                  {
                     continue;
                  }
                  else if (j < col)
                  {
                     cofactorMatrixArray[i][j] = this.matrix[i][j];
                  }
                  else // if (j > col)
                  {
                     cofactorMatrixArray[i][j - 1] = this.matrix[i][j];
                  }
               }
            }
            else // if (i > row)
            {
               for (int j = 0; j < 4; j++)
               {
                  if (j == col)
                  {
                     continue;
                  }
                  else if (j < col)
                  {
                     cofactorMatrixArray[i - 1][j] = this.matrix[i][j];
                  }
                  else // if (j > col)
                  {
                     cofactorMatrixArray[i - 1][j - 1] =
                           this.matrix[i][j];
                  }
               }
            }
         }
         
         cofactorMatrix = new Matrix(cofactorMatrixArray);
         
         cofactor = /*this.matrix[row][col] **/ negativeOneResult *
               cofactorMatrix.determinant3();
      }
      else if (this.m > 4)
      {
         int dimension = this.m;
         
         cofactorMatrixArray = 
               new double[dimension - 1][dimension - 1];
         
         for (int i = 0; i < dimension - 1; i++)
         {
            if (i == row)
            {
               continue;
            }
            else if (i < row)
            {
               for (int j = 0; j < dimension - 1; j++)
               {
                  if (j == col)
                  {
                     continue;
                  }
                  else if (j < col)
                  {
                     cofactorMatrixArray[i][j] = this.matrix[i][j];
                  }
                  else // if (j > col)
                  {
                     cofactorMatrixArray[i][j] = this.matrix[i][j - 1];
                  }
               }
            }
            else // if (i > row)
            {
               for (int j = 0; j < 2; j++)
               {
                  if (j == col)
                  {
                     continue;
                  }
                  else if (j < col)
                  {
                     cofactorMatrixArray[i][j] = this.matrix[i - 1][j];
                  }
                  else // if (j > col)
                  {
                     cofactorMatrixArray[i][j] =
                           this.matrix[i - 1][j - 1];
                  }
               }
            }
         }
         
         cofactorMatrix = new Matrix(cofactorMatrixArray);
         if (this.numberOfZeros == 0)
         {
            cofactorMatrix.setNumberOfZeros(0);
         }
         
         cofactor = /*this.matrix[row][col] **/ negativeOneResult *
               cofactorMatrix.determinantFourPlus();
      }
      else
      {
         throw new Exception("Passed matrix is not large enough.");
      }
      
      return cofactor;
   }
   
   private Matrix findAdjugate() throws Exception
   {
      double[][] adjugateArray = new double[this.m][this.n];
      
      for (int i = 0; i < this.m; i++)
      {
         for (int j = 0; j < this.n; j++)
         {
            adjugateArray[i][j] = this.findCofactor(i, j);
         }
      }
      
      Matrix adjugateMatrix = new Matrix(adjugateArray);
      adjugateMatrix = adjugateMatrix.transpose();
      return adjugateMatrix;
   }
   
   // finds a matrix's inverse through cofactors and the adjugate
   public Matrix inverseCramer() throws Exception
   {
      if (!this.squareMatrix())
      {
         throw new Exception("A non-square matrix is not invertible.");
      }
      
      Matrix adjugate = this.findAdjugate();
      Matrix inverseMatrix;
      
      if (this.m == 2)
      {
         inverseMatrix = adjugate.multiplyByScalar(1 / this.determinant2());
      }
      else if (this.m == 3)
      {
         inverseMatrix = adjugate.multiplyByScalar(
               1 / this.determinant3());
      }
      else
      {
         inverseMatrix = adjugate.multiplyByScalar( 
               (1 / this.determinantFourPlus()));
      }
      return inverseMatrix;
   }
   
   // finds the determinant of a passed 2x2 matrix - BUILD EXCEPTIONS!
   public double determinant2() throws Exception
   {
      if (!this.squareMatrix())
      {
         throw new Exception("A non-square matrix is not invertible.");
      }
      
      double determinant =
            (this.matrix[0][0] * this.matrix[1][1]) -
            (this.matrix[0][1] * this.matrix[1][0]);
      return determinant;
   }
   
   // finds the determinant of a passed 3x3 matrix - BUILD EXCEPTIONS!
   public double determinant3() throws Exception
   {
      if (!this.squareMatrix())
      {
         throw new Exception("A non-square matrix is not invertible.");
      }
      
      double determinant =
            -(this.matrix[2][0] * this.matrix[1][1]
                  * this.matrix[0][2]) -
            (this.matrix[2][1] * this.matrix[1][2]
                  * this.matrix[0][0]) -
            (this.matrix[2][2] * this.matrix[1][0]
                  * this.matrix[0][1]) +
            (this.matrix[0][0] * this.matrix[1][1]
                  * this.matrix[2][2]) +
            (this.matrix[0][1] * this.matrix[1][2]
                  * this.matrix[2][0]) +
            (this.matrix[0][2] * this.matrix[1][0]
                  * this.matrix[2][1]);
      
      return determinant;
   }
   
   // finds the determinant of a matrix with dimensions of 4 or greater
   public double determinantFourPlus() throws Exception
   {
      if (!this.squareMatrix())
      {
         throw new Exception("A non-square matrix is not invertible.");
      }
      if (this.m < 4)
      {
         throw new Exception("Passed matrix is too small.");
      }
      
      double determinant = 0;
      
      /*
      if (this.isTriangular())
      {
         //double determinant = 0;
         determinant = this.matrix[0][0];
         for (int d = 1; d < this.m; d++)
         {
            determinant = determinant * this.matrix[d][d];
         }
         return determinant;
      }
      */
      
      int dimension = this.m;
      //double determinant = 0;
      
      if (dimension == 4)
      {
         int negativeOneResult;
         
         for (int i = 0; i < dimension; i++)
         {
            if (this.matrix[0][i] == 0)
            {
               continue;
            }
            if ((i + 2) % 2 == 0)
            {
               negativeOneResult = 1;
            }
            else
            {
               negativeOneResult = -1;
            }
            double[][] newMatrixArray =
                  new double[this.m - 1][this.n - 1];
            newMatrixArray = 
                  this.getCofactorMatrix(0, true, i);
            Matrix newMatrix = new Matrix(newMatrixArray);
            
            determinant = determinant + (this.matrix[0][i] *
                  negativeOneResult * newMatrix.determinant3());
         }
      }
      
      if (dimension > 4)
      {
         //List<MatrixZeros> zeroMap = new ArrayList<MatrixZeros>();
         int[] rowZeros = new int[dimension];
         int[] colZeros = new int[dimension];
         int rowWithMost = 0;
         int colWithMost = 0;
         int currentRowValue = 0;
         int currentColValue = 0;
         int mostZerosLine;
         boolean isRow;
         //List<Matrix> cofactorMatrices = new ArrayList<Matrix>();
         
         // make sure the search for zeros doesn't run if already set to 0!
         
         // NEED to store each step for the purpose of explanation/solution
         
         
         if (this.numberOfZeros != 0)
         {
            for (int h = 0; h < dimension; h++)
            {
               rowZeros[h] = 0;
               colZeros[h] = 0;
            }
         
            for (int i = 0; i < dimension; i++)
            {
               for (int j = 0; j < dimension; j++)
               {
                  if (this.matrix[i][j] == 0)
                  {
                     this.zeroPositionMap.add(new MatrixZeros(i, j));
                     rowZeros[i]++;
                     colZeros[j]++;
                  }
               }
            }
         }
         
         /* TESTING
         for (MatrixZeros zero : zeroMap)
         {
            System.out.println(Integer.toString(zero.getRow()));
            System.out.println(Integer.toString(zero.getCol()));
            System.out.println("\n");
         }
         
         System.out.println(zeroMap.size());
         */
         
         // TESTING
         //for (int p = 0; p < rowZeros.length; p++)
         //{
         //   System.out.println(Integer.toString(rowZeros[p]));
         //}
         //System.out.println("\n");
         
         //for (int q = 0; q < colZeros.length; q++)
         //{
         //   System.out.println(Integer.toString(colZeros[q]));
         //}
         //System.out.println("\n");
         
         if (this.zeroPositionMap.size() == 0)
         {
            mostZerosLine = 0;
            isRow = true; // just use row 0 by default if there are no zeros
         
            int negativeOneResult;
            
            for (int i = 0; i < dimension; i++)
            {
               if (this.matrix[0][i] == 0)
               {
                  continue;
               }
               if ((i + 2) % 2 == 0)
               {
                  negativeOneResult = 1;
               }
               else
               {
                  negativeOneResult = -1;
               }
               double[][] newMatrixArray =
                     new double[this.m - 1][this.n - 1];
               newMatrixArray = 
                     this.getCofactorMatrix(0, true, i);
               Matrix newMatrix = new Matrix(newMatrixArray);
               
               determinant = determinant + (this.matrix[0][i] *
                     negativeOneResult * newMatrix.determinantFourPlus());
            }
            
         }
         else
         {
            // TESTING
            //for (int t = 0; t < dimension; t++)
            //{
            //   System.out.println(Integer.toString(colZeros[t]));
            //}
            //System.out.println("\n");
            
            for (int k = 0; k < dimension; k++)
            {
               if (rowZeros[k] > currentRowValue)
               {
                  currentRowValue = rowZeros[k];
                  rowWithMost = k;
               }
            }
         
            for (int m = 0; m < dimension; m++)
            {
               if (colZeros[m] > currentColValue)
               {
                  currentColValue = colZeros[m];
                  colWithMost = m;
               }
            }
            
            if (currentRowValue >= currentColValue)
            {
               mostZerosLine = rowWithMost;
               isRow = true;
            }
            else
            {
               mostZerosLine = colWithMost;
               isRow = false;
            }
            
            int negativeOneResult;
            
            if (isRow)
            {
               for (int i = 0; i < this.n; i++)
               {
                  if (this.matrix[mostZerosLine][i] == 0)
                  {
                     continue;
                  }
                  if (((mostZerosLine + 1) + (i + 1)) % 2 == 0)
                  {
                     negativeOneResult = 1;
                  }
                  else
                  {
                     negativeOneResult = -1;
                  }
                  double[][] newMatrixArray =
                        new double[this.m - 1][this.n - 1];
                  newMatrixArray =
                        this.getCofactorMatrix(mostZerosLine, isRow, i);
                  Matrix newMatrix = new Matrix(newMatrixArray);
                  this.setZeroPositions(newMatrix, mostZerosLine,
                        true);
                  determinant =
                        determinant + this.matrix[mostZerosLine][i] *
                        negativeOneResult * newMatrix.determinantFourPlus();
               }
            }
            else
            {
               for (int i = 0; i < this.m; i++)
               {
                  if (this.matrix[i][mostZerosLine] == 0)
                  {
                     continue;
                  }
                  if (((mostZerosLine + 1) + (i + 1)) % 2 == 0)
                  {
                     negativeOneResult = 1;
                  }
                  else
                  {
                     negativeOneResult = -1;
                  }
                  double[][] newMatrixArray =
                        new double[this.m - 1][this.n - 1];
                  newMatrixArray =
                        this.getCofactorMatrix(mostZerosLine, isRow, i);
                  Matrix newMatrix = new Matrix(newMatrixArray);
                  this.setZeroPositions(newMatrix, mostZerosLine,
                        false);
                  determinant =
                        determinant + this.matrix[i][mostZerosLine] *
                        negativeOneResult * newMatrix.determinantFourPlus();
               }
            }
         }
      }
      return determinant;
   }
   
   
   // if there are no zeros, don't go to this method! fix this in D4P method!
   private void setZeroPositions(Matrix cofactorMatrix, int mostZerosLine,
         boolean isRow)
   {
      // potentially also need to set numberOfZeros attribute?
      // determinantFourPlus checks if numberOfZeros == 0
      // but it's never set unless the number is actually 0
      
      int row = -1;
      int col = -1;
      
      if (isRow)
      {
         for (int i = 0; i < this.m; i++)
         {
            for (MatrixZeros zero : this.zeroPositionMap)
            {
               if (zero.getRow() > mostZerosLine)
               {
                  row = zero.getRow() - 1;
               }
               else
               {
                  row = zero.getRow();
               }
               
               if (zero.getCol() > i)
               {
                  col = zero.getCol() - 1;
               }
               else
               {
                  col = zero.getCol();
               }
               cofactorMatrix.zeroPositionMap.add(new MatrixZeros(row, col));
            }
         }
      }
      else
      {
         for (int i = 0; i < this.m; i++)
         {
            for (MatrixZeros zero : this.zeroPositionMap)
            {
               if (zero.getCol() > mostZerosLine)
               {
                  col = zero.getCol() - 1;
               }
               else
               {
                  col = zero.getCol();
               }
               
               if (zero.getRow() > i)
               {
                  row = zero.getRow() - 1;
               }
               else
               {
                  row = zero.getRow();
               }
               cofactorMatrix.zeroPositionMap.add(new MatrixZeros(row, col));
            }
         }
      }
   }
   
   public String determinantToString() throws Exception
   {
      String determinantString = "The determinant of the given matrix is: ";
      df.setRoundingMode(RoundingMode.HALF_UP);
      
      if (this.m < 2)
      {
         throw new Exception("Given matrix does not have a determinant.");
      }
      else if (this.m == 2)
      {
         determinantString =
               determinantString + df.format(this.determinant2());
      }
      else if (this.m == 3)
      {
         determinantString =
               determinantString + df.format(this.determinant3());
      }
      else
      {
         determinantString = determinantString +
               df.format(this.determinantFourPlus());
      }
      return determinantString;
   }
   
   public void findEigenvalues() throws Exception
   {
      if (!this.squareMatrix())
      {
         throw new Exception("A non-square matrix is not invertible.");
      }
      
      if (this.isTriangular())
      {
         for (int i = 0; i < this.m; i++)
         {
            this.realEigenvalueList.add(this.matrix[i][i]);
         }
         Collections.sort(this.realEigenvalueList);
         return;
      }
      
      if (this.m == 2)
      {
         if (findTypeEigenvaluesQuadratic(this.calculateDiscriminant()))
         {
            this.findRealEigenvaluesQuadratic();
         }
         else
         {
            this.findComplexEigenvaluesQuadratic();
         }
      }
      else if (this.m == 3)
      {
         if (findTypeEigenvaluesCubic(this.calculateHCubic()))
         {
            this.findRealEigenvaluesCubic();
         }
         else
         {
            this.findComplexEigenvaluesCubic();
         }
      }
      else if (this.m == 4)
      {
         if (this.findTypeEigenvaluesQuartic())
         {
            this.findRealEigenvaluesQuartic();
         }
         else
         {
            this.findComplexEigenvaluesQuartic();
         }
      }
      Collections.sort(this.realEigenvalueList);
   }
   
   public String eigenvaluesToString()
   {
      String eigenOutput = "";
      String realEigenOutput = 
            "The real eigenvalues for the given matrix are: ";
      String complexEigenOutput =
            "The complex eigenvalues for the given matrix are: ";
      if (this.realEigenvalueList.size() != 0)
      {
         eigenOutput = eigenOutput + realEigenOutput;
         for (double eV : this.realEigenvalueList)
         {
            eigenOutput = eigenOutput + "\n" + df.format(eV);
         }
      }
      if (this.complexEigenvalueList.size() != 0)
      {
         eigenOutput = eigenOutput + complexEigenOutput;
         for (ComplexNumber eV : this.complexEigenvalueList)
         {
            eigenOutput = eigenOutput + "\n" + eV.getNumberRep();
         }
      }
      eigenOutput = eigenOutput + "\n";
      
      return eigenOutput;
   }
   
   private double calculateDiscriminant()
   {
      double a = this.coefficientArray[0];
      double b = this.coefficientArray[1];
      double c = this.coefficientArray[2];
      
      double discriminant = (b * b) - (4 * a * c);
      return discriminant;
   }
   
   // finds roots for a quadratic equation (for eigenvalues for a 2x2 matrix)
   private boolean findTypeEigenvaluesQuadratic(double discriminant)
   {  
      // real roots --> true, complex roots --> false
      if (discriminant >= 0)
      {
         return true;
      }
      else
      {
         return false;
      }
   }
   
   private void findRealEigenvaluesQuadratic()
   {
      double rootRealPositive;
      double rootRealNegative;
      double onlyRealRoot;
      List<Double> realRootList = new ArrayList<Double>();
      
      double a = this.coefficientArray[0];
      double b = this.coefficientArray[1];
      double discriminant = this.calculateDiscriminant();
      
      if (discriminant < 0)
      {
         this.findComplexEigenvaluesQuadratic();
         return;
      }
      if (discriminant > 0)
      {
         rootRealPositive =
               ((-b + Math.sqrt(discriminant)) / (2 * a));
         rootRealNegative =
               ((-b - Math.sqrt(discriminant)) / (2 * a));
         
         realRootList.add(rootRealPositive);
         realRootList.add(rootRealNegative);
      }
      else if (discriminant == 0)
      {
         onlyRealRoot = (-b / (2 * a));
         
         realRootList.add(onlyRealRoot);
         realRootList.add(onlyRealRoot);
      }
   }
   
   private void findComplexEigenvaluesQuadratic()
   {
      double discriminant = this.calculateDiscriminant();
      
      if (discriminant >= 0)
      {
         this.findRealEigenvaluesQuadratic();
         return;
      }
      else
      {
         double realCoefficient =
               -(this.coefficientArray[1]) / 
                     (2 * this.coefficientArray[0]);
         double imagCoefficient =
               Math.sqrt(-(discriminant)) / 
                     (2 * this.coefficientArray[0]);
         
         ComplexNumber complex1 =
               new ComplexNumber(realCoefficient, imagCoefficient, 1);
         ComplexNumber complex2 =
               new ComplexNumber(realCoefficient, -imagCoefficient, 1);
         
         this.complexEigenvalueList.add(complex1);
         this.complexEigenvalueList.add(complex2);
      }
   }
   
   private double calculateHCubic()
   {
      double a = this.coefficientArray[0];
      double b = this.coefficientArray[1];
      double c = this.coefficientArray[2];
      double d = this.coefficientArray[3];
      
      double f = (((3 * c) / a) - ((b * b) / (a * a))) / 3;
      double g =
            (((2 * b * b * b) / (a * a * a)) - ((9 * b * c) / (a * a)) +
                  ((27 * d) / a)) / 27;
      double h = ((g * g) / 4) + ((f * f * f) / 27);
      
      if (h < 0.00001 || h > -0.00001)
      {
         h = 0;
      }
      
      return h;
   }
   
   private static double calculateHCubic(double a, double b,
         double c, double d)
   {
      double f = (((3 * c) / a) - ((b * b) / (a * a))) / 3;
      double g =
            (((2 * b * b * b) / (a * a * a)) - ((9 * b * c) / (a * a)) +
                  ((27 * d) / a)) / 27;
      double h = ((g * g) / 4) + ((f * f * f) / 27);
      
      if (h < 0.00001 || h > -0.00001)
      {
         h = 0;
      }
      
      return h;
   }
   
   private static boolean findTypeEigenvaluesCubic(double h)
   {
      // three real eigenvalues
      if (h <= 0)
      {
         return true;
      }
      // one real, two complex
      else
      {
         return false;
      }
   }
   
   private void findRealEigenvaluesCubic()
   {
      double a = this.coefficientArray[0];
      double b = this.coefficientArray[1];
      double c = this.coefficientArray[2];
      double d = this.coefficientArray[3];
      
      double f = (((3 * c) / a) - ((b * b) / (a * a))) / 3;
      double g =
            (((2 * b * b * b) / (a * a * a)) - ((9 * b * c) / (a * a)) +
                  ((27 * d) / a)) / 27;
      double h = ((g * g) / 4) + ((f * f * f) / 27);
      
      if (h < 0.00001 || h > -0.00001)
      {
         h = 0;
      }
      
      if (h > 0)
      {
         this.findComplexEigenvaluesCubic();
         return;
      }
      else if (h == 0 && f == 0 && g == 0) // all 3 roots real AND equal
      {
         double onlyRealRoot = -(Math.cbrt(d / a));
         this.realEigenvalueList.add(onlyRealRoot);
      }
      else
      {
         double i = Math.sqrt(((g * g) / 4) - h);
         double j = Math.cbrt(i);
         double k = Math.acos(-(g / (2 * i)));
         double l = -j;
         double m = Math.cos(k / 3);
         double n = (Math.sqrt(3) * Math.sin(k / 3));
         double p = (-(b / (3 * a)));
         
         double realRootOne = (2 * j * m) + p;
         double realRootTwo = (l * (m + n)) + p;
         double realRootThree = (l * (m - n)) + p;
         
         this.realEigenvalueList.add(realRootOne);
         this.realEigenvalueList.add(realRootTwo);
         this.realEigenvalueList.add(realRootThree);
      }
   }
   
   private static List<Double> findRealEigenvaluesCubic(double a, double b,
         double c, double d, double h)
   {
      double f = (((3 * c) / a) - ((b * b) / (a * a))) / 3;
      double g =
            (((2 * b * b * b) / (a * a * a)) - ((9 * b * c) / (a * a)) +
                  ((27 * d) / a)) / 27;
      List<Double> auxRoots = new ArrayList<Double>();
      
      if (h > 0)
      {
         findComplexEigenvaluesCubic(a, b, c, d, h);
         return null;
      }
      else if (h == 0 && f == 0 && g == 0) // all 3 roots real AND equal
      {
         double onlyRealRoot = -(Math.cbrt(d / a));
         auxRoots.add(onlyRealRoot);
         auxRoots.add(onlyRealRoot);
         auxRoots.add(onlyRealRoot);
      }
      else
      {
         double i = Math.sqrt(((g * g) / 4) - h);
         double j = Math.cbrt(i);
         double k = Math.acos(-(g / (2 * i)));
         double l = -j;
         double m = Math.cos(k / 3);
         double n = (Math.sqrt(3) * Math.sin(k / 3));
         double p = (-(b / (3 * a)));
         
         double realRootOne = (2 * j * m) + p;
         double realRootTwo = (l * (m + n)) + p;
         double realRootThree = (l * (m - n)) + p;
         
         auxRoots.add(realRootOne);
         auxRoots.add(realRootTwo);
         auxRoots.add(realRootThree);
      }
      return auxRoots;
   }
   
   private void findComplexEigenvaluesCubic()
   {
      /* this *****NEEDS***** include a conversion of the real complex root
       * to a double, and elimination from the complex list and addition
       * to the real list, so the output is formatted/attributed properly.
       */
      
      double a = this.coefficientArray[0];
      double b = this.coefficientArray[1];
      double c = this.coefficientArray[2];
      double d = this.coefficientArray[3];
      
      double f = (((3 * c) / a) - ((b * b) / (a * a))) / 3;
      double g =
            (((2 * b * b * b) / (a * a * a)) - ((9 * b * c) / (a * a)) +
                  ((27 * d) / a)) / 27;
      double h = ((g * g) / 4) + ((f * f * f) / 27);
      
      if (h < 0.00001 || h > -0.00001)
      {
         h = 0;
      }
      
      if (h <= 0)
      {
         this.findRealEigenvaluesCubic();
         return;
      }
      else
      {
         double r = -(g / 2) + Math.sqrt(h);
         double s = Math.cbrt(r);
         double t = -(g / 2) - Math.sqrt(h);
         double u = Math.cbrt(t);
         
         double realRoot = (s + u) - (b / (3 * a));
         this.realEigenvalueList.add(realRoot);
         
         double realCoefficientComplex =
               -((s + u) / 2) - (b / (3 * a));
         double imagCoefficientComplex =
               ((s - u) * Math.sqrt(3)) / 2;
         
         ComplexNumber complexRoot1 =
               new ComplexNumber(realCoefficientComplex,
                     (imagCoefficientComplex), 1);
         ComplexNumber complexRoot2 =
               new ComplexNumber(realCoefficientComplex,
                     -(imagCoefficientComplex), 1);
         
         this.complexEigenvalueList.add(complexRoot1);
         this.complexEigenvalueList.add(complexRoot2);
      }
   }
   
   private static List<ComplexNumber> findComplexEigenvaluesCubic(double a,
         double b, double c, double d, double h)
   {
      //double f = (((3 * c) / a) - ((b * b) / (a * a))) / 3;
      double g =
            (((2 * b * b * b) / (a * a * a)) - ((9 * b * c) / (a * a)) +
                  ((27 * d) / a)) / 27;
      List<ComplexNumber> auxRoots = new ArrayList<ComplexNumber>();
      
      if (h <= 0)
      {
         findRealEigenvaluesCubic(a, b, c, d, h);
         return null;
      }
      else
      {
         double r = -(g / 2) + Math.sqrt(h);
         double s = Math.cbrt(r);
         double t = -(g / 2) - Math.sqrt(h);
         double u = Math.cbrt(t);
         
         double realRoot = (s + u) - (b / (3 * a));
         
         System.out.println(Double.toString(realRoot));
         
         ComplexNumber complexVersionOfReal =
               new ComplexNumber(realRoot, 0, 0);
         auxRoots.add(complexVersionOfReal);
         
         double realCoefficientComplex =
               -((s + u) / 2) - (b / (3 * a));
         double imagCoefficientComplex =
               ((s - u) * Math.sqrt(3)) / 2;
         
         ComplexNumber complexRoot1 =
               new ComplexNumber(realCoefficientComplex,
                     (imagCoefficientComplex), 1);
         ComplexNumber complexRoot2 =
               new ComplexNumber(realCoefficientComplex,
                     -(imagCoefficientComplex), 1);
         
         System.out.println(complexRoot1.getNumberRep());
         System.out.println(complexRoot2.getNumberRep());
         //System.out.println("\n");
         
         auxRoots.add(complexRoot1);
         auxRoots.add(complexRoot2);
      }
      return auxRoots;
   }
   
   /* Possible results for roots of a quartic:
    * (X = tested and working)
    *    X-Two distinct reals, two complex conjugates (R, C)
    *    X-Four distinct reals (R)
    *    -Two pairs of complex conjugates (C)
    *    X-Double real, two simple reals (R)
    *    -Double real, two complex conjugates (R, C)
    *    X-Triple real, one simple real (R)
    *    -Two double reals (R)
    *    X-Two complex conjugate double roots (C)
    *    -Four real roots equal to -b / 4a (R)
    */
   
   private boolean findTypeEigenvaluesQuartic()
   {
      double b = this.coefficientArray[1];
      double c = this.coefficientArray[2];
      double d = this.coefficientArray[3];
      double e = this.coefficientArray[4];
      
      double f = c - ((3 * b * b) / 8);
      double g = d + ((b * b * b) / 8) - ((b * c) / 2);
      double h =
            e - ((3 * b * b * b * b) / 256) +
            ((b * b * c) / 16) - ((b * d) / 4);
      
      double aCubic = 1;
      double bCubic = f / 2;
      double cCubic = ((f * f) - (4 * h)) / 16;
      double dCubic = -(g * g) / 64;
      
      this.internalCubicCoeffs.add(aCubic);
      this.internalCubicCoeffs.add(bCubic);
      this.internalCubicCoeffs.add(cCubic);
      this.internalCubicCoeffs.add(dCubic);
      
      double hCubic = calculateHCubic(aCubic, bCubic, cCubic, dCubic);
      boolean rootType = findTypeEigenvaluesCubic(hCubic);
      
      if (rootType)
      {
         List<Double> realRoots = 
               findRealEigenvaluesCubic(aCubic, bCubic, cCubic, dCubic,
                     hCubic);
         if (realRoots.get(0) < 0 || realRoots.get(1) < 0 ||
               realRoots.get(2) < 0)
         {
            ComplexNumber fakeComplex1 =
                  new ComplexNumber(realRoots.get(0), 0, 0);
            ComplexNumber fakeComplex2 =
                  new ComplexNumber(realRoots.get(1), 0, 0);
            ComplexNumber fakeComplex3 =
                  new ComplexNumber(realRoots.get(2), 0, 0);
            this.internalCubicRootsComplex.add(fakeComplex1);
            this.internalCubicRootsComplex.add(fakeComplex2);
            this.internalCubicRootsComplex.add(fakeComplex3);
            return false;
         }
         else
         {
            this.internalCubicRootsReal.add(realRoots.get(0));
            this.internalCubicRootsReal.add(realRoots.get(1));
            this.internalCubicRootsReal.add(realRoots.get(2));
            return true;
         }
      }
      else
      {
         List<ComplexNumber> complexRoots =
               findComplexEigenvaluesCubic(aCubic, bCubic, cCubic, dCubic,
                     hCubic);
         this.internalCubicRootsComplex.add(complexRoots.get(0));
         this.internalCubicRootsComplex.add(complexRoots.get(1));
         this.internalCubicRootsComplex.add(complexRoots.get(2));
         return false;
      }
   }
   
   private void findRealEigenvaluesQuartic()
   {
      double a = this.coefficientArray[0];
      double b = this.coefficientArray[1];
      double c = this.coefficientArray[2];
      double d = this.coefficientArray[3];
      
      double g = d + ((b * b * b) / 8) - ((b * c) / 2);
      
      double cubicRootOne = this.internalCubicRootsReal.get(0);
      double cubicRootTwo = this.internalCubicRootsReal.get(1);
      double cubicRootThree = this.internalCubicRootsReal.get(2);
      
      double p;
      double q;
      
      // roots used in p and q must be nonzero! add check for this!
      // what if more than one cubic root is zero?!?!?
      if (cubicRootOne == 0)
      {
         p = Math.sqrt(cubicRootTwo);
         q = Math.sqrt(cubicRootThree);
      }
      else if (cubicRootTwo == 0)
      {
         p = Math.sqrt(cubicRootOne);
         q = Math.sqrt(cubicRootThree);
      }
      else
      {
         p = Math.sqrt(cubicRootOne);
         q = Math.sqrt(cubicRootTwo);
      }

      double r;
      
      if (p * q == 0)
      {
         r = 0;
      }
      else
      {
         r = -g / (8 * p * q);
      }
      
      double s = b / (4 * a);
      
      double rootOne = p + q + r - s;
      double rootTwo = p - q - r - s;
      double rootThree = -p + q - r - s;
      double rootFour = -p - q + r - s;
      
      this.realEigenvalueList.add(rootOne);
      this.realEigenvalueList.add(rootTwo);
      this.realEigenvalueList.add(rootThree);
      this.realEigenvalueList.add(rootFour);
   }
   
   private void findComplexEigenvaluesQuartic()
   {
      double a = this.coefficientArray[0];
      double b = this.coefficientArray[1];
      double c = this.coefficientArray[2];
      double d = this.coefficientArray[3];
      
      double g = d + ((b * b * b) / 8) - ((b * c) / 2);
      
      // findTypeEigenvaluesQuartic must be called for these to already be set
      
      //ComplexNumber cubicRootOne =
      //      passedMatrix.internalCubicRootsComplex.get(0);
      ComplexNumber cubicRootTwo =
            this.internalCubicRootsComplex.get(1);
      ComplexNumber cubicRootThree =
            this.internalCubicRootsComplex.get(2);
      
      // List<ComplexNumber> complexQuarticRoots =
      //       new ArrayList<ComplexNumber>();
      ComplexNumber p;
      ComplexNumber q;
      
      // the complex roots of the cubic equation have to be used in p and q!
      
      // p and q are two conjugates! which one to pick?
      if ((cubicRootTwo.complexSqrt().get(0))
            .getRealCoeff() > 0)
      {
         p = cubicRootTwo.complexSqrt().get(0);
      }
      else
      {
         if (b == 0)
         {
            p = cubicRootTwo.complexSqrt().get(0);
         }
         else
         {
            p = cubicRootTwo.complexSqrt().get(1);
         }
      }
      
      if ((cubicRootThree.complexSqrt().get(0))
            .getRealCoeff() > 0)
      {
         q = cubicRootThree.complexSqrt().get(0);
      }
      else
      {
         if (b == 0)
         {
            q = cubicRootThree.complexSqrt().get(0);
         }
         else
         {
            q = cubicRootThree.complexSqrt().get(1);
         }
      }
      
      System.out.println(p.getNumberRep());
      System.out.println(q.getNumberRep());
      
      ComplexNumber pq = p.multiply(q);
      ComplexNumber r;
      
      if (pq.getRealCoeff() == 0 && pq.getImagCoeff() == 0)
      {
         r = new ComplexNumber(0, 0, 0);
      }
      else
      {
         r = (new ComplexNumber(-g, 0, 0)).divide(pq.multiplyByScalar(8));
      }
      ComplexNumber s = (new ComplexNumber(b, 0, 0)).divide(
                  ((new ComplexNumber(a, 0, 0)).multiplyByScalar(4)));
      
      System.out.println(r.getNumberRep());
      System.out.println(s.getNumberRep());
      System.out.println("\n");
      
      System.out.println((p.multiply(q)).getNumberRep());
      System.out.println("\n");
      
      // p + q + r - s
      // p - q - r - s
      // -p + q - r - s
      // -p - q + r - s
      ComplexNumber rootOne = (p.add(q)).add(r.subtract(s));
      ComplexNumber rootTwo = (p.subtract(q)).subtract(r.add(s));
      ComplexNumber rootThree = (q.subtract(p)).subtract(r.add(s));
      ComplexNumber rootFour = (r.subtract(s)).subtract(p.add(q));
      
      this.complexEigenvalueList.add(rootOne);
      this.complexEigenvalueList.add(rootTwo);
      this.complexEigenvalueList.add(rootThree);
      this.complexEigenvalueList.add(rootFour);
   }
   
   /*
   // finds the characteristic equation for a 2x2 matrix
   private String eigenDeterminant2() throws Exception
   {
      // this method also technically should have already validated squareness
      if (!this.squareMatrix())
      {
         throw new Exception("A non-square matrix is not invertible.");
      }
      
      
      // comes out as "ad - (ad)LAMBDA + LAMBDA*LAMBDA"
      // make sure to subtract the quantity (bc)!!!
      
      // finally ends up as "(ad - bc) - (ad)LAMBDA + LAMBDA*LAMBDA"
      
      String eigenDeterminant =
            Double.toString((this.matrix[0][0] * this.matrix[1][1]) -
                  (this.matrix[0][1] * this.matrix[1][0])) + " - " +
                  Double.toString(this.matrix[0][0] * this.matrix[1][1]) +
                  LAMBDA_STRING + " + " + LAMBDA_STRING + "*" + LAMBDA_STRING;
      
      return eigenDeterminant;
   }
   */
   
   /*
   // finds the characteristic equation for a 3x3 matrix
   private String eigenDeterminant3() throws Exception
   {
      // this method also technically should have already validated squareness
      if (!this.squareMatrix())
      {
         throw new Exception("A non-square matrix is not invertible.");
      }
      
      
      /* DETERMINANT =
       * aei + bfg + cdh - afh - bdi - ceg
       *    + (cg)LAMBDA + (fh)LAMBDA + (bd)LAMBDA
       *    - (ae)LAMBDA - (ai)LAMBDA - (ei)LAMBDA
       *    + (a)LAMBDA*LAMBDA + (e)LAMBDA*LAMBDA + (i)LAMBDA*LAMBDA
       *    - LAMBDA*LAMBDA*LAMBDA
       */
      
   /*
      String eigenDeterminant;
      double exponentZero;
      String exponentOne;
      String exponentTwo;
      String exponentThree;
      
      exponentZero =
            (this.matrix[0][0] * this.matrix[1][1] * this.matrix[2][2]) +
            (this.matrix[0][1] * this.matrix[1][2] * this.matrix[2][0]) +
            (this.matrix[0][2] * this.matrix[1][0] * this.matrix[2][1]) -
            (this.matrix[0][0] * this.matrix[1][2] * this.matrix[2][1]) -
            (this.matrix[0][1] * this.matrix[1][0] * this.matrix[2][2]) -
            (this.matrix[0][2] * this.matrix[1][1] * this.matrix[2][0]);
      
      exponentOne =
            Double.toString(this.matrix[0][2] * this.matrix[2][0]) +
               LAMBDA_STRING + " + " + Double.toString(this.matrix[1][2] *
                  this.matrix[2][1]) + LAMBDA_STRING + " + " +
                  Double.toString(this.matrix[0][1] * this.matrix[1][0]) +
                  LAMBDA_STRING + " - " + Double.toString(this.matrix[0][0] *
                  this.matrix[1][1]) + LAMBDA_STRING + " - " +
                  Double.toString(this.matrix[0][0] * this.matrix[2][2]) +
                  LAMBDA_STRING + " - " + Double.toString(this.matrix[1][1] *
                  this.matrix[2][2]) + LAMBDA_STRING;
      
      exponentTwo =
            Double.toString(this.matrix[0][0]) + LAMBDA_STRING + "*" +
                  LAMBDA_STRING + " + " +
                  Double.toString(this.matrix[1][1]) + LAMBDA_STRING +
                  "*" + LAMBDA_STRING + " + " +
                  Double.toString(this.matrix[2][2]) + LAMBDA_STRING +
                  "*" + LAMBDA_STRING;
      
      exponentThree =
            "-" + LAMBDA_STRING + "*" + LAMBDA_STRING + "*" + LAMBDA_STRING;
      
      eigenDeterminant = Double.toString(exponentZero) + exponentOne +
            exponentTwo + exponentThree;
            
      return eigenDeterminant;
   }
   
   /*
   // don't even finish this, it doesn't matter. refer to actual coefficients
   private String eigenDeterminant4Plus() throws Exception
   {
      // this method also technically should have already validated squareness
      if (!this.squareMatrix())
      {
         throw new Exception("A non-square matrix is not invertible.");
      }
      
      // ALSO NEED OPTIONS FOR LARGER EXPONENTS
      
      /* DETERMINANT =
       */
      
   /*
      String eigenDeterminant;
      double exponentZero;
      String exponentOne;
      String exponentTwo;
      String exponentThree;
      String exponentFour;
      
      
      coefficientArray[2] =
            (a * f) + (a * k) + (a * p) + (f * k) + (f * p) + (k * p)
            - (b * e) - (c * i) - (d * m) - (g * j) - (h * n) - (l * o);
      coefficientArray[3] =
            (a * g * j) + (a * h * n) + (a * l * o) + (b * e * k)
            + (b * e * p) + (c * f * i) + (c * i * p) + (d * f * m)
            + (d * k * m) + (f * l * o) + (g * j * p) + (h * k * n)
            - (a * f * k) - (a * f * p) - (a * k * p) - (b * g * i)
            - (b * h * m) - (c * e * j) - (c * l * m) - (d * e * n)
            - (d * i * o) - (f * k * p) - (g * l * n) - (h * j * o);
      coefficientArray[4] =
            (a * f * k * p) + (a * g * l * n) + (a * h * j * o)
            + (b * e * l * o) + (b * g * i * p) + (b * h * k * m)
            + (c * e * j * p) + (c * f * l * m) + (c * h * i * n)
            + (d * e * k * n) + (d * f * i * o) + (d * g * j * m)
            - (a * f * l * o) - (a * g * j * p) - (a * h * k * n)
            - (b * e * k * p) - (b * g * l * m) - (b * h * i * o)
            - (c * e * l * n) - (c * f * i * p) - (c * h * j * m)
            - (d * e * j * o) - (d * f * k * m) - (d * g * i * n);
      
      
      //exponentZero =;
      //exponentOne =;
      //exponentTwo =;
      exponentThree = "-" + Double.toString(this.matrix[0][0]) +
            LAMBDA_STRING + "*" + LAMBDA_STRING + "*" + LAMBDA_STRING + " -" +
            Double.toString(this.matrix[1][1]) + LAMBDA_STRING +
            "*" + LAMBDA_STRING + "*" + LAMBDA_STRING + " -" +
            Double.toString(this.matrix[2][2]) + LAMBDA_STRING +
            "*" + LAMBDA_STRING + "*" + LAMBDA_STRING + " -" +
            Double.toString(this.matrix[3][3]) + LAMBDA_STRING +
            "*" + LAMBDA_STRING + "*" + LAMBDA_STRING;
      
      exponentFour = LAMBDA_STRING + "*" + LAMBDA_STRING +
            "*" + LAMBDA_STRING + "*" + LAMBDA_STRING;
      
      eigenDeterminant = Double.toString(exponentZero) + exponentOne +
            exponentTwo + exponentThree + exponentFour;
      return eigenDeterminant;
      // figure out later
   }
   */
   
   // finds sample eigenvectors for a passed matrix
   public List<Vector> findEigenvectors() throws Exception
   {
      // this method also technically should have already validated squareness
      if (!this.squareMatrix())
      {
         throw new Exception("A non-square matrix is not invertible.");
      }
      
      int dimension = this.m;
      
      if (this.realEigenvalueList.size() == 0)
      {
         this.findEigenvalues();
      }
      
      // FOR MAKING SURE DUPLICATES AREN'T REHASHED
      boolean[] uniqueEVArray = new boolean[this.realEigenvalueList.size()];
      uniqueEVArray[0] = true;
      
      for (int i = 1; i < this.realEigenvalueList.size(); i++)
      {
         if (this.realEigenvalueList.get(i).equals(
               this.realEigenvalueList.get((i - 1))))
         {
            uniqueEVArray[i] = false;
         }
         else
         {
            uniqueEVArray[i] = true;
         }
      }
      
      // make sure to include complex list too
      List<Vector> vectorList = new ArrayList<Vector>();
      
      for (int j = 0; j < uniqueEVArray.length; j++)
      {
         if (uniqueEVArray[j])
         {
            Matrix copyMatrix = this.copy();
            Matrix eigenIdentity =
                  identityMatrix(dimension).multiplyByScalar(
                        this.realEigenvalueList.get(j) * -1);
            copyMatrix = copyMatrix.add(eigenIdentity);
            
            copyMatrix = copyMatrix.rref();
            
            int numberOfVectors = dimension - copyMatrix.numPivots;
            
            List<Vector> tempEigenvectors = new ArrayList<Vector>();
            // create array of zeros here
            double[] zeroArray = new double[dimension];
            for (int f = 0; f < dimension; f++)
            {
               zeroArray[f] = 0;
            }
            
            for (int g = numberOfVectors - 1; g >= 0; g--)
            {
               Vector tempVector = new Vector(zeroArray);
               tempVector.setValue(1, dimension - 1 - g);
               tempEigenvectors.add(tempVector.copy());
            }
            
            boolean pivotFound;
            
            for (int h = 0; h < copyMatrix.numPivots; h++)
            {
               pivotFound = false;
               for (int i = 0; i < dimension; i++)
               {
                  if (copyMatrix.matrix[h][i] != 0 && !pivotFound)
                  {
                     pivotFound = true;
                     continue;
                  }
                  else if (copyMatrix.matrix[h][i] != 0 && pivotFound)
                  {
                     tempEigenvectors.get(i - copyMatrix.numPivots).setValue
                           (-(copyMatrix.matrix[h][i]), h);
                  }
               }
            }
            
            for (Vector tempEV : tempEigenvectors)
            {
               vectorList.add(tempEV.copy());
            }
         }
      }
      return vectorList;
   }
   
   public String eigenvectorsToString() throws Exception
   {
      // find a way to tie eigenvectors to eigenvalues (multiplicities, etc.)
      // do it based on eigenvalues, since they're in sequential order already
      
      String eigenVectorStringOutput = "";
      List<Vector> eVectorList = this.findEigenvectors();
      boolean[] uniqueEVArray = new boolean[this.realEigenvalueList.size()];
      uniqueEVArray[0] = true;
      
      for (int i = 1; i < this.realEigenvalueList.size(); i++)
      {
         if (this.realEigenvalueList.get(i).equals(
               this.realEigenvalueList.get((i - 1))))
         {
            uniqueEVArray[i] = false;
         }
         else
         {
            uniqueEVArray[i] = true;
         }
      }
      
      for (int j = 0; j < uniqueEVArray.length; j++)
      {
         if (uniqueEVArray[j])
         {
            eigenVectorStringOutput =
                  eigenVectorStringOutput +
                  "The eigenvectors corresponding to the eigenvalue " +
                  LAMBDA_STRING + " = " + 
                  df.format(this.realEigenvalueList.get(j)) +
                  //Double.toString(this.realEigenvalueList.get(j)) +
                  " are as follows: \n";
         }
         eigenVectorStringOutput = eigenVectorStringOutput +
               eVectorList.get(j).vectorToString() + "\n";
      }
      return eigenVectorStringOutput;
   }
   
   // determines an A = PDP^-1 diagonalization of a matrix, if possible
   public List<Matrix> diagonalize(List<Double> eigenvalues,
         List<Vector> eigenvectors) throws Exception
   {
      // this method also should technically already have squareness validated
      if (!this.squareMatrix())
      {
         throw new Exception("A non-square matrix is not invertible.");
      }
      
      
      // make sure it's actually diagonalizable, too!
      if (eigenvalues.size() != eigenvectors.size())
      {
         return null;
      }
      
      List<Matrix> pdpInverse = new ArrayList<Matrix>();
      
      Matrix tempIdentity = identityMatrix(this.m);
      for (int i = 0; i < this.m; i++)
      {
         tempIdentity.matrix[i][i] =
               tempIdentity.matrix[i][i] * eigenvalues.get(i);
      }
      
      matrixD = tempIdentity;
      
      for (int j = 0; j < this.n; j++)
      {
         for (int k = 0; k < this.m; k++)
         {
            matrixP.matrix[k][j] = eigenvectors.get(j).getValue(k);
         }
      }
      
      matrixPInverse = matrixP.inverseGaussJordan();
      
      pdpInverse.add(matrixP);
      pdpInverse.add(matrixD);
      pdpInverse.add(matrixPInverse);
      return pdpInverse;
   }
   
   private void diagonalizeBothReal()
   {
      
   }
   
   private void diagonalizeValuesComplex()
   {
      
   }
   
   private void diagonalizeVectorsComplex()
   {
      
   }
   
   private void diagonalizeBothComplex()
   {
      
   }
   
   public String diagonalizeToString()
   {
      
   }
   
}