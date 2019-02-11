/****************************************
 * Linear Algebra Calculator
 * 
 * Supports representation of complex numbers, matrix multiplication
 *    and addition, row reduction, transposing, methods for finding
 *    determinants, several methods for inverses, eigenvalues, eigenvectors,
 *    diagonalization, and operations related to orthogonality.
 * 
 * 
 * ComplexMatrix.java
 * 
 * API:
 *    public ComplexMatrix(int m, int n)
 *    public ComplexMatrix(ComplexNumber[][] matrixEntries)
 *    private ComplexMatrix complexCopy()
 *    private ComplexMatrix multiplyByScalar(ComplexNumber complexScalar)
 *    private ComplexMatrix add(ComplexMatrix matrixOne, ComplexMatrix matrixTwo)
 *    public ComplexMatrix multiplyMatrices(ComplexMatrix matrixTwo)
 *    public String multiplyToString()
 *    private boolean validateSameSize(ComplexMatrix matrixOne,
 *                          ComplexMatrix matrixTwo)
 *    private boolean squareMatrix()
 *    private boolean equalMatrices(ComplexMatrix matrixOne,
 *                         ComplexMatrix matrixTwo)
 *    private boolean isTriangular() throws Exception
 *    private boolean isUpperTriangular() throws Exception
 *    private boolean isLowerTriangular() throws Exception
 *    public int getCurrentRow()
 *    public int getNumPivots()
 *    public ComplexMatrix rref()
 *    public List<ComplexMatrix> getRREFSteps()
 *    public String rRefToString()
 *    private void findPivotColumn()
 *    private ComplexMatrix reduceUpward()
 *    public ComplexMatrix ref()
 *    public List<ComplexMatrix> getRefSteps()
 *    public String refToString()
 *    private void findPivotRow()
 *    private ComplexMatrix reduceDownward()
 *    private ComplexMatrix interchange(int pivotRow, int otherRow)
 *    private ComplexMatrix complexScale()
 *    public ComplexMatrix transpose()
 *    public String transposeToString()
 *    public ComplexMatrix inverse2()
 *    public String inverseToString()
 *    public ComplexMatrix inverseGaussJordan()
 *    private ComplexMatrix identityMatrix(int dimension)
 *    public List<ComplexNumber> cramersRule(ComplexVector equalsVector)
 *    public String cramersRuleToString(List<ComplexNumber> variables)
 *    private ComplexNumber findCofactor(int row, int col)
 *    private ComplexMatrix findAdjugate()
 *    public ComplexMatrix inverseCramer()
 *    public ComplexNumber determinant2()
 *    public ComplexNumber determinant3()
 *    public String determinantToString() throws Exception
 *    public void findEigenvalues()
 *    public List<ComplexNumber> complexMatrixEValues()
 *    public String eigenvaluesToString()
 *    private void findCMEigenvaluesQuadratic()
 *    
 *    
 * Dependencies:
 *    -ComplexNumber.java
 *    -ComplexVector.java
 *    -MatrixZeros.java
 * 
 ****************************************/

import java.util.List;
import java.math.RoundingMode;
import java.text.DecimalFormat;
import java.util.ArrayList;
import java.util.Collections;

// class for matrices containing any complex values
public class ComplexMatrix
{
   public static final ComplexNumber COMPLEX_ZERO = new ComplexNumber(0, 0, 0);
   public static final String MATRIX_SPACING = "                   ";
   
   // private member data
   private final int m;
   private final int n;
   private final ComplexNumber[][] matrix;
   private int                   numberOfZeros;
   
   private ComplexNumber[]       coefficientArray;
   private List<MatrixZeros>     zeroPositionMap;
   private boolean               zerosSet;
   private List<Double>          realEigenvalueList;
   private List<ComplexNumber>   complexEigenvalueList;
   private List<Double>          internalCubicCoeffs;
   private List<Double>          internalCubicRootsReal;
   private List<ComplexNumber>   internalCubicRootsComplex;
   private List<String>          listOfRefSteps;
   private List<String>          listOfRREFSteps;
   
   // private data for use in multiple methods
   private int                   pivotColumn;
   private int                   currentRow;
   private int                   numPivots;
   private int                   lastRowWithPivot;
   private boolean               refDoneFlag;
   private List<ComplexMatrix>   refStepsList;
   private List<ComplexMatrix>   rRefStepsList;
   private List<String>          refExplanationList;
   private List<String>          rRefExplanationList;
   private Matrix                matrixP;
   private Matrix                matrixD;
   private Matrix                matrixPInverse;
   private static DecimalFormat  df = new DecimalFormat("#.#####");
   
   
   // RENAME ALL METHODS TO INCLUDE THE WORD "COMPLEX".
   
   /*
   // constructor with dimensions
   public ComplexMatrix(int m, int n)
   {
      this.m = m;
      this.n = n;
      matrix = new ComplexNumber[m][n];
   }
   */
   
   // constructor with actual matrix entries
   public ComplexMatrix(ComplexNumber[][] matrixEntries)
   {
      m = matrixEntries.length;
      n = matrixEntries[0].length;
      this.matrix = new ComplexNumber[m][n];
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
      
      refStepsList = new ArrayList<ComplexMatrix>();
      rRefStepsList = new ArrayList<ComplexMatrix>();
      
      refExplanationList = new ArrayList<String>();
      rRefExplanationList = new ArrayList<String>();
      
      listOfRefSteps = new ArrayList<String>();
      listOfRREFSteps = new ArrayList<String>();
      
      if (m == n)
      {
         zeroPositionMap = new ArrayList<MatrixZeros>();
         zerosSet = false;
         realEigenvalueList = new ArrayList<Double>();
         complexEigenvalueList = new ArrayList<ComplexNumber>();
         
         coefficientArray = new ComplexNumber[m + 1];
         
         if (m == 2)
         {
            coefficientArray[0] = new ComplexNumber(1, 0, 0);
            coefficientArray[1] = // -(a + d)
                  ((this.matrix[0][0]).add(this.matrix[1][1])).negate();
                  //-(this.matrix[0][0] + this.matrix[1][1]);
            coefficientArray[2] = // (ad - bc)
                  ((this.matrix[0][0]).multiply(this.matrix[1][1])).subtract(
                        ((this.matrix[0][1].multiply(this.matrix[1][0]))));
                  //((this.matrix[0][0] * this.matrix[1][1])
                  //- (this.matrix[0][1] * this.matrix[1][0]));
         }
         else if (m == 3)
         {
            coefficientArray[0] = new ComplexNumber(-1, 0, 0);
            coefficientArray[1] = //a + e + i
                  this.matrix[0][0].add((this.matrix[1][1]).add(
                        this.matrix[2][2]));
                  
                  //(this.matrix[0][0] + this.matrix[1][1] +
                  //      this.matrix[2][2]);
            coefficientArray[2] = // cg + fh + bd - ae - ai - ei
                  (((this.matrix[0][2].multiply(this.matrix[2][0])).add(
                        (this.matrix[1][2].multiply(this.matrix[2][1])))).add(
                        ((this.matrix[0][1].multiply(this.matrix[1][0]))
                        .subtract((this.matrix[0][0].multiply(
                        this.matrix[1][1]))))).add(
                        (((this.matrix[0][0].multiply(this.matrix[2][2]))
                        .negate()).subtract((this.matrix[1][1].multiply(
                        this.matrix[2][2]))))));
                  
                  //(this.matrix[0][2] * this.matrix[2][0]) +
                  //(this.matrix[1][2] * this.matrix[2][1]) +
                  //(this.matrix[0][1] * this.matrix[1][0]) -
                  //(this.matrix[0][0] * this.matrix[1][1]) -
                  //(this.matrix[0][0] * this.matrix[2][2]) -
                  //(this.matrix[1][1] * this.matrix[2][2]);
            coefficientArray[3] = // aei + bfg + cdh - afh - bdi - ceg
                  ((((this.matrix[0][0].multiply(this.matrix[1][1])).multiply(
                        this.matrix[2][2])).add((this.matrix[0][1].multiply(
                        this.matrix[1][2])).multiply(this.matrix[2][0]))).add(
                        (((this.matrix[0][2].multiply(this.matrix[1][0]))
                        .multiply(this.matrix[2][1])).subtract((
                        this.matrix[0][0].multiply(this.matrix[1][2]))
                        .multiply(this.matrix[2][1]))))).add(
                        ((this.matrix[0][1].multiply(this.matrix[1][0]))
                        .multiply(this.matrix[2][2])).negate().subtract(
                        (this.matrix[0][2].multiply(
                        this.matrix[1][1])).multiply(this.matrix[2][0])));
                  
                  
                  //(this.matrix[0][0] * this.matrix[1][1] *
                  //      this.matrix[2][2]) +
                  //(this.matrix[0][1] * this.matrix[1][2] *
                  //      this.matrix[2][0]) +
                  //(this.matrix[0][2] * this.matrix[1][0] *
                  //      this.matrix[2][1]) -
                  //(this.matrix[0][0] * this.matrix[1][2] *
                  //      this.matrix[2][1]) -
                  //(this.matrix[0][1] * this.matrix[1][0] *
                  //      this.matrix[2][2]) -
                  //(this.matrix[0][2] * this.matrix[1][1] *
                  //      this.matrix[2][0]);
         }
         else if (m == 4)
         {
            internalCubicCoeffs = new ArrayList<Double>();
            internalCubicRootsReal = new ArrayList<Double>();
            internalCubicRootsComplex = new ArrayList<ComplexNumber>();
            
            ComplexNumber a = this.matrix[0][0];
            ComplexNumber b = this.matrix[0][1];
            ComplexNumber c = this.matrix[0][2];
            ComplexNumber d = this.matrix[0][3];
            ComplexNumber e = this.matrix[1][0];
            ComplexNumber f = this.matrix[1][1];
            ComplexNumber g = this.matrix[1][2];
            ComplexNumber h = this.matrix[1][3];
            ComplexNumber i = this.matrix[2][0];
            ComplexNumber j = this.matrix[2][1];
            ComplexNumber k = this.matrix[2][2];
            ComplexNumber l = this.matrix[2][3];
            ComplexNumber m = this.matrix[3][0];
            ComplexNumber n = this.matrix[3][1];
            ComplexNumber o = this.matrix[3][2];
            ComplexNumber p = this.matrix[3][3];
            
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
            
            coefficientArray[0] = new ComplexNumber(1, 0, 0);
            coefficientArray[1] = 
                  ((a.add(f)).add(k.add(p))).negate();
                  //- (a + f + k + p);
            coefficientArray[2] =
                  (((a.multiply(f)).add((a.multiply(k)))).add((a.multiply(
                        p)).add(f.multiply(k)))).add(((f.multiply(p)).add(
                        k.multiply(p)))).add(((((b.multiply(e)).add(
                        (c.multiply(i)))).add((d.multiply(m)).add(g.multiply(
                        j)))).add(((h.multiply(n)).add(l.multiply(o))))
                        .negate()));
                  
                  //(a * f) + (a * k) + (a * p) + (f * k) + (f * p) + (k * p) -
                  //(b * e) - (c * i) - (d * m) - (g * j) - (h * n) - (l * o);
            coefficientArray[3] =
                  ((((((a.multiply(g)).multiply(j)).add(
                        (a.multiply(h)).multiply(n))).add(
                        ((a.multiply(l)).multiply(o)).add(
                        (b.multiply(e)).multiply(k)))).add(
                        (((b.multiply(e)).multiply(p)).add(
                        (c.multiply(f)).multiply(i))).add(
                        ((c.multiply(i)).multiply(p)).add(
                        (d.multiply(f)).multiply(m))))).add(
                        (((d.multiply(k)).multiply(m)).add(
                        (f.multiply(l)).multiply(o))).add(
                        ((g.multiply(j)).multiply(p)).add(
                        (h.multiply(k)).multiply(n))))).add(
                        ((((((a.multiply(f)).multiply(k)).add(
                        (a.multiply(f)).multiply(p))).add(
                        ((a.multiply(k)).multiply(p)).add(
                        (b.multiply(g)).multiply(i)))).add(
                        (((b.multiply(h)).multiply(m)).add(
                        (c.multiply(e)).multiply(j))).add(
                        ((c.multiply(l)).multiply(m)).add(
                        (d.multiply(e)).multiply(n))))).add(
                        (((d.multiply(i)).multiply(o)).add(
                        (f.multiply(k)).multiply(p))).add(
                        ((g.multiply(l)).multiply(n)).add(
                        (h.multiply(j)).multiply(o))))).negate());
                  
                  /*
                  (a * g * j) + (a * h * n) + (a * l * o) + (b * e * k)
                  + (b * e * p) + (c * f * i) + (c * i * p) + (d * f * m)
                  + (d * k * m) + (f * l * o) + (g * j * p) + (h * k * n)
                  - (a * f * k) - (a * f * p) - (a * k * p) - (b * g * i)
                  - (b * h * m) - (c * e * j) - (c * l * m) - (d * e * n)
                  - (d * i * o) - (f * k * p) - (g * l * n) - (h * j * o);
                  */
            coefficientArray[4] =
                  ((((((a.multiply(f)).multiply(k.multiply(p))).add(
                        (a.multiply(g)).multiply(l.multiply(n)))).add(
                        ((a.multiply(h)).multiply(j.multiply(o))).add(
                        (b.multiply(e)).multiply(l.multiply(o))))).add(
                        (((b.multiply(g)).multiply(i.multiply(p))).add(
                        (b.multiply(h)).multiply(k.multiply(m)))).add(
                        ((c.multiply(e)).multiply(j.multiply(p))).add(
                        (c.multiply(f)).multiply(l.multiply(m)))))).add(
                        (((c.multiply(h)).multiply(i.multiply(n))).add(
                        (d.multiply(e)).multiply(k.multiply(n)))).add(
                        ((d.multiply(f)).multiply(i.multiply(o))).add(
                        (d.multiply(g)).multiply(j.multiply(m)))))).add(
                        ((((((a.multiply(f)).multiply(l.multiply(o))).add(
                        (a.multiply(g)).multiply(j.multiply(p)))).add(
                        ((a.multiply(h)).multiply(k.multiply(n))).add(
                        (b.multiply(e)).multiply(k.multiply(p))))).add(
                        (((b.multiply(g)).multiply(l.multiply(m))).add(
                        (b.multiply(h)).multiply(i.multiply(o)))).add(
                        ((c.multiply(e)).multiply(l.multiply(n))).add(
                        (c.multiply(f)).multiply(i.multiply(p)))))).add(
                        (((c.multiply(h)).multiply(j.multiply(m))).add(
                        (d.multiply(e)).multiply(j.multiply(o)))).add(
                        ((d.multiply(f)).multiply(k.multiply(m))).add(
                        (d.multiply(g)).multiply(i.multiply(n)))))).negate());
                  
                  
                  /*
                  (a * f * k * p) + (a * g * l * n) + (a * h * j * o)
                  + (b * e * l * o) + (b * g * i * p) + (b * h * k * m)
                  + (c * e * j * p) + (c * f * l * m) + (c * h * i * n)
                  + (d * e * k * n) + (d * f * i * o) + (d * g * j * m)
                  - (a * f * l * o) - (a * g * j * p) - (a * h * k * n)
                  - (b * e * k * p) - (b * g * l * m) - (b * h * i * o)
                  - (c * e * l * n) - (c * f * i * p) - (c * h * j * m)
                  - (d * e * j * o) - (d * f * k * m) - (d * g * i * n);
                  */
         }
      }
   }
   
   private ComplexMatrix complexCopy()
   {
      ComplexNumber[][] complexCopyArray = new ComplexNumber[this.m][this.n];
      for (int i = 0; i < this.m; i++)
      {
         for (int j = 0; j < this.n; j++)
         {
            complexCopyArray[i][j] = this.matrix[i][j];
         }
      }
      ComplexMatrix matrixCopy = new ComplexMatrix(complexCopyArray);
      return matrixCopy;
   }
   
   // multiplies a matrix or vector by a scalar quantity
   private ComplexMatrix multiplyByScalar(ComplexNumber complexScalar)
   {
      for (int i = 0; i < this.m; i++)
      {
         for (int j = 0; j < this.n; j++)
         {
            this.matrix[i][j] = this.matrix[i][j].multiply(complexScalar);
         }
      }
      return this;
   }
   
   // adds two matrices or vectors together
   private ComplexMatrix add(ComplexMatrix matrixOne, ComplexMatrix matrixTwo)
   {
      if (validateSameSize(matrixOne, matrixTwo))
      {
         // modify max for particular matrices
         for (int i = 0; i < matrixOne.m; i++)
         {
            for (int j = 0; j < matrixOne.n; j++)
            {
               matrixOne.matrix[i][j] = 
                     matrixOne.matrix[i][j].add(matrixTwo.matrix[i][j]);
            }
         }
      }
      return matrixOne;
   }
   
   public ComplexMatrix multiplyMatrices(ComplexMatrix matrixTwo)
         throws Exception
   {
      if (this.n != matrixTwo.m)
      {
         throw new Exception("Matrices cannot be multiplied.");
      }
      
      ComplexNumber[][] multipliedArray =
            new ComplexNumber[this.m][matrixTwo.n];
      
      // need to initialize this array, here and in the regular matrix class!!!
      
      for (int i = 0; i < this.m; i++)
      {
         for (int j = 0; j < matrixTwo.n; j++)
         {
            multipliedArray[i][j] = new ComplexNumber(0, 0, 0);
            for (int k = 0; k < this.n; k++)
            {
               // check to see if this is even right.
               multipliedArray[i][j] = (multipliedArray[i][j]).add(
                     ((this.matrix[i][k]).multiply(matrixTwo.matrix[k][j])));
            }
         }
      }
      ComplexMatrix multResult = new ComplexMatrix(multipliedArray);
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
            multipliedString =
                  multipliedString + (this.matrix[i][j]).getNumberRep()
            + MATRIX_SPACING.substring(0, (MATRIX_SPACING.length() - 
                  (((this.matrix[i][j]).getNumberRep()).length()))) + "   ";
         }
         multipliedString = multipliedString + "\n";
      }
      return multipliedString;
   }
   
   // makes sure that two matrices/vectors are same size and can be added
   private boolean validateSameSize(ComplexMatrix matrixOne,
         ComplexMatrix matrixTwo)
   {
      if (matrixOne.m == matrixTwo.m && matrixOne.n == matrixTwo.n)
      {
         return true;
      }
      return false;
   }
   
   // validates whether or not a matrix is square
   // needed for inverses, eigenstuff, Cramer's, ... 
   private boolean squareMatrix()
   {
      if (this.m == this.n)
      {
         return true;
      }
      return false;
   }
   
   // checks if two matrices/vectors are perfectly equal
   private boolean equalMatrices(ComplexMatrix matrixOne,
         ComplexMatrix matrixTwo)
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
   
   // FIX THIS METHOD.
   private boolean isTriangular() throws Exception
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
   
   private boolean isUpperTriangular() throws Exception
   {
      if (!this.squareMatrix())
      {
         throw new Exception("Matrix must be square to be triangular.");
      }
      
      boolean upperTriangular = true;
      
      for (int i = 1; i < this.m; i++)
      {
         for (int j = i - 1; j < this.n - 1; j++)
         {
            // upper triangular
            if (this.matrix[i][j].getRealCoeff() != 0 ||
                  this.matrix[i][j].getImagCoeff() != 0 ||
                  this.matrix[i][j].getIExponent() != 0)
            //if (this.matrix[i][j] != 0)
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
            if (this.matrix[k][l].getRealCoeff() != 0 ||
                  this.matrix[k][l].getImagCoeff() != 0 ||
                  this.matrix[k][l].getIExponent() != 0)
            //if (this.matrix[k][l] != 0)
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
   
   private int getCurrentRow()
   {
      return this.currentRow;
   }
   
   private int getNumPivots()
   {
      return this.numPivots;
   }
   
   // determines the row-reduced echelon form of the passed matrix
   public ComplexMatrix rref()
   {    
      ComplexMatrix refPart = this.ref();
      refPart.rRefStepsList = refPart.refStepsList;
      refPart.rRefExplanationList = refPart.refExplanationList;
      
      refPart.currentRow = refPart.numPivots - 1; // what if there's one pivot?
      
      while (refPart.currentRow > 0)
      {
         refPart.findPivotColumn();
         /*
         if (refPart.matrix[refPart.currentRow][refPart.pivotColumn] != 1)
         {
            refPart = refPart.scale();
            //rRefStepsList.add(refPart.copy());
         }
         */
         refPart = refPart.reduceUpward();
         //rRefStepsList.add(refPart.copy());
      }
      
      //refPart.findPivotColumn();
      /*
      if (refPart.matrix[refPart.currentRow][refPart.pivotColumn] != 1)
      {
         refPart = refPart.scale();
      }
      */
      return refPart;
   }
   
   public List<ComplexMatrix> getRREFSteps()
   {
      return this.rRefStepsList;
   }
   
   public String rRefToString()
   {
      String rRefOutput = "The steps leading to a *reduced* row echelon form "
            + "of the given matrix are: \n";
      // this runs in cubic time!
      for (ComplexMatrix rS : this.rRefStepsList)
      {
         rRefOutput = rRefOutput + "\n" +
               this.refExplanationList.get(this.refStepsList.indexOf(rS)) +
               "\n";
         for (int i = 0; i < rS.m; i++)
         {
            for (int j = 0; j < rS.n; j++)
            {
               rRefOutput = rRefOutput + rS.matrix[i][j].getNumberRep()
                     + MATRIX_SPACING.substring(0, (MATRIX_SPACING.length() - 
                           (rS.matrix[i][j].getNumberRep().length()))) + "   ";
               
               /*
               if (rS.matrix[i][j] >= 0)
               {
                  rRefOutput = rRefOutput + " " +
                        df.format(rS.matrix[i][j]) + "   ";
                        //Double.toString(rS.matrix[i][j]) + "   ";
               }
               else
               {
                  rRefOutput = rRefOutput +
                     df.format(rS.matrix[i][j]) + "   ";
                     //Double.toString(rS.matrix[i][j]) + "   ";
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
   
   private void findPivotColumn()
   {
      // need to find pivot position when reducing upward
      // many methods in RREF rely on pivotColumn value for calculations
      
      for (int j = 0; j < this.n; j++)
      {
         if (!(this.matrix[this.currentRow][j].getRealCoeff() == 0.0 &&
               this.matrix[this.currentRow][j].getImagCoeff() == 0.0 &&
               this.matrix[this.currentRow][j].getIExponent() == 0.0))
         //if (!(this.matrix[this.currentRow][j].equals(COMPLEX_ZERO)));
         {
            this.pivotColumn = j;
            return;
         }
      }
   }
   
   private ComplexMatrix reduceUpward()
   {
      for (int i = 0; i < this.currentRow; i++)
      {
         ComplexNumber currMultFactor =
               (new ComplexNumber(
                     this.matrix[i][this.pivotColumn].getRealCoeff(),
                     this.matrix[i][this.pivotColumn].getImagCoeff(), 1)
                     ).negate();
         
         for (int j = 0; j < this.n; j++)
         {
            this.matrix[i][j] = this.matrix[i][j].add(
                  (this.matrix[this.currentRow][j]).multiply(currMultFactor));
            
            if (this.matrix[i][j].getRealCoeff() == -0.0)
            {
               this.matrix[i][j].setRealCoeff(0.0);
            }
            if (this.matrix[i][j].getImagCoeff() == -0.0)
            {
               this.matrix[i][j].setImagCoeff(0.0);
            }
            if (this.matrix[i][j].getRealCoeff() < 0.00001 &&
                  this.matrix[i][j].getRealCoeff() > -0.00001)
            {
               this.matrix[i][j].setRealCoeff(0.0);
            }
            if (this.matrix[i][j].getImagCoeff() < 0.00001 &&
                  this.matrix[i][j].getImagCoeff() > -0.00001)
            {
               this.matrix[i][j].setImagCoeff(0.0);
            }
         }
      }
      // so pivot finding and reduction continues upward in matrix
      this.currentRow--;
      rRefExplanationList.add("Reducing the matrix upward from row " +
            Integer.toString(currentRow + 1) + " yields: ");
      rRefStepsList.add(this.complexCopy());
      return this;
   }
   
   // determines one (of many) reduced echelon form of the passed matrix
   public ComplexMatrix ref()
   {
      ComplexMatrix reductionStep = this.complexCopy();
      reductionStep.refStepsList.add(reductionStep.complexCopy());
      reductionStep.refExplanationList.add("Original matrix: ");
      
      while (!reductionStep.refDoneFlag)
      {
         reductionStep.findPivotRow();
         reductionStep.complexScale();
         reductionStep = reductionStep.reduceDownward();
         
         if (reductionStep.currentRow == reductionStep.m)
         {
            reductionStep.refDoneFlag = true;
         }
      }
      return reductionStep;
   }
   
   public List<ComplexMatrix> getRefSteps()
   {
      return this.refStepsList;
   }
   
   public String refToString()
   {
      String refOutput = "The steps leading to a row echelon form of the "
            + "given matrix are: \n";
      df.setRoundingMode(RoundingMode.HALF_UP);
      
      //System.out.println(Integer.toString(this.refStepsList.size()));
      
      // this runs in cubic time!
      for (ComplexMatrix rS : this.refStepsList)
      {
         //System.out.println(Integer.toString(this.refStepsList.indexOf(rS)));
         refOutput = refOutput + "\n" +
               this.refExplanationList.get(this.refStepsList.indexOf(rS)) +
               "\n";
         for (int i = 0; i < rS.m; i++)
         {
            for (int j = 0; j < rS.n; j++)
            {
               refOutput = refOutput + rS.matrix[i][j].getNumberRep()
                     + MATRIX_SPACING.substring(0, (MATRIX_SPACING.length() - 
                           (rS.matrix[i][j].getNumberRep().length()))) + "   ";
               
               
               /*
               if (rS.matrix[i][j].getRealCoeff() >= 0)
               {
                  refOutput = refOutput + " " +
                        rS.matrix[i][j].getNumberRep() + "   ";
                        //df.format(rS.matrix[i][j]) + "   ";
                        //Double.toString(rS.matrix[i][j]) + "   ";
               }
               else
               {
                  refOutput = refOutput +
                        rS.matrix[i][j].getNumberRep() + "   ";
                        //df.format(rS.matrix[i][j]) + "   ";
                        //Double.toString(rS.matrix[i][j]) + "   ";
               }
               */
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
            if (!(this.matrix[i][j].getRealCoeff() == 0.0 &&
                  this.matrix[i][j].getImagCoeff() == 0.0 &&
                  this.matrix[i][j].getIExponent() == 0.0))
            
            //if (!this.matrix[i][j].equals(COMPLEX_ZERO))
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
   private ComplexMatrix reduceDownward()
   {      
      // SHOULD ONLY ACTUALLY BE USED IF A ROW BELOW ALSO CONTAINS
      // NONZERO VALUE IN SAME COLUMN. WASTE OF TIME OTHERWISE.
      
      boolean reductionOccurred = false;
      
      for (int i = this.currentRow + 1; i < this.m; i++)
      {
         // only if matrixArray[i][pivotColumn] != 0
         // only reduces a row if the value in the pivot column is nonzero
         if (this.matrix[i][this.pivotColumn].equals(COMPLEX_ZERO))
         {
            continue;
         }
         else
         {
            reductionOccurred = true;
            // changes with each following row
            ComplexNumber currMultFactor =
                  (new ComplexNumber(
                        this.matrix[i][this.pivotColumn].getRealCoeff(),
                        this.matrix[i][this.pivotColumn].getImagCoeff(), 1)
                        ).negate();
                  //(this.matrix[i][this.pivotColumn]).negate();
            
            for (int j = this.pivotColumn; j < this.n; j++)
            {
               this.matrix[i][j] = // is this going to change internal values?
                     (this.matrix[i][j]).add(
                           (this.matrix[this.currentRow][j]).multiply(
                                 currMultFactor));
               
               if (this.matrix[this.currentRow][j].getRealCoeff() == -0.0)
               {
                  this.matrix[this.currentRow][j].setRealCoeff(0.0);
               }
               if (this.matrix[this.currentRow][j].getImagCoeff() == -0.0)
               {
                  this.matrix[this.currentRow][j].setImagCoeff(0.0);
               }
               if (this.matrix[this.currentRow][j].getRealCoeff() < 0.00001 &&
                     this.matrix[this.currentRow][j].getRealCoeff() > -0.00001)
               {
                  this.matrix[this.currentRow][j].setRealCoeff(0.0);
               }
               if (this.matrix[this.currentRow][j].getImagCoeff() < 0.00001 &&
                     this.matrix[this.currentRow][j].getImagCoeff() > -0.00001)
               {
                  this.matrix[this.currentRow][j].setImagCoeff(0.0);
               }
            }
         }
      }
      
      if (reductionOccurred)
      {
         this.refStepsList.add(this.complexCopy());
      }
      
      refExplanationList.add("Performing row replacement based on row "
            + Integer.toString(this.currentRow + 1) + " yields: ");
      // so pivot finding and reduction continues downward in matrix
      this.currentRow++;
      return this;
   }
   
   // moves row with earlier pivot above other rows
   // pivot row first, less significant row second
   private ComplexMatrix interchange(int pivotRow, int otherRow)
   {
      ComplexNumber[] tempArray = new ComplexNumber[this.n];
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
      this.refStepsList.add(this.complexCopy());
      return this;
   }
   
   
   // multiplies a row of a matrix by a scalar quantity to make its pivot 1
   private ComplexMatrix complexScale()
   {
      // for loop finds pivot
      ComplexNumber scaleBase = this.matrix[this.currentRow][this.pivotColumn];
      if (this.matrix[this.currentRow][this.pivotColumn].getRealCoeff() != 1 &&
            this.matrix[this.currentRow][this.pivotColumn].getImagCoeff() != 0)
      {
         this.matrix[this.currentRow][this.pivotColumn] =
               new ComplexNumber(1, 0, 0);
         for (int j = this.pivotColumn + 1; j < this.n; j++)
         {  
            this.matrix[this.currentRow][j] =
                  this.matrix[this.currentRow][j].divide(scaleBase);
            
            if (this.matrix[this.currentRow][j].getRealCoeff() == -0.0)
            {
               this.matrix[this.currentRow][j].setRealCoeff(0.0);
            }
            if (this.matrix[this.currentRow][j].getImagCoeff() == -0.0)
            {
               this.matrix[this.currentRow][j].setImagCoeff(0.0);
            }
            if (this.matrix[this.currentRow][j].getRealCoeff() < 0.00001 &&
                  this.matrix[this.currentRow][j].getRealCoeff() > -0.00001)
            {
               this.matrix[this.currentRow][j].setRealCoeff(0.0);
            }
            if (this.matrix[this.currentRow][j].getImagCoeff() < 0.00001 &&
                  this.matrix[this.currentRow][j].getImagCoeff() > -0.00001)
            {
               this.matrix[this.currentRow][j].setImagCoeff(0.0);
            }
         }
      }
      refExplanationList.add("Scaling row " + 
            Integer.toString(currentRow + 1) + " by a factor of " +
            (scaleBase.getNumberRep()) + " yields: ");
      refStepsList.add(this.complexCopy());
      return this;
   }
   
   
   // *** FULLY TESTED ***
   // determines the transpose of a passed matrix
   // should run in n^2 - n time for square matrices
   // n^2 - n array accesses, BUT n^2 comparisons
   public ComplexMatrix transpose()
   {  
      ComplexNumber[][] auxArray = new ComplexNumber[this.n][this.m];
      //ComplexMatrix auxMatrix = new ComplexMatrix(this.n, this.m);
      ComplexNumber tempValue;
      
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
         ComplexMatrix auxMatrix = new ComplexMatrix(auxArray);
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
      String transposeString = "The transpose of the given matrix is: ";
      df.setRoundingMode(RoundingMode.HALF_UP);
      
      
      transposeString = transposeString + "\n\n";
      for (int i = 0; i < this.m; i++)
      {
         for (int j = 0; j < this.n; j++)
         {
            transposeString =
                  transposeString + this.matrix[i][j].getNumberRep()
                  + MATRIX_SPACING.substring(0, (MATRIX_SPACING.length() - 
                  (this.matrix[i][j].getNumberRep().length()))) + "   ";
            
            /*
            if (this.matrix[i][j] >= 0)
            {
               transposeString = transposeString + " " +
                     df.format(this.matrix[i][j]) + "   ";
                     //Double.toString(transpose.matrix[i][j]) + "   ";
            }
            else
            {
               transposeString = transposeString +
                     df.format(this.matrix[i][j]) + "   ";
                     //Double.toString(transpose.matrix[i][j]) + "   ";
            }
            */
         }
         transposeString = transposeString + "\n";
      }
      return transposeString;
   }
   
   // returns the inverse of a 2x2 matrix, if it is invertible
   public ComplexMatrix inverse2() throws Exception
   {
      if (!this.squareMatrix())
      {
         throw new Exception("A non-square matrix is not invertible.");
      }
      
      ComplexNumber determinant = this.determinant2();
      if (determinant.getRealCoeff() == 0 && determinant.getImagCoeff() == 0 &&
            determinant.getIExponent() == 0)
      {
         // fix this later to make easier for user, but works as placeholder
         throw new Exception("This matrix is not invertible.");
      }
      
      ComplexMatrix inverseMatrix = this.complexCopy();
      inverseMatrix.matrix[0][1] = (inverseMatrix.matrix[0][1]).negate();
      inverseMatrix.matrix[1][0] = (inverseMatrix.matrix[1][0]).negate();
      ComplexNumber temp = inverseMatrix.matrix[0][0];
      inverseMatrix.matrix[0][0] = inverseMatrix.matrix[1][1];
      inverseMatrix.matrix[1][1] = temp;
            
      inverseMatrix = inverseMatrix.multiplyByScalar(
            (new ComplexNumber(1, 0, 0)).divide(determinant));
      return inverseMatrix;
   }
   
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
            inverseString = inverseString + this.matrix[i][j].getNumberRep()
            + MATRIX_SPACING.substring(0, (MATRIX_SPACING.length() - 
                  (this.matrix[i][j].getNumberRep().length()))) + "   ";
            
            /*
            if (this.matrix[i][j] >= 0)
            {
               inverseString = inverseString + " " +
                     df.format(this.matrix[i][j]) + "   ";
                     //Double.toString(inverse.matrix[i][j]) + "   ";
            }
            else
            {
               inverseString = inverseString +
                     df.format(this.matrix[i][j]) + "   ";
                     //Double.toString(inverse.matrix[i][j]) + "   ";
            }
            */
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
    *    (b) A is row equivalent to the n x n identity matrix.
    *    (c) A has n pivot positions.
    *    (d) The equation Ax = 0 has only the trivial solution.
    *    (e) The columns of A form a linearly independent set.
    *    (f) The linear transformation x |-> Ax is one-to-one.
    *    (g) The equation Ax = b has at least one solution for each b in R^n.
    *    (h) The columns of A span R^n.
    *    (i) The linear transformation x |-> Ax maps R^n onto R^n.
    *    (j) There is an n x n matrix C such that CA = I.
    *    (k) There is an n x n matrix D such that AD = I.
    *    (l) A^T is an invertible matrix.
    *    (m) The columns of A form a basis of R^n.
    *    (n) Col A = R^n.
    *    (o) dim Col A = n.
    *    (p) rank A = n.
    *    (q) Nul A = {0}.
    *    (r) dim Nul A = 0.
    * 
    */
   public ComplexMatrix inverseGaussJordan() throws Exception
   {
      if (!this.squareMatrix())
      {
         throw new Exception("A non-square matrix is not invertible.");
      }
      
      int length = matrix.length;
      ComplexNumber[][] newMatrixArray = new ComplexNumber[length][length * 2];
      //Matrix newMatrix = new Matrix(length, length * 2);
      ComplexMatrix sizedIdentity = identityMatrix(length);
      
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
      
      ComplexMatrix newMatrix = new ComplexMatrix(newMatrixArray);
      
      // somehow continue with row reduction.
      // left side has to be "equal" to the original identity matrix
      // therefore create method to test if a matrix is "equal" to another
      
      
      ComplexMatrix newMatrixInverse = newMatrix.rref();
      
      for (int j = 0; j < length; j++)
      {
         if (newMatrixInverse.matrix[j][j].getRealCoeff() == 0 &&
               newMatrixInverse.matrix[j][j].getImagCoeff() == 0 &&
               newMatrixInverse.matrix[j][j].getIExponent() == 0)
         {
            return null;
         }
      }
      
      // this algorithm is really just finding the rref of the doubled matrix,
      // then extracting the non-I portion (the right side) and processing
      // and returning it as the inverse of the original matrix.
      
      //Matrix newMatrixInverse = newMatrix.rref();
      
      // make sure to check here that the rref is actually invertible!
      // this is where it could get complicated if a pivot is lost on the
      //    left-hand side of the matrix!
      
      ComplexNumber[][] inverseGJArray = new ComplexNumber[length][length];
      
      for (int p = 0; p < length; p++)
      {
         for (int q = 0; q < length; q++)
         {
            inverseGJArray[p][q] =
                  newMatrixInverse.matrix[p][q + length];
         }
      }
      ComplexMatrix inverseGJ = new ComplexMatrix(inverseGJArray);
      
      return inverseGJ;
      
      /*
      // enclose the part below in a loop of row operations
      
      boolean fullyReduced = false;
      boolean onesOnDiagonal = true;
      double[][] checkArray = new double[length][length];
      Matrix checkMatrix;
      
      // row reduction needs to occur within this functionality!
      while (!fullyReduced)
      {
         // make sure the left side of the matrix is scaled on the diagonal
         for (int i = 0; i < length; i++)
         {
            if (newMatrix.matrix[i][i] != 1)
            {
               onesOnDiagonal = false;
               break;
            }
         }
         if (!onesOnDiagonal)
         {
            continue;
         }
         
         for (int j = 0; j < length; j++)
         {
            for (int k = 0; k < length; k++)
            {
               checkArray[j][k] = newMatrix.matrix[j][k];
            }
         }
         checkMatrix = new Matrix(checkArray);
         if (checkMatrix.isTriangular())
         {
            fullyReduced = true;
         }
      }
      
      double[][] inverseGJArray = new double[length][length];
      
      for (int p = 0; p < length; p++)
      {
         for (int q = 0; q < length; q++)
         {
            inverseGJArray[p][q] = newMatrix.matrix[p + length][q + length];
         }
      }
      Matrix inverseGJ = new Matrix(inverseGJArray);
      
      return inverseGJ;
      */
   }
   
   
   private ComplexMatrix identityMatrix(int dimension)
   {
      ComplexNumber[][] identityArray = new ComplexNumber[dimension][dimension];
      //ComplexMatrix identityInstance = new ComplexMatrix(dimension, dimension);
      for (int i = 0; i < dimension; i++)
      {
         for (int j = 0; j < dimension; j++)
         {
            if (i == j)
            {
               identityArray[i][j] = new ComplexNumber(1, 0, 0);
            }
            else
            {
               identityArray[i][j] = new ComplexNumber(0, 0, 0);
            }
         }
      }
      ComplexMatrix identityMatrix = new ComplexMatrix(identityArray);
      return identityMatrix;
   }
   
   // Cramer's is really inefficient - rarely used and only on small matrices
   // solves systems of equations (passedMatrix = equalsVector) w/Cramer's rule
   // must be a square matrix - build in exceptions!
   public List<ComplexNumber> cramersRule(ComplexVector equalsVector)
         throws Exception
   {
      if (!this.squareMatrix())
      {
         throw new Exception("A non-square matrix is not invertible.");
      }
      
      int numberOfVariables = this.n; // number of variables
      //ComplexNumber[][] newMatrixArray =
      //      new ComplexNumber[numberOfVariables][numberOfVariables];
      ComplexNumber determinantPassed;
      ComplexNumber determinantCramer;
      //ComplexNumber currentVariable;
      List<ComplexNumber> variableValues = new ArrayList<ComplexNumber>();
      
      if (numberOfVariables == 2)
      {
         determinantPassed = this.determinant2();
      }
      else if (numberOfVariables == 3)
      {
         determinantPassed = this.determinant3();
      }
      else // just for now so testing will work
      {
         determinantPassed = new ComplexNumber(1, 0, 0);
      }
      //else
      //{
      //   determinantPassed = this.determinantFourPlus();
      //}
      
      for (int j = 0; j < numberOfVariables; j++)
      {
         ComplexMatrix newMatrix = this.complexCopy();
         
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
         else // just to get the method to work for testing for now
         {
            determinantCramer = new ComplexNumber(0, 0, 0);
         }
         //else
         //{
         //   determinantCramer = newMatrix.determinantFourPlus();
         //}
         
         ComplexNumber currentVariable =
               (determinantCramer.divide(determinantPassed));
         variableValues.add(currentVariable);
      }
      return variableValues;
   }
   
   public String cramersRuleToString(List<ComplexNumber> variables)
   {
      String cramerOutput = "By Cramer's Rule, the values of the variables "
            + "in the original system are: \n\n";
      
      for (int i = 0; i < variables.size(); i++)
      {
         cramerOutput = cramerOutput + "x" + Integer.toString(i + 1) + " = "
               + (variables.get(i)).getNumberRep() + "\n";
      }
      return cramerOutput;
   }
   
   
   private ComplexNumber findCofactor(int row, int col) throws Exception
   {  
      // technically methods using this should already check for squareness
      if (!this.squareMatrix())
      {
         throw new Exception("Matrix must be square for this operation.");
      }
      
      ComplexNumber cofactor;
      ComplexNumber negativeOneResult;
      ComplexNumber[][] cofactorMatrixArray;
      ComplexMatrix cofactorMatrix;
      
      if (((row + 1) + (col + 1)) % 2 == 0)
      {
         negativeOneResult = new ComplexNumber(1, 0, 0);
      }
      else
      {
         negativeOneResult = new ComplexNumber(-1, 0, 0);
      }
      
      if (this.m == 3)
      {
         cofactorMatrixArray = new ComplexNumber[2][2];
         
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
         
         cofactorMatrix = new ComplexMatrix(cofactorMatrixArray);
         
         cofactor = negativeOneResult.multiply(cofactorMatrix.determinant2());
      }
      else if (this.m == 4)
      {
         cofactorMatrixArray = new ComplexNumber[3][3];
         
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
         
         cofactorMatrix = new ComplexMatrix(cofactorMatrixArray);
         
         cofactor = negativeOneResult.multiply(cofactorMatrix.determinant3());
         
         /*
         for (int p = 0; p < cofactorMatrix.m; p++)
         {
            for (int q = 0; q < cofactorMatrix.n; q++)
            {
               System.out.println(Double.toString(cofactorMatrix.matrix[p][q]));
            }
         }
         
         System.out.println("\n");
         */
         
         /*
         System.out.println(Double.toString(passedMatrix.matrix[row][col]));
         System.out.println(Integer.toString(negativeOneResult));
         System.out.println(Double.toString(determinant3(cofactorMatrix)));
         System.out.println("\n");
         */
         
         /*
         System.out.println("Cofactor: " + Double.toString(cofactor));
         System.out.println("\n");
         */
         
      }
      else if (this.m > 4) // JUST FOR TESTING FOR NOW!
      {
         cofactor = new ComplexNumber(1, 0, 0);
      }
      /*
      else if (this.m > 4)
      {
         int dimension = this.m;
         
         cofactorMatrixArray = 
               new ComplexNumber[dimension - 1][dimension - 1];
         
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
         
         cofactorMatrix = new ComplexMatrix(cofactorMatrixArray);
         //if (this.numberOfZeros == 0)
         //{
         //   cofactorMatrix.setNumberOfZeros(0);
         //}
         
         cofactor = negativeOneResult.multiply(
               cofactorMatrix.determinantFourPlus());
      }
      */
      else
      {
         throw new Exception("Passed matrix is not large enough.");
      }
      
      return cofactor;
   }
   
   private ComplexMatrix findAdjugate() throws Exception
   {
      ComplexNumber[][] adjugateArray = new ComplexNumber[this.m][this.n];
      
      for (int i = 0; i < this.m; i++)
      {
         for (int j = 0; j < this.n; j++)
         {
            adjugateArray[i][j] = this.findCofactor(i, j);
         }
      }
      
      ComplexMatrix adjugateMatrix = new ComplexMatrix(adjugateArray);
      adjugateMatrix = adjugateMatrix.transpose();
      return adjugateMatrix;
   }
   
   // finds a matrix's inverse through cofactors and the adjugate
   public ComplexMatrix inverseCramer() throws Exception
   {
      if (!this.squareMatrix())
      {
         throw new Exception("A non-square matrix is not invertible.");
      }

      ComplexMatrix adjugate = this.findAdjugate();
      ComplexMatrix inverseMatrix;
      
      if (this.m == 2)
      {
         inverseMatrix = adjugate.multiplyByScalar(
               (new ComplexNumber(1, 0, 0)).divide(this.determinant2()));
      }
      else if (this.m == 3)
      {
         inverseMatrix = adjugate.multiplyByScalar(
               (new ComplexNumber(1, 0, 0)).divide(this.determinant3()));
      }
      else // JUST FOR TESTING FOR NOW!
      {
         inverseMatrix = null;
      }
      /*
      else
      {
         inverseMatrix = adjugate.multiplyByScalar( 
               (new ComplexNumber(1, 0, 0)).divide(this.determinantFourPlus()));
      }
      */
      return inverseMatrix;
   }
   
   
   // finds the determinant of a passed 2x2 matrix - BUILD EXCEPTIONS!
   public ComplexNumber determinant2() throws Exception
   {
      if (!this.squareMatrix())
      {
         throw new Exception("A non-square matrix is not invertible.");
      }
      
      ComplexNumber determinant = //this.coefficientArray[2];
            ((this.matrix[0][0]).multiply(this.matrix[1][1])).subtract(
                  ((this.matrix[0][1].multiply(this.matrix[1][0]))));
      return determinant;
   }
   
   // finds the determinant of a passed 3x3 matrix - BUILD EXCEPTIONS!
   public ComplexNumber determinant3() throws Exception
   {
      if (!this.squareMatrix())
      {
         throw new Exception("A non-square matrix is not invertible.");
      }
      
      ComplexNumber determinant = //this.coefficientArray[3];
            ((((this.matrix[0][0].multiply(this.matrix[1][1])).multiply(
                  this.matrix[2][2])).add((this.matrix[0][1].multiply(
                  this.matrix[1][2])).multiply(this.matrix[2][0]))).add(
                  (((this.matrix[0][2].multiply(this.matrix[1][0]))
                  .multiply(this.matrix[2][1])).subtract((
                  this.matrix[0][0].multiply(this.matrix[1][2]))
                  .multiply(this.matrix[2][1]))))).add(
                  ((this.matrix[0][1].multiply(this.matrix[1][0]))
                  .multiply(this.matrix[2][2])).negate().subtract(
                  (this.matrix[0][2].multiply(
                  this.matrix[1][1])).multiply(this.matrix[2][0])));
      
      return determinant;
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
               determinantString + (this.determinant2()).getNumberRep();
      }
      else if (this.m == 3)
      {
         determinantString =
               determinantString + (this.determinant3()).getNumberRep();
      }
      //else
      //{
      //   determinantString = determinantString +
      //         (this.determinantFourPlus()).getNumberRep();
      //}
      return determinantString;
   }
   
   
   public void findEigenvalues() throws Exception
   {
      if (!this.squareMatrix())
      {
         throw new Exception("A non-square matrix is not invertible.");
      }
      
      /*
      if (this.isTriangular())
      {
         for (int i = 0; i < this.m; i++)
         {
            this.complexEigenvalueList.add(this.matrix[i][i]);
         }
         Collections.sort(this.complexEigenvalueList);
         // WRITE A COMPARETO METHOD FOR COMPLEXNUMBER
         return;
      }
      */
      
      if (this.m == 2)
      {
         findCMEigenvaluesQuadratic();
         /*
         if (findTypeEigenvaluesQuadratic(this.calculateDiscriminant()))
         {
            this.findRealEigenvaluesQuadratic();
         }
         else
         {
            this.findComplexEigenvaluesQuadratic();
         }
         */
      }
      /*
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
      */
      //Collections.sort(this.complexEigenvalueList);
   }
   
   
   
   /*
   public List<ComplexNumber> complexMatrixEValues()
   {
      if (!this.squareMatrix())
      {
         throw new Exception("A non-square matrix is not invertible.");
      }
      
      //String eigenDeterminant;
      
      ComplexMatrix eigenCMatrix = 
            new ComplexMatrix(this.m, this.n);
      
      // creates Ax - LAMBDAx matrix
      for (int i = 0; i < this.m; i++)
      {
         for (int j = 0; j < this.n; j++)
         {
            eigenCMatrix.matrix[i][j] = this.matrix[i][j];
         }
      }
      
      ComplexNumber[] coefficientArray = new ComplexNumber[this.m + 1];
      List<ComplexNumber> eigenvalueList = new ArrayList<ComplexNumber>();
      
      if (this.m == 2)
      {
         //eigenDeterminant = eigenDeterminant2(eigenMatrix);
         coefficientArray[0] = new ComplexNumber(1, 0, 0);
         coefficientArray[1] = // -(a + d)
               ((eigenCMatrix.matrix[0][0]).add(
                     eigenCMatrix.matrix[1][1])).negate();
         coefficientArray[2] = // (ad - bc)
               ComplexNumber.subtract(
                     (ComplexNumber.multiply(eigenCMatrix.matrix[0][0],
                           eigenCMatrix.matrix[1][1])),
                     (ComplexNumber.multiply(eigenCMatrix.matrix[0][1],
                           eigenCMatrix.matrix[1][0])));
         eigenvalueList = findComplexRootsQuadratic(coefficientArray);
         // return eigenvalueList;
      }
      else if (this.m == 3)
      {
         //eigenDeterminant = eigenDeterminant3(eigenMatrix);
         coefficientArray[0] = new ComplexNumber(-1, 0, 0);
         coefficientArray[1] = // a + e + i
               ComplexNumber.add(eigenCMatrix.matrix[0][0],
                     (ComplexNumber.add(eigenCMatrix.matrix[1][1],
                           eigenCMatrix.matrix[2][2])));
         
         // cg + fh + bd - ae - ai - ei        
         ComplexNumber cgfhbd =
               ComplexNumber.add(ComplexNumber.add(
                     (ComplexNumber.multiply(eigenCMatrix.matrix[0][2],
                           eigenCMatrix.matrix[2][0])), // cg
                     (ComplexNumber.multiply(eigenCMatrix.matrix[1][2],
                           eigenCMatrix.matrix[2][1]))), // fh
                           (ComplexNumber.multiply(eigenCMatrix.matrix[0][1],
                           eigenCMatrix.matrix[1][0]))); // bd
         
         ComplexNumber aeaiei = 
               ComplexNumber.subtract(ComplexNumber.subtract(
                     (ComplexNumber.multiply(eigenCMatrix.matrix[0][0],
                           eigenCMatrix.matrix[1][1])), // ae
                     (ComplexNumber.multiply(eigenCMatrix.matrix[0][0],
                           eigenCMatrix.matrix[2][2]))), // ai
                     (ComplexNumber.multiply(eigenCMatrix.matrix[1][1],
                           eigenCMatrix.matrix[2][2]))); // ei
         
         coefficientArray[2] = 
               ComplexNumber.subtract(cgfhbd, aeaiei);
         
         coefficientArray[3] = // aei + bfg + cdh - afh - bdi - ceg
         
               
               
               
               
               
               
               
               
               
               
               (passedMatrix.matrix[0][0] * passedMatrix.matrix[1][1] *
                     passedMatrix.matrix[2][2]) +
               (passedMatrix.matrix[0][1] * passedMatrix.matrix[1][2] *
                     passedMatrix.matrix[2][0]) +
               (passedMatrix.matrix[0][2] * passedMatrix.matrix[1][0] *
                     passedMatrix.matrix[2][1]) -
               (passedMatrix.matrix[0][0] * passedMatrix.matrix[1][2] *
                     passedMatrix.matrix[2][1]) -
               (passedMatrix.matrix[0][1] * passedMatrix.matrix[1][0] *
                     passedMatrix.matrix[2][2]) -
               (passedMatrix.matrix[0][2] * passedMatrix.matrix[1][1] *
                     passedMatrix.matrix[2][0]);
         
         
         
         
         eigenvalueList = findComplexRootsCubic(coefficientArray);
         //return eigenvalueList;
      }
      return eigenvalueList;
      
      //else
      //{
      //   // figure out these coefficients later - requires cofactors
      //   eigenDeterminant = eigenDeterminant4Plus(eigenMatrix);
      //}
      
   }
   */
   
   
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
         //Double.toString(eV);
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
   
   
   // finds complex roots for a quadratic (for eigenvalues for a 2x2 matrix)
   private void findCMEigenvaluesQuadratic()
   {  
      ComplexNumber a = this.coefficientArray[0];
      ComplexNumber b = this.coefficientArray[1];
      ComplexNumber c = this.coefficientArray[2];
      
      ComplexNumber discriminant =
            (b.multiply(b)).subtract((a.multiply(c.multiply(
                        new ComplexNumber(4, 0, 0)))));
      
      ComplexNumber rootComplexPositive;
      ComplexNumber rootComplexNegative;
      //List<ComplexNumber> rootList = new ArrayList<ComplexNumber>();
      
      rootComplexPositive =
            (b.negate()).add(discriminant.complexSqrt().get(0)).divide(
                  (a.multiply(new ComplexNumber(2, 0, 0))));
      
      rootComplexNegative =
            (b.negate()).add(discriminant.complexSqrt().get(1)).divide(
                  (a.multiply(new ComplexNumber(2, 0, 0))));
      
      this.complexEigenvalueList.add(rootComplexPositive);
      this.complexEigenvalueList.add(rootComplexNegative);
   }
   
   public static void main(String[] args) throws Exception
   {
      /*
      ComplexNumber[][] complexExample = new ComplexNumber[2][2];
      complexExample[0][0] = new ComplexNumber(1, -2, 1);
      complexExample[0][1] = new ComplexNumber(2, 0, 1); // both double and complex?
      complexExample[1][0] = new ComplexNumber(-1, 0, 1);
      complexExample[1][1] = new ComplexNumber(3, 1, 1);
      ComplexMatrix complexMatrix = new ComplexMatrix(complexExample);
      */
      
      /*
      ComplexNumber[][] complexRefArray = new ComplexNumber[3][3];
      complexRefArray[0][0] = new ComplexNumber(4, 1, 1);
      complexRefArray[0][1] = new ComplexNumber(6, 0, 0);
      complexRefArray[0][2] = new ComplexNumber(-2, 0, 0);
      complexRefArray[1][0] = new ComplexNumber(1, 1, 1);
      complexRefArray[1][1] = new ComplexNumber(3, -1, 1);
      complexRefArray[1][2] = new ComplexNumber(-1, 0, 0);
      complexRefArray[2][0] = new ComplexNumber(-2, 0, 0);
      complexRefArray[2][1] = new ComplexNumber(1, -3, 1);
      complexRefArray[2][2] = new ComplexNumber(2, -1, 1);
      ComplexMatrix complexRefMatrix = new ComplexMatrix(complexRefArray);
      
      //ComplexMatrix testRREF = complexRefMatrix.rref();
      //System.out.println(testRREF.rRefToString());
      
      ComplexMatrix refIGJMatrix = complexRefMatrix.inverseGaussJordan();
      System.out.println(refIGJMatrix.inverseToString());
      //*/
      
      /*
      ComplexNumber[][] complexDeterminant2Array = new ComplexNumber[2][2];
      complexDeterminant2Array[0][0] = new ComplexNumber(3, -1, 1);
      complexDeterminant2Array[0][1] = new ComplexNumber(2, 2, 1);
      complexDeterminant2Array[1][0] = new ComplexNumber(4, -6, 1);
      complexDeterminant2Array[1][1] = new ComplexNumber(5, 1, 1);
      ComplexMatrix complexDeterminant2Matrix = new ComplexMatrix(complexDeterminant2Array);
      
      //double testDeterminant = complexDeterminant2Matrix.determinant2();
      System.out.println(complexDeterminant2Matrix.determinantToString());
      */
      
      /*
      ComplexNumber[][] complexRefArray = new ComplexNumber[3][3];
      complexRefArray[0][0] = new ComplexNumber(4, 1, 1);
      complexRefArray[0][1] = new ComplexNumber(6, 0, 0);
      complexRefArray[0][2] = new ComplexNumber(-2, 0, 0);
      complexRefArray[1][0] = new ComplexNumber(1, 1, 1);
      complexRefArray[1][1] = new ComplexNumber(3, -1, 1);
      complexRefArray[1][2] = new ComplexNumber(-1, 0, 0);
      complexRefArray[2][0] = new ComplexNumber(-2, 0, 0);
      complexRefArray[2][1] = new ComplexNumber(1, -3, 1);
      complexRefArray[2][2] = new ComplexNumber(2, -1, 1);
      ComplexMatrix complexRefMatrix = new ComplexMatrix(complexRefArray);
      
      System.out.println(complexRefMatrix.determinantToString());
      */
      
      /*
      ComplexNumber[][] complexDeterminant2Array = new ComplexNumber[2][2];
      complexDeterminant2Array[0][0] = new ComplexNumber(3, -1, 1);
      complexDeterminant2Array[0][1] = new ComplexNumber(2, 2, 1);
      complexDeterminant2Array[1][0] = new ComplexNumber(4, -6, 1);
      complexDeterminant2Array[1][1] = new ComplexNumber(5, 1, 1);
      ComplexMatrix complexDeterminant2Matrix = new ComplexMatrix(complexDeterminant2Array);
      
      complexDeterminant2Matrix.findEigenvalues();
      System.out.println(complexDeterminant2Matrix.eigenvaluesToString());
      */
      
      /*
      ComplexNumber[][] complexRefArray = new ComplexNumber[3][3];
      complexRefArray[0][0] = new ComplexNumber(4, 1, 1);
      complexRefArray[0][1] = new ComplexNumber(6, 0, 0);
      complexRefArray[0][2] = new ComplexNumber(-2, 0, 0);
      complexRefArray[1][0] = new ComplexNumber(1, 1, 1);
      complexRefArray[1][1] = new ComplexNumber(3, -1, 1);
      complexRefArray[1][2] = new ComplexNumber(-1, 0, 0);
      complexRefArray[2][0] = new ComplexNumber(-2, 0, 0);
      complexRefArray[2][1] = new ComplexNumber(1, -3, 1);
      complexRefArray[2][2] = new ComplexNumber(2, -1, 1);
      ComplexMatrix complexRefMatrix = new ComplexMatrix(complexRefArray);
      
      ComplexMatrix complexTranspose = complexRefMatrix.transpose();
      System.out.println(complexTranspose.transposeToString());
      */
      
      /*
      ComplexNumber[][] complexDeterminant2Array = new ComplexNumber[2][2];
      complexDeterminant2Array[0][0] = new ComplexNumber(3, -1, 1);
      complexDeterminant2Array[0][1] = new ComplexNumber(2, 2, 1);
      complexDeterminant2Array[1][0] = new ComplexNumber(4, -6, 1);
      complexDeterminant2Array[1][1] = new ComplexNumber(5, 1, 1);
      ComplexMatrix complexDeterminant2Matrix = new ComplexMatrix(complexDeterminant2Array);
      
      ComplexMatrix complexInverse2 = complexDeterminant2Matrix.inverse2();
      System.out.println(complexInverse2.inverseToString());
      */
      
      /*
      ComplexNumber[][] complexDeterminant2Array = new ComplexNumber[2][2];
      complexDeterminant2Array[0][0] = new ComplexNumber(3, -1, 1);
      complexDeterminant2Array[0][1] = new ComplexNumber(2, 2, 1);
      complexDeterminant2Array[1][0] = new ComplexNumber(4, -6, 1);
      complexDeterminant2Array[1][1] = new ComplexNumber(5, 1, 1);
      ComplexMatrix complexDeterminant2Matrix = new ComplexMatrix(complexDeterminant2Array);
      ComplexNumber[] vectorArray = new ComplexNumber[2];
      vectorArray[0] = new ComplexNumber(1, 1, 1);
      vectorArray[1] = new ComplexNumber(-2, -1, 1);
      ComplexVector testVector = new ComplexVector(vectorArray);
      
      List<ComplexNumber> testList = complexDeterminant2Matrix.cramersRule(testVector);
      System.out.println(complexDeterminant2Matrix.cramersRuleToString(testList));
      */
      
      /*
      ComplexNumber[][] complexRefArray = new ComplexNumber[3][3];
      complexRefArray[0][0] = new ComplexNumber(4, 1, 1);
      complexRefArray[0][1] = new ComplexNumber(6, 0, 0);
      complexRefArray[0][2] = new ComplexNumber(-2, 0, 0);
      complexRefArray[1][0] = new ComplexNumber(1, 1, 1);
      complexRefArray[1][1] = new ComplexNumber(3, -1, 1);
      complexRefArray[1][2] = new ComplexNumber(-1, 0, 0);
      complexRefArray[2][0] = new ComplexNumber(-2, 0, 0);
      complexRefArray[2][1] = new ComplexNumber(1, -3, 1);
      complexRefArray[2][2] = new ComplexNumber(2, -1, 1);
      ComplexMatrix complexRefMatrix = new ComplexMatrix(complexRefArray);
      
      ComplexMatrix refCramerMatrix = complexRefMatrix.inverseCramer();
      System.out.println(refCramerMatrix.inverseToString());
      */
      
      // M > N
      ComplexNumber[][] firstMatrixArray = new ComplexNumber[3][3];
      ComplexNumber[][] secondMatrixArray = new ComplexNumber[3][2];
      firstMatrixArray[0][0] = new ComplexNumber(4, -1, 1);
      firstMatrixArray[0][1] = new ComplexNumber(-2, 0, 0);
      firstMatrixArray[0][2] = new ComplexNumber(-1, 0, 0);
      firstMatrixArray[1][0] = new ComplexNumber(3, -6, 1);
      firstMatrixArray[1][1] = new ComplexNumber(8, 0, 0);
      firstMatrixArray[1][2] = new ComplexNumber(0, -1, 1);
      firstMatrixArray[2][0] = new ComplexNumber(4, 0, 0);
      firstMatrixArray[2][1] = new ComplexNumber(3, 1, 1);
      firstMatrixArray[2][2] = new ComplexNumber(1, 2, 1);
      secondMatrixArray[0][0] = new ComplexNumber(7, 0, 0);
      secondMatrixArray[0][1] = new ComplexNumber(8, -2, 1);
      secondMatrixArray[1][0] = new ComplexNumber(4, 1, 1);
      secondMatrixArray[1][1] = new ComplexNumber(3, 0, 0);
      secondMatrixArray[2][0] = new ComplexNumber(1, -1, 1);
      secondMatrixArray[2][1] = new ComplexNumber(-2, 0, 0);
      
      ComplexMatrix firstMatrix = new ComplexMatrix(firstMatrixArray);
      ComplexMatrix secondMatrix = new ComplexMatrix(secondMatrixArray);
      
      ComplexMatrix newMatrix = firstMatrix.multiplyMatrices(secondMatrix);
      System.out.println(newMatrix.multiplyToString());
   }
}