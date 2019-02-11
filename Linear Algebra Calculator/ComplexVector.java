/****************************************
 * Linear Algebra Calculator
 * 
 * Supports representation of complex numbers, matrix multiplication
 *    and addition, row reduction, transposing, methods for finding
 *    determinants, several methods for inverses, eigenvalues, eigenvectors,
 *    diagonalization, and operations related to orthogonality.
 * 
 * 
 * ComplexVector.java
 * 
 * API:
 *    public ComplexVector(int rows)
 *    public ComplexVector(ComplexNumber[] values)
 *    public ComplexVector copy()
 *    public ComplexNumber getValue(int row)
 *
 *    
 * Dependencies:
 *    -ComplexNumber.java
 * 
 ****************************************/


public class ComplexVector
{
   private ComplexNumber[] vector;
   
   public ComplexVector(int rows)
   {
      vector = new ComplexNumber[rows];
   }
      
   public ComplexVector(ComplexNumber[] values)
   {
      this.vector = new ComplexNumber[values.length];
      
      for (int i = 0; i < values.length; i++)
      {
         vector[i] = values[i];
      }
   }
   
   public ComplexVector copy()
   {
      ComplexNumber[] copyArray = new ComplexNumber[this.vector.length];
      for (int i = 0; i < this.vector.length; i++)
      {
         copyArray[i] = this.vector[i];
      }
      ComplexVector vectorCopy = new ComplexVector(copyArray);
      return vectorCopy;
   }
   
   public ComplexNumber getValue(int row)
   {
      return vector[row];
   }
   
}