/****************************************
 * Linear Algebra Calculator
 * 
 * Supports representation of complex numbers, matrix multiplication
 *    and addition, row reduction, transposing, methods for finding
 *    determinants, several methods for inverses, eigenvalues, eigenvectors,
 *    diagonalization, and operations related to orthogonality.
 * 
 * 
 * MatrixZeros.java
 * 
 * API:
 *    public MatrixZeros(int row, int col)
 *    public int getRow()
 *    public void setRow(int row)
 *    public int getCol()
 *    public void setCol(int row)
 *
 *    
 * Dependencies: None
 * 
 ****************************************/

public class MatrixZeros
{
   private int row;
   private int col;
   
   public MatrixZeros(int row, int col)
   {
      this.row = row;
      this.col = col;
   }
   
   public int getRow()
   {
      return row;
   }
   
   public void setRow(int row)
   {
      this.row = row;
   }
   
   public int getCol()
   {
      return col;
   }
   
   public void setCol(int col)
   {
      this.col = col;
   }
}