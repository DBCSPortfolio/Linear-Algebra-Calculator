/****************************************
 * Linear Algebra Calculator
 * 
 * Supports representation of complex numbers, matrix multiplication
 *    and addition, row reduction, transposing, methods for finding
 *    determinants, several methods for inverses, eigenvalues, eigenvectors,
 *    diagonalization, and operations related to orthogonality.
 * 
 * 
 * Vector.java
 * 
 * API:
 *    public Vector(int rows)
 *    public Vector(double[] values)
 *    public Vector copy()
 *    public Vector multiplyByScalar(double scalar)
 *    public double getValue(int row)
 *    public int getDimension()
 *    public double getLength()
 *    public Vector normalize()
 *    public void setValue(double value, int entry)
 *    public String vectorToString()
 *
 *    
 * Dependencies: None
 * 
 ****************************************/

import java.text.DecimalFormat;

// special subclass for vectors (matrices with one column)
public class Vector
{
   private double[]              vector;
   
   private static DecimalFormat  df = new DecimalFormat("#.#####");
   
   public Vector(int rows)
   {
      vector = new double[rows];
   }
      
   public Vector(double[] values)
   {
      this.vector = new double[values.length];
      
      for (int i = 0; i < values.length; i++)
      {
         vector[i] = values[i];
      }
   }
   
   public Vector copy()
   {
      double[] copyArray = new double[this.vector.length];
      for (int i = 0; i < this.vector.length; i++)
      {
         copyArray[i] = this.vector[i];
      }
      Vector vectorCopy = new Vector(copyArray);
      return vectorCopy;
   }
   
   public Vector multiplyByScalar(double scalar)
   {
      for (int i = 0; i < this.getDimension(); i++)
      {
         this.vector[i] = this.vector[i] * scalar;
      }
      return this;
   }
   
   public double getValue(int row)
   {
      return vector[row];
   }
   
   public int getDimension()
   {
      return vector.length;
   }
   
   // returns the length of a passed vector
   public double getLength()
   {
      double length = 0;
      for (int i = 0; i < this.getDimension(); i++)
      {
         length = length + (this.getValue(i) * this.getValue(i));
      }
      length = Math.sqrt(length);
      return length;
   }
   
   // normalizes a passed vector (makes it a unit vector - length 1)
   public Vector normalize()
   {
      double length = this.getLength();
      Vector normalizedVector =
            this.multiplyByScalar((1 / length));
      return normalizedVector;
   }
   
   public void setValue(double value, int entry)// throws Exception
   {
      //if (!validateValue(value))
      //{
      //   throw new Exception("Invalid value passed.");
      //}
      this.vector[entry] = value;
   }
   
   //private boolean validateValue(double value)
   //{
   //   
   //}
   
   public String vectorToString()
   {
      String vectorString = "";
      int dimension = this.getDimension();
      
      for (int i = 0; i < dimension; i++)
      {
         vectorString = vectorString + df.format(this.vector[i]) + "\n";
      }
      return vectorString;
   }
}