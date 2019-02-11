/****************************************
 * Linear Algebra Calculator
 * 
 * Supports representation of complex numbers, matrix multiplication
 *    and addition, row reduction, transposing, methods for finding
 *    determinants, several methods for inverses, eigenvalues, eigenvectors,
 *    diagonalization, and operations related to orthogonality.
 * 
 * 
 * ComplexNumber.java
 * 
 * API:
 *    public ComplexNumber(double real, double imaginary, double exponent)
 *    public int compareTo(ComplexNumber complex2)
 *    public boolean equals(Object y)
 *    public ComplexNumber multiplyByScalar(double scalar)
 *    public ComplexNumber multiply(ComplexNumber complex2)
 *    public ComplexNumber divide(ComplexNumber complex2)
 *    public ComplexNumber add(ComplexNumber complex2)
 *    public ComplexNumber subtract(ComplexNumber complex2)
 *    public ComplexNumber negate()
 *    public ComplexNumber getConjugate()
 *    public List<ComplexNumber> complexSqrt()
 *    public double getRealCoeff()
 *    public double getImagCoeff()
 *    public double getIExponent()
 *    public void setRealCoeff(double realCoeff)
 *    public void setImagCoeff(double imagCoeff)
 *    public String getNumberRep()
 *    private double convertComplexToDouble()
 *
 *    
 * Dependencies: None
 * 
 ****************************************/


import java.util.List;
import java.util.ArrayList;
import java.text.DecimalFormat;
import java.math.RoundingMode;

public final class ComplexNumber
{
   private double realCoeff;
   private double imagCoeff;
   private double iExponent;
   private static DecimalFormat df = new DecimalFormat("#.#####");
   
   public ComplexNumber(double real, double imaginary, double exponent)
   {
      realCoeff = real;
      imagCoeff = imaginary;
      iExponent = exponent;
   }
   
   public int compareTo(ComplexNumber complex2)
   {
      if (this.realCoeff > complex2.realCoeff)
      {
         return 1;
      }
      if (this.realCoeff < complex2.realCoeff)
      {
         return -1;
      }
      if (this.imagCoeff > complex2.imagCoeff)
      {
         return 1;
      }
      if (this.imagCoeff < complex2.imagCoeff)
      {
         return -1;
      }
      return 0;
   }
   
   /*
   public boolean equals(Object y)
   {
      if (y == this)
      {
         return true;
      }
      if (y == null)
      {
         return false;
      }
      if (y.getClass() != this.getClass())
      {
         return false;
      }
      
      ComplexNumber that = (ComplexNumber) y;
      if (Double.compare(this.realCoeff, that.realCoeff) != 0)
      {
         return false;
      }
      if (Double.compare(this.imagCoeff, that.imagCoeff) != 0)
      {
         return false;
      }
      if (Double.compare(this.iExponent, that.iExponent) != 0)
      {
         return false;
      }
      return true;
   }
   */
   
   public ComplexNumber multiplyByScalar(double scalar)
   {
      ComplexNumber complexScalar = new ComplexNumber(scalar, 0, 0);
      ComplexNumber result = this.multiply(complexScalar);
      return result;
   }
   
   public ComplexNumber multiply(ComplexNumber complex2)
   {
      double newReal;
      double newImag;
      
      newReal =
            (this.realCoeff * complex2.realCoeff) -
            (this.imagCoeff * complex2.imagCoeff);
      newImag =
            ((this.realCoeff * complex2.imagCoeff) +
            (this.imagCoeff * complex2.realCoeff));
      
      if (newReal == 0 && newImag == 0)
      {
         return new ComplexNumber(0, 0, 0);
      }
      else
      {
         return new ComplexNumber(newReal, newImag, 1);
      }
      //ComplexNumber newComplex = new ComplexNumber(newReal, newImag, 1);
      //return newComplex;
   }
   
   public ComplexNumber divide(ComplexNumber complex2)
   {
      ComplexNumber quotient;
      ComplexNumber numerator = this.multiply(complex2.getConjugate());
      ComplexNumber denominator = complex2.multiply(complex2.getConjugate());
      
      /*
      if (getRealCoeff(numerator) == -0.0 && )
      {
         numerator = new ComplexNumber(-(numerator.realCoeff), numerator.imagCoeff, 1);
      }
      */
      
      if (numerator.getRealCoeff() == 0 && numerator.getImagCoeff() == 0)
      {
         quotient = new ComplexNumber(0, 0, 0);
      }
      else
      {
         // need to include control for when denominator is zero!!!
         double newReal = numerator.realCoeff / denominator.realCoeff;
         double newImag = numerator.imagCoeff / denominator.realCoeff;
         quotient = new ComplexNumber(newReal, newImag, 1);
      }
      return quotient;
   }
   
   public ComplexNumber add(ComplexNumber complex2)
   {
      // parse first number and find real component
      // parse second number and find real component
      // add real components
      // parse first number and find imagCoefficient
      // parse second number and find imagCoefficient
      // add imaginary components
      // maybe add a catch for if the imagCoefficient becomes 0?
      
      double newReal =
            (this.realCoeff + complex2.realCoeff);
      double newImag =
            (this.imagCoeff + complex2.imagCoeff);
      
      if (newReal == 0 && newImag == 0)
      {
         return new ComplexNumber(0, 0, 0);
      }
      else
      {
         return new ComplexNumber(newReal, newImag, 1);
      }
      //ComplexNumber newComplex = new ComplexNumber(newReal, newImag, 1);
      //return newComplex;
   }
   
   public ComplexNumber subtract(ComplexNumber complex2)
   {
      double newReal =
            (this.realCoeff - complex2.realCoeff);
      double newImag =
            (this.imagCoeff - complex2.imagCoeff);
      
      if (newReal == 0 && newImag == 0)
      {
         return new ComplexNumber(0, 0, 0);
      }
      else
      {
         return new ComplexNumber(newReal, newImag, 1);
      }
      
      //ComplexNumber newComplex = new ComplexNumber(newReal, newImag, 1);
      //return newComplex;
   }
   
   public ComplexNumber negate()
   {
      if (this.realCoeff != 0.0)
      {
         this.realCoeff = -(this.realCoeff);
      }
      if (this.imagCoeff != 0.0)
      {
         this.imagCoeff = -(this.imagCoeff);
      }
      return this;
   }
   
   public ComplexNumber getConjugate()
   {
      if (this.imagCoeff == 0)
      {
         return this;
      }
      else
      {
         double newReal = this.realCoeff;
         double newImag = -(this.imagCoeff);
         ComplexNumber conjugate = new ComplexNumber(newReal, newImag, 1);
         return conjugate;
      }
   }
   
   public List<ComplexNumber> complexSqrt()
   {
      List<ComplexNumber> sqrt = new ArrayList<ComplexNumber>();
      
      double a = this.realCoeff;
      double b = this.imagCoeff;
      
      double xPositive;
      double xNegative;
      double yPositive;
      double yNegative;
      
      // x^2 - y^2 = a
      // 2xy = b
      // x^2 + y^2 = modulus = sqrt(a^2 + b^2)
      
      double modulus = Math.sqrt((a * a) + (b * b));
      yPositive = Math.sqrt((modulus - a) / 2);
      yNegative = -yPositive;
      
      if (b > 0)
      {
         xPositive = b / (2 * yPositive);
         xNegative = -xPositive;
         
         sqrt.add(new ComplexNumber(xPositive, yPositive, 1));
         sqrt.add(new ComplexNumber(xNegative, yNegative, 1));
      }
      else if (b < 0)
      {
         xPositive = b / (2 * yNegative);
         xNegative = -xPositive;
         
         sqrt.add(new ComplexNumber(xPositive, yNegative, 1));
         sqrt.add(new ComplexNumber(xNegative, yPositive, 1));
      }
      else
      {
         if (a >= 0)
         {
            double oneRealCoeff = Math.sqrt(a);
            sqrt.add(new ComplexNumber(oneRealCoeff, 0, 0));
         }
         else
         {
            double oneImagCoeff = Math.sqrt(-a);
            sqrt.add(new ComplexNumber(0, oneImagCoeff, 1));
         }
      }
      
      return sqrt;
   }
   
   public double getRealCoeff()
   {
      return this.realCoeff;
   }
   
   public double getImagCoeff()
   {
      return this.imagCoeff;
   }
   
   public double getIExponent()
   {
      return this.iExponent;
   }
   
   // accessibility of mutators?
   public void setRealCoeff(double realCoeff)
   {
      this.realCoeff = realCoeff;
   }
   
   public void setImagCoeff(double imagCoeff)
   {
      this.imagCoeff = imagCoeff;
   }
   
   public String getNumberRep()
   {
      df.setRoundingMode(RoundingMode.HALF_UP);
      
      String roundedReal = df.format(this.realCoeff);
      String roundedImag = df.format(this.imagCoeff);
      String numberRepresentation;
      
      if (roundedReal.equals("0") && roundedImag.equals("0"))
      {
         numberRepresentation = "0";
      }
      else if (roundedReal.equals("0") && !roundedImag.equals("0"))
      {
         if (roundedImag.equals("1"))
         {
            numberRepresentation = "i";
         }
         else if (roundedImag.equals("-1"))
         {
            numberRepresentation = "-i";
         }
         else
         {
            numberRepresentation = roundedImag + "i";
         }
      }
      else if (roundedImag.equals("0"))
      {
         numberRepresentation = roundedReal;
      }
      else if (roundedImag.equals("1") || roundedImag.equals("-1"))
      {
         if (roundedImag.equals("1"))
         {
            numberRepresentation = roundedReal + " + " + "i";
         }
         else
         {
            numberRepresentation = roundedReal + " - " + "i";
         }
      }
      else
      {
         if (roundedImag.substring(0, 1).equals("-"))
         {
            numberRepresentation = roundedReal + " - " + df.format(-this.imagCoeff) + "i";
         }
         else
         {
            numberRepresentation = roundedReal + " + " + roundedImag + "i";
         }
      }
      
      return numberRepresentation;
   }
   
   private double convertComplexToDouble()
         throws Exception
   {
      if (this.imagCoeff == 0)
      {
         return this.realCoeff;
      }
      else
      {
         throw new Exception("Passed number is still complex.");
      }
   }
}