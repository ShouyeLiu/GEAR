package score;

import weka.core.matrix.Matrix;
/**
 * This class implements AbstractScore.  It is utilized to calculate score with linear regression when without adjustment
 * @author Guo-Bo Chen, Zhejiang University, gc5k@zju.edu.cn
 *
 */
public class LinearMeanScore implements AbstractScore {
	private double[] Residual;//residual score
	private Matrix Response;//Y

	/**
	 * Construct the object
	 * @param Y Response
	 */
	public LinearMeanScore (double[][] Y)
	{
		Response = new Matrix(Y);
		Residual = new double[Response.getRowDimension()];
	}

	/**
	 * Get the coefficients of the linear regression model constructed.
	 */
	public double[] getCoefficients()
	{
		double mean[]=new double[Response.getColumnDimension()];
		double mu[] = new double[1];
		for( int i=0; i< Response.getColumnDimension(); i++)
		{
			for( int j=0; j < Response.getRowDimension(); j++ )
			{
				mean[i] += Response.get(j, i);
			}
			mean[i]/= Response.getRowDimension();
		}
		mu[0]=mean[0];
		return mu;
	}

	/**
	 * Calculate score by a mean model.
	 */
	public void CalculateScore() {
		double mean[]=new double[Response.getColumnDimension()];
		for( int i=0; i< Response.getColumnDimension(); i++)
		{
			for( int j=0; j < Response.getRowDimension(); j++ )
			{
				mean[i] += Response.get(j, i);
			}
			mean[i]/= Response.getRowDimension();
			for( int j=0; j < Response.getRowDimension(); j++ )
			{
				Residual[j] = Response.get(j, i)-mean[i];
			}
		}
	}

	/**
	 * Get score calculated by the mean model.
	 */
	public double[] GetScore() {
		return Residual;
	}

	public void SetRidge(double r) {
	}
}