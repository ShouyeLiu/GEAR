package gear.subcommands.eigengwas;

import gear.util.Logger;
import gear.util.stat.PrecisePvalue;

import org.apache.commons.math.MathException;
import org.apache.commons.math.distribution.NormalDistributionImpl;

public class EigenGWASResult
{
	public EigenGWASResult(String snp, String chr, long BP, char A1, char A2, double n1, double n2, double N, double Freq, double freq1, double freq2, double b, double se, double p)
	{
		this.SNP = snp;
		this.chr = chr;
		this.BP = BP;
		this.A1 = A1;
		this.A2 = A2;

		this.n1 = n1;
		this.n2 = n2;
		this.N = N;
		
		this.freq = Freq;
		this.freq1 = freq1;
		this.freq2 = freq2;
		
		this.b = b;
		this.se = se;
		this.z = b/se;
		this.p = p;
		this.chi = this.z * this.z;
		CalFst();
	}

	public void SetGC(double gc)
	{
		if(gc < 1) 
		{
			pgc = p;
		}
		else
		{
			NormalDistributionImpl unitNormal = new NormalDistributionImpl(0, 1);
			double z1 = b / (Math.sqrt(gc) * se);
			try
			{
				if (Math.abs(z1) < 8)
				{
					pgc = (1-unitNormal.cumulativeProbability(Math.abs(z1)))*2;
				}
				else
				{
					pgc = PrecisePvalue.TwoTailZcumulativeProbability(Math.abs(z1));
				}
			}
			catch (MathException e)
			{
				Logger.printUserError(e.toString());
			}
			
		}

	}

	private void CalFst()
	{
		Fst = 2*(n1/N*(freq1-freq)*(freq1-freq) + n2/N * (freq2-freq) * (freq2-freq))/(freq * (1-freq));
	}

	public String toString()
	{
		StringBuffer sb = new StringBuffer();
		sb.append(SNP + "\t");
		sb.append(chr + "\t");
		sb.append(BP + "\t");
		sb.append(A1 + "\t");
		sb.append(A2 + "\t");
		sb.append(freq + "\t");
		sb.append(b + "\t");
		sb.append(se + "\t");
		sb.append(p + "\t");
		sb.append(pgc + "\t");
		sb.append(chi + "\t");
		sb.append(n1 + "\t");
		sb.append(freq1 + "\t");
		sb.append(n2 + "\t");
		sb.append(freq2 + "\t");
		sb.append(Fst);
		return sb.toString();
	}
	private String SNP;
	private String chr;
	private long BP;
	private char A1;
	private char A2;
	private double n1;
	private double n2;
	private double N;
	
	private double freq;
	private double freq1;
	private double freq2;
	
	private double b;
	private double se;
	private double z;
	private double p;
	private double pgc;
	private double chi;
	private double Fst;

}
