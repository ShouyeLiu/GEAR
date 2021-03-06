package gear.subcommands.labpop;

import java.io.BufferedReader;
import java.io.BufferedWriter;
import java.io.DataOutputStream;
import java.io.FileOutputStream;
import java.io.FileWriter;
import java.io.IOException;
import java.io.PrintWriter;
import java.util.Arrays;

import org.apache.commons.math.random.RandomDataImpl;
import org.apache.commons.math.stat.StatUtils;

import gear.ConstValues;
import gear.subcommands.CommandArguments;
import gear.subcommands.CommandImpl;
import gear.util.FileUtil;
import gear.util.Logger;
import gear.util.Sample;

public class LabPopCommandImpl extends CommandImpl
{

	@Override
	public void execute(CommandArguments cmdArgs) 
	{
		lpArgs = (LabPopCommandArguments) cmdArgs;

		M = lpArgs.getNumberOfMarkers();
		nullM = lpArgs.getNullMarkerNum();
		rec = lpArgs.getRec();
		gm = new int[lpArgs.getSampleSize()][lpArgs.getNumberOfMarkers()];
		bv = new double[lpArgs.getSampleSize()];
		phe = new double[lpArgs.getSampleSize()][lpArgs.getReplication()];

		idx = Sample.SampleIndex(0, M-1, M-nullM);
		Arrays.sort(idx);

		if (lpArgs.is1234mode())
		{
			A[0] = "1";
			A[1] = "2";
		}

		if (lpArgs.isATGCmode())
		{
			A[0] = "A";
			A[1] = "C";
		}

		rnd = new RandomDataImpl();
		rnd.reSeed(lpArgs.getSeed());

		getAddEffect();
		getDomEffect();

		generateGenotype();
		generatePhenotype();
		
		if (lpArgs.isVerbose())
		{
			writeEffFile();			
		}
		writePheno();

		if (lpArgs.getMakeBed())
		{
			writeBFile();
		}
		else
		{
			writeFile();
		}
	}

	private void generatePhenotype()
	{
		if (lpArgs.getHsq() == 0)
		{
			Arrays.fill(polyEffect, 0);
		}

		for (int i = 0; i < lpArgs.getSampleSize(); i++)
		{
			for (int j = 0; j < lpArgs.getNumberOfMarkers(); j++)
			{
				bv[i] += gm[i][j] * polyEffect[j];
				if (gm[i][j] == 1)
				{
					bv[i] += deffect[j];
				}
			}
		}

		double vg = StatUtils.variance(bv);
		Logger.printUserLog("Genetic variation is "+vg);
		double ve = 0;

		double H2=lpArgs.getHsq() + lpArgs.getHsqDom();
		if (lpArgs.getHSQB() > 0)
		{
			H2 = lpArgs.getHSQB();
		}

		if (H2 == 0)
		{
			Arrays.fill(bv, 0);
			ve = 1;
		}
		else
		{
			ve = vg * (1 - H2) / H2;
		}

		for (int i = 0; i < lpArgs.getSampleSize(); i++)
		{
			for (int j = 0; j < lpArgs.getReplication(); j++)
			{
				phe[i][j] = bv[i] + rnd.nextGaussian(0, Math.sqrt(ve));
			}
		}
	}

	private void getAddEffect()
	{
		polyEffect = new double[M];

		Sample.setSeed(lpArgs.getSeed()+1);

		for (int i = 0; i < idx.length; i++)
		{
			polyEffect[idx[i]] = rnd.nextGaussian(0, 1);
		}
		
		if (lpArgs.isPolyEffectFile())
		{
			BufferedReader reader = FileUtil.FileOpen(lpArgs.getPolyEffectFile());
			int c = 0;
			String line = null;
			try
			{
				while ((line = reader.readLine()) != null)
				{
					if(c >= M)
					{
						Logger.printUserLog("Have already read " + M + " allelic effects. Ignore the rest of the content in '" + lpArgs.getPolyEffectFile() + "'.");
						break;
					}

					line.trim();
					String[] l = line.split(ConstValues.WHITESPACE_DELIMITER);
					if (l.length < 1) continue;
					polyEffect[c++] = Double.parseDouble(l[0]);
				}
				reader.close();
			}
			catch (IOException e)
			{
				Logger.handleException(e,
						"An exception occurred when reading the frequency file '"
								+ lpArgs.getPolyEffectFile() + "'.");
			}
		}
		else if (lpArgs.getHsq() == 0)
		{
			Arrays.fill(polyEffect, 0);
		}
		return;
	}
	
	private void getDomEffect()
	{
		deffect = new double[M];

		Sample.setSeed(lpArgs.getSeed()+1);

		for (int i = 0; i < idx.length; i++)
		{
			deffect[idx[i]] = rnd.nextGaussian(0, 1);
		}

		double t = 0;
		double adj=1;
		if (lpArgs.getHsq() > 0)
		{
			t = lpArgs.getHsq()/lpArgs.getHsqDom();
			adj = Math.sqrt(1/(0.5*t));
			for (int i = 0; i < idx.length; i++)
			{
				deffect[idx[i]] = rnd.nextGaussian(0, adj);
			}
		}
		else if (lpArgs.isPolyDomEffectFile())
		{
			BufferedReader reader = FileUtil.FileOpen(lpArgs.getPolyDomEffectFile());
			int c = 0;
			String line = null;
			try
			{
				while ((line = reader.readLine()) != null)
				{
					if(c >= M)
					{
						Logger.printUserLog("Have already read " + M + " allelic effects. Ignore the rest of the content in '" + lpArgs.getPolyDomEffectFile() + "'.");
						break;
					}

					line.trim();
					String[] l = line.split(ConstValues.WHITESPACE_DELIMITER);
					if (l.length < 1) continue;
					deffect[c++] = Double.parseDouble(l[0]);
				}
				reader.close();
			}
			catch (IOException e)
			{
				Logger.handleException(e,
						"An exception occurred when reading the frequency file '"
								+ lpArgs.getPolyDomEffectFile() + "'.");
			}
		}
		else if (lpArgs.getHsqDom() == 0)
		{
			Arrays.fill(deffect, 0);
		}
		return;
	}

	private void generateGenotype()
	{
		if (lpArgs.isBC())
		{
			for (int i = 0; i < gm.length; i++)
			{
				gm[i] = sampleHaploid();
			}
		}
		else if (lpArgs.isDH())
		{
			for (int i = 0; i < gm.length; i++)
			{
				int[] hp = sampleHaploid();
				System.arraycopy(hp, 0, gm[i], 0, hp.length);

				for (int j = 0; j < lpArgs.getNumberOfMarkers(); j++)
				{
					gm[i][j] *= 2; 
				}
			}
		}
		else if (lpArgs.isF2())
		{
			for (int i = 0; i < gm.length; i++)
			{
				int[] hp1 = sampleHaploid();
				int[] hp2 = sampleHaploid();

				for (int j = 0; j < lpArgs.getNumberOfMarkers(); j++)
				{
					gm[i][j] = hp1[j] + hp2[j];
				}
			}
		}
		else if (lpArgs.isRIL()) 
		{
			for (int i = 0; i < gm.length; i++)
			{
				int[] hp = new int[lpArgs.getNumberOfMarkers()];

				for (int j = 0; j < lpArgs.getNumberOfMarkers(); j++)
				{
					double r = rnd.nextUniform(0, 1);

					if (j == 0)
					{
						hp[j] = r < rec[j] ? 0 : 1;
					}
					else
					{
						hp[j] = r < ( 1 / (1 + 2*rec[j]) ) ? hp[j-1] : (1 - hp[j-1]);
					}
				}

				for (int j = 0; j < lpArgs.getNumberOfMarkers(); j++)
				{
					gm[i][j] = hp[j] * 2;
				}
			}
		}
		else if (lpArgs.isIF2())
		{
			for (int i = 0; i < gm.length; i++)
			{
				int[][] hp2 = new int[2][lpArgs.getNumberOfMarkers()];

				for (int h = 0; h < 2; h++)
				{
					for (int j = 0; j < lpArgs.getNumberOfMarkers(); j++)
					{
						double r = rnd.nextUniform(0, 1);

						if (j == 0)
						{
							hp2[h][j] = r < rec[j] ? 0 : 1;
						}
						else
						{
							hp2[h][j] = r < ( 1 / (1 + 2*rec[j]) ) ? hp2[h][j-1] : (1 - hp2[h][j-1]);
						}
					}
				}

				for (int j = 0; j < lpArgs.getNumberOfMarkers(); j++)
				{
					gm[i][j] = hp2[0][j] + hp2[1][j];
				}
			}
		}
	}

	private int[] sampleHaploid()
	{
		int[] hp = new int[lpArgs.getNumberOfMarkers()];
		for (int j = 0; j < lpArgs.getNumberOfMarkers(); j++)
		{
			double r = rnd.nextUniform(0, 1);
			if (j == 0)
			{
				hp[j] = r < rec[j] ? 0 : 1;
			}
			else
			{
				hp[j] = r < rec[j] ? (1-hp[j-1]) : hp[j-1];
			}
		}
		return hp;
	}

	private void writeEffFile()
	{
		PrintWriter eff = null;
		try
		{
			eff = new PrintWriter(new BufferedWriter(new FileWriter(lpArgs.getOutRoot() + ".rnd")));
		}
		catch (IOException e)
		{
			Logger.handleException(e,
					"An exception occurred when writing files.");
		}

		for (int i = 0; i < M; i++)
		{
			eff.println("rs" + i + " " + A[0] + " " + polyEffect[i] + " " + deffect[i]);
		}
		eff.close();
		
		PrintWriter nullMFile = null;
		PrintWriter QTLFile = null;

		try
		{
			nullMFile = new PrintWriter(new BufferedWriter(new FileWriter(lpArgs.getOutRoot() + ".nqtl")));
			QTLFile = new PrintWriter(new BufferedWriter(new FileWriter(lpArgs.getOutRoot() + ".qtl")));
		}
		catch (IOException e)
		{
			Logger.handleException(e,
					"An exception occurred when writing files.");
		}

		int cnt=0;
		for (int i = 0; i < M; i++)
		{
			if (cnt < idx.length && i == idx[cnt])
			{
				QTLFile.println("rs" + i);
				cnt++;
				continue;
			}
			nullMFile.println("rs" + i);
		}
		nullMFile.close();
		QTLFile.close();
	}

	public void writePheno()
	{
		PrintWriter pheno = null;
		PrintWriter breed = null;
		PrintWriter add = null;
		try
		{
			pheno = new PrintWriter(new BufferedWriter(new FileWriter(lpArgs.getOutRoot() + ".phe")));
			breed = new PrintWriter(new BufferedWriter(new FileWriter(lpArgs.getOutRoot() + ".breed")));
			add = new PrintWriter(new BufferedWriter(new FileWriter(lpArgs.getOutRoot() + ".add")));
		}
		catch (IOException e)
		{
			Logger.handleException(e, "An I/O exception occurred when creating the .phe file.");
		}

		for (int i = 0; i < lpArgs.getSampleSize(); i++)
		{
			int fid = (i + 1) * 10000;

			pheno.print(fid + " ");
			pheno.print(fid + " ");
			for (int j = 0; j < lpArgs.getReplication(); j++)
			{
				if (j != (lpArgs.getReplication() - 1))
				{
					pheno.print(phe[i][j] + " ");			
				}
				else
				{
					pheno.print(phe[i][j]);					
				}
			}
			pheno.println();
		}
		pheno.close();
		
		for (int i = 0; i < lpArgs.getSampleSize(); i++)
		{
			int fid = (i + 1) * 10000;
			breed.println(fid + " " + fid + " " + bv[i]);
		}
		breed.close();
		
		for (int i = 0; i < lpArgs.getNumberOfMarkers(); i++)
		{
			add.println("rs" + i + " " + polyEffect[i]);
		}
		add.close();
	}

	public void writeFile()
	{
		PrintWriter ped = null;
		PrintWriter map = null;

		try
		{
			ped = new PrintWriter(new BufferedWriter(new FileWriter(lpArgs.getOutRoot() + ".ped")));
			map = new PrintWriter(new BufferedWriter(new FileWriter(lpArgs.getOutRoot() + ".map")));
		}
		catch (IOException e)
		{
			Logger.handleException(e, "An I/O exception occurred when creating the .ped and .map files.");
		}

		for (int h = 0; h < lpArgs.getSampleSize(); h++)
		{
			int fid = (h + 1) * 10000;

			ped.print(fid + " ");
			ped.print(fid + " ");
			ped.print(0 + " ");
			ped.print(0 + " ");
			ped.print(1 + " ");
			ped.print(1 + " ");

			for (int i = 0; i < lpArgs.getNumberOfMarkers(); i++)
			{
				StringBuilder sb = new StringBuilder();
				switch (gm[h][i])
				{
				case 0:
					sb.append(A[0] + " " + A[0]);
					break;
				case 1:
					sb.append(A[0] + " " + A[1]);
					break;
				case 2:
					sb.append(A[1] + " " + A[1]);
					break;
				default:
					break;
				}
				if (i == (lpArgs.getNumberOfMarkers() - 1))
				{
					ped.print(sb.toString());			
				}
				else
				{
					ped.print(sb.toString() + " ");
				}
			}
			ped.print("\n");
		}

		for (int i = 0; i < lpArgs.getNumberOfMarkers(); i++)
		{
			map.print(1 + " ");
			map.print("rs" + i + " ");
			map.print(i / (lpArgs.getNumberOfMarkers() * 1.0) + " ");
			map.println(i * 100);
		}

		ped.close();
		map.close();
	}

	public void writeBFile()
	{
		DataOutputStream bedout = null;
		PrintWriter fam = null;
		PrintWriter bim = null;

		try
		{
			bedout = new DataOutputStream(new FileOutputStream(lpArgs.getOutRoot() + ".bed"));
			fam = new PrintWriter(new BufferedWriter(new FileWriter(lpArgs.getOutRoot() + ".fam")));
			bim = new PrintWriter(new BufferedWriter(new FileWriter(lpArgs.getOutRoot() + ".bim")));

		} 
		catch (IOException e)
		{
			Logger.handleException(e, "An I/O exception occurred when creating the .bed, .fam and .bim files.");
		}

		int fid = 0;
		for (int i = 0; i < gm.length; i++)
		{
			fid = (i+1) * 10000;
			fam.print(fid + " ");
			fam.print(fid + " ");
			fam.print(0 + " ");
			fam.print(0 + " ");
			fam.print(0 + " ");
			fam.print(-9 + "\n");
		}
		fam.close();

		try
		{
			bedout.writeByte(ConstValues.PLINK_BED_BYTE1);
			bedout.writeByte(ConstValues.PLINK_BED_BYTE2);
			bedout.writeByte(ConstValues.PLINK_BED_BYTE3);
			for (int i = 0; i < lpArgs.getNumberOfMarkers(); i++)
			{
				byte gbyte = 0;
				int idx = 0;
				for (int j = 0; j < gm.length; j++)
				{
					int g = (int) gm[j][i];
					switch (g)
					{
					case 0:
						g = 0;
						break;
					case 1:
						g = 2;
						break;
					case 2:
						g = 3;
						break;
					default:
						g = 1;
						break; // missing
					}

					g <<= 2 * idx;
					gbyte |= g;
					idx++;

					if (j != (gm.length - 1))
					{
						if (idx == 4)
						{
							bedout.writeByte(gbyte);
							gbyte = 0;
							idx = 0;
						}
					}
					else
					{
						bedout.writeByte(gbyte);
					}
				}
			}
			bedout.close();
		}
		catch (IOException e)
		{
			Logger.handleException(e, "An I/O exception occurred when writing the .bed file.");
		}

		for (int i = 0; i < lpArgs.getNumberOfMarkers(); i++)
		{
			bim.print(1 + " ");
			bim.print("rs" + i + " ");
			bim.print(i / (lpArgs.getNumberOfMarkers() * 1.0) + " ");
			bim.print(i * 100 + " ");
			bim.println(A[0] + " " + A[1]);
		}
		bim.close();

	}

	private int M;
	private int nullM;
	private int[] idx;
	private LabPopCommandArguments lpArgs;
	private RandomDataImpl rnd = new RandomDataImpl();

	private String[] A = { "A", "C" };

	private int[][] gm = null;
	private double[] bv = null;
	private double[][] phe = null;

	private double[] polyEffect;
	private double[] deffect;
	private double[] rec;

	PrintWriter ibdF = null;
	PrintWriter ibdSibF = null;

}
