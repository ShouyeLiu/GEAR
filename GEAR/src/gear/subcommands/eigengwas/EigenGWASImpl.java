package gear.subcommands.eigengwas;

import java.io.PrintStream;
import java.util.ArrayList;

import org.apache.commons.math.MathException;
import org.apache.commons.math.distribution.NormalDistributionImpl;
import org.apache.commons.math.stat.regression.SimpleRegression;

import gear.data.InputDataSet;
import gear.family.pedigree.file.MapFile;
import gear.family.pedigree.file.SNP;
import gear.family.plink.PLINKParser;
import gear.family.popstat.GenotypeMatrix;
import gear.family.qc.rowqc.SampleFilter;
import gear.subcommands.CommandArguments;
import gear.subcommands.CommandImpl;
import gear.sumstat.qc.rowqc.SumStatQC;
import gear.util.FileUtil;
import gear.util.Logger;
import gear.util.NewIt;
import gear.util.stat.PrecisePvalue;

public class EigenGWASImpl extends CommandImpl
{

	@Override
	public void execute(CommandArguments cmdArgs)
	{
		eigenArgs = (EigenGWASArguments)cmdArgs; 

		traitIdx = eigenArgs.getMphneo();
		data.readSubjectIDFile(eigenArgs.getFam());
		data.readPhenotypeFile(eigenArgs.getPhenotypeFile());

		PLINKParser pp = PLINKParser.parse(eigenArgs);
		sf = new SampleFilter(pp.getPedigreeData(), pp.getMapData());
		ssQC = new SumStatQC(pp.getPedigreeData(), pp.getMapData(), sf);
		mapFile = ssQC.getMapFile();
		gm = new GenotypeMatrix(ssQC.getSample());
		
		eigenGWAS();

	}
	
	private void eigenGWAS()
	{

		NormalDistributionImpl unitNormal = new NormalDistributionImpl(0, 1);
		
		double[] Y = new double[data.getNumberOfSubjects()];
		ArrayList<Integer> pheIdx = NewIt.newArrayList();
		ArrayList<Double> pArray = NewIt.newArrayList();
		double threshold = 0;

		for(int subjectIdx = 0; subjectIdx < Y.length; subjectIdx++)
		{
			if(data.isPhenotypeMissing(subjectIdx, traitIdx))
			{
				continue;
			}
			pheIdx.add(subjectIdx);
			Y[subjectIdx] = data.getPhenotype(subjectIdx, traitIdx);
			threshold += Y[subjectIdx];
		}
		threshold /= pheIdx.size();

		PrintStream eGWAS = FileUtil.CreatePrintStream(eigenArgs.getOutRoot() + ".egwas");
		eGWAS.println("SNP\tCHR\tBP\tRefAllele\tAltAllele\tfreq\tBeta\tSE\tP\tChi\tn1\tfreq1\tn2\tfreq2\tFst");

		ArrayList<SNP> snpList = mapFile.getMarkerList();

		for(int i = 0; i < gm.getNumMarker(); i++)
		{
			SNP snp = snpList.get(i);

			if (eigenArgs.isChrFlagOn() && Integer.parseInt(snp.getChromosome()) != eigenArgs.getChr()) continue;

			SimpleRegression sReg = new SimpleRegression();
			double n1=0, n2=0, N=0, freq1=0, freq2=0, freq=0;
			for (int j = 0; j < pheIdx.size(); j++)
			{
				int idx = pheIdx.get(j).intValue();
				int g = gm.getAdditiveScoreOnFirstAllele(idx, i);
				if (g != 3)
				{
					sReg.addData(g, Y[idx]);
					if (Y[idx] < threshold)
					{
						n1 = n1 + 1;
						freq1 = freq1 + g;
					}
					else
					{
						n2 = n2 + 1;
						freq2 = freq2 + g;
					}
					N = N + 1;
					freq = freq + g;
				}
			}

			freq1 = freq1/(2*n1);
			freq2 = freq2/(2*n2);
			freq = freq/(2*N);

			double fst = (n1/N*(freq1-freq)*(freq1-freq) + n2/N * (freq2-freq) * (freq2-freq))/(freq * (1-freq));

			eGWAS.print(snp.getName() + "\t" + snp.getChromosome() + "\t" + snp.getPosition() + "\t" + snp.getFirstAllele() + "\t" + snp.getSecAllele()+ "\t" + freq + "\t");
			double b = sReg.getSlope();
			double b_se = sReg.getSlopeStdErr();
			double z = b/b_se;
			double p = 1;

			try
			{
				if (Math.abs(z) < 8)
				{
					p = (1-unitNormal.cumulativeProbability(Math.abs(z)))*2;
				}
				else
				{
					p = PrecisePvalue.TwoTailZcumulativeProbability(Math.abs(z));
				}
			}
			catch (MathException e)
			{
				Logger.printUserError(e.toString());
			}
			pArray.add(p);
			eGWAS.println(b + "\t" + b_se + "\t" + "\t" + p + "\t" + z * z + "\t" + n1 + "\t" + freq1 + "\t" + n2 + "\t" + freq2 + "\t" + fst);
		}
		eGWAS.close();
		
		double lambdaGC = PrecisePvalue.getRealGC(pArray);
		Logger.printUserLog("Lambda GC is : " + lambdaGC);

	}

	private EigenGWASArguments eigenArgs;
	private MapFile mapFile;
	private SampleFilter sf;
	private SumStatQC ssQC;
	private GenotypeMatrix gm;
	private int traitIdx;
	private InputDataSet data = new InputDataSet();
}
