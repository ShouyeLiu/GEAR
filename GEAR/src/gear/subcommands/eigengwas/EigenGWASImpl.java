package gear.subcommands.eigengwas;

import java.io.PrintStream;
import java.util.ArrayList;

import org.apache.commons.math.MathException;
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
import gear.util.NewIt;

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
		double[] Y = new double[data.getNumberOfSubjects()];
		ArrayList<Integer> pheIdx = NewIt.newArrayList();
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
		eGWAS.println("SNP\tCHR\tBP\tRefAllele\tAltAllele\tBeta\tSE\tP\tChi\tFst");

		ArrayList<SNP> snpList = mapFile.getMarkerList();

		for(int i = 0; i < gm.getNumMarker(); i++)
		{
			SimpleRegression sReg = new SimpleRegression();
			SNP snp = snpList.get(i);
			double w1=0, w2=0, w=0, p1=0, p2=0, p=0;
			for(int j = 0; j < pheIdx.size(); j++)
			{
				int idx = pheIdx.get(j).intValue();
				int g = gm.getAdditiveScoreOnFirstAllele(idx, i);
				if (g != 3)
				{
					sReg.addData(g, Y[idx]);
					if (Y[idx] < threshold)
					{
						w1 = w1 + 1;
						p1 = p1 + g;
					}
					else
					{
						w2 = w2 + 1;
						p2 = p2 + g;
					}
					w = w + 1;
					p = p + g;
				}
			}
			
			p1 = p1/(2*w1);
			p2 = p2/(2*w2);
			p = p/(2*w);
			w1 = w1/w;
			w2 = w2/w;

			double fst = (w1*(p1-p)*(p1-p) + w2 * (p2-p) * (p2-p))/(p * (1-p));
			
			eGWAS.print(snp.getName() + "\t" + snp.getChromosome() + "\t" + snp.getPosition() + "\t" + snp.getFirstAllele() + "\t" + snp.getSecAllele()+"\t");
			
			try
			{
				eGWAS.println(sReg.getSlope() + " " + sReg.getSlopeStdErr() + " " + sReg.getSignificance() + " " + p + " " + (sReg.getSlope() / sReg.getSlopeStdErr()) * (sReg.getSlope() / sReg.getSlopeStdErr()) + " " + fst);
			}
			catch (MathException e)
			{
				// TODO Auto-generated catch block
				e.printStackTrace();
			}
		}
		
		eGWAS.close();
	}

	private EigenGWASArguments eigenArgs;
	private MapFile mapFile;
	private SampleFilter sf;
	private SumStatQC ssQC;
	private GenotypeMatrix gm;
	private int traitIdx;
	private InputDataSet data = new InputDataSet();
}
