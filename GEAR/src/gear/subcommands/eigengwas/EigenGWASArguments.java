package gear.subcommands.eigengwas;

import gear.subcommands.CommandArguments;
import gear.util.Logger;

public class EigenGWASArguments extends CommandArguments
{
	public String getPhenotypeFile()
	{
		return pheFile;
	}

	public void setPhenotypeFile(String pheFile)
	{
		this.pheFile = pheFile;
	}

	public void setPhentypeIndex(int i)
	{
		mPheno = i - 1;
	}
	
	public int getMphneo()
	{
		return mPheno;
	}

	public void setChr(String c)
	{
		chr = Integer.parseInt(c);
		if (chr < 1)
		{
			Logger.printUserLog("Chromosome should be greater than 0.\n GEAR quitted");
			System.exit(1);
		}
		chrFlag = true;
	}
	
	public int getChr()
	{
		return chr;
	}

	public boolean isChrFlagOn()
	{
		return chrFlag;
	}

	private String pheFile;
	private int chr;
	private boolean chrFlag = false;
	private int mPheno;
}
