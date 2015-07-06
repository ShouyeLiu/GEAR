package gear.subcommands.eigengwas;

import gear.subcommands.CommandArguments;

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
		mPheno = i;
	}
	
	public int getMphneo()
	{
		return mPheno;
	}

	private String pheFile;

	private int mPheno;
}
