package gear.subcommands.eigengwas;

import org.apache.commons.cli.CommandLine;
import org.apache.commons.cli.OptionBuilder;
import org.apache.commons.cli.Options;

import gear.subcommands.Command;
import gear.subcommands.CommandArgumentException;
import gear.subcommands.CommandArguments;
import gear.subcommands.CommandImpl;

public class EigenGWASCommand extends Command
{

	public EigenGWASCommand()
	{
		addAlias("egwas");
	}

	@Override
	public String getName()
	{
		return "eigengwas";
	}

	@Override
	public String getDescription()
	{
		return "Eigen GWAS";
	}

	@SuppressWarnings("static-access")
	@Override
	public void prepareOptions(Options options)
	{
		options.addOption(OptionBuilder.withDescription(OPT_BFILE_DESC).withLongOpt(OPT_BFILE_LONG).hasArg().isRequired().create());
		options.addOption(OptionBuilder.withDescription(OPT_PHE_DESC).withLongOpt(OPT_PHE_LONG).hasArg().isRequired().create());
		options.addOption(OptionBuilder.withDescription(OPT_MPHE_DESC).withLongOpt(OPT_MPHE_LONG).hasArg().create());
		options.addOption(OptionBuilder.withDescription(OPT_CHR_DESC).hasArg().create(OPT_CHR));
	}

	@Override
	public CommandArguments parse(CommandLine cmdLine) throws CommandArgumentException
	{
		EigenGWASArguments eigenArgs = new EigenGWASArguments();
		parseFileArguments(eigenArgs, cmdLine);
		eigenArgs.setPhentypeIndex(parseIntOptionValue(cmdLine, OPT_MPHE_LONG, "1"));
		eigenArgs.setPhenotypeFile(cmdLine.getOptionValue(OPT_PHE_LONG));
		return eigenArgs;
	}

	private void parseFileArguments(EigenGWASArguments eigenArgs, CommandLine cmdLine) throws CommandArgumentException
	{
		String bfile = cmdLine.getOptionValue(OPT_BFILE_LONG);
		String file = cmdLine.getOptionValue(OPT_FILE_LONG);

		if (bfile == null && file == null)
		{
			throw new CommandArgumentException("No genotypes are provided. Either --" + OPT_BFILE_LONG + " or --" + OPT_FILE_LONG + " must be set.");
		}

		if (bfile != null && file != null)
		{
			throw new CommandArgumentException("--" + OPT_BFILE_LONG + " and --" + OPT_FILE_LONG + " cannot be set together.");
		}

		eigenArgs.setBFile(bfile);
		eigenArgs.setFile(file);
		
		if (cmdLine.hasOption(OPT_CHR))
		{
			eigenArgs.setChr(cmdLine.getOptionValue(OPT_CHR));
		}
	}

	@Override
	protected CommandImpl createCommandImpl()
	{
		return new EigenGWASImpl();
	}

	private final static String OPT_PHE_LONG = "pheno";
	private final static String OPT_PHE_DESC = "Specify the phenotype file (individual eigenvector)";

	private final static String OPT_MPHE_LONG = "mpheno";
	private final static String OPT_MPHE_DESC = "Specify the phenotype index";
	
	private final static String OPT_CHR = "chr";
	private final static String OPT_CHR_DESC = "Specify the chromosomes for analysis";

}
