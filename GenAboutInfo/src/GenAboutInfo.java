import java.io.*;
import java.util.Calendar;

public class GenAboutInfo {
	private static String genWelcomeMessage(String svnRev) 
	{
		int width = 66; // No. of characters of each line of the welcome message
		String line;
		String welcomeMessage = "\\n" + "******************************************************************\\n"
				+ "| GEAR [GEnetic Analysis Repository] version 0.7.7               |\\n"
				+ "| (C) 2013 Guo-Bo Chen, Zhi-Xiang Zhu                            |\\n";
		line = "| Git r" + svnRev + ", built on ";
		line += Calendar.getInstance().getTime();
		for (int i = width - line.length() - 1; i > 0; --i) {
			line += " ";
		}
		line += "|\\n";
		welcomeMessage += line;
		welcomeMessage += "| GNU General Public License, v2                                 |\\n"
				+ "******************************************************************\\n";
		return welcomeMessage;
	}

	/**
	 * @param args
	 *            [package-name]
	 */
	public static void main(String[] args) throws Exception 
	{
		String packageName = null;
		if (args.length >= 1) {
			packageName = args[0];
		}
		String packageDeclaration = "", folder = "src/";
		if (packageName != null) {
			packageDeclaration = "package" + ";\n";
			folder += packageName.replace('.', '/') + "/";
		}
		// Read the SVN revision number
		// Process proc = Runtime.getRuntime().exec("svnversion -version");
		Process proc = Runtime.getRuntime().exec("git rev-list head");
		// Process proc = Runtime.getRuntime().exec("which svnversion");

		// Process proc = Runtime.getRuntime().exec("date");

		BufferedReader reader = new BufferedReader(new InputStreamReader(proc.getInputStream()));
		int cnt = 0;
		while (reader.readLine() != null) 
		{
			cnt++;
		}
		// Welcome message
		String svnRev = (new Integer(cnt)).toString();
		String welcomeMessage = genWelcomeMessage(svnRev);
		// Write the revision number to a java source file
		String writeContent = packageDeclaration + "\n" + "public class AboutInfo {\n"
				+ "    public static final String WELCOME_MESSAGE = \"" + welcomeMessage + "\";\n"
				+ "    public static final String SVN_REVISION = \"" + svnRev + "\";\n" + "}";
		FileWriter writer = new FileWriter(folder + "AboutInfo.java");
		writer.write(writeContent);
		writer.close();
		System.out.println("Git Revision: " + svnRev);
	}
}