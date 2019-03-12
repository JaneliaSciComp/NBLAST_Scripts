import ij.*;
import ij.io.*;
import ij.process.*;
import ij.gui.*;
import ij.plugin.PlugIn;

import java.awt.*;
import java.io.*;
import java.util.*;
import java.util.regex.Pattern;
import java.math.BigInteger;

import org.apache.commons.io.*;
import org.apache.commons.lang3.SystemUtils;

public class nblast_search_batch implements PlugIn {

	static final String RSCRIPT = "NBSH.RScript.string";
	static final String RESAMPLE = "NBSH.resample.double";
	static final String KVAL = "NBSH.kval.int";
	static final String NORMALIZATION = "NBSH.normalization.string";
	static final String RESULTNUM = "NBSH.resultnum.int";

	static final String SEARCH_SCRIPT = IJ.getDirectory("plugins") + File.separator + "nblast_search.R";

	static final Pattern NUMBERS = Pattern.compile("(?<=\\D)(?=\\d)|(?<=\\d)(?=\\D)");

	static final String[] Methods = {"forward", "mean"};

	String m_defaultRpath = "RScript";
	String m_rscript = Prefs.get(RSCRIPT, m_defaultRpath);
	double m_rsmp = (double)Prefs.get(RESAMPLE, 4.0);
	int m_k = (int)Prefs.get(KVAL, 3);
	String m_mtd = Prefs.get(NORMALIZATION, "mean");
	int m_rnum = (int)Prefs.get(RESULTNUM, 300);
	

	public String getDefaultRPath() {
		String defaultRpath = "RScript";
		IJ.log("Searching Rscript...");
		try {
			if (SystemUtils.IS_OS_WINDOWS) {
				String r64dir = System.getenv("ProgramFiles") + File.separator + "R";
				String r86dir = System.getenv("ProgramFiles(X86)") + File.separator + "R";
				String rdir = "";
				if (new File(r64dir).isDirectory())
					rdir = r64dir;
				else if (new File(r86dir).isDirectory())
					rdir = r86dir;
				IJ.log("r64dir: " + r64dir);
				IJ.log("r86dir: " + r86dir);
				IJ.log("R Root: " + rdir);
				if (!rdir.isEmpty()) {
					ArrayList<String> matches = new ArrayList<String>();
					Pattern p = Pattern.compile("R-*");
					final File folder = new File(rdir);
					final File[] filelist = folder.listFiles();
					if (filelist.length > 0) {
						File latest = filelist[0];
						for (int i = 1; i < filelist.length; i++) {
							String Rpathtest1 = latest.getAbsolutePath() + File.separator + "bin" + File.separator + "RScript.exe";
							String Rpathtest2 = filelist[i].getAbsolutePath() + File.separator + "bin" + File.separator + "RScript.exe";
							if (!new File(Rpathtest1).isFile()) {
								latest = filelist[i];
								continue;
							}
							if (!new File(Rpathtest2).isFile())
								continue;
							
							String[] split1 = NUMBERS.split(latest.getAbsolutePath());
							String[] split2 = NUMBERS.split(filelist[i].getAbsolutePath());

							int length = Math.min(split1.length, split2.length);
							int cmp = 0;
							for (int j = 0; j < length; j++) {
 								char c1 = split1[j].charAt(0);
 								char c2 = split2[j].charAt(0);
 								if (c1 >= '0' && c1 <= '9' && c2 >= 0 && c2 <= '9')
 									cmp = new BigInteger(split1[j]).compareTo(new BigInteger(split2[j]));
 								if (cmp == 0)
 									cmp = split1[j].compareTo(split2[j]);
 								if (cmp != 0)
 									break;
 							}
 							if (cmp == 0)
 								cmp = split1.length - split2.length;

							if (cmp < 0)
								latest = filelist[i];
						}
						String Rpathtest = latest.getAbsolutePath() + File.separator + "bin" + File.separator + "RScript.exe";
						if (new File(Rpathtest).isFile())
							defaultRpath = Rpathtest;
					}
				}
			} else if (SystemUtils.IS_OS_MAC_OSX) {
				defaultRpath = "usr/local/bin/Rscript";
			} else if (SystemUtils.IS_OS_LINUX) {
				defaultRpath = "usr/local/bin/Rscript";
			}
			IJ.log("RScript: " + defaultRpath);
		} catch (Exception e) {
			e.printStackTrace();
		}

		return defaultRpath;
	}

	public boolean showDialog() {
        GenericDialog gd = new GenericDialog("Build NBLAST Database");

        if ( !new File(m_rscript).isFile() )
        	m_rscript = getDefaultRPath();
		
		gd.addStringField("RScript",  m_rscript);
		gd.addNumericField("Resample",  m_rsmp, 2);
		gd.addNumericField("K",  m_k, 0);
		gd.addChoice("Scoring method", Methods, m_mtd);
		gd.addNumericField("Number of Results",  m_rnum, 0);
		
		gd.showDialog();
		if(gd.wasCanceled()){
			return false;
		}
		
		m_rscript = gd.getNextString();
		m_rsmp = (double)gd.getNextNumber();
		m_k = (int)gd.getNextNumber();
		m_mtd = gd.getNextChoice();
		m_rnum = (int)gd.getNextNumber();
		
		Prefs.set(RSCRIPT, m_rscript);
		Prefs.set(RESAMPLE, m_rsmp);
		Prefs.set(KVAL, m_k);
		Prefs.set(NORMALIZATION, m_mtd);
		Prefs.set(RESULTNUM, m_rnum);
		
        return true;
    }

	public void run(String arg) {

		IJ.log("Choose an input directory");
		DirectoryChooser dcin = new DirectoryChooser("input directory");
		String indir = dcin.getDirectory();
		if (indir == null) return;

		IJ.log("Choose a NBLAST database");
		OpenDialog indbdlg = new OpenDialog("NBLAST database");
		if (indbdlg.getFileName() == null) return;
		String indbdir = indbdlg.getDirectory();
		String indbpath = indbdir + indbdlg.getFileName();
		
		IJ.log("Choose an output directory");
		DirectoryChooser dcout = new DirectoryChooser("output directory");
		String outdir = dcout.getDirectory();
		if (outdir == null) return;

		String dbname = FilenameUtils.getBaseName(indbdlg.getFileName());
		
		IJ.log("INPUT: " + indir);
		IJ.log("DBPATH: " + indbpath);
		IJ.log("OUTDIR: " + outdir);
		
		if (!showDialog())
			return;

		int count = 0;

		try {
			IJ.log("Running NBLAST Search...");
			String[] listCommands = {
						m_rscript,
						SEARCH_SCRIPT,
						indir,
						indbpath,
						"none",
						outdir,
						String.valueOf(m_rnum),
						dbname,
						m_mtd,
						String.valueOf(m_rsmp),
						String.valueOf(m_k)
					};
			RThread rth = new RThread(listCommands);
			rth.start();
			rth.join();
			IJ.log(rth.getStdOut());
			IJ.log(rth.getStdErr());
			IJ.log("DONE");
			
		} catch (Exception e) {
			e.printStackTrace();
		}
	}

	class RThread extends Thread{
        String[] command;
        StringBuffer stdout_sb = new StringBuffer();
        StringBuffer stderr_sb = new StringBuffer();
        RThread(String[] command){        
            this.command=command;            
        }
        public void run(){      
			try{          
				String s = null;
                Process process = new ProcessBuilder(command).start();                
                BufferedReader stdInput = new BufferedReader(new InputStreamReader(process.getInputStream()));
                BufferedReader stdError = new BufferedReader(new InputStreamReader(process.getErrorStream()));
                // read the output from the command                    
                while ((s = stdInput.readLine()) != null)
                    stdout_sb.append(s+"\n");
                stdInput.close();
                // read any errors from the attempted command                
                while ((s = stdError.readLine()) != null)
                    stderr_sb.append(s+"\n");    
            }catch(Exception ex){
                System.out.println(ex.toString());
            }
        }
        public String getStdOut() { return stdout_sb.toString(); }
        public String getStdErr() { return stderr_sb.toString(); }
    }
}
