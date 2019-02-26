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

public class nblast_clustering implements PlugIn {

	static final String RSCRIPT = "NBC.RScript.string";
	static final String METHOD = "NBC.method.string";
	static final String KVAL = "NBC.kval.int";
	static final String HVAL = "NBC.hval.double";
	static final String USEH = "NBC.korh.boolean";

	static final String CLUSTERING_SCRIPT = IJ.getDirectory("plugins") + File.separator + "cluster.R";

	static final String[] Methods = {"ward.D", "ward.D2", "single", "complete", "average", "mcquitty", "median", "centroid"};

	static final Pattern NUMBERS = Pattern.compile("(?<=\\D)(?=\\d)|(?<=\\d)(?=\\D)");

	String m_defaultRpath = "RScript";
	String m_rscript = Prefs.get(RSCRIPT, m_defaultRpath);
	String m_mtd = Prefs.get(METHOD, "ward.D2");
	int m_k = (int)Prefs.get(KVAL, 50);
	double m_h = (double)Prefs.get(HVAL, 0.5);
	boolean m_useh = Prefs.get(USEH, false);

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
        GenericDialog gd = new GenericDialog("NBLAST Clustering");

        if ( !new File(m_rscript).isFile() )
        	m_rscript = getDefaultRPath();
		
		gd.addStringField("RScript", m_rscript);
		gd.addChoice("method", Methods, m_mtd);
		gd.addNumericField("N_clusters", m_k, 0);
		gd.addNumericField("height cutoff value", m_h, 2);
		String[] items = {"Set number of clusters", "Use height cutoff value"};
		int id = m_useh ? 1 : 0;
		gd.addRadioButtonGroup(null, items, 2, 1, items[id]);
		
		gd.showDialog();
		if(gd.wasCanceled()){
			return false;
		}
		
		m_rscript = gd.getNextString();
		m_mtd = gd.getNextChoice();
		m_k = (int)gd.getNextNumber();
		m_h = (double)gd.getNextNumber();
		m_useh = gd.getNextRadioButton().equals(items[1]);

		if (m_useh) m_k = -1;
		else m_h = -1.0;
		
		Prefs.set(RSCRIPT, m_rscript);
		Prefs.set(METHOD, m_mtd);
		Prefs.set(KVAL, m_k);
		Prefs.set(HVAL, m_h);
		Prefs.set(USEH, m_useh);
		
        return true;
    }

	public void run(String arg) {

		IJ.log("input score matrix");
		OpenDialog ind = new OpenDialog("input score matrix");
		if (ind.getFileName() == null) return;
		String indir = ind.getDirectory();
		String inpath = indir + ind.getFileName();
		String label = FilenameUtils.getBaseName(ind.getFileName()) + "_";
		
		IJ.log("Choose an output directory");
		DirectoryChooser dcout = new DirectoryChooser("output directory");
		String outdir = dcout.getDirectory();
		if (outdir == null) return;

		indir.replace(File.separator, "/");
		inpath.replace(File.separator, "/");
		outdir.replace(File.separator, "/");

		String swcdir = indir + label + "swc" + File.separator;
		String mipdir = indir + label + "mip" + File.separator;
		
		IJ.log("inpput dir: " + indir);
		IJ.log("score matrix: " + inpath);
		IJ.log("swc dir: " + swcdir);
		IJ.log("mip dir: " + mipdir);
		IJ.log("output dir: " + outdir);
		
		if (!showDialog())
			return;

		int count = 0;

		try {
			IJ.log("-------- NBLAST Clustering --------");
			String[] listCommands = {
						m_rscript,
						CLUSTERING_SCRIPT,
						inpath,
						swcdir,
						mipdir,
						outdir,
						m_mtd,
						String.valueOf(m_k),
						String.valueOf(m_h)
					};
			RThread rth = new RThread(listCommands);
			rth.start();
			rth.join();
			IJ.log(rth.getStdOut());
			IJ.log(rth.getStdErr());
			IJ.log("--------         (DONE)          --------\n");

			
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
