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

public class build_nblastdb implements PlugIn {

	static final String RSCRIPT = "NBDB.RScript.string";
	static final String RESAMPLE = "NBDB.resample.double";
	static final String KVAL = "NBDB.kval.int";

	static final String INSTALL_SCRIPT = IJ.getDirectory("plugins") + File.separator + "install.R";
	static final String NEW_SCRIPT = IJ.getDirectory("plugins") + File.separator + "new.R";
	static final String APPEND_SCRIPT = IJ.getDirectory("plugins") + File.separator + "append.R";

	static final Pattern NUMBERS = Pattern.compile("(?<=\\D)(?=\\d)|(?<=\\d)(?=\\D)");

	String m_defaultRpath = "RScript";
	String m_rscript = Prefs.get(RSCRIPT, m_defaultRpath);
	double m_rsmp = (double)Prefs.get(RESAMPLE, 3.0);
	int m_k = (int)Prefs.get(KVAL, 3);

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
				defaultRpath = "/usr/local/bin/RScript";
			} else if (SystemUtils.IS_OS_LINUX) {
				defaultRpath = "/usr/local/bin/RScript";
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
		
		gd.showDialog();
		if(gd.wasCanceled()){
			return false;
		}
		
		m_rscript = gd.getNextString();
		m_rsmp = (double)gd.getNextNumber();
		m_k = (int)gd.getNextNumber();
		
		Prefs.set(RSCRIPT, m_rscript);
		Prefs.set(RESAMPLE, m_rsmp);
		Prefs.set(KVAL, m_k);
		
        return true;
    }

	public void run(String arg) {

		IJ.log("Choose an input directory");
		DirectoryChooser dcin = new DirectoryChooser("input directory");
		String indir = dcin.getDirectory();
		if (indir == null) return;
		IJ.log("Save NBLAST database");
		SaveDialog outdb = new SaveDialog("save nblast database", "database", ".rds");
		if (outdb.getFileName() == null) return;
		String outdir = outdb.getDirectory();
		String outpath = outdir + outdb.getFileName();

		boolean append = false;
		if ( new File(outpath).isFile() )
			append = true;

		indir.replace(File.separator, "/");
		outdir.replace(File.separator, "/");
		outpath.replace(File.separator, "/");
		
		IJ.log("INPUT: " + indir);
		IJ.log("OUTDIR: " + outdir);
		IJ.log("OUTPUT: " + outpath);
		
		if (!showDialog())
			return;

		int count = 0;

		try {
			IJ.log("-------- Checking R Libraries --------");
			String[] install_listCommands = {
						m_rscript,
						INSTALL_SCRIPT
					};
			RThread rth = new RThread(install_listCommands);
			rth.start();
			rth.join();
			IJ.log(rth.getStdOut());
			IJ.log(rth.getStdErr());
			IJ.log("--------        (DONE)        --------\n");


			if (append) {
				IJ.log("-------- Adding Neurons into Database --------");
				String[] append_listCommands = {
							m_rscript,
							APPEND_SCRIPT,
							indir,
							outpath,
							String.valueOf(m_rsmp),
							String.valueOf(m_k)
						};
				rth = new RThread(append_listCommands);
				
			} else {
				IJ.log("-------- Creating New Database --------");
				String[] new_listCommands = {
							m_rscript,
							NEW_SCRIPT,
							indir,
							outpath,
							String.valueOf(m_rsmp),
							String.valueOf(m_k)
						};
				rth = new RThread(new_listCommands);
			}
			rth.start();
			rth.join();
			IJ.log(rth.getStdOut());
			IJ.log(rth.getStdErr());
			IJ.log("--------           (DONE)            --------\n");

			IJ.log("-------- Generating MIPs --------");
			count = 0;
			String swcdirpath = outdb.getDirectory() + "swc" + File.separator;
			String mipdirpath = outdb.getDirectory() + "swc_prev" + File.separator;
			File mipdir = new File(mipdirpath);
			if (!mipdir.exists())
				mipdir.mkdirs();
			final File folder = new File(dcin.getDirectory());
			for (final File fileEntry : folder.listFiles()) {
				String ext = FilenameUtils.getExtension(fileEntry.getName());
				if (!ext.equals("swc") && !ext.equals("nrrd"))
					continue;
				String swcname = FilenameUtils.getBaseName(fileEntry.getName()) + ".swc";
				String mipname = FilenameUtils.getBaseName(fileEntry.getName()) + ".png";
				String srcswcpath = swcdirpath + swcname;
				String dstmippath = mipdirpath + mipname;
				if ( new File(srcswcpath).isFile() ) {
					IJ.run( "swc draw single", 
							"input=" + srcswcpath + " " +
							"output=" + dstmippath + " " +
							"width=1210 " +
							"height=566 " +
							"depth=174 " +
							"voxel_w=0.5189161 " +
							"voxel_h=0.5189161 " +
							"voxel_d=1.0000000 " +
							"radius=1");
					count++;
				}
			}
			IJ.log("Files: " + folder.listFiles().length);
			IJ.log("Generated MIPs: " + count);
			IJ.log("--------     (DONE)      --------\n");
/*
			IJ.log("-------- Generating MIPs --------");
			String mipdirpath = outdb.getDirectory() + "swc_prev" + File.separator;
			File mipdir = new File(mipdirpath);
			if (!mipdir.exists())
				mipdir.mkdirs();
			IJ.run("swc draw2", 
				"input=" + indir + " " +
				"output=" + mipdirpath + " " +
				"width=1210 " +
				"height=566 " +
				"depth=174 " +
				"voxel_w=0.5189161 " +
				"voxel_h=0.5189161 " +
				"voxel_d=1.0000000 " +
				"radius=1");
			IJ.log("--------     (DONE)      --------\n");
*/

			IJ.log("-------- Generating Thumbnails --------");
			count = 0;
			String thumbdirpath = outdb.getDirectory() + "swc_thumb" + File.separator;
			File thumbdir = new File(thumbdirpath);
			if (!thumbdir.exists())
				thumbdir.mkdirs();

			//final File folder = new File(dcin.getDirectory());
			for (final File fileEntry : folder.listFiles()) {
				String ext = FilenameUtils.getExtension(fileEntry.getName());
				if (!ext.equals("swc") && !ext.equals("nrrd"))
					continue;
				String mipname = FilenameUtils.getBaseName(fileEntry.getName()) + ".png";
				String srcmippath = mipdirpath + mipname;
				String dstmippath = thumbdirpath + mipname;
				if ( new File(srcmippath).isFile() ) {
					ImagePlus imp = IJ.openImage(srcmippath);
					IJ.run(imp, "Scale...", "x=- y=- width=153 height=76 interpolation=Bicubic average");
					IJ.run(imp, "Canvas Size...", "width=153 height=76 position=Center");
					IJ.saveAs(imp , "png", dstmippath);
					imp.close();
					count++;
				}
			}
			IJ.log("Files: " + folder.listFiles().length);
			IJ.log("Generated Thumbnails: " + count);
			IJ.log("--------        (DONE)        --------\n");

/*
			IJ.log("-------- Copying SWCs --------");
			String dstswcdirpath = outdb.getDirectory() + "swc" + File.separator;
			File dstswcdir = new File(dstswcdirpath);
			if (!dstswcdir.exists())
				dstswcdir.mkdirs();
			for (final File fileEntry : folder.listFiles()) {
				if (!FilenameUtils.getExtension(fileEntry.getName()).equals("swc"))
					continue;
				String dstpath = dstswcdirpath + fileEntry.getName();
				FileUtils.copyFile(fileEntry, new File(dstpath));
			}
			IJ.log("--------    (DONE)    --------\n");
*/
			
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
