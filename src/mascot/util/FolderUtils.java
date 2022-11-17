package mascot.util;

import java.io.File;

public class FolderUtils {
	private static File lastDir = null;

//	    public static JFileChooser getFileChooser() {
//	        if(lastDir != null) {
//	            JFileChooser fc = new JFileChooser(lastDir);
//	            return fc;
//	        } else {
//	            JFileChooser fc = new JFileChooser();
//	            return fc;
//	        }
//	    }

	public static void setLastDirToParentOf(File file) {
		lastDir = file.getParentFile();
	}

	public static File getLastDir() {
		return lastDir;
	}
}
	
