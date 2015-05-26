/*
 * Copyright (C) 2014 Brian L. Browning
 *
 * This file is part of Beagle
 *
 * Beagle is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 *
 * Beagle is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with this program.  If not, see <http://www.gnu.org/licenses/>.
 */
package blbutil;

import java.io.BufferedInputStream;
import java.io.BufferedOutputStream;
import java.io.BufferedWriter;
import java.io.DataInputStream;
import java.io.DataOutputStream;
import java.io.File;
import java.io.FileInputStream;
import java.io.FileNotFoundException;
import java.io.FileOutputStream;
import java.io.FileWriter;
import java.io.IOException;
import java.io.OutputStream;
import java.io.PrintWriter;
import java.util.zip.GZIPOutputStream;
import net.sf.samtools.util.BlockCompressedOutputStream;

/**
 * Class {@code FileUtil} contains static methods for working with files.
 *
 * @author Brian L. Browning {@code <browning@uw.edu>}
 */
public class FileUtil {

    /**
     * The default number of bytes in the buffer for buffered file writers.
     * The value of this field is {@code 10 * (1 << 20)}.
     */
    public static final int DEFAULT_BUFFER_SIZE = 10 * (1 << 20);

    private FileUtil() {
        // private constructor prevents instantiation
    }

    /**
     * Returns a {@code java.io.PrintWriter} writing to
     * the specified file and having the specified buffer size.
     * The output will be compressed using the BGZIP compression algorithm.
     * Any existing file corresponding to the specified file will be deleted.
     * If the file cannot be opened, an error message will be printed and
     * the java interpreter will exit.
     *
     * @param file the file to be opened for output.
     * @param size the buffer size in bytes.
     * @return a {@code java.io.PrintWriter} writing to
     * the specified file.
     * @throws IllegalArgumentException if {@code size<=0}
     */
    public static PrintWriter bgzipPrintWriter(File file, int size) {
        PrintWriter out = null;
        try {
            OutputStream fout = new FileOutputStream(file);
            out = new PrintWriter(new BlockCompressedOutputStream(
                    new BufferedOutputStream(fout, size), file));
        } catch (FileNotFoundException e) {
            Utilities.exit("Error opening " + file, e);
        }
        return out;
    }

    /**
     * Returns a {@code java.io.PrintWriter} writing to
     * the specified file.  The resulting file will be compressed using
     * the BGZIP compression algorithm.  Any existing file corresponding
     * to the specified file will be deleted.  If the file
     * cannot be opened, an error message will be printed and the
     * java interpreter will exit.
     *
     * @param file the file to be opened for output.
     * @return a {@code java.io.PrintWriter} writing to
     * the specified file.
     */
    public static PrintWriter bgzipPrintWriter(File file) {
        return bgzipPrintWriter(file, DEFAULT_BUFFER_SIZE);
    }

    /**
     * Returns a data input stream reading from the specified file.  If the
     * input stream cannot be opened, an error message will be printed and the
     * java interpreter will exit.
     * @param file a binary file.
     * @return a data input stream reading from the specified file.
     */
    public static DataInputStream dataInputStream(File file) {
        DataInputStream dis = null;
        try {
            dis = new DataInputStream(new BufferedInputStream(
                    new FileInputStream(file), DEFAULT_BUFFER_SIZE));
        } catch (FileNotFoundException e) {
            Utilities.exit("Error opening " + file, e);
        }
        return dis;
    }

    /**
     * Returns a {@code java.io.DataOutputStream} writing to
     * the specified file and having the specified buffer size.
     * Any existing file will be overwritten. If the file cannot be opened,
     * an error message will be printed and the java interpreter will exit.
     * @param file the file to be opened for output.
     * @param size size of the buffer in bytes.
     * @return a {@code java.io.DataOutputStream} writing to
     * the specified file.
     */
    public static DataOutputStream dataOutputStream(File file, int size) {
        OutputStream dos = null;
        try {
            dos = new FileOutputStream(file);
        } catch (FileNotFoundException e) {
            Utilities.exit("Error opening " + file, e);
        }
        DataOutputStream out = new DataOutputStream(
                new BufferedOutputStream(dos, size));
        return out;
    }

    /**
     * Returns a {@code java.io.DataOutputStream} writing to
     * the specified file.  Any existing file corresponding to the
     * {@code File} object will be deleted.   If the file cannot be opened,
     * an error message will be printed and the java interpreter will exit.
     * @param file the file to be opened for output.
     * @return a {@code java.io.DataOutputStream} writing to
     * the specified file.
     */
    public static DataOutputStream dataOutputStream(File file) {
        return dataOutputStream(file, DEFAULT_BUFFER_SIZE);
    }

    /**
     * Returns a {@code java.io.PrintWriter} writing to
     * the specified file.  The resulting file will be compressed using
     * the GZIP compression algorithm.  Any existing file corresponding
     * to the specified file will be deleted.  If the file
     * cannot be opened, an error message will be printed and the
     * java interpreter will exit.
     * @param file the file to be opened for output.
     * @return a {@code java.io.PrintWriter} writing to
     * the specified file.
     */
    public static PrintWriter gzipPrintWriter(File file) {
        return gzipPrintWriter(file, DEFAULT_BUFFER_SIZE);
    }

    /**
     * Returns a {@code java.io.PrintWriter} writing to
     * the specified file and having the specified buffer size.
     * The resulting file will be compressed using the GZIP compression
     * algorithm.  Any existing file corresponding to the specified file
     * will be deleted.  If the file cannot be opened, an error message
     * will be printed and the java interpreter will exit.
     *
     * @param file the file to be opened for output.
     * @param size the buffer size in bytes.
     * @return a {@code java.io.PrintWriter} writing to
     * the specified file.
     * @throws IllegalArgumentException if {@code size<=0}
     */
    public static PrintWriter gzipPrintWriter(File file, int size) {
        PrintWriter out = null;
        try {
            out = new PrintWriter(
                    new GZIPOutputStream(new FileOutputStream(file), size));
        } catch (IOException e) {
            Utilities.exit("Error opening " + file, e);
        }
        return out;
    }

    /**
     * Returns a {@code java.io.PrintWriter} writing to
     * the specified file.  Any existing file corresponding
     * to the specified filename will be deleted.  If the file
     * cannot be opened, an error message will be printed and the
     * java interpreter will exit.
     * @param file the file to be opened for output.
     * @return a {@code java.io.PrintWriter} writing to
     * the specified file.
     */
    public static PrintWriter printWriter(File file) {
        return printWriter(file, false);
    }

    /**
     * Returns a {@code java.io.PrintWriter} that writes
     * to standard out.
     *
     * @return a {@code java.io.PrintWriter} that writes
     * to standard out.
     */
    public static PrintWriter stdOutPrintWriter() {
        return new PrintWriter(
                new BufferedOutputStream(System.out, DEFAULT_BUFFER_SIZE));
    }

    /**
     * Returns a {@code java.io.PrintWriter} writing to
     * the specified file.  If {@code append==false}
     * any existing file corresponding to the specified file will be deleted.
     * If the file cannot be opened, an error message will be printed and the
     * java interpreter will exit.
     *
     * @param file the file to be opened for output.
     * @param append if {@code true}, then data will be appended
     * to the end of any existing file.
     * @return a {@code java.io.PrintWriter} writing to
     * the specified file.
     */
    public static PrintWriter printWriter(File file, boolean append) {
        return printWriter(file, append, DEFAULT_BUFFER_SIZE);
    }

    /**
     * Returns a {@code java.io.PrintWriter} writing to
     * the specified file and having the specified buffer size.
     * If {@code append==false} any existing file corresponding
     * to the specified file will be deleted.  If {@code size==0}
     * there will be no buffering.  If the file cannot be opened,
     * an error message will be printed and the java interpreter will exit.
     *
     * @param file the file to be opened for output.
     * @param append if {@code true}, then data will be appended
     * to the end of any existing file.
     * @param size the buffer size in bytes.
     * @return a {@code java.io.PrintWriter} writing to
     * the specified file.
     * @throws IllegalArgumentException if {@code size < 0}.
     */
    public static PrintWriter printWriter(File file, boolean append, int size) {
        if (size==0) {
            return nonBufferedFileWriter(file, append);
        }
        PrintWriter out = null;
        try {
            out = new PrintWriter(
                    new BufferedWriter(new FileWriter(file, append), size));
        } catch (IOException e) {
            Utilities.exit("Error opening " + file, e);
        }
        return out;
    }

    private static PrintWriter nonBufferedFileWriter(File file, boolean append) {
        boolean autoflush = true;
        PrintWriter pw = null;
        try {
            pw = new PrintWriter(new FileWriter(file, append), autoflush);
        } catch (IOException e) {
            Utilities.exit("Error opening " + file, e);
        }
        return pw;
    }

    /**
     * Returns a temporary {@code File} that will be deleted when
     * the Java virtual machine exits.
     *
     * @param prefix the prefix string to be used in generating the file's
     * name.  The prefix must be at least three characters long.
     *
     * @return a {@code File} with an abstract pathname denoting a
     * newly created empty file.
     *
     * @throws IllegalArgumentException if the {@code prefix}
     * argument contains fewer than three characters.
     */
    public static File tempFile(String prefix) {
        File tempFile = null;
        try {
            tempFile = File.createTempFile(prefix, null);
            tempFile.deleteOnExit();
        } catch (IOException e) {
            Utilities.exit("Exception thrown by createTempFile: ", e);
        }
        return tempFile;
    }
}
