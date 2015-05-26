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
import java.io.BufferedReader;
import java.io.File;
import java.io.FileInputStream;
import java.io.FileNotFoundException;
import java.io.IOException;
import java.io.InputStream;
import java.io.InputStreamReader;
import java.util.NoSuchElementException;
import java.util.zip.GZIPInputStream;
import net.sf.samtools.util.BlockCompressedInputStream;

/**
 * <p>Class {@code InputIterator} is an iterator whose {@code next()}
 * method returns lines of a text input stream.  Class
 * {@code InputIterator} uses buffering when reading from the text
 * input stream.
 * </p>
 * If an {@code IOException} is thrown when an {@code InputIterator}
 * instance reads from the text input stream, the {@code IOException}
 * is trapped, an error message is written to standard out, and the
 * Java Virtual Machine is terminated.
 *
 * @author Brian L. Browning {@code <browning@uw.edu>}
 */
public class InputIterator implements FileIterator<String> {

    /**
     * The default buffer size, which is 8,388,608 bytes.
     */
    public static final int DEFAULT_BUFFER_SIZE = 1<<22;

    private final BufferedReader in;
    private String next = null;

    /**
     * Constructs a new {@code InputStreamIterator} with default buffer
     * size that will iterate through lines of the specified input stream.
     *
     * @param is input stream of text data.
     *
     * @see #DEFAULT_BUFFER_SIZE
     */
    public InputIterator(InputStream is) {
        this(is, DEFAULT_BUFFER_SIZE);
    }

    /**
     * Constructs a new {@code InputStreamIterator} that will iterate through
     * the lines of the specified input stream.
     *
     * @param is input stream of text data.
     * @param bufferSize the buffer size in bytes.
     *
     * @throws IllegalArgumentException if {@code bufferSize<0}.
     */
    public InputIterator(InputStream is, int bufferSize) {
        BufferedReader br = null;
        try {
            InputStreamReader isr = new InputStreamReader(is);
            br = new BufferedReader(isr, bufferSize);
            next = br.readLine();
        }
        catch(IOException e) {
            Utilities.exit("Error reading " + is, e);
        }
        this.in = br;
    }

    @Override
    public File file() {
        return null;
    }

    /**
     * Returns {@code true} if the iteration has more elements.
     * @return {@code true} if the iteration has more elements.
     */
    @Override
    public boolean hasNext() {
        return (next != null);
    }

    /**
     * Returns the next element in the iteration.
     * @return the next element in the iteration.
     * @throws NoSuchElementException if the iteration has no more elements.
     */
    @Override
    public String next() {
        if (!hasNext()) {
            throw new NoSuchElementException();
        }
        String current = next;
        try {
            next = in.readLine();
        }
        catch (IOException e) {
            Utilities.exit("Error reading " + in, e);
        }
        return current;
    }

    /**
     * The {@code remove} method is not supported by this iterator.
     * @throws UnsupportedOperationException if this method is invoked.
     */
    @Override
    public void remove() {
        String s = "remove() is not supported by LineIterator";
        throw new UnsupportedOperationException(s);
    }

    @Override
    public void close() {
        try {
            in.close();
        }
        catch (IOException e) {
            Utilities.exit("Error closing " + in, e);
        }
        next=null;
    }

    /**
     * Returns a string representation of this iterator.  The exact details
     * of the representation are unspecified and subject to change.
     * @return a string representation of this iterator.
     */
    @Override
    public String toString() {
        StringBuilder sb = new StringBuilder(200);
        sb.append("[next = ");
        sb.append(next);
        sb.append("\nBufferedReader = ");
        sb.append(in);
        return sb.toString();
    }

    /**
     * Returns a buffered {@code InputIterator} that iterates through
     * the lines of the specified compressed or uncompressed text file.
     * The buffer size is {@code InputIterator.DEFAULT_BUFFER_SIZE}.
     * If the filename ends in ".gz", the file must be either BGZIP-compressed
     * or GZIP-compressed.
     *
     * @param file a compressed or uncompressed text file.
     * @return a buffered {@code FileIterator<String>} that iterates through
     * the lines of the specified compressed or uncompressed text file.
     *
     * @throws NullPointerException if {@code file==null}.
     */
    public static InputIterator fromGzipFile(File file) {
        try {
            InputStream is = new FileInputStream(file);
            if (file.getName().endsWith(".gz")) {
                if (isBGZipFile(file)) {
                    return new InputIterator(
                            new BlockCompressedInputStream(is));
                }
                else {
                    return new InputIterator(new GZIPInputStream(is));
                }
            }
            else {
                return new InputIterator(is);
            }
        }
        catch(FileNotFoundException e) {
            Utilities.exit("Error opening " + file, e);
        }
        catch(IOException e) {
            Utilities.exit("Error reading " + file, e);
        }
        assert false;
        return null;
    }

    private static boolean isBGZipFile(File file) throws IOException {
        boolean result = false;
        try (InputStream is=new BufferedInputStream(new FileInputStream(file))) {
            result=BlockCompressedInputStream.isValidFile(is);
        }
        return result;
     }

    /**
     * Returns a buffered {@code InputIterator} that iterates through
     * the lines of the specified compressed or uncompressed text file.
     * The buffer size is {@code InputIterator.DEFAULT_BUFFER_SIZE}.
     * If the filename ends in ".gz", the file must be either BGZIP-compressed
     * or GZIP-compressed.
     *
     * @param filename a file name of a compressed or uncompressed text file.
     * @return a buffered {@code FileIterator<String>} that iterates through
     * the lines of the specified compressed or uncompressed text file.
     *
     * @throws NullPointerException if {@code file==null}.
     */
    public static InputIterator fromGzipFile(String filename) {
        return fromGzipFile(new File(filename));
    }

     /**
     * Returns a buffered {@code InputIterator} that iterates through
     * the lines of the specified text file.  The buffer size is
     * {@code InputIterator.DEFAULT_BUFFER_SIZE}.
     *
     * @param file a text file.
     * @return a buffered {@code FileIterator<String>} that iterates through
     * the lines of the specified text file.
     *
     * @throws NullPointerException if {@code filename==null}.
     */
    public static InputIterator fromTextFile(File file) {
        try {
            return new InputIterator(new FileInputStream(file));
        }
        catch(FileNotFoundException e) {
            Utilities.exit("Error opening " + file, e);
        }
        assert false;
        return null;
    }

    /**
     * Returns a buffered {@code InputIterator} that iterates through
     * the lines of the specified text file.  The buffer size is
     * {@code InputIterator.DEFAULT_BUFFER_SIZE}.
     *
     * @param filename name of a text file.
     * @return a buffered {@code FileIterator<String>} that iterates through
     * the lines of the specified text file.
     *
     * @throws NullPointerException if {@code filename==null}.
     */
    public static InputIterator fromTextFile(String filename) {
        return fromTextFile(new File(filename));
    }
}
