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

import java.io.File;
import java.io.PrintWriter;
import java.text.SimpleDateFormat;
import java.util.Collections;
import java.util.Date;
import java.util.HashSet;
import java.util.Set;

/**
 * Class {@code Utilities} contains miscellaneous static utility methods.
 *
 * @author Brian L. Browning {@code <browning@uw.edu>}
 */
public class Utilities {

    private Utilities() {
        // private constructor to prevent instantiation
    }

    /**
     * Prints a summary of memory use at the time of method invocation
     * to standard output.
     * @param msg a string a message to be printed with the summary
     * of memory use
     */
    public static void printMemoryUse(String msg) {
        long Mb = 1024*1024;
        Runtime rt = Runtime.getRuntime();
        System.out.println(Const.nl + msg
                + Const.tab + "maxMb=" + (rt.maxMemory()/Mb)
                + Const.tab + "totalMb=" + (rt.totalMemory()/Mb)
                + Const.tab + "freeMb=" + (rt.freeMemory()/Mb)
                + Const.tab + "usedMb=" + ((rt.totalMemory() - rt.freeMemory())/Mb));
    }

    /**
     * Returns the current local time as a string.  The
     * exact details of the string representation
     * are unspecified and subject to change.
     *
     * @return the current local time as a string.
     */
    public static String timeStamp() {
        Date now = new Date();
        SimpleDateFormat sdf =
                new SimpleDateFormat("hh:mm a z 'on' dd MMM yyyy");
        return sdf.format(now);
    }

    /**
     * <p>Returns a set of identifiers found in a text file that has
     * one identifier per line.  The empty set is returned if
     * {@code file == null}. Blank lines are ignored, and white-space that
     * begins or ends a line is ignored.
     * </p>
     * If an {@code IOException} is thrown, an error message is printed
     * to standard error and the Java virtual machine is forced to terminate.
     *
     * @param file a text file with one identifier per line
     * @return a set of identifiers
     *
     * @throws IllegalArgumentException if the specified file does not exist
     * @throws IllegalArgumentException if the specified file is a directory
     * @throws IllegalArgumentException if any line of the specified
     * file contains two non-white-space characters separated by one or
     * more white-space characters
     */
    public static Set<String> idSet(File file) {
        if (file==null) {
            return Collections.emptySet();
        }
        else {
            if (file.exists()==false) {
                String s = "file does not exist: " + file;
                throw new IllegalArgumentException(s);
            }
            if (file.isDirectory()) {
                String s = "file is a directory: " + file;
                throw new IllegalArgumentException(s);
            }
            Set<String> idSet = new HashSet<>();
            try (FileIt<String> it = InputIt.fromGzipFile(file)) {
                while (it.hasNext()) {
                    String line = it.next().trim();
                    if (line.length() > 0) {
                        if (StringUtil.countFields(line) > 1) {
                            String s = "line has >1 white-space delimited fields: "
                                    + line;
                            throw new IllegalArgumentException(s);
                        }
                        idSet.add(line);
                    }
                }
            }
            return idSet;
        }
    }

    /**
     * Prints the specified string to the specified {@code PrintWriter} and
     * to standard out.  The line separator string is not appended to the
     * specified string before printing.
     *
     * @param out a print writer
     * @param s a string to be printed
     *
     * @throws NullPointerException if {@code out == null}
     */
    public static void duoPrint(PrintWriter out, CharSequence s) {
        System.out.print(s);
        out.print(s);
    }

   /**
     * Prints the specified string to the specified {@code PrintWriter} and
     * to standard out.  The line separator string is appended to the
     * specified string before printing.
     *
     * @param out a print writer
     * @param s a string to be printed
     *
     * @throws NullPointerException if {@code out == null}
     */
    public static void duoPrintln(PrintWriter out, CharSequence s) {
        System.out.println(s);
        out.println(s);
    }

     /**
     * Returns a string representation of the specified elapsed time
     * in the format "H hours M minutes S seconds".
     *
     * @param nanoseconds the elapsed time in nanoseconds
     *
     * @return a string representation of the specified elapsed time
     */
    public static String elapsedNanos(long nanoseconds) {
        long seconds = Math.round(nanoseconds /1000000000.0);
        StringBuilder sb = new StringBuilder(80);
        if (seconds >= 3600) {
            long hours = seconds / 3600;
            sb.append(hours);
            sb.append(hours==1 ? " hour " : " hours ");
            seconds %= 3600;

        }
        if (seconds >= 60) {
            long minutes = seconds / 60;
            sb.append(minutes);
            sb.append(minutes==1 ? " minute " : " minutes ");
            seconds %= 60;
        }
        sb.append(seconds);
        sb.append(seconds==1 ? " second" : " seconds");
        return sb.toString();
    }

    /**
     * Prints the specified exception, its stack trace, and
     * the specified string to standard out and then terminates the
     * Java virtual machine.
     *
     * @param s a string to be printed to standard err
     * @param e an exception or error to be printed to standard err
     *
     * @throws NullPointerException if {@code e == null}
     */
    public static void exit(String s, Throwable e) {
        e.printStackTrace(System.out);
        System.out.println(e);
        System.out.println(s);
        System.out.println("terminating program.");
        System.exit(1);
    }

    /**
     * Prints the specified string to standard out and then terminates the
     * Java virtual machine.
     *
     * @param s a string to be written to standard output
     */
    public static void exit(String s) {
        System.out.println(s);
        System.out.flush();
        System.exit(0);
    }
}

