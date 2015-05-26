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
 * Class {@code Utilities} contains static utility methods.
 *
 * @author Brian L. Browning {@code <browning@uw.edu>}
 */
public class Utilities {

    private Utilities() {
        // private constructor to prevent instantiation
    }

   /**
    * Returns a string representation of the command line arguments.
    * The exact details of the representation are unspecified and
    * subject to change.
    *
    * @param jarFile the name of the program's jar file.
    * @param args command line arguments.
    * @return a string representation of the command line arguments.
    */
    public static String commandLine(String jarFile, String[] args) {
        StringBuilder sb = new StringBuilder(args.length*20);
        long maxMemory = Runtime.getRuntime().maxMemory();
        sb.append(Const.nl);
        sb.append("Command line: java");
        if (maxMemory!=Long.MAX_VALUE) {
            long maxMb = maxMemory / (1024*1024);
            sb.append(" -Xmx");
            sb.append(maxMb);
            sb.append("m");
        }
        sb.append(" -jar ");
        sb.append(jarFile);
        sb.append(Const.nl);
        for (int j = 0; j < args.length; ++j) {
            sb.append("  ");
            sb.append(args[j]);
            sb.append(Const.nl);
        }
        return sb.toString();
    }

    /**
     * Prints a summary of memory use at the time of method invocation
     * to standard output.
     * @param msg a string a message to be printed with the summary
     * of memory use.
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
     * Returns a string representing the current local time.  The
     * exact details of the representation are unspecified and
     * subject to change.
     *
     * @return a string representing the current local time.
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
     * {@code file==null}. Blank lines are ignored, and white-space that
     * begins or ends a line is ignored.
     * </p>
     * If an {@code IOException} is thrown, an error message is printed
     * to standard error and the Java virtual machine is forced to terminate.
     *
     * @param file a text file with one identifier per line.
     * @return a set of identifiers.
     *
     * @throws IllegalArgumentException if the specified file does not exist,
     * if the specified file is a directory, or if any line of the specified
     * file contains two non-white-space characters separated by one or
     * more white-space characters.
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
            try (FileIterator<String> it=InputIterator.fromGzipFile(file)) {
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
     * to standard out.  A line separator character is not appended to the
     * specified string.
     *
     * @param out a character-output stream.
     * @param s a string to be printed.
     *
     * @throws NullPointerException if {@code out==null}.
     */
    public static void duoPrint(PrintWriter out, String s) {
        System.out.print(s);
        out.print(s);
    }

   /**
     * Prints the specified string to the specified {@code PrintWriter} and
     * to standard out.  A line separator character is appended to the
     * specified string.
     *
     * @param out a character-output stream.
     * @param s a string to be printed.
     *
     * @throws NullPointerException if {@code out==null}.
     */
    public static void duoPrintln(PrintWriter out, String s) {
        System.out.println(s);
        out.println(s);
    }

     /**
     * Returns a string representation of the specified elapsed time
     * in the format "H hours M minutes S seconds".
     *
     * @param milliseconds the elapsed time in milliseconds
     *
     * @return a string representation of the specified elapsed time.
     */
    public static String elapsedMillis(long milliseconds) {
        return elapsedNanos(1000000*milliseconds);
    }

     /**
     * Returns a string representation of the specified elapsed time
     * in the format "H hours M minutes S seconds".
     *
     * @param nanoseconds the elapsed time in nanoseconds
     *
     * @return a string representation of the specified elapsed time.
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
     * Returns a string with specified message following by the elapsed time
     * (in hours, minutes, and seconds).
     *
     * @param message a message preceding the elapsed time.
     * @param millis the number of elapsed milliseconds.
     *
     * @return a string with specified message following by the elapsed time
     * (in hours, minutes, and seconds).
     *
     * @throws NullPointerException if {@code message==null}.
     */
    public static String printElapsedTime(String message, long millis) {
        StringBuilder sb = new StringBuilder(message.length() + 30);
        sb.append(message);
        sb.append(Utilities.elapsedMillis(millis));
        sb.append(Const.nl);
        return sb.toString();
    }

    /**
     * Prints the specified exception (including stack trace) and
     * the specified string to standard error, and then terminates the
     * Java virtual machine.
     *
     * @param s a message to be printed to standard err.
     * @param e an exception or error to be printed to standard err.
     *
     * @throws NullPointerException if {@code e==null}.
     */
    public static void exit(String s, Throwable e) {
        e.printStackTrace(System.err);
        System.err.println(e);
        System.err.println(s);
        System.err.println("terminating program.");
        System.exit(1);
    }

    /**
     * Prints the specified message to standard out and then terminates the
     * Java virtual machine.
     *
     * @param s a message to be written to standard output.
     */
    public static void exit(String s) {
        System.out.println(s);
        System.out.println();
        System.out.flush();
        System.exit(0);
    }
}

