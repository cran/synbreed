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
package main;

import blbutil.FileUtil;
import java.io.File;
import java.io.PrintWriter;

/**
 * <p>Class {@code Logger} is a light-weight singleton class for writing
 * log message.
 * </p>
 * Class {@code Logger} is thread-safe.
 *
 * @author Brian L. Browning {@code <browning@uw.edu>}
 */
public class Logger {

    private static final Logger instance = new Logger();

    private File file = null;
    private boolean append = true;
    private PrintWriter out = open(file, append);
    private int logCnt = 0;

    private Logger() {
        // private constructor prevents instantiation
    }

    /**
     * Returns the singleton Logger instance.  If the {@code setTarget()}
     * method has not yet been invoked, log messages will be written to
     * standard output.
     * @return the singleton Logger instance.
     */
    public static Logger getInstance() {
        return instance;
    }

    /**
     * Sets the file to which log messages will be written.
     *
     * @param target a target file to which log messages will be written,
     * or {@code null} if log messages will be written to standard out.
     *
     * @param append {@code true} if log messages will be written
     * to the end of the specified  file, and {@code false} if log messages
     * will overwrite any existing file.  The target file is not created
     * or modified unless one or more log messages are written.
     * The {@code append} parameter has no effect if {@code target==null}.
     */
    public synchronized void setTarget(File target, boolean append) {
        this.close();
        this.file = target;
        this.append = append;
    }

    /**
     * Sets the file to which log messages will be written.
     * If {@code target!=null}, log messages will overwrite any
     * existing file.
     *
     * @param target a target file to which log messages will be written,
     * or {@code null} if log messages will be written to standard out.
     */
    public synchronized void setTarget(File target) {
        boolean appendFile = false;
        setTarget(target, appendFile);
    }

    /**
     * Returns the file to which log messages will be written, or returns
     * {@code null} if log messages will be written to standard out.
     *
     * @return the file to which log messages will be written, or
     * {@code null} if log messages will be written to standard out.
     */
    public synchronized File target() {
        return file;
    }

    private static PrintWriter open(File file, boolean append) {
        if (file==null) {
            return FileUtil.stdOutPrintWriter();
        }
        else {
            return FileUtil.printWriter(file, append);
        }
    }

    /**
     * Closes the current log target, and sets the log target to standard out.
     */
    public synchronized void close() {
        if (logCnt > 0) {
            if (file!=null) {
                out.close();
                file=null;
            }
            else {
                out.flush();
            }
        }
        out = null;
        logCnt = 0;
    }

    /**
     * Writes the specified log message and appends the new line character.
     * @param msg a log message.
     * @throws NullPointerException if {@code msg==null}.
     */
    public synchronized void println(String msg) {
        if (logCnt==0) {
            this.out = open(file, append);
        }
        ++logCnt;
        if (out!=null) {
            out.println(msg);
        }
    }

    /**
     * Writes the specified log message without appending a new line character.
     * @param msg a log message.
     * @throws NullPointerException if {@code msg==null}.
     */
    public synchronized void print(String msg) {
        if (logCnt==0) {
            this.out = open(file, append);
        }
        ++logCnt;
        if (out!=null) {
            out.print(msg);
        }
    }

    /**
     * Returns the number of {@code Logger.print()} and
     * {@code Logger.println()} method calls for the current log target.
     *
     * @return the number of {@code Logger.print()} and
     * {@code Logger.println()} method calls for the current log target.
     *
     * @see #setTarget(java.io.File)
     * @see #close()
     */
    public synchronized int logCnt() {
        return logCnt;
    }
}
