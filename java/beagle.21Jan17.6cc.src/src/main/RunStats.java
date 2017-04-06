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

import vcf.Markers;
import vcf.Marker;
import vcf.Data;
import blbutil.Const;
import blbutil.FileUtil;
import blbutil.Utilities;
import dag.Dag;
import dag.DagUtil;
import java.io.File;
import java.io.PrintWriter;
import java.text.DecimalFormat;

/**
 * Class {@code RunStats} contains methods for storing and printing
 * statistics describing a Beagle analysis.
 *
 * @author Brian L. Browning {@code <browning@uw.edu>}
 */
public class RunStats {

    private static final DecimalFormat df2 = new DecimalFormat("0.00");

    private final Par par;
    private final PrintWriter log;
    private final long startNanos;

    private long buildNanos = 0;
    private long lastBuildNanos = 0;

    private long sampleNanos = 0;
    private long totalSampleNanos = 0;

    private long imputeNanos = 0;
    private long totalImputeNanos = 0;

    private long totalIbdNanos = 0;

    private String dagStats = null;

    /**
     * Constructs a new {@code RunStats} instance.
     * @param par the analysis parameters
     * @throws NullPointerException if {@code par == null}
     */
    RunStats(Par par) {
        this.startNanos = System.nanoTime();
        this.par = par;
        this.log = log(par.out());
    }

    private static PrintWriter log(String outPrefix) {
        File logFile = new File(outPrefix + ".log");
        boolean append = false;
        return FileUtil.nonBufferedPrintWriter(logFile, append);
    }

    /**
     * Prints initial information about the analysis to a log
     * file and to standard output.
     */
    public void printStartInfo() {
        Utilities.duoPrint(log, Main.shortHelp + Const.nl);
        Utilities.duoPrintln(log, "Start time: " + Utilities.timeStamp());
        Utilities.duoPrint(log, commandLine("beagle.jar", par.args()));
        if (par.ped() != null) {
            String s = Const.nl + "WARNING: This version will not model"
                    + " duos or trios in the pedigree file";
            Utilities.duoPrintln(log, s);
        }
        if (par.map() == null) {
            String s = Const.nl + "No genetic map is specified: using 1 cM = 1 Mb";
            Utilities.duoPrintln(log, s);
        }
        if (par.gt()==null && par.ref()!=null && par.impute()==true) {
            assert par.gl()!=null || par.gtgl()!=null;
            String s = Const.nl + "WARNING: Imputation of ungenotyped markers will not be performed."
                     + Const.nl + "         Imputation requires the \"gt=\" argument and called genotypes.";
            Utilities.duoPrintln(log, s);
        }
        if (par.gt()==null && par.ibd()) {
            assert par.gl()!=null || par.gtgl()!=null;
            String s = Const.nl + "WARNING: IBD segment detection will not be performed."
                     + Const.nl + "         IBD analysis requires the \"gt=\" argument and called genotypes.";
            Utilities.duoPrintln(log, s);
        }
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
    private static String commandLine(String jarFile, String[] args) {
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
     * Prints information about the complete analysis to a log
     * file and to standard output, and closes the log file.
     * @param nTargetMarkers the total number of target markers analyzed
     * @param nMarkers the total number of markers analyzed
     */
    public void printSummaryAndClose(int nTargetMarkers, int nMarkers) {
        long totalTime = System.nanoTime() - startNanos;
        if (nTargetMarkers == nMarkers) {
            Utilities.duoPrint(log, Const.nl);
            Utilities.duoPrint(log, "Number of markers:             ");
            Utilities.duoPrintln(log, String.format("%7d", nMarkers));
        }
        else {
            Utilities.duoPrint(log, Const.nl);
            Utilities.duoPrint(log, "Number of reference markers:   ");
            duoPrintln7d(nMarkers);
            Utilities.duoPrint(log, "Number of target markers:      ");
            duoPrintln7d(nTargetMarkers);
        }
        if (buildNanos > 0) {
            duoPrintNanos("Total time for building model: ", buildNanos);
        }
        if (totalSampleNanos > 1000) {
            duoPrintNanos("Total time for sampling:       ", totalSampleNanos);
        }
        if (par.ibd()==true) {
            duoPrintNanos("Total time for IBD detection:  ", totalIbdNanos);
        }
        if (totalImputeNanos > 0) {
            duoPrintNanos("Total time for imputation:     ", totalImputeNanos);
        }
        duoPrintNanos("Total run time:                ", totalTime);
        Utilities.duoPrintln(log, Const.nl + "End time: "
                + Utilities.timeStamp());
        Utilities.duoPrintln(log, Main.program + " finished");
        log.close();
    }

    /**
     * Increases the cumulative time to build the DAG models by the
     * specified number of nanoseconds.
     * @param nanos the nanoseconds required to build an instance
     * of the DAG model
     */
    public void buildNanos(long nanos) {
        buildNanos += nanos;
    }

    /**
     * Stores the time for sampling new haplotypes and increases the
     * cumulative sampling time by the specified number of nanoseconds.
     * @param nanos the nanoseconds required to sample new haplotypes
     */
    public void sampleNanos(long nanos) {
        sampleNanos = nanos;
        totalSampleNanos += nanos;
    }

    /**
     * Stores the time for imputing ungenotyped marker and increases
     * the cumulative imputation time by the specified number
     * of nanoseconds.
     * @param nanos the nanoseconds required to impute ungenotyped
     * markers
     */
    public void imputationNanos(long nanos) {
        imputeNanos = nanos;
        totalImputeNanos += nanos;
    }

   /**
     * Increases the cumulative time for detecting identity-by-descent
     * by the specified number of nanoseconds.
     * @param nanos the nanoseconds required to perform IBD detection
     */
    public void ibdNanos(long nanos) {
        totalIbdNanos += nanos;
    }

    /**
     * Stores statistics for the DAG model used to sample single individuals.
     * @param dag the DAG model used to sample individuals
     */
    public void setDagStats(Dag dag) {
        dagStats = (dag==null) ? null : DagUtil.dagStats(dag);
    }

    /**
     * Prints information about the Refined IBD analysis to a log
     * file and to standard output.
     * @param ibdScale the value used to multiplicatively scales the node
     * similarity threshold when building the DAG model
     * @param ibdDag the DAG model used for IBD detection
     * @param nanos the nanoseconds required for IBD detection
     * @throws NullPointerException if {@code ibdDag == null}
     */
    public void printRefinedIbdUpdate(float ibdScale, Dag ibdDag, long nanos) {
        Utilities.duoPrintln(log, Const.nl + "Refined IBD");
        Utilities.duoPrintln(log, "model scale: " + df2.format(ibdScale));
        duoPrintNanos("run time:    ", nanos);
        Utilities.duoPrint(log, Const.nl + DagUtil.dagStats(ibdDag));
    }

    /**
     * Prints run time for most recent imputation to a log file
     * and to standard output.
     */
    public void printImputationUpdate() {
        Utilities.duoPrint(log, Const.nl);
        duoPrintNanos("Imputation time (this window): ", imputeNanos);
    }

   /**
     * Prints information about the samples to a log
     * file and to standard output.
     * @param fam the parent-offspring relationships
     * @param data the input genotype data
     */
    public void printSampleSummary(NuclearFamilies fam, Data data) {
        Utilities.duoPrint(log, Const.nl);
        Utilities.duoPrint(log, "reference samples: ");
        duoPrintln7d(data.nRefSamples());
        Utilities.duoPrint(log, "target samples:    ");
        duoPrintln7d(data.nTargetSamples());
        if (par.ped() != null) {
            Utilities.duoPrint(log, "  ");
            Utilities.duoPrint(log, String.valueOf(fam.nSingles()));
            Utilities.duoPrintln(log, " singles");
            Utilities.duoPrint(log, "  ");
            Utilities.duoPrint(log, String.valueOf(fam.nDuos()));
            Utilities.duoPrintln(log, " duos");
            Utilities.duoPrint(log, "  ");
            Utilities.duoPrint(log, String.valueOf(fam.nTrios()));
            Utilities.duoPrintln(log, " trios");
        }
    }

   /**
     * Prints information about the marker window to a log
     * file and to standard output.
     * @param data the input genotype data
     */
    public void printWindowUpdate(Data data) {
        Markers markers = data.markers();
        Marker first = markers.marker(0);
        Marker last = markers.marker(markers.nMarkers() - 1);
        StringBuilder sb = new StringBuilder(30);
        sb.append(Const.nl);
        sb.append("Window ");
        sb.append(data.window());
        sb.append(" [ ");
        String chr = first.chrom();
        if (chr.equals(Const.MISSING_DATA_STRING)==false) {
            sb.append(chr);
            sb.append(Const.colon);
        }
        sb.append(first.pos());
        sb.append(Const.hyphen);
        if (chr.equals(last.chrom())==false) {
            sb.append(last.chrom());
            sb.append(Const.colon);
        }
        sb.append(last.pos());
        sb.append(" ]");
        sb.append(Const.nl);
        if (data.nRefSamples()>0) {
            sb.append("reference markers: ");
            sb.append(String.format("%7d", data.nMarkers()));
            sb.append(Const.nl);
        }
        sb.append("target markers:    ");
        sb.append(String.format("%7d", data.nTargetMarkers()));
        Utilities.duoPrintln(log, sb.toString());
    }

    /**
     * Prints the specified string to the log file and to standard out.
     * @param msg the message to be printed
     */
    public void println(String msg) {
        Utilities.duoPrintln(log, msg);
    }

    /**
     * Prints information about the specified iteration.
     * @param window the window
     * @param iter the iteration
     */
    public void printIterationUpdate(int window, int iter) {
        long buildTime = buildNanos - lastBuildNanos;
        lastBuildNanos = buildNanos;
        Utilities.duoPrint(log, Const.nl + "Window=" + window
                + " Iteration=" + iter + Const.nl);
        duoPrintNanos("Time for building model:         ", buildTime);
        if (dagStats != null) {
            duoPrintNanos("Time for sampling (singles):     ", sampleNanos);
            sampleNanos = 0;
        }
        if (dagStats != null) {
            Utilities.duoPrint(log, "DAG statistics" + Const.nl);
            Utilities.duoPrint(log, dagStats);
        }
        log.flush();
    }

    /**
     * Returns a string with specified message following by the elapsed time
     * (in hours, minutes, and seconds).
     *
     * @param message a message preceding the elapsed time
     * @param nanos the number of elapsed nanoseconds
     *
     * @return a string with specified message following by the elapsed time
     * (in hours, minutes, and seconds)
     */
    private static String elapsedNanos(String message, long nanos) {
        StringBuilder sb = new StringBuilder(message.length() + 30);
        sb.append(message);
        sb.append(Utilities.elapsedNanos(nanos));
        sb.append(Const.nl);
        return sb.toString();
    }

    private void duoPrintNanos(String message, long nanos) {
        Utilities.duoPrint(log, elapsedNanos(message, nanos));
    }

    private void duoPrintln7d(int i) {
        Utilities.duoPrintln(log, String.format("%7d", i));
    }
}
