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
import dag.DagUtils;
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

    private static final String shortHelp = Main.version + Const.nl
            + Main.copyright + Const.nl
            + "Enter \"java -jar beagle.jar\" for a summary of command line "
            + "arguments.";
    private static final DecimalFormat df2 = new DecimalFormat("0.00");

    private final Parameters par;
    private final PrintWriter log;
    private final long startMillis;

    private long buildMillis = 0;
    private long lastBuildMillis = 0;

    private long singleSampleMillis = 0;
    private long duoSampleMillis = 0;
    private long trioSampleMillis = 0;
    private long totalSampleMillis = 0;

    private long totalIbdMillis = 0;

    private String singleDagStats = null;
    private String duoDagStats = null;
    private String trioDagStats = null;

    /**
     * Constructs a new {@code RunStats} instance.
     * @param par the analysis parameters.
     * @throws NullPointerException if {@code par==null}
     */
    RunStats(Parameters par) {
        this.startMillis = System.currentTimeMillis();
        this.par = par;
        this.log = log(par.out());
    }

    private static PrintWriter log(String outPrefix) {
        File logFile = new File(outPrefix + ".log");
        boolean append = false;
        int printBuffer = 0;
        return FileUtil.printWriter(logFile, append, printBuffer);
    }

    /**
     * Prints information about the start of the analysis to a log
     * file and to standard output.
     */
    void printStartInfo() {
        Utilities.duoPrint(log, shortHelp + Const.nl);
        Utilities.duoPrintln(log, "Start time: " + Utilities.timeStamp());
        Utilities.duoPrint(log, Utilities.commandLine("beagle.jar", par.args()));
    }

    /**
     * Prints information about the complete analysis to a log
     * file and to standard output, and closes the log file.
     * @param cumMarkerCnt the total number of markers analyzed.
     */
    void printSummaryAndClose(int cumMarkerCnt) {
        long totalTime = System.currentTimeMillis() - startMillis;
        Utilities.duoPrint(log, Const.nl + "Number of markers:             "
                + cumMarkerCnt + Const.nl);
        Utilities.duoPrint(log,
                Utilities.printElapsedTime(    "Total time for building model: ",
                        buildMillis));
        Utilities.duoPrint(log,
                Utilities.printElapsedTime(    "Total time for sampling:       ",
                        totalSampleMillis));
        if (par.ibd()==true) {
            Utilities.duoPrint(log,
                Utilities.printElapsedTime(    "Total time for IBD detection:  ",
                        totalIbdMillis));
        }
        Utilities.duoPrint(log,
                Utilities.printElapsedTime(    "Total run time:                ",
                        totalTime));
        Utilities.duoPrintln(log, Const.nl + "End time: "
                + Utilities.timeStamp());
        Utilities.duoPrintln(log, Main.version + " finished");
        log.close();
    }

    /**
     * Increases the cumulative time to build the DAG models by the
     * specified amount.
     * @param milliseconds the milliseconds required to build an instance
     * of the DAG model.
     */
    void buildMillis(long milliseconds) {
        buildMillis += milliseconds;
    }

    /**
     * Stores the time for sampling new haplotypes for single individuals
     * and increases the cumulative sampling time by this amount.
     * @param milliseconds the milliseconds required to sample new haplotypes
     * for single individuals.
     */
    void singleSampleMillis(long milliseconds) {
        singleSampleMillis = milliseconds;
        totalSampleMillis += milliseconds;
    }

   /**
     * Stores the time for sampling new haplotypes for parent-offspring
     * duos and increases the cumulative sampling time by this amount.
     * @param milliseconds the milliseconds required to sample new haplotypes
     * for parent-offspring duos.
     */
    void duoSampleMillis(long milliseconds) {
        duoSampleMillis = milliseconds;
        totalSampleMillis += milliseconds;
    }

   /**
     * Stores the time for sampling new haplotypes for parent-offspring
     * trio and increases the cumulative sampling time by the this amount.
     * @param milliseconds the milliseconds required to sample new haplotypes
     * for parent-offspring trios.
     */
    void trioSampleMillis(long milliseconds) {
        trioSampleMillis = milliseconds;
        totalSampleMillis += milliseconds;
    }

   /**
     * Increases the cumulative time for detecting identity-by-descent
     * by the specified amount.
     * @param milliseconds the milliseconds required to sample new haplotypes
     * for parent-offspring duos.
     */
    void ibdMillis(long milliseconds) {
        totalIbdMillis += milliseconds;
    }

    /**
     * Stores statistics for the DAG model used to sample single individuals.
     * @param dag the DAG model used to sample individuals.
     */
    void setSingleDagStats(Dag dag) {
        singleDagStats = (dag==null) ? null : DagUtils.dagStats(dag);
    }

    /**
     * Stores statistics for the DAG model used to sample parent-offspring duos.
     * @param dag the DAG model used to sample parent-offspring duos.
     */
    void setDuoDagStats(Dag dag) {
        duoDagStats = (dag==null) ? null : DagUtils.dagStats(dag);
    }

    /**
     * Stores statistics for the DAG model used to sample parent-offspring trios.
     * @param dag the DAG model used to sample parent-offspring trios.
     */
    void setTrioDagStats(Dag dag) {
        trioDagStats = (dag==null) ? null : DagUtils.dagStats(dag);
    }

    /**
     * Prints information about the Refined IBD analysis to a log
     * file and to standard output.
     * @param ibdScale the value used to multiplicatively scales the node
     * similarity threshold when building the DAG model.
     * @param ibdDag the DAG model used for IBD detection.
     * @param ibdMillis the milliseconds required for IBD detection.
     */
    void printRefinedIbdUpdate(float ibdScale, Dag ibdDag, long ibdMillis) {
        Utilities.duoPrint(log, Const.nl + "Refined IBD");
        Utilities.duoPrint(log, Const.nl + "model scale: " + df2.format(ibdScale));
        Utilities.duoPrint(log, Const.nl
                + Utilities.printElapsedTime("run time:    ", ibdMillis));
        Utilities.duoPrint(log, Const.nl + DagUtils.dagStats(ibdDag));
    }

   /**
     * Prints information about the samples to a log
     * file and to standard output.
     * @param the parent-offspring relationships.
     * @param data the input genotype data.
     */
    void printSampleSummary(NuclearFamilies fam, Data data) {
        Utilities.duoPrint(log, Const.nl + "reference samples: "
                + data.nRefSamples());
        Utilities.duoPrint(log, Const.nl + "target samples:    "
                + data.nNonRefSamples());
        Utilities.duoPrint(log, Const.nl + "  singles: "
                + fam.nSingles());
        Utilities.duoPrint(log, Const.nl + "  duos:    "
                + fam.nDuos());
        Utilities.duoPrintln(log, Const.nl +"  trios:   "
                + fam.nTrios());
    }

   /**
     * Prints information about the marker window to a log
     * file and to standard output.
     * @param data the input genotype data.
     */
    void printWindowUpdate(Data data) {
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
            sb.append(data.nMarkers());
            sb.append(Const.nl);
        }
        sb.append("target markers:    ");
        sb.append(data.nNonRefMarkers());
        Utilities.duoPrintln(log, sb.toString());
    }

    /**
     * Prints the specified string to the log file and to standard out.
     * @param msg the message to be printed.
     */
    void println(String msg) {
        Utilities.duoPrintln(log, msg);
    }

    /**
     * Prints information about the specified iteration.
     * @param windowIndex the marker window index.
     * @param iter the iteration
     */
    void printIterationUpdate(int window, int iter) {
        long buildTime = buildMillis - lastBuildMillis;
        lastBuildMillis = buildMillis;
        Utilities.duoPrint(log, Const.nl + "Window=" + window
                + " Iteration=" + iter + Const.nl);
        Utilities.duoPrint(log,
                    Utilities.printElapsedTime("Time for building model:         ", buildTime ));
        if (singleDagStats != null) {
            Utilities.duoPrint(log,
                    Utilities.printElapsedTime("Time for sampling (singles):     ", singleSampleMillis));
        }
        if (duoDagStats != null) {
            Utilities.duoPrint(log,
                    Utilities.printElapsedTime("Time for sampling (duos):        ", duoSampleMillis));
        }
        if (trioDagStats != null) {
            Utilities.duoPrint(log,
                    Utilities.printElapsedTime("Time for sampling (trios):       ", trioSampleMillis));
        }
        if (singleDagStats != null) {
            Utilities.duoPrint(log, "Singles model" + Const.nl);
            Utilities.duoPrint(log, singleDagStats);
        }
        if (duoDagStats != null) {
            Utilities.duoPrint(log, "Duos model" + Const.nl);
            Utilities.duoPrint(log, duoDagStats);
        }
        if (trioDagStats != null) {
            Utilities.duoPrint(log, "Trios model" + Const.nl);
            Utilities.duoPrint(log, trioDagStats);
        }
        log.flush();
    }
}
