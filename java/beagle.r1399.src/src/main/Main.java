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

import beagleutil.ChromInterval;
import beagleutil.Samples;
import blbutil.Const;
import blbutil.Filter;
import blbutil.FilterUtils;
import blbutil.IntPair;
import blbutil.Utilities;
import haplotype.HapPair;
import haplotype.SampleHapPairs;
import haplotype.SampleHapPairsSplicer;
import haplotype.Weights;
import ibd.IbdSegment;
import java.io.File;
import java.util.List;
import java.util.Locale;
import java.util.Map;
import java.util.Random;
import java.util.Set;
import vcf.AllData;
import vcf.Data;
import vcf.GL;
import vcf.Marker;
import vcf.MarkerFilterUtils;
import vcf.Markers;
import vcf.NonRefData;
import vcf.VcfRecord;

/**
 * See {@code Parameters.usage()} for usage instructions.
 *
 * @author Brian L. Browning {@code <browning@uw.edu>}
 */
public class Main {

    /**
     * The program name and version.
     */
    public static final String version = "beagle.jar (r1399)";

    /**
     * The copyright string.
     */
    public static final String copyright = "Copyright (C) 2014 Brian L. Browning";

    private final Parameters par;
    private final GeneticMap genMap;
    private final Data data;
    private final RunStats runStats;
    private final WindowWriter windowOut;

    /**
     * Entry point to Beagle program.  See .pdf documentation for
     * information on program arguments.
     *
     * @param args command line arguments
     */
    public static void main(String[] args) {
	Locale.setDefault(Locale.US);
        if (args.length==0) {
            System.out.println(version);
            System.out.println(copyright);
            System.out.println(Parameters.usage());
            System.exit(0);
        }

        Parameters par = parameters(args);
        setLoggingTarget(par.out());
        RunStats runStats = new RunStats(par);
        runStats.printStartInfo();
        Data data = (par.ref()==null) ? nonRefData(par) : allData(par);
        GeneticMap genMap = geneticMap(par);
        WindowWriter windowOut = new WindowWriter(data.nonRefSamples(), par.out());

        Main main = new Main(par, data, genMap, windowOut, runStats);
        main.phaseData();

        data.close();
        windowOut.close();
        closeLogger(par.out(), runStats);
        runStats.printSummaryAndClose(data.cumMarkerCnt());
    }

    private Main(Parameters par, Data data, GeneticMap genMap,
            WindowWriter windowWriter, RunStats runStats) {
        assert par!=null;
        assert data!=null;
        assert windowWriter!=null;
        assert runStats!=null;
        this.par = par;
        this.genMap = genMap;
        this.data = data;
        this.runStats = runStats;
        this.windowOut = windowWriter;
    }

    /*
     * Determines parameters for phasing the data from command line
     * arguments and phases the data.
     */
    private void phaseData() {
        NuclearFamilies fam = new NuclearFamilies(data.nonRefSamples(), par.ped());
        Weights weights = new Weights(fam);
        runStats.printSampleSummary(fam, data);
        SampleHapPairs prevNonRefHaps = null;
        Random random = new Random(par.seed());
        MainHelper mh = new MainHelper(par, genMap, fam, weights, runStats,
                random);
        while (data.canAdvanceWindow()) {
            advanceWindow();
            int lastSplice = data.overlap()/2;
            int nextOverlap = nextOverlap(data, par.overlap());
            int nextSplice = nextSplice(data, par.overlap());

            Markers markers = data.markers();
            List<HapPair> restrictedRefHaps = data.restrictedRefHaps();
            GL gl = data.nonRefEmissions();
            GenotypeValues gv = gv(markers, gl.samples());
            GenotypeValues restrictedGV = restrictAndInitializeGV(gv, gl);

            SampleHapPairs nonRefHaps = mh.sample(data, restrictedGV);

            if (prevNonRefHaps != null) {
                int nextStart = targetIndex(lastSplice);
                nonRefHaps = SampleHapPairsSplicer.spliceNext(prevNonRefHaps,
                        nonRefHaps, nextStart, data.nonRefOverlap());
            }
            Map<IntPair, List<IbdSegment>> ibd = mh.refinedIbd(restrictedRefHaps,
                    gl, nonRefHaps, weights);
            prevNonRefHaps = nonRefHaps;
            SampleHapPairs impHaps = mh.impute(data, nonRefHaps, gv);
            windowOut.print(impHaps, gv, ibd, lastSplice, nextOverlap,
                    nextSplice);
        }
    }

    /* returns the first index in the next overlap */
    static int nextOverlap(Data data, int overlap) {
        if (data.canAdvanceWindow() && data.lastWindowOnChrom()==false) {
            return data.nMarkers() - overlap;
        }
        else {
            return data.nMarkers();
        }
    }

    /* returns the first index after the next splice point */
    static int nextSplice(Data data, int overlap) {
        if (data.canAdvanceWindow() && data.lastWindowOnChrom()==false) {
            return data.nMarkers() - overlap + (overlap/2);
        }
        else {
            return data.nMarkers();
        }
    }

    private GenotypeValues gv(Markers markers, Samples samples) {
        GenotypeValues gv = null;
        if (par.gprobs()) {
            gv = new BasicGenotypeValues(markers, samples);
        }
        return gv;
    }

    /* first target index on or after specified ref index */
    private int targetIndex(int refIndex) {
        int i=0;
        while (i<data.nNonRefMarkers() && data.markerIndex(i)<refIndex) {
            ++i;
        }
        return i;
    }

    /*
     * Insures known genotypes have correct GL values
     */
    private GenotypeValues restrictAndInitializeGV(
            GenotypeValues gv, GL gl) {
        if (par.gprobs()==false) {
            return gv;
        }
        GenotypeValues restrictGV
                =  new RestrictedGenotypeValues(gv, gl.markers());
        if (par.gprobs()==true) {
            int nMarkers = gl.nMarkers();
            int nSamples = gl.nSamples();
            for (int m=0; m<nMarkers; ++m) {
                for (int s=0; s<nSamples; ++s) {
                    byte a1 = gl.allele1(m, s);
                    byte a2 = gl.allele2(m, s);
                    if (a1>=0 && a2>=0) {
                        int gt = VcfRecord.gtIndex(a1, a2);
                        restrictGV.add(m, s, gt, 1.0);
                    }
                }
            }
        }
        return restrictGV;
    }

    private static void setLoggingTarget(String prefix) {
        File logFile = new File(prefix + ".warnings");
        Logger.getInstance().setTarget(logFile);
    }

    private static void closeLogger(String prefix, RunStats runStats) {
        if (Logger.getInstance().logCnt() > 0) {
            String s = "Warnings printed: see " + prefix + ".warnings file";
            runStats.println(s);
        }
        Logger.getInstance().close();
    }

    private static Filter<Marker> markerFilter(Parameters par) {
        Filter<Marker> markerFilter = null;
        if (par.excludemarkers()!=null) {
            Set<String> excludedIds = Utilities.idSet(par.excludemarkers());
            markerFilter = MarkerFilterUtils.excludeIdFilter(excludedIds);
        }
        return markerFilter;
    }

    private static Filter<String> sampleFilter(Parameters par) {
        Filter<String> sampleFilter = null;
        if (par.excludesamples()!=null) {
            Set<String> exclude = Utilities.idSet(par.excludesamples());
            sampleFilter = FilterUtils.excludeFilter(exclude);
        }
        return sampleFilter;
    }

    private static ChromInterval chromInterval(Parameters par) {
        ChromInterval chromInterval = null;
        if (par.chrom()!=null) {
            chromInterval = ChromInterval.parse(par.chrom());
        }
        return chromInterval;
    }

    private static Data nonRefData(Parameters par) {
        Filter<String> sampleFilter = sampleFilter(par);
        Filter<Marker> markerFilter = markerFilter(par);
        ChromInterval chromInterval = chromInterval(par);
        if (par.gt()!=null) {
            assert par.gl()==null && par.gtgl()==null;
            return NonRefData.gt(par.gt(), sampleFilter, markerFilter,
                    chromInterval, par.ped(), par.usephase());
        }
        else if (par.gl()!=null) {
            assert par.gt()==null && par.gtgl()==null;
            return NonRefData.gl(par.gl(), sampleFilter, markerFilter,
                    chromInterval, par.ped(), par.maxlr());
        }
        else {
            assert par.gt()==null && par.gl()==null;
            return NonRefData.gtgl(par.gtgl(), par.ped(), par.usephase(),
                    par.maxlr(), sampleFilter, markerFilter, chromInterval);
        }
    }

    private static Data allData(Parameters par) {
        Filter<String> sampleFilter = sampleFilter(par);
        Filter<Marker> markerFilter = markerFilter(par);
        ChromInterval chromInterval = chromInterval(par);
        if (par.gt()!=null) {
            assert par.gl()==null && par.gtgl()==null;
            return AllData.gt(par.ref(), par.gt(), sampleFilter, markerFilter,
                    chromInterval, par.ped(), par.usephase(), par.impute());
        }
        else if (par.gl()!=null) {
            assert par.gt()==null && par.gtgl()==null;
            return AllData.gl(par.ref(), par.gl(), sampleFilter,
                    markerFilter, chromInterval, par.ped(), par.maxlr(),
                    par.impute());
        }
        else {
            assert par.gt()==null && par.gl()==null && par.gtgl()!=null;
            return AllData.gtgl(par.ref(), par.gtgl(), sampleFilter,
                    markerFilter, chromInterval, par.ped(), par.usephase(),
                    par.maxlr(), par.impute());
        }
    }

    private static GeneticMap geneticMap(Parameters par) {
        if (par.map()==null) {
            return null;
        }
        else {
            String chrom = extractChrom(par.chrom());
            if (chrom==null) {
                return GeneticMap.fromPlinkMapFile(par.map());
            }
            else {
                return GeneticMap.fromPlinkMapFile(par.map(), chrom);
            }
        }
    }

    private static String extractChrom(String chromParameter) {
        if (chromParameter!=null && chromParameter.length()>0) {
            int index = chromParameter.indexOf(Const.colon);
            if (index == -1) {
                return chromParameter;
            }
            else {
                return new String(chromParameter.substring(0, index));
            }
        }
        else {
            return null;
        }
    }

    /*
     * Checks that certain parameters are consistent, and prints error
     * message and exits if parameters are inconsistent.
     *
     * @param args the command line arguments.
     */
    private static Parameters parameters(String[] args) {
        Parameters par = new Parameters(args);
        checkForOneInputFile(par);
        checkOutputPrefix(par);
        if (par.overlap() >= par.window()) {
            String s = "ERROR: \"overlap\" parameter (" + par.overlap() +
                    " must less than \"windowsize\" parameter ("
                    + par.window() + ")";
            Utilities.exit(Parameters.usage() + s);
        }
        if (par.chrom()!=null && ChromInterval.parse(par.chrom())==null) {
            String s = "ERROR: invalid \"chrom\" parameter: \""
                    + par.chrom() + "\"";
            Utilities.exit(Parameters.usage() + s);
        }
        if (par.ibd()==true && par.ref()!=null && par.impute()==true) {
            String s = "ERROR: The \"impute=false\" parameter is required "
                    + "when a reference panel is specified and \"ibd=true\"";
            Utilities.exit(Parameters.usage() + s);
        }
        if (par.burnin_its()==0 && par.phase_its()==0 && par.ped()!=null) {
            String s = "ERROR: The \"ped\" parameter cannot be used when there"
                    + "are no burnin or phasing iterations ("
                    + "\"burnin-its=0 phase-its=0\")";
            Utilities.exit(Parameters.usage() + s);
        }
        return par;
    }

    private static void checkForOneInputFile(Parameters par) {
        int cnt = 0;
        if (par.gt()!=null) {
            ++cnt;
        }
        if (par.gl()!=null) {
            ++cnt;
        }
        if (par.gtgl()!=null) {
            ++cnt;
        }
        if (cnt != 1) {
            String s = "ERROR: exactly one \"gt\", \"gl\" or \"gtgl\" "
                    + "parameter is required.";
            Utilities.exit(Parameters.usage() + s);
        }
    }

    private static void checkOutputPrefix(Parameters par) {
        File outPrefix = new File(par.out());
        if (outPrefix.isDirectory()) {
            String s = "ERROR: \"out\" parameter cannot be a directory: \""
                    + par.out() + "\"";
            Utilities.exit(Parameters.usage() + s);
        }

        File vcfOut = new File(par.out() + ".vcf.gz");
        if (vcfOut.equals(par.ref())) {
            String s = "ERROR: VCF output file equals input file: " + par.ref();
            Utilities.exit(Parameters.usage() + s);
        }
        if (vcfOut.equals(par.gt())) {
            String s = "ERROR: VCF output file equals input file: " + par.gt();
            Utilities.exit(Parameters.usage() + s);
        }
        if (vcfOut.equals(par.gl())) {
            String s = "ERROR: VCF output file equals input file: " + par.gl();
            Utilities.exit(Parameters.usage() + s);
        }
        if (vcfOut.equals(par.gtgl())) {
            String s = "ERROR: VCF output file equals input file: " + par.gtgl();
            Utilities.exit(Parameters.usage() + s);
        }
    }

    private void advanceWindow() {
        if (data.canAdvanceWindow()) {
            data.advanceWindow(par.overlap(), par.window());
        }
        runStats.printWindowUpdate(data);
    }
}
