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
import blbutil.FileIt;
import blbutil.Filter;
import blbutil.InputIt;
import blbutil.IntPair;
import blbutil.SampleFileIt;
import blbutil.Utilities;
import haplotype.BasicSampleHapPairs;
import haplotype.BitHapPair;
import haplotype.HapPair;
import haplotype.SampleHapPairs;
import ibd.IbdSegment;
import java.io.File;
import java.util.ArrayList;
import java.util.HashSet;
import java.util.List;
import java.util.Locale;
import java.util.Map;
import java.util.Set;
import vcf.AllData;
import vcf.VcfIt;
import vcf.BrefIt;
import vcf.Data;
import vcf.IntervalVcfIt;
import vcf.Marker;
import vcf.FilterUtil;
import vcf.GL;
import vcf.Markers;
import vcf.TargetData;
import vcf.RefIt;
import vcf.VcfEmission;
import vcf.VcfRecord;

/**
 * Class {@code Main} is the entry class for the Beagle program.
 * See {@code Par.usage()} and online program documentation for usage
 * instructions.
 *
 * @author Brian L. Browning {@code <browning@uw.edu>}
 */
public class Main {

    /**
     * The program name and version.
     */
    public static final String program = "beagle.21Jan17.6cc.jar (version 4.1)";
    public static final String command = "java -jar beagle.21Jan17.6cc.jar";

    /**
     * The copyright string.
     */
    public static final String copyright = "Copyright (C) 2014-2015 Brian L. Browning";

    /**
     * The program name and a brief help message.
     */
    public static final String shortHelp = Main.program
            + Const.nl + Main.copyright
            + Const.nl + "Enter \"java -jar beagle.21Jan17.6cc.jar\" for a "
            + "summary of command line " + "arguments.";

    private final Par par;
    private final GeneticMap genMap;
    private final Data data;
    private final RunStats runStats;
    private final WindowWriter windowOut;

    /**
     * Entry point to Beagle program.  See {@code Parameters.usage()} and
     * online program documentation for usage instructions.
     *
     * @param args command line arguments
     */
    public static void main(String[] args) {
	Locale.setDefault(Locale.US);
        if (args.length==0) {
            System.out.println(program);
            System.out.println(copyright);
            System.out.println(Par.usage());
            System.exit(0);
        }
        Par par = parameters(args);
        System.setProperty("java.util.concurrent.ForkJoinPool.common.parallelism",
                String.valueOf(par.nthreads()));
        RunStats runStats = new RunStats(par);
        runStats.printStartInfo();
        GeneticMap genMap = geneticMap(par);

        try (Data data = (par.ref()==null) ? nonRefData(par) : allData(par);
                WindowWriter winOut = new WindowWriter(
                        data.targetSamples(), par.out())) {
            Main main = new Main(par, data, genMap, winOut, runStats);
            main.phaseData();
            runStats.printSummaryAndClose(data.nTargetMarkersSoFar(),
                    data.nMarkersSoFar());
        }
    }

    private Main(Par par, Data data, GeneticMap genMap,
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
     * Phases the data, imputes ungenotyped markers, and performed IBD segment
     * detection.
     */
    private void phaseData() {
        NuclearFamilies fam = new NuclearFamilies(data.targetSamples(), par.ped());
        runStats.printSampleSummary(fam, data);
        MainHelper mh = new MainHelper(par, genMap, runStats);
        SampleHapPairs overlapHaps = null;
        int overlap = 0;
        while (data.canAdvanceWindow()) {
            advanceWindow(overlap, par.window());
            CurrentData cd = new CurrentData(par, genMap, data, overlapHaps, fam);
            GenotypeValues gv = gv(par, cd);
            SampleHapPairs targetHapPairs = mh.phase(cd, gv);
            // targetHapPairs required to be aligned, GT-consistent with input data

            if (gv!=null) {
                windowOut.printGV(cd, gv);
            }
            else {
                Map<IntPair, List<IbdSegment>> ibd = mh.refinedIbd(cd, targetHapPairs);
                AlleleProbs alProbs = mh.LSImpute(cd, targetHapPairs);
                printOutput(cd, targetHapPairs, alProbs, ibd);
            }
            overlapHaps = overlapHaps(cd, targetHapPairs);
            overlap = cd.nMarkers() - cd.nextOverlapStart();
        }
    }

    private static GenotypeValues gv(Par par, CurrentData cd) {
        GenotypeValues gv = null;
        if (par.gt()==null) {
            gv = new BasicGenotypeValues(cd.targetMarkers(), cd.targetSamples());
            if (par.gtgl() != null) {
                initializeGV(gv, cd.targetGL());
            }
        }
        return gv;
    }

    /*
     * Initialize GenotypeValues to have values 1 at known genotypes.
     */
    private static void initializeGV(GenotypeValues gv, GL gl) {
        assert gv.markers().equals(gl.markers());
        assert gv.samples().equals(gl.samples());
        int nMarkers = gl.nMarkers();
        int nSamples = gl.nSamples();
        for (int m=0; m<nMarkers; ++m) {
            for (int s=0; s<nSamples; ++s) {
                int a1 = gl.allele1(m, s);
                int a2 = gl.allele2(m, s);
                if (a1>=0 && a2>=0) {
                    int gt = VcfRecord.gtIndex(a1, a2);
                    gv.add(m, s, gt, 1.0);
                }
            }
        }
    }

    private void printOutput(CurrentData cd, SampleHapPairs targetHapPairs,
            AlleleProbs alProbs, Map<IntPair, List<IbdSegment>> ibd) {
        assert par.gt()!=null;
        boolean markersAreImputed = cd.nTargetMarkers() < cd.nMarkers();
        boolean[] isImputed = isImputed(cd);
        int start = cd.prevSpliceStart();
        int end = cd.nextSpliceStart();
        boolean dose = markersAreImputed;
        boolean gprobs = markersAreImputed && par.gprobs();
        int nThreads = par.nthreads();
        if (markersAreImputed){
            alProbs = new ConstrainedAlleleProbs(targetHapPairs, alProbs,
                    cd.targetMarkerIndices());
        }
        windowOut.print(alProbs, isImputed, start, end, dose, gprobs, nThreads);
        if (par.ibd()) {
            windowOut.printIbd(cd, ibd);
        }
    }

    private static boolean[] isImputed(CurrentData cd) {
        boolean[] ba = new boolean[cd.nMarkers()];
        if (cd.nTargetMarkers()<ba.length) {
            for (int j=0; j<ba.length; ++j) {
                if (cd.targetMarkerIndex(j) == -1) {
                    ba[j] = true;
                }
            }
        }
        return ba;
    }

    private SampleHapPairs overlapHaps(CurrentData cd,
            SampleHapPairs targetHapPairs) {
        int nextOverlap = cd.nextTargetOverlapStart();
        int nextSplice = cd.nextTargetSpliceStart();
        if (cd.nextOverlapStart() == cd.nextSpliceStart()) {
            return null;
        }
        int nSamples = targetHapPairs.nSamples();
        int nMarkers = nextSplice - nextOverlap;
        Markers markers = targetHapPairs.markers().restrict(nextOverlap, nextSplice);
        Samples samples = targetHapPairs.samples();
        List<HapPair> list = new ArrayList<>(nSamples);
        int[] a1 = new int[nMarkers];
        int[] a2 = new int[nMarkers];
        for (int s = 0; s < nSamples; ++s) {
            for (int m = 0; m < nMarkers; ++m) {
                a1[m] = targetHapPairs.allele1(nextOverlap + m, s);
                a2[m] = targetHapPairs.allele2(nextOverlap + m, s);
            }
            list.add(new BitHapPair(markers, samples, s, a1, a2));
        }
        return new BasicSampleHapPairs(targetHapPairs.samples(), list);
    }

    private static Data nonRefData(Par par) {
        Filter<String> sampleFilter = FilterUtil.sampleFilter(
                par.excludesamples());
        Filter<Marker> markerFilter = FilterUtil.markerFilter(
                par.excludemarkers());
        ChromInterval chromInterval = ChromInterval.parse(par.chrom());
        SampleFileIt<? extends VcfEmission> targIt;
        // NB: originally used NuclearFamilies object in VcfEmission constructors.
        //     If this is still required, will need to read nonRefFile
        //     to get samples required to construct NuclearFamilies object.
        if (par.gt()!=null) {
            assert par.gl()==null && par.gtgl()==null;
            FileIt<String> it = InputIt.fromGzipFile(par.gt());
            targIt = VcfIt.create(it, sampleFilter, markerFilter, VcfIt.toBitSetGT);
        }
        else if (par.gl()!=null) {
            assert par.gt()==null && par.gtgl()==null;
            FileIt<String> it = InputIt.fromGzipFile(par.gl());
            targIt = VcfIt.create(it, sampleFilter, markerFilter, VcfIt.toGLRec);
        }
        else {
            assert par.gt()==null && par.gl()==null;
            FileIt<String> it = InputIt.fromGzipFile(par.gtgl());
            targIt = VcfIt.create(it, sampleFilter, markerFilter, VcfIt.toGTGLRec);
        }
        if (chromInterval!=null) {
            targIt = new IntervalVcfIt<>(targIt, chromInterval);
        }
        return TargetData.targetData(targIt);
    }

    private static Data allData(Par par) {
        Filter<String> sampleFilter = FilterUtil.sampleFilter(
                par.excludesamples());
        Filter<Marker> markerFilter = FilterUtil.markerFilter(
                par.excludemarkers());
        ChromInterval chromInterval = ChromInterval.parse(par.chrom());

        File targFile;
        SampleFileIt<? extends VcfEmission> targIt;
        SampleFileIt<VcfEmission> refIt;

        if (par.gt()!=null) {
            assert par.gl()==null && par.gtgl()==null;
            targFile = par.gt();
            FileIt<String> it = InputIt.fromGzipFile(targFile);
            targIt = VcfIt.create(it, sampleFilter, markerFilter, VcfIt.toBitSetGT);
        }
        else if (par.gl()!=null) {
            assert par.gt()==null && par.gtgl()==null;
            targFile = par.gl();
            FileIt<String> it = InputIt.fromGzipFile(targFile);
            targIt = VcfIt.create(it, sampleFilter, markerFilter, VcfIt.toGLRec);
        }
        else {
            assert par.gt()==null && par.gl()==null && par.gtgl()!=null;
            targFile = par.gtgl();
            FileIt<String> it = InputIt.fromGzipFile(targFile);
            targIt = VcfIt.create(it, sampleFilter, markerFilter, VcfIt.toGTGLRec);
        }

        if (par.impute()==false || par.gt()==null) {
            markerFilter = restrictToVcfMarkers(targFile, markerFilter,
                    chromInterval);
        }
        if (par.ref().toString().endsWith(".bref")) {
            refIt = new BrefIt(par.ref(), markerFilter);
        }
        else {
            FileIt<String> it = InputIt.fromGzipFile(par.ref());
            refIt = RefIt.create(it, sampleFilter, markerFilter,
                    RefIt.DEFAULT_EM_BUFFER_SIZE);
        }
        if (chromInterval!=null) {
            targIt = new IntervalVcfIt<>(targIt, chromInterval);
            refIt = new IntervalVcfIt<>(refIt, chromInterval);

        }
        return AllData.allData(refIt, targIt);
    }

    private static Filter<Marker> restrictToVcfMarkers(File vcfFile,
            Filter<Marker> markerFilter, ChromInterval chromInterval) {
        Set<Marker> includedMarkers = new HashSet<>(50000);
        try (FileIt<String> it = InputIt.fromGzipFile(vcfFile)) {
            Filter<String> sampleFilter = null;
            SampleFileIt<VcfRecord> vcfIt = VcfIt.create(it, sampleFilter,
                    markerFilter, VcfIt.toGTGLRec);
            if (chromInterval != null) {
                vcfIt = new IntervalVcfIt<>(vcfIt, chromInterval);
            }
            while (vcfIt.hasNext()) {
                includedMarkers.add(vcfIt.next().marker());
            }
            vcfIt.close();
        }
        return Filter.includeFilter(includedMarkers);
    }

    private static GeneticMap geneticMap(Par par) {
        if (par.map()==null) {
            return null;
        }
        else {
            String chrom = extractChrom(par.chrom());
            if (chrom==null) {
                return PlinkGeneticMap.fromPlinkMapFile(par.map());
            }
            else {
                return PlinkGeneticMap.fromPlinkMapFile(par.map(), chrom);
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
    private static Par parameters(String[] args) {
        // warnings are printed in RunStats.startInfo() method
        Par par = new Par(args);
        checkForOneInputFile(par);
        checkOutputPrefix(par);
        if (par.overlap() >= par.window()/2) {
            String s = shortHelp + Const.nl
                    + Const.nl + "ERROR: The \"window\" parameter must be at least "
                    + "two times the \"overlap\" parameter"
                    + Const.nl + "Exiting program.";
            Utilities.exit(s);
        }
        if (par.chrom()!=null && ChromInterval.parse(par.chrom())==null) {
            String s = shortHelp + Const.nl
                    + Const.nl + "ERROR: invalid \"chrom\" parameter: \""
                    + par.chrom() + "\""
                    + Const.nl + "Exiting program.";
            Utilities.exit(s);
        }
        if (par.ibd()==true && par.ref()!=null && par.impute()) {
            String s = shortHelp + Const.nl
                    + Const.nl + "ERROR: The \"impute=false\" argument is "
                    + "required when a reference panel"
                    + Const.nl + "       is specified and \"ibd=true\"."
                    + Const.nl + "Exiting program.";
            Utilities.exit(s);
        }
        return par;
    }

    private static void checkForOneInputFile(Par par) {
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
            Utilities.exit(Par.usage() + s);
        }
    }

    private static void checkOutputPrefix(Par par) {
        File outPrefix = new File(par.out());
        if (outPrefix.isDirectory()) {
            String s = "ERROR: \"out\" parameter cannot be a directory: \""
                    + par.out() + "\"";
            Utilities.exit(Par.usage() + s);
        }

        File vcfOut = new File(par.out() + ".vcf.gz");
        if (vcfOut.equals(par.ref())) {
            String s = "ERROR: VCF output file equals input file: " + par.ref();
            Utilities.exit(Par.usage() + s);
        }
        if (vcfOut.equals(par.gt())) {
            String s = "ERROR: VCF output file equals input file: " + par.gt();
            Utilities.exit(Par.usage() + s);
        }
        if (vcfOut.equals(par.gl())) {
            String s = "ERROR: VCF output file equals input file: " + par.gl();
            Utilities.exit(Par.usage() + s);
        }
        if (vcfOut.equals(par.gtgl())) {
            String s = "ERROR: VCF output file equals input file: " + par.gtgl();
            Utilities.exit(Par.usage() + s);
        }
    }

    private void advanceWindow(int overlap, int window) {
        if (data.canAdvanceWindow()) {
            data.advanceWindow(overlap, window);
        }
        runStats.printWindowUpdate(data);
    }
}
