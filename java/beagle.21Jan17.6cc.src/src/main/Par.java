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

import blbutil.Const;
import blbutil.Validate;
import java.io.File;
import java.util.Map;

/**
 * <p>Class {@code Parameters} represents the parameters for a Beagle analysis.
 * </p>
 * <p>Instances of class {@code Parameters} are immutable.
 * </p>
 *
 * @author Brian L. Browning {@code <browning@uw.edu>}
 */
public final class Par {

    private final String[] args;

    // data input/output parameters
    private final File gt;
    private final File gl;
    private final File gtgl;
    private final File ref;
    private final File dag;
    private final String out;
    private final File excludesamples;
    private final File excludemarkers;
    private final File ped;
    private final File map;
    private final String chrom;
    private final float maxlr;

    // algorithm parameters
    private final int nthreads;
    private final boolean lowmem;
    private final int window;
    private final int overlap;
    private final boolean impute;
    private final boolean gprobs;
    private final int niterations;
    private final float mapscale;
    private final float ne;
    private final float err;
    private final float cluster;
    private final long seed;

    // ibd parameters
    private final boolean ibd;
    private final float ibdlod;
    private final float ibdcm;
    private final float ibdscale;
    private final int ibdtrim;

    // expert parameters
    private final float modelscale;

    // undocumented parameters
    private final int burnin_its;
    private final int phase_its;
    private final int nsamples;
    private final float ibdlength;
    private final float ibdextend;

    /**
     * Constructs a new {@code Parameters} instance from the specified
     * command line arguments.
     * @param args the Beagle command line arguments
     * @throws IllegalArgumentException if a command line argument
     * is incorrectly specified
     * @throws NumberFormatException if a numeric value for a parameter
     * is incorrectly specified
     */
    public Par(String[] args) {

        int IMAX = Integer.MAX_VALUE;
        long LMIN = Long.MIN_VALUE;
        long LMAX = Long.MAX_VALUE;
        float FMIN = Float.MIN_VALUE;
        float FMAX = Float.MAX_VALUE;

        this.args = args.clone();
        Map<String, String> argsMap = Validate.argsToMap(args, '=');

        // data input/output parameters
        gt = Validate.getFile(
                Validate.stringArg("gt", argsMap, false, null, null));
        gl = Validate.getFile(
                Validate.stringArg("gl", argsMap, false, null, null));
        gtgl = Validate.getFile(
                Validate.stringArg("gtgl", argsMap, false, null, null));
        ref = Validate.getFile(
                Validate.stringArg("ref", argsMap, false, null, null));
        dag = Validate.getFile(
                Validate.stringArg("dag", argsMap, false, null, null));
        out = Validate.stringArg("out", argsMap, true, null, null);
        excludesamples = Validate.getFile(
                Validate.stringArg("excludesamples", argsMap, false, null, null));
        excludemarkers = Validate.getFile(
                Validate.stringArg("excludemarkers", argsMap, false, null, null));
        ped = Validate.getFile(
                Validate.stringArg("ped", argsMap, false, null, null));
        map = Validate.getFile(Validate.stringArg("map", argsMap, false, null, null));
        chrom = Validate.stringArg("chrom", argsMap, false, null, null);
        maxlr = Validate.floatArg("maxlr", argsMap, false, 5000.0f, 1.1f, FMAX);

        // algorithm parameters
        nthreads = modNthreads(Validate.intArg("nthreads", argsMap, false, IMAX, 0, IMAX));
        lowmem = Validate.booleanArg("lowmem", argsMap, false, true);
        window = Validate.intArg("window", argsMap, false, 50000, 1, IMAX);
        overlap = Validate.intArg("overlap", argsMap, false, 3000, 0, IMAX);
        niterations = Validate.intArg("niterations", argsMap, false, 5, 0, IMAX);
        impute = Validate.booleanArg("impute", argsMap, false, true);
        gprobs = Validate.booleanArg("gprobs", argsMap, false, false);
        ne = Validate.floatArg("ne", argsMap, false, 1_000_000f, FMIN, FMAX);
        err = Validate.floatArg("err", argsMap, false, 0.0001f, 0.0f, FMAX);
        cluster = Validate.floatArg("cluster", argsMap, false, 0.005f, 0.0f, FMAX);
        seed = Validate.longArg("seed", argsMap, false, -99999, LMIN, LMAX);

        // ibd parameters
        ibd = Validate.booleanArg("ibd", argsMap, false, false);
        ibdlod = Validate.floatArg("ibdlod", argsMap, false, 3.0f, FMIN, FMAX);
        ibdcm = Validate.floatArg("ibdcm", argsMap, false, 1.5f, FMIN, FMAX);
        ibdscale = Validate.floatArg("ibdscale", argsMap, false, 0.0f, 0.0f, FMAX);
        ibdtrim = Validate.intArg("ibdtrim", argsMap, false, 40, 0, IMAX);

        // expert parameters
        modelscale = Validate.floatArg("modelscale", argsMap, false, 0.8f, FMIN, FMAX);

        // undocumented parameters
        burnin_its = 5;
        phase_its = 5;
        nsamples = 4;
        mapscale = Validate.floatArg("mapscale", argsMap, false, 1.0f, FMIN, FMAX);
        ibdlength = Validate.floatArg("ibdlength", argsMap, false, 0.07f, FMIN, FMAX);
        ibdextend = Validate.floatArg("ibdextend", argsMap, false, 0.13f, 0.0f, FMAX);
        Validate.confirmEmptyMap(argsMap);
    }

    /**
     * Returns the Beagle command line arguments.
     * @return the Beagle command line arguments
     */
    public String[] args() {
        return args.clone();
    }

    /**
     * Returns a description of the possible Beagle command line arguments.
     * @return a description of the possible Beagle command line arguments
     */
    public static String usage() {
        String nl = Const.nl;
        return  "Command line syntax: " + Main.command + " [arguments]" + nl
                + nl
                + "data input/output parameters ..." + nl
                + "  gt=<VCF file: use GT field>                        (optional)" + nl
                + "  gl=<VCF file: use GL/PL field>                     (optional)" + nl
                + "  gtgl=<VCF file: use GT (preferred) or GL/PL field> (optional)" + nl
                + "  ref=<VCF file with phased genotypes>               (optional)" + nl
                + "  out=<output file prefix>                           (required)" + nl
                + "  excludesamples=<file with 1 sample ID per line>    (optional)" + nl
                + "  excludemarkers=<file with 1 marker ID per line>    (optional)" + nl
//                + "  ped=<linkage format pedigree file>                 (optional)" + nl
                + "  map=<PLINK map file with cM units>                 (optional)" + nl
                + "  chrom=<[chrom] or [chrom]:[start]-[end]>           (optional)" + nl
                + "  maxlr=<max GL/PL likelihood ratio>                 (default=5000)" + nl + nl

                + "general parameters ..." + nl
                + "  nthreads=<number of threads>                       (default: machine-dependent)" + nl
                + "  lowmem=<use low-memory algorithm (true/false)>     (default=false)" + nl
                + "  window=<markers per window>                        (default=50000)" + nl
                + "  overlap=<overlap between windows>                  (default=3000)" + nl
                + "  seed=<random seed>                                 (default=-99999)" + nl + nl

                + "phasing and imputation parameters ..." + nl
                + "  niterations=<number of phasing iterations>         (default=5)" + nl
                + "  impute=<impute ungenotyped markers (true/false)>   (default=true)" + nl
                + "  gprobs=<print GP field for imputed markers>        (default=false)" + nl
                + "  ne=<effective population size>                     (default=1000000)" + nl
                + "  err=<allele miscall rate>                          (default=0.0001)" + nl
                + "  cluster=<max cM in a marker cluster>               (default=0.005)" + nl + nl

                + "IBD parameters ..." + nl
                + "  ibd=<perform IBD detection (true/false)>           (default=false)" + nl
                + "  ibdlod=<min LOD score of reported IBD segments>    (default=3.0)" + nl
                + "  ibdcm=<min cM length of reported IBD segments>     (default=1.5)" + nl
                + "  ibdscale=<model scale factor for Refined IBD>      (default: data-dependent)" + nl
                + "  ibdtrim=<markers at each segment end>              (default=40)" + nl;
    }

    /**
     * Returns a sample-size-adjusted IBD scale parameter equal to
     * {@code Math.max(2.0f, (float) Math.sqrt(nSamples/100.0))} if
     * {@code this.ibdscale() == 0f}, and returns
     * {@code this.ibdscale()} otherwise.
     *
     * @param nSamples the number of samples
     * @return a sample-size-adjusted IBD scale parameter if
     * {@code this.ibdscale() == 0f}, and {@code this.ibdscale()} otherwise
     * @throws IllegalArgumentException if {@code nSamples < 0}
     */
    public float adjustedIbdScale(int nSamples) {
        if (nSamples <= 0) {
            throw new IllegalArgumentException(String.valueOf(nSamples));
        }
        if (ibdscale==0) {
            return Math.max(2.0f, (float) Math.sqrt(nSamples/100.0));
        }
        else {
            return ibdscale;
        }
    }

    /**
     * Returns the nthreads parameter, which is equal to
     * {@code Runtime.getRuntime().availableProcessors()} if
     * {@code nthreads == Integer.MAX_VALUE}.
     * @return the nthreads parameter
     */
    private static int modNthreads(int nthreads) {
        if (nthreads==Integer.MAX_VALUE) {
            return Runtime.getRuntime().availableProcessors();
        }
        else {
            return nthreads;
        }
    }

    // data input/output parameters

    /**
     * Returns the gt parameter or {@code null} if no gt parameter was
     * specified.
     * @return the gt parameter or {@code null} if no gt parameter was
     * specified
     */
    public File gt() {
        return gt;
    }

    /**
     * Returns the gl parameter or {@code null} if no gl parameter was
     * specified.
     * @return the gl parameter or {@code null} if no gl parameter was
     * specified
     */
    public File gl() {
        return gl;
    }

    /**
     * Returns the gtgl parameter or {@code null} if no gtgl parameter was
     * specified.
     * @return the gtgl parameter or {@code null} if no gtgl parameter was
     * specified.
     */
    public File gtgl() {
        return gtgl;
    }

    /**
     * Returns the ref parameter or {@code null} if no ref parameter was
     * specified.
     * @return the ref parameter or {@code null} if no ref parameter was
     * specified
     */
    public File ref() {
        return ref;
    }

    /**
     * Returns the dag parameter or {@code null} if no ref parameter was
     * specified.
     * @return the dag parameter or {@code null} if no ref parameter was
     * specified
     */
    public File dag() {
        return dag;
    }

    /**
     * Returns the out parameter.
     * @return the out parameter
     */
    public String out() {
        return out;
    }

    /**
     * Returns the excludesamples parameter or {@code null}
     * if no excludesamples parameter was specified.
     *
     * @return the excludesamples parameter or {@code null}
     * if no excludesamples parameter was specified
     */
    public File excludesamples() {
        return excludesamples;
    }

    /**
     * Returns the excludemarkers parameter or {@code null}
     * if no excludemarkers parameter was specified.
     *
     * @return the excludemarkers parameter or {@code null}
     * if no excludemarkers parameter was specified
     */
    public File excludemarkers() {
        return excludemarkers;
    }

    /**
     * Returns the ped parameter or {@code null}
     * if no ped parameter was specified.
     *
     * @return the ped parameter or {@code null}
     * if no ped parameter was specified
     */
    public File ped() {
        return ped;
    }

    /**
     * Returns the map parameter.
     * @return the map parameter
     */
    public File map() {
        return map;
    }

    /**
     * Returns the chrom parameter or {@code null}
     * if no chrom parameter was specified.
     *
     * @return the chrom parameter or {@code null}
     * if no chrom parameter was specified
     */
    public String chrom() {
        return chrom;
    }

    /**
     * Returns the maxlr parameter.
     * @return the maxlr parameter
     */
    public float maxlr() {
        return maxlr;
    }

    // general parameters

    /**
     * Returns the nthreads parameter.
     * @return the nthreads parameter
     */
    public int nthreads() {
        return nthreads;
    }

    /**
     * Returns the lowmem parameter.
     * @return the lowmem parameter
     */
    public boolean lowmem() {
        return lowmem;
    }

    /**
     * Returns the window parameter.
     * @return the window parameter
     */
    public int window() {
        return window;
    }

    /**
     * Return the overlap parameter.
     * @return the overlap parameter.
     */
    public int overlap() {
        return overlap;
    }

    /**
     * Returns the seed parameter.
     * @return the seed parameter
     */
    public long seed() {
        return seed;
    }

    // phasing and imputation parameters

    /**
     * Returns the niterations parameter.
     * @return the niterations parameter
     */
    public int niterations() {
        return niterations;
    }

    /**
     * Returns the impute parameter.
     * @return the impute parameter
     */
    public boolean impute() {
        return impute;
    }

    /**
     * Returns the gprobs parameter.
     * @return the gprobs parameter
     */
    public boolean gprobs() {
        return gprobs;
    }

    /**
     * Returns the ne parameter
     * @return the ne parameter
     */
    public float ne() {
        return ne;
    }

    /**
     * Returns the err parameter.
     * @return the err parameter
     */
    public float err() {
        return err;
    }

    /**
     * Returns the cluster parameter.
     * @return the cluster parameter
     */
    public float cluster() {
        return cluster;
    }

    // ibd parameters

    /**
     * Returns the ibd parameter.
     * @return the ibd parameter
     */
    public boolean ibd() {
        return ibd;
    }

    /**
     * Returns the ibdlod parameter.
     * @return the ibdlod parameter
     */
    public float ibdlod() {
        return ibdlod;
    }

    /**
     * Returns the ibdcm parameter.
     * @return the ibdcm parameter
     */
    public float ibdcm() {
        return ibdcm;
    }

    /**
     * Returns the ibdscale parameter.
     * @return the ibdscale parameter
     */
    public float ibdscale() {
        return ibdscale;
    }

    /**
     * Returns the ibdtrim parameter.
     * @return the ibdtrim parameter
     */
    public int ibdtrim() {
        return ibdtrim;
    }

    // expert parameters

    /**
     * Returns the modelscale parameter.
     * @return the modelscale parameter
     */
    public float modelscale() {
        return modelscale;
    }

    // undocumented parameters

    /**
     * Returns the burnin-its parameter.
     * @return the burnin-its parameter
     */
    public int burnin_its() {
        return burnin_its;
    }

    /**
     * Returns the phase-its parameter.
     * @return the phase-its parameter
     */
    public int phase_its() {
        return phase_its;
    }

    /**
     * Return the nsamples parameter.
     * @return the nsamples parameter
     */
    public int nsamples() {
        return nsamples;
    }

    /**
     * Returns the mapscale parameter.
     * @return the mapscale parameter
     */
    public float mapscale() {
        return mapscale;
    }

    /**
     * Returns the ibdlength parameter.
     * @return the ibdlength parameter
     */
    public float ibdlength() {
        return ibdlength;
    }

    /**
     * Returns the ibdextend parameter.
     * @return the ibdextend parameter
     */
    public float ibdextend() {
        return ibdextend;
    }
}
