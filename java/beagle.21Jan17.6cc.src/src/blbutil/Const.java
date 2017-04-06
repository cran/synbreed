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

/**
 * Class {@code Const} provides public static final fields with
 * string and character constants.
 *
 * @author Brian L Browning
 */
public class Const {

    private Const() {
        // private constructor to prevent instantiation.
    }

    /**
     * The system-dependent string representing a new line-line:
     * {@code System.getProperty("line.separator")}
     */
    public static final String nl = System.getProperty("line.separator");

    /**
     * The VCF missing-data symbol as a string: {@code "."}
     */
    public static final String MISSING_DATA_STRING = ".";

    /**
     * The VCF missing-data symbol as a character: {@code '.'}
     */
    public static final char MISSING_DATA_CHAR = '.';

    /**
     * The colon character: {@code ':'}
     */
    public static final char colon = ':';

    /**
     * The hyphen character: {@code '-'}
     */
    public static final char hyphen = '-';

    /**
     * The tab character: {@code '\t'}
     */
    public static final char tab = '\t';

    /**
     * The semicolon character: {@code ';'}
     */
    public static final char semicolon = ';';

    /**
     * The comma character: {@code ','}
     */
    public static final char comma = ',';

    /**
     * The phased allele separator: {@code '|'}
     */
    public static final char phasedSep = '|';

    /**
     *  The unphased allele separator: {@code '/'}
     */
    public static final char unphasedSep = '/';

    /**
     * The value 1,000,000,000
     */
    public static final int giga = 1000000000;

    /**
     * The value 1,000,000
     */
    public static final int mega = 1000000;
}
