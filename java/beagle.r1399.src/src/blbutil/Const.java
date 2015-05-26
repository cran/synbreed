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
     * A constant holding {@code System.getProperty("line.separator")}.
     */
    public static final String nl = System.getProperty("line.separator");

    /**
     * A constant holding {@code ""}, the empty string.
     */
    public static final String EMPTY_STRING = "";

    /**
     * A string constant holding {@code "."}, the VCF missing-data symbol.
     */
    public static final String MISSING_DATA_STRING = ".";

    /**
     * A character constant holding {@code '.'}, the VCF missing-data
     * symbol.
     */
    public static final char MISSING_DATA_CHAR = '.';

    /**
     * A constant holding {@code ':'}.
     */
    public static final char colon = ':';

    /**
     * A constant holding {@code '-'}.
     */
    public static final char hyphen = '-';

    /**
     * A constant holding {@code '\t'}.
     */
    public static final char tab = '\t';

    /**
     * A constant holding {@code ';'}.
     */
    public static final char semicolon = ';';

    /**
     * A constant holding {@code ','}.
     */
    public static final char comma = ',';

    /**
     * A constant holding {@code '|'}, the phased allele separator.
     */
    public static final char phasedSep = '|';

    /**
     * A constant holding  {@code '/'}, the unphased allele separator.
     */
    public static final char unphasedSep = '/';

    /**
     * A constant holding 1,000,000,000.
     */
    public static final int giga = 1000000000;

    /**
     * A constant holding 1,000,000.
     */
    public static final int mega = 1000000;
}
