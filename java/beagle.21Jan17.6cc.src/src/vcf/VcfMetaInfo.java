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
package vcf;

/**
 * <p>Class {@code VcfMetaInfo} represents a VCF meta-information line.
 * </p>
 * <p>Instances of class {@code VcfMetaInfo} are immutable.
 * </p>
 *
 * @author Brian L. Browning {@code <browning@uw.edu>}
 */
public final class VcfMetaInfo {

    /**
     * The VCF meta-information line prefix: "##"
     */
    public static final String PREFIX = "##";

    /**
     * The VCF meta-information line key-value delimiter: '='
     */
    public static final char DELIMITER = '=';

    private final String line;
    private final String key;
    private final String value;

    /**
     * Constructs a {@code VcfMetaInfo} instance representing
     * the specified VCF meta-information line.
     *
     * @param line a VCF meta-information line
     *
     * @throws IllegalArgumentException if the specified information line,
     * after removing any beginning and ending white-space, does not begin with
     * {@code VcfMetaInfo.PREFIX}, and does not contain non-empty key and
     * value strings separated by the  {@code VcfMetaInfo.DELIMITER} character
     *
     * @throws NullPointerException if {@code line == null}
     */
    public VcfMetaInfo(String line) {
        line = line.trim();
        if (line.startsWith(PREFIX)==false) {
            String s = "VCF meta-information line: missing starting \""
                    + PREFIX + "\": " + line;
            throw new IllegalArgumentException(s);
        }
        int index = line.indexOf(DELIMITER);
        if (index <=0 || index == line.length() - 1) {
            String s = "VCF meta-information line: missing \""
                    + DELIMITER + "\"";
            throw new IllegalArgumentException(s);
        }
        this.line = line;
        this.key = line.substring(2, index);
        this.value = line.substring(index+1);
    }

    /**
     * Returns the VCF meta-information line key.
     * @return the VCF meta-information line key
     */
    public String key() {
        return key;
    }

    /**
     * Returns the VCF meta-information line value.
     * @return the VCF meta-information line value
     */
    public String value() {
        return value;
    }

    /**
     * Returns the VCF meta-information line represented by {@code this}.
     *
     * @return the VCF meta-information line represented by {@code this}
     */
    @Override
    public String toString() {
        return line;
    }
}
