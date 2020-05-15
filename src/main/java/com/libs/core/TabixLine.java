package com.libs.core;

import htsjdk.tribble.Feature;

/**
 * Created by wei_qiang on 1/18/19.
 */
public class TabixLine implements Feature
{
    private String chrom;
    private int start;
    private int end;

    TabixLine(String line, int chromIndex, int startIndex, int endIndex)
    {

        String[] tokens = line.split("\t");
        try {
            this.chrom = tokens[chromIndex];
            this.start = Integer.parseInt(tokens[startIndex]);
            this.end = (tokens.length<3 ? start  : Integer.parseInt(tokens[endIndex]));
        }
        catch(final NumberFormatException ignored) {}
    }


    public String getContig() {
        return chrom;
    }

    public int getStart() {
        return start;
    }

    public int getEnd() {
        return end;
    }

}
