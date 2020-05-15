package com.libs.core;

import htsjdk.samtools.util.BlockCompressedInputStream;
import htsjdk.tribble.index.Index;
import htsjdk.tribble.index.tabix.TabixFormat;
import htsjdk.tribble.index.tabix.TabixIndexCreator;

import java.io.File;
import java.io.FileInputStream;
import java.io.IOException;

/**
 * Created by wei_qiang on 1/19/19.
 */

public class MakeTabixCompressedIndex
{
    private final String input;
    private int sequence = 2;
    private int start = 14;
    private int end = 4;

    MakeTabixCompressedIndex(String input)
    {
        this.input = input;
    }

    public MakeTabixCompressedIndex(String input, int sequence, int start, int end)
    {
        this(input);
        this.sequence = sequence;
        this.start = start;
        this.end = end;
    }

    void run() throws IOException
    {
        TabixFormat tabixFormat = new TabixFormat(0,sequence,start,end,'#',0);
        TabixIndexCreator indexCreator = new TabixIndexCreator(tabixFormat);

        BlockCompressedInputStream inputStream = new BlockCompressedInputStream(new FileInputStream(input));

        long filePosition = inputStream.getFilePointer();
        String line;

        while ((line = inputStream.readLine())!= null)
        {
            if (line.startsWith("##"))
            {
                continue;
            }

            TabixLine bed = new TabixLine(line, sequence-1,start-1, end-1);

            if (!line.startsWith("#"))
            {
                indexCreator.addFeature(bed, filePosition);
            }
            filePosition = inputStream.getFilePointer();
        }

        Index index = indexCreator.finalizeIndex(filePosition);
        index.writeBasedOnFeatureFile(new File(input));

        inputStream.close();
    }
}
