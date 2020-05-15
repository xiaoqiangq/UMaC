package com.libs.core;

import com.libs.file.DeClassLoader;
import htsjdk.tribble.readers.TabixReader;

import java.io.BufferedReader;
import java.io.BufferedWriter;
import java.io.ByteArrayInputStream;
import java.io.File;
import java.io.FileWriter;
import java.io.IOException;
import java.io.InputStreamReader;
import java.util.Properties;

/**
 * Created by wei_qiang on 8/23/19.
 */
public class Segment
{
    private final String input;
    private final String output;
    private final File windows;
    private final int cutoff = 2;

    public Segment(Properties properties)
    {
        this.input = properties.getProperty("input");
        this.output = properties.getProperty("output");

        this.windows = new File(properties.getProperty("window"));
    }

    public void run() throws IOException {

        FileWriter fw = new FileWriter(output);
        BufferedWriter bw = new BufferedWriter(fw);

        DeClassLoader loaderWindows = new DeClassLoader();
        ByteArrayInputStream windowsByte = new ByteArrayInputStream(loaderWindows.file2Bytes(windows));

        BufferedReader bfReader = new BufferedReader(new InputStreamReader(windowsByte));

        String window;

        String[] str_split;

        TabixReader tabixReader= new TabixReader(input);

        while((window = bfReader.readLine()) != null)
        {
            if (window.startsWith("#") || window.startsWith("fixedStep") )
            {
                continue;
            }
            str_split = window.split("\t");

            String chrom = str_split[0];
            int start = Integer.parseInt(str_split[1]);
            int end = Integer.parseInt(str_split[2]);

            TabixReader.Iterator iterator = tabixReader.query(chrom, start, end);

            String line;
            String[] segment;
            double tumorReads = 0.0;
            double normalReads = 0.0;

            while ((line = iterator.next()) != null)
            {
                if (line.startsWith("#"))
                {
                    continue;
                }

                segment = line.split("\t");

                if (Integer.parseInt(segment[10]) >= cutoff)
                {
                    double tumortmp = Double.parseDouble(segment[5]);
                    double normaltmp = Double.parseDouble(segment[6]);

                    tumorReads += tumortmp;
                    normalReads += normaltmp;
                }
            }

            if(tumorReads+normalReads == 0)
            {
                bw.write(chrom+"\t"+start+"\t"+end+"\t"+str_split[3]+"\t"+0+"\t"+0+"\t"+0+"\n");
            }else
            {
                int d1 = (int) tumorReads;
                int d2 = (int) normalReads;
                double ratio = tumorReads/(tumorReads+normalReads);

                bw.write(chrom+"\t"+start+"\t"+end+"\t"+str_split[3]+"\t"+d1+"\t"+d2+"\t"+ratio+"\n");
            }
        }

        bw.close();
        bfReader.close();
        windowsByte.close();
    }
}
