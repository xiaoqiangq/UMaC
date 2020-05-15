package com.libs.file;

import htsjdk.tribble.readers.TabixReader;

import java.io.File;
import java.io.IOException;
import java.util.HashMap;

/**
 * Created by wei_qiang on 11/15/18.
 */

public class ReadFile
{
    private final File path;
    private final String cmd;
    private final String interval;

    private int min = 50;
    private int max = 800;

    public ReadFile(File path, String cmd, String interval) {
        this.path = path;
        this.cmd = cmd;
        this.interval  = interval;
    }

    public ReadFile(File path, String cmd, String interval, int min, int max) {
        this(path,cmd,interval);
        this.min = min;
        this.max = max;
    }

    public HashMap<String, float[]> run() throws IOException
    {

        HashMap<String,float[]> hashMethylationRate = new HashMap<String, float[]>();

        TabixReader tbx = new TabixReader(path.getAbsolutePath());
        TabixReader.Iterator iterator = tbx.query(interval);

        String line;
        String[] str_split;

        String regex = "\\s+";
        if(cmd.equals("saveRef"))
        {
            while(true)
            {
                line = iterator.next();

                if (line == null) {
                    break;
                }

                str_split = line.split(regex);

                if (str_split.length <4) {
                    continue;
                }

                String[] str_split3 = str_split[3].split(";"); //[0].split(":");


                String[] str_split2 = str_split3[0].split(":");

                float totalDepth = Float.parseFloat(str_split2[2]);

                if (totalDepth > min && totalDepth < max)
                {
                    float[] depth = new float[4];
                    depth[0] = Float.parseFloat(str_split2[3])/totalDepth;
                    depth[1] = Float.parseFloat(str_split2[4])/totalDepth;
                    depth[2] = Float.parseFloat(str_split2[5])/totalDepth;
                    depth[3] = Float.parseFloat(str_split2[6])/totalDepth;

                    hashMethylationRate.put(str_split[0]+":"+str_split[1], depth);

                }
            }
            tbx.close();
        }
        else if(cmd.equals("saveTumor"))
        {
            while(true)
            {
                line = iterator.next();

                if (line == null) {
                    break;
                }

                str_split = line.split(regex);

                if (str_split.length <4) {
                    continue;
                }

                String[] str_split3 = str_split[3].split(";"); //[0].split(":");


                if (str_split3.length <2) {
                   continue;
                }

                String[] str_split2 = str_split3[0].split(":");

                float totalDepth = Float.parseFloat(str_split2[2]);

                if (totalDepth > min && totalDepth < max)
                {
                    float[] depth = new float[4];
                    depth[0] = Float.parseFloat(str_split2[3])/totalDepth;
                    depth[1] = Float.parseFloat(str_split2[4])/totalDepth;
                    depth[2] = Float.parseFloat(str_split2[5])/totalDepth;
                    depth[3] = Float.parseFloat(str_split2[6])/totalDepth;

                    hashMethylationRate.put(str_split[0]+":"+str_split[1], depth);

                }
            }
            tbx.close();
        }

        return hashMethylationRate;
    }

}
