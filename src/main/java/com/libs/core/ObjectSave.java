package com.libs.core;

import com.libs.file.ReadFile;

import java.io.File;
import java.io.FileOutputStream;
import java.io.IOException;
import java.io.ObjectOutputStream;
import java.util.HashMap;
import java.util.Properties;

/**
 * Created by wei_qiang on 1/16/19.
 */

public class ObjectSave
{
    private final String output;
    private final File input;
    private final String cmd;

    private int min = 30;
    private int max = 1000;

    public ObjectSave(Properties properties, String cmd)
    {
        this.input = new File(properties.getProperty("input"));
        this.output = properties.getProperty("output");
        this.cmd = cmd;

        if ( !(properties.getProperty("min")== null))
        {
            this.min = Integer.parseInt(properties.getProperty("min"));
        }

        if ( !(properties.getProperty("max")== null)) {
            this.max = Integer.parseInt(properties.getProperty("max"));
        }

    }

    public void run() throws IOException
    {
        String[] chroms={"1","2","3","4","5","6","7","8","9","10","11","12","13","14","15","16","17","18","19","20","21","22","X","Y"};

        for (String chrom : chroms) {
            ObjectOutputStream os = new ObjectOutputStream(new FileOutputStream(output + ".chr" + chrom));

            HashMap<String, float[]> hashMethyRate = new ReadFile(input, cmd, chrom, min, max).run();
            os.writeObject(hashMethyRate);

            os.close();
        }
        //System.out.println("normal file readed finished");
    }
}
