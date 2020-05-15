package com.libs.core;

import com.libs.file.DeClassLoader;

import java.io.BufferedReader;
import java.io.BufferedWriter;
import java.io.ByteArrayInputStream;
import java.io.File;
import java.io.FileWriter;
import java.io.IOException;
import java.io.InputStream;
import java.io.InputStreamReader;
import java.util.Properties;

/**
 * Created by wei_qiang on 8/24/19.
 */

public class IchorCNAMethylation
{
    private final String input;

    public IchorCNAMethylation(Properties properties)
    {
        this.input = properties.getProperty("output");
    }

    public void run() throws IOException, InterruptedException {

        String userDir = System.getProperty("user.dir");
        String fileName = new File(System.getProperty("java.class.path")).getParent();

        String id = input+".cna";
        String rscript = userDir+"/"+input+".temp";


        FileWriter fw = new FileWriter(rscript);
        BufferedWriter bw = new BufferedWriter(fw);

        DeClassLoader loaderRscript = new DeClassLoader();
        ByteArrayInputStream RscriptByte = new ByteArrayInputStream(loaderRscript.file2Bytes(new File(fileName+"/resource/ichorCNA.check4.R.en")));

        BufferedReader bfReader = new BufferedReader(new InputStreamReader(RscriptByte));

        String line;
        while((line = bfReader.readLine()) != null)
        {
            if (line.startsWith("#"))
            {
                continue;
            }

            bw.write(line+"\n");
        }

        bw.close();
        bfReader.close();
        RscriptByte.close();

        String cmds = "Rscript "+rscript+" --id "+id+" --BED "+input;

        System.out.println(cmds);

        Process pro = Runtime.getRuntime().exec(cmds);
        pro.waitFor();
        InputStream in = pro.getInputStream();
        BufferedReader read = new BufferedReader(new InputStreamReader(in));

        while((line = read.readLine())!=null){
            //System.out.println(line);
        }

        read.close();
        in.close();

        File file = new File(rscript);
        //file.delete();

        //file = new File(id);
        //file.delete();
    }

}
