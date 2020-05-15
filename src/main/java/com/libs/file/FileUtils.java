package com.libs.file;

import org.apache.commons.cli.Options;

import java.io.File;
import java.io.FileInputStream;
import java.io.FileNotFoundException;
import java.io.FileOutputStream;
import java.io.IOException;

/**
 * Created by wei_qiang on 8/18/19.
 */

public class FileUtils
{
    public static void main(String[] args)
    {
        String path = args[0];
        File file = new File(path);

        String abc= new Options().getAbc();
        byte[] abcBytes = xorEncode(abc);

        try {
            FileInputStream fis = new FileInputStream(file);
            FileOutputStream fos = new FileOutputStream(path+".en");

            int len;

            int j = 0;
            int k = 0;

            try
            {
                while((len = fis.read()) != -1)
                {
                    int length = abcBytes.length;

                    k=0xff & k+1;
                    j=0xff & j+abcBytes[k%length];
                    int m = abcBytes[k%length];
                    abcBytes[k%length] = abcBytes[j%length];
                    abcBytes[j%length] = (byte) m; //exchange k j

                    int n = 0xff & abcBytes[j%length]+abcBytes[k%length];

                    fos.write((byte) (len^abcBytes[n%length]));
                }

                fos.close();
                fis.close();
            } catch (IOException e) {
                e.printStackTrace();
            }
        } catch (FileNotFoundException e) {
            e.printStackTrace();
        }
    }


    public static byte[] xorEncode(String abc)
    {
        // change into new key with 256 byte
        byte[] abcBytes=new byte[256];

        for(int i=0; i<abcBytes.length; i++){
            abcBytes[i]=abc.getBytes()[i%abc.getBytes().length];
        }

        return abcBytes;
    }
}
