package com.libs.file;

import org.apache.commons.cli.Options;

import java.io.ByteArrayOutputStream;
import java.io.File;
import java.io.FileInputStream;
import java.io.IOException;
import java.io.InputStream;

/**
 * Created by wei_qiang on 8/18/19.
 */
public class DeClassLoader extends ClassLoader
{

    @Override
    protected Class<?> findClass(String name) throws ClassNotFoundException
    {
        String fileName = getName(name);

        //System.out.println("filename:"+fileName);
        //System.out.println("name:"+name);
        //URL fileURL=this.getClass().getResource(fileName);
        //System.out.println(fileURL.getFile());

        InputStream is = this.getClass().getResourceAsStream(fileName);

        byte[] data = file2Bytes(is);

        if(data == null)
        {
            return super.findClass(name);
        }
        else
        {
            return defineClass(name,data,0,data.length);
        }
    }

    private String getName(String name)
    {
        name = "/"+name.replace(".","/");

        return name+".class.en";
    }


    public byte[] file2Bytes(InputStream is)
    {
        String abc= new Options().getAbc();
        XorEncode encode = new XorEncode(abc);

        byte[] data = null;

        try
        {
            ByteArrayOutputStream bos = new ByteArrayOutputStream();

            int len;
            byte b;

            try {
                while ((len = is.read()) != -1)
                {
                    b = encode.run((byte) len);
                    //b = (byte) len;
                    bos.write(b);
                }
            } catch (IOException e) {
                e.printStackTrace();
            }

            data = bos.toByteArray();
            is.close();
            bos.close();

        } catch (IOException e) {
            e.printStackTrace();
        }
        return data;
    }


    public byte[] file2Bytes(File file)
    {
        String abc= new Options().getAbc();
        XorEncode encode = new XorEncode(abc);

        byte[] data = null;

        try {
            FileInputStream is = new FileInputStream(file);
            ByteArrayOutputStream bos = new ByteArrayOutputStream();

            int len;
            byte b;

            try {
                while ((len = is.read()) != -1)
                {
                    b = encode.run((byte) len);
                    bos.write(b);
                }
            } catch (IOException e) {
                e.printStackTrace();
            }

            data = bos.toByteArray();
            is.close();
            bos.close();

        } catch (IOException e) {
            e.printStackTrace();
        }
        return data;
    }
}


