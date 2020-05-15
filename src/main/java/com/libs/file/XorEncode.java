package com.libs.file;

/**
 * Created by wei_qiang on 8/18/19.
 */

public class XorEncode
{
    private final byte[] abcBytes;
    private int j = 0;
    private int k = 0;

    public XorEncode(String abc)
    {
        // change into new key with 256 byte
        abcBytes=new byte[256];

        for(int i=0; i<abcBytes.length; i++){
            abcBytes[i]=abc.getBytes()[i%abc.getBytes().length];
        }
    }

    public byte[] run(byte[] data)
    {
        byte[] encodeBytes=new byte[data.length];

        for(int i=0; i<data.length; i++)
        {
            encodeBytes[i]=(this.run(data[i]));
        }
        return encodeBytes;
    }


    byte run(byte data)
    {
        int length = abcBytes.length;

        k=0xff & k+1;
        j=0xff & j+abcBytes[k%length];
        int m = abcBytes[k%length];
        abcBytes[k%length] = abcBytes[j%length];
        abcBytes[j%length] = (byte) m; //exchange k j

        int n = 0xff & abcBytes[j%length]+abcBytes[k%length];

        return (byte) (data^abcBytes[n%length]);
    }
}
