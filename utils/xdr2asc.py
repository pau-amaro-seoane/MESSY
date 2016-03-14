#!/usr/bin/python

"""xdr2asc.py: read a portable binary XDR file with the format
used by ME(SSY)**2, and output ascii data to
standard output. If all arrays have same dimension, a flat, multi-column
output is produced (with a 2-line header). Otherwise, each array is
output after the other, separated by a header (not tested yet).

Calling sequence:

   xdr2asc.py [array 1 array2 ...] xdr_file
   
The last argument must be the name of the xdr file.
Optionally, other arguments can be given to specify which arrays to output.

For example, 'xdr2asc.py MyRun%0001200000%AMAS.xdr' will output all
the basic particle quantities: T_SE, J_SE, M_SE, R_SE, i.e. specific
kinetic energy, specific angular momentum, mass and radius, all in
N-body units. And 'xdr2asc.py iTet_SE MyRun%0001200000%ETOILES.xdr' will
output the stellar types only (array 'iTet_SE').

--- Marc Dewi Freitag, IoA, 10/03/08"""


def my_zip(array_of_array) :
    res=[]
    array_len=len(array_of_array[0])
    n_arrays=len(array_of_array)
    for i in range(0,array_len) :
        qq=[]
        for j in range(0,n_arrays) :
            qq.append(array_of_array[j][i])
        res.append(tuple(qq))
    return res
    
import sys, xdrlib

from numpy import array

if len(sys.argv)==1:
    print >> sys.stderr, "Are you joking? I need a file name to work!"
    sys.exit(1)

arrays_to_ouput=[]
if len(sys.argv)>2: # a list of arrays to output can be provided
    arrays_to_ouput=sys.argv[1:-1]
    sys.stderr.write('> looking for arrays: %s\n' % ', '.join(arrays_to_ouput))
    
FileName = sys.argv[-1] # last argument on command line should be file name


types_dic={
    '_i_':'integer', '_ti_':'integer',
    '_r_':'real*4', '_tr_':'real*4',
    '_d_':'real*8', '_td_':'real*8',
    '_s_':'string'
       }
is_array_dic={
    '_i_':False, '_ti_':True,
    '_r_':False, '_tr_':True,
    '_d_':False, '_td_':True,
    '_s_':False
       }

output_formats={ 'integer':'%i', 'real*4':'%13.7e', 'real*8':'%22.15e', 'string':'%s' }

# Open and read file

# Will produce xdr data one "line" at a time from the opened file
input = open(FileName, mode='rb').read() # 'r' for 'read', 'b' for 'binary'

unpacker = xdrlib.Unpacker(input) # class to unpack xdr data from the stream "data"

FileType = unpacker.unpack_string()

arrays=[]
names=[]
types=[]

DimPrec=-1
Narrays=0
SameDim = True

while 1:
    try:
        VarName = unpacker.unpack_string()
        VarType = unpacker.unpack_string()
        VarIsArray = is_array_dic[VarType]
        VarType = types_dic[VarType]
        VarDim = unpacker.unpack_int()
        if VarIsArray :
            sys.stderr.write('> Reading array %s with %s elements of type %s' %(VarName,VarDim,VarType))
        else :
            sys.stderr.write('> Reading single value %s of type %s' %(VarName,VarType))
            
        if VarType=='real*8' :
            method=unpacker.unpack_double
        elif VarType=='real*4' :
            method=unpacker.unpack_float
        elif VarType=='integer' :
            method=unpacker.unpack_int
        elif VarType=='string' :
            method=unpacker.unpack_string
        try: 
            if VarIsArray:
                arrays.append(array(unpacker.unpack_array(method)))
            else:
                arrays.append(array([method()]))
            
        except xdrlib.ConversionError, instance:
            sys.stderr.write('\n!!! Error xdrlib.ConversionError: ',instance.msg)
            sys.exit(2)
            
        names.append(VarName)
        types.append(VarType)
        sys.stderr.write(' <\n')
        if VarName=='R_SE' or VarName=='M_SE' :
            sys.stderr.write('  > Removing dummy element 0 of %s\n' % VarName)
            arrays[-1] = arrays[-1][1:]
            VarDim-=1
        SameDim = DimPrec<0 or (SameDim and DimPrec==VarDim)
        DimPrec = VarDim
        Narrays=Narrays+1
    except EOFError:
        break

if arrays_to_ouput : # only output a sub-selection of arrays
    
    DimPrec=-1
    Narrays=0
    SameDim = True
    
    all_arrays=arrays
    arrays=[]
    all_names=names
    names=[]
    all_types=types
    types=[]
    for name in arrays_to_ouput :
        if name in all_names:
            pos=all_names.index(name)
            arrays.append(all_arrays[pos])
            names.append(all_names[pos])
            types.append(all_types[pos])
            VarDim=len(arrays[-1])
            SameDim = DimPrec<0 or (SameDim and DimPrec==VarDim)
            DimPrec = VarDim
            Narrays=Narrays+1
    del all_arrays
    del all_names
    del all_types
    
if Narrays<1 :
    sys.stderr.write('!!! No data found !!!\n')
    sys.exit(1)
    
if SameDim and VarDim>1:
    sys.stderr.write('\n> Got %i arrays of same dimension (%i) => Can write a flat, multi-column file\n' %(Narrays,VarDim))
    header='# Content of XDR file\n# '
    format=''
    
    for i in range(0,Narrays) :
        name=names[i]
        type=types[i]
        header+=r' %i: %s' %(i+1,name)
        format+=' %s' % output_formats[type]
        
    print header

    for i in range(0,VarDim) :
        my_list=[]
        for quantity in arrays :
            my_list.append(quantity[i])
                    
        print format % tuple(my_list)

else :
    if SameDim :
        sys.stderr.write('\n> Only got single values. Will write one after the other\n')

    else :
        sys.stderr.write('\n> Got arrays of different length. Will write one after the other\n')

    print '# Content of XDR file\n#'
    for i in range(0,Narrays) :
        VarDim=len(arrays[i])
                
        print '# Name : %s' % names[i]
        print '# Type : %s' % types[i]
        print '# Dim  : %i' % VarDim
        format=output_formats[types[i]]
        for j in range(0,VarDim) :
            print format % arrays[i][j]

sys.exit(0)
