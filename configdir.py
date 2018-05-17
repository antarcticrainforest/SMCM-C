'''
Created on August 13, 2013

@author: Martin Bergemann

@institution: Monash University School of Mathematical

@description: Collect information from a given config file and return the
              content as a dict-object
'''
import os, sys, string, glob
class Config(dict):
    """
    Config-class : Config(filename,**kwargs)
    Read a configuration file and return it's content stored in a dictionary 
    object.

    Parameters
    ----------
    filename : (str-obj) 
         the name of the configuration file to be read
    Keywords
    --------
    maketuple : (bool) Default: True
        if maketuple is True the method tries to interprete values with , 
        as seperaters for different tuple values
    skipwhitespace: (bool) Default: True
        whitspaces wont be considered if set True
    split : (str) Default =
        the seperator to seperate key and value in the configuration file

    Example
    -------
    Consider the following simple configuration file:
    #Filename of the test data
    filename = 'foo.nc' #
    variable = bar # The variable to be considered
     x1 = 9.0 # First index
    x2 =10  # Last index
    update = true
    times = 1,2,3 # Some time steps

    >>> C = Config('files.conf')
    >>> print(C)
    >>> Keys     | Values
    >>> -------------------------
    >>> update   | True
    >>> times    | (1.0, 2.0, 3.0)
    >>> x2       | 10
    >>> filename | foo.nc
    >>> variable | bar
    >>> x1       | 9.0

    """
    def __setattr__(self,k,v):
        if k in self.keys():
            self[k] = v
        elif not hasattr(self,k):
            self[k] = v
        else:
            raise AttributeError( "Cannot set '%s', class attribute already \
                    exists" % (k, ))
    def __getattr__(self, k):
        if k in self.keys():
            return self[k]
        raise AttributeError("Attribute '%s', deos not exist, available\
                attirbutes are: %s" %(k,self.keys().__repr__().strip(']')\
                .strip('[')))
    def __repr__(self):
        a=8
        b=8
        for v,k in self.items():
            if len(str(v)) > a:
                a = len(str(v))
            if len(str(k)) > b:
                b = len(str(k))
        a = a + 1
        b = b + 1
        s="Keys"+(a-4)*' '+"| Values\n"
        s=s+(a+b)*'-'+'\n'
        for v,k in self.items():
            s=s+str(v)+(a-len(str(v)))*' '+'| '+str(k)+'\n'
        return s
    def __init__(self,filename,maketuple=True,skipwhitespace=True,split='='):
        """ Try to open the file filename, read it and create a class instance 
            for every key entry 
            NOTE: Every key is an instance of Config but Config itself is of type dict
                  therefore the intances can be accesses both ways
        Example:
        -------
        from configdir import Config
        C = Config('Path/to/run.conf')
        >>> T = C['times']
        >>> t = C.times
        >>> T == t
        >>> Ture


        """

        try:
            if isinstance(filename,str):
                f=open(filename)
            file = f.readlines()
            f.close()
        except IOError:
            print( "Could not open %s: No such file or directory" %filename)
            return
        self.__get(file,maketuple,skipwhitespace,split=split)
        for key,value in self.items():
            try:
                if "$" in value:
                    var=value.split('/')[0]
                    try:
                        path=os.environ[var.replace('$','')]
                        self[key]=value.replace(var,path)
                    except KeyError:
                        raise KeyError('Environmet variable %s not set'%var)
            except TypeError:
                pass
                


    def __get(self,file,maketuple,skipwhitespace,split='='):
        """Explanation: This function takes a variable-name and extracts
        the according value form the config-file
        varname : the name of the variable that should be looked up"""
        
        #Define a blacklist of characters
        blacklist=[']','[','{','}','@','#','"',"'"]
        strip={'\t':'=','=':'\t'}[split]
        for i,line in enumerate(file):
            try:
                if not line.strip('\n')[0] in blacklist and split in line:
                    var=line.strip('\n').strip(strip).strip()\
                            .split(split)
                    if skipwhitespace :
                        value=var[1].replace(' ','')
                    else:
                        value=var[1]
                    for pos,v in enumerate(value):
                        if v.startswith('#'):
                            value=value[:pos+1]
                            break

                    for b in blacklist:
                        value=value.replace(b,'')
                    try:
                        value = int(value)
                    except ValueError:
                        try :
                            value = float(value)
                        except ValueError:
                            if value.lower() == 'false':
                                value = False
                            elif value.lower() == 'true':
                                value = True
                            elif value.lower() == 'none':
                                value = None
                    if isinstance(value,str):
                            #value=value.strip()
                        if ',' in value and maketuple:
                            value2=[]
                            for v in value.split(','):
                                try:
                                    value2.append(float(v.strip(')').strip('(')))
                                except ValueError:
                                    value2.append(os.path.expanduser(v).strip(')').strip('('))
                            value=tuple(value2)
                        else:
                            value = os.path.expanduser(value)
                    self[var[0].replace(' ','').strip()]=value
            except IndexError:
                pass
