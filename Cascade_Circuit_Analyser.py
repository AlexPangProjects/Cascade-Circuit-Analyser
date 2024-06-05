import sys
import os
import numpy as np
import csv

# T1, T2, T3 tests to check for input file blocks
# T4 read component values
# main program

class Component:
    
    #creates a class called Component that contains all the values that the component would have
    #e.g. start and end nodes, impedance etc
    
    def __init__(self, n1, n2, R, C, L, G, Z):
        self.n1 = n1 #node 1
        self.n2 = n2 #node 2
        self.R = R #resistance
        self.C = C #capcitance
        self.L = L #inductance
        self.G = G #conductance
        self.Z = Z #impedance
    
    def T4(self):    #function that prints the component values to aid in troubleshooting and testing
        sn1 = str(self.n1)
        sn2 = str(self.n2)
        sR = str(self.R)
        sC = str(self.C)
        sL = str(self.L)
        sG = str(self.G)
        sZ = str(self.Z)
        print("Component values: start node=", sn1, " end node=", sn2, " R=", sR, " C=", sC, " L=", sL, " G=", sG, " Z=", sZ)
        
class Terms:
    
    #creates a class called Terms that contains all the values found in the <Terms> block
    
    def __init__(self, VT, RS, RL, Fstart, Fend, Nfreqs):
        self.VT = VT #thevenin voltage
        self.RS = RS #source resistance
        self.RL = RL #load resistance
        self.Fstart = Fstart #start frequency
        self.Fend = Fend #end frequency
        self.Nfreqs = Nfreqs #n of frequency steps

    def T5(self): #function that prints the terms values to aid in troubleshooting and testing
        sVT = str(self.VT)
        sRS = str(self.RS)
        sRL = str(self.RL)
        sFstart = str(self.Fstart)
        sFend = str(self.Fend)
        sNfreqs = str(self.Nfreqs)
        print("Terms: VT=", sVT, " RS=", sRS, " RL=", sRL, " Fstart=", sFstart, " Fend=", sFend, " Nfreqs=", sNfreqs)

def output_empty_csv(outname):
    # Create an empty CSV file with the given filename
    with open(outname, 'w', newline='') as file:
        writer = csv.writer(file)
        # Optionally write headers or leave completely empty
        
def blockcheck(block,file):
    
    #checks if block delimiters are present
    
    T1 = file.find(block) #returns -1 if a block delimited cannot be found in the file
    if T1 == -1:
        raise Exception(block + " BLOCK NOT FOUND")
    
        output_empty_csv(sys.argv[2])
        
        sys.exit(3)
           
def commentremover(infile):
    
    #removes comments from the input file and returns a string without the comments
    #:param infile: input file
    #:return: string containg cleaned version of input file
    
    cleanstring = "" #generates empty string
    for line in open(infile):
        li=line.strip() #removes leading and trailing whitespaces
        if not li.startswith("#"): #checks if the line starts with #
            cleanstring = cleanstring +  line.rstrip() + "\n" #if the line does not start with # it appends this to the empty string with a newline character
    return cleanstring
           
def extract(name, string1):
    
    #extracts value from a string and returns it as a float
    #:param name: unit value that needs to be extracted
    #:param string1: string to be searched for unit value
    #:return: unit value
    
    EX3 = 0  # Initialize with a default value
    EX = string1.find(name) #searches for unit, will return -1 if not found   
    if EX > -1: #if unit value is found
        EX1 = string1.index(name) #returns index of the unit
        EX2 = EX1 + len(name) #finds character next to unit, this is the value desired
        EX5 = string1.find(" ", EX2) 
        EX6 = string1.find("\n", EX2)
        #finds next closest spaces or newlines
        if EX5 > -1 and EX6 > -1:
            EX4 = min(EX5, EX6) 
        elif EX5 > -1:
            EX4 = EX5
        elif EX6 > -1:
            EX4 = EX6
        else:
            # If no space or newline is found, take the substring till the end
            EX4 = len(string1)
        EX3 = string1[EX2:EX4] #gives characters between unit and next space or newline which is the value desired
 
        EX3 = float(EX3) #converts string to float
        
            
    return EX3

def extractint(name, string2):
   
    #extracts value from a string and returns it as a float
    #:param name: unit value that needs to be extracted
    #:param string1: string to be searched for unit value
    #:return: unit value
    
    EX3 = 0  # Initialize with a default value
    EX = string2.find(name) #searches for unit, will return -1 if not found   
    if EX > -1: #if unit value is found
        EX1 = string2.index(name) #returns index of the unit
        EX2 = EX1 + len(name) #finds character next to unit, this is the value desired
        EX5 = string2.find(" ", EX2) 
        EX6 = string2.find("\n", EX2)
        #finds next closest spaces or newlines
        if EX5 > -1 and EX6 > -1:
            EX4 = min(EX5, EX6) 
        elif EX5 > -1:
            EX4 = EX5
        elif EX6 > -1:
            EX4 = EX6
        else:
            # If no space or newline is found, take the substring till the end
            EX4 = len(string2)
        EX3 = string2[EX2:EX4] #gives characters between unit and next space or newline which is the value desired
 
        EX3 = int(EX3) #converts string to integer
        
            
    return EX3
     
def outputscraper(name1, name2, string):
    
    #returns the output file format
    #:param name1: delimiter the function should start copying from
    #:param name2: delimiter the function should stop copying at
    #:param string: string that the function should copy from
    #:return: copied text from string (output file format)
    
    outputflag = 0
    outlist = []

    for line in string.split('\n'):
        outon = line.find(name1)
        outoff = line.find(name2)
        if outoff>-1:
            break
        if outputflag == 1:
            outlist.append(line)
        if outon>-1:
            outputflag = 1
    return outlist
            
def genfreaklist(Fe,Fs,Fn):
    
    #returns the list of frequencies that the circuit needs to be tested at
    #:param Fe: frequency end
    #:param Fs: frequency start
    #:param Fn: number of frequency steps
    #:return: list of frequences
    
    
    Fend1 = int(Fe)
    Fstart1 = int(Fs)
    Nfreqs1 = int(Fn)

    decrement = (Fend1 - Fstart1) / (Nfreqs1 - 1)
    frequencies = [Fstart1] 

    for i in range(1, Nfreqs1 - 1):
        next_freq = Fstart1 + i * decrement
        frequencies.append(next_freq)
    frequencies.append(Fend1)
    return frequencies  

def formatlistnumerical(list,numformat):
    
    #formats the values in a list
    #:param list: list to be formatted
    #:param numformat: format that the items in the list should follow
    #:return: list with desired formatting
    
    formattedlist = []
    for i in range(len(list)):
        formattedlist.append(numformat.format(list[i]))
    return formattedlist

def custom_sort_key(obj, attributes):
    return tuple(getattr(obj, attr) for attr in attributes)       
    
def capacitorimpedance(C,angfrq):
    
    #calculates capacitor impedance
    #:param C: capacitance
    #:param angfrq: angular frequency
    #:return: impedance
    
    impedance = -1j*(1/(C*angfrq))
    return impedance

def inductorimpedance(L,angfrq):
    
    #calculates inductor impedance
    #:param L: inductance
    #:param angfrq: angular frequency
    #:return: impedance
    
    impedance = 1j * angfrq * L
    return impedance

def conductanceimpedance(G):
    
    #calculates conductance impedance
    #:param G: conductance
    #:return: impedance
    
    impedance = 1/G
    return impedance

def attrcheck(obj, attribute):
    
    #Check if an object attribute has a non-zero value.
    #:param obj: The object whose attribute needs to be checked.
    #:param attribute: The name of the attribute to check.
    #:returns: True if the attribute has a non-zero value, False otherwise.
    
    try:
        attr_value = getattr(obj, attribute)
        return attr_value != (0)
    except AttributeError:
        return False
         
def impedancecalc(obj,angfrq):
    
    #Calculates the impedance of an object
    #:param obj: object to calculate impedance for
    #:param angfrq: angular frequency at which to calculate impedance
    #:return: calculated impedance
    
    Rflag = attrcheck(obj, "R")
    Cflag = attrcheck(obj, "C")
    Lflag = attrcheck(obj, "L")
    Gflag = attrcheck(obj, "G")
    if Rflag == 1:
        return obj.R
        
    if Cflag == 1:
        Z = capacitorimpedance(obj.C,angfrq)
        return Z
    
    if Lflag == 1:
        Z = inductorimpedance(obj.L,angfrq)
        return Z
    
    if Gflag == 1:
        Z = (conductanceimpedance(obj.G))
        return Z

def seriesmatrixcalc(Z):
    
    #creates abcd matrix for a series element
    
    matrixA = np.array([[1 , Z],[0 , 1]])
    return matrixA
            
def shuntmatrixcalc(Z):
    
    #creates abcd matrix for a shunt element
    
    Z = 1/Z
    matrixB = np.array([[1 , 0],[Z , 1]])
    return matrixB    
    
def matrixcalc(obj):
    
    #checks if element is shunt or series then calculates the appropriate abcd matrix
    
    seriesflag = attrcheck(obj, "n2")
    if seriesflag == 1:
        matrix = seriesmatrixcalc(obj.Z)
        return matrix
    else:
        matrix = shuntmatrixcalc(obj.Z)
        return matrix
        
def cascadematrix(matrixlist):
    
    #calculates the overall cascade abcd matrix 
    
    hold = matrixlist[0]
    count = (len(matrixlist)) - 1
    for i in range(count):
        hold = np.dot(hold, (matrixlist[i+1]))
    return hold

def zincalc(cascmatrix,rl):
    
    #calculates zin
    
    zin = (cascmatrix[0][0]*rl + cascmatrix[0][1])/(cascmatrix[1][0]*rl + cascmatrix[1][1])
    return zin

def zoutcalc(cascmatrix,rs):
    
    #calculates zout
    
    zout = (cascmatrix[1][1]*rs + cascmatrix[0][1])/(cascmatrix[1][0]*rs + cascmatrix[0,0])
    return zout

def vincalc(cascmatrix,rl,vt,rs):
    
    #calculates vin
    
    vin = vt*(zincalc(cascmatrix,rl)/(rs+zincalc(cascmatrix,rl)))
    return vin

def voutcalc(cascmatrix,rs,rl,vt):
    
    #calculates vout
    
    vout = rl*ioutcalc(cascmatrix,vt,rl,rs)
    return vout

def iincalc(cascmatrix,rl,vt,rs):
    
    #calculates Iin
    
    iin = vincalc(cascmatrix,rl,vt,rs)/zincalc(cascmatrix,rl)
    return iin

def ioutcalc(cascmatrix,vt,rl,rs):
    
    #calculates Iout
    
    iout = iincalc(cascmatrix,rl,vt,rs)/(cascmatrix[1][0]*rl + cascmatrix[1][1])
    return iout

def avcalc(cascmatrix,rs,rl,vt):
    
    #calculates av
    
    av = voutcalc(cascmatrix,rs,rl,vt)/vincalc(cascmatrix,rl,vt,rs)
    return av

def aicalc(cascmatrix,vt,rl,rs):
    
    #calculates ai
    
    ai = ioutcalc(cascmatrix,vt,rl,rs)/iincalc(cascmatrix,rl,vt,rs)
    return ai

def aiconj(cascmatrix,vt,rl,rs):
    
    #calculates the conjugate of ai
    
    ai = ioutcalc(cascmatrix,vt,rl,rs)/iincalc(cascmatrix,rl,vt,rs)
    aiconjugate = ai.conjugate()
    return aiconjugate 

def pincalc(cascmatrix,rl,vt,rs):
    
    #calculates pin
    
    pin = (vincalc(cascmatrix,rl,vt,rs)*np.conj(iincalc(cascmatrix,rl,vt,rs)))
    return pin

def apcalc(cascmatrix,rs,rl,vt,):
    
    #calculates ap
    
    av = avcalc(cascmatrix,rs,rl,vt)
    ap = av * aiconj(cascmatrix,vt,rl,rs)
    return ap

def poutcalc(cascmatrix,rs,rl,vt):
    
    #calculates pout
    
    pout = pincalc(cascmatrix,rl,vt,rs) * apcalc(cascmatrix,rs,rl,vt)
    return pout
    
def outputvaluelist(outlist,termsobj,cascmatrix):
    
    #returns a list of output values according to the outputfile format
    #:param outlist: output file format
    #:param termsobj: Terms class/object
    #:param cascmatrix: cascade abcd matrix
    #:return: list of output values 
    
    outvalueslist = []
    for i in range(len(outlist)):
        if outlist[i] == "Vin V":
            outvalueslist.append(vincalc(cascmatrix,termsobj.RL,termsobj.VT,termsobj.RS))
        elif outlist[i] == "Vout V":
            outvalueslist.append(voutcalc(cascmatrix,termsobj.RS,termsobj.RL,termsobj.VT))
        elif outlist[i] == "Iin A":
            outvalueslist.append(iincalc(cascmatrix,termsobj.RL,termsobj.VT,termsobj.RS))
        elif outlist[i] == "Iout A":
            outvalueslist.append(ioutcalc(cascmatrix,termsobj.VT,termsobj.RL,termsobj.RS))
        elif outlist[i] == "Pin W":
            outvalueslist.append(pincalc(cascmatrix,termsobj.RL,termsobj.VT,termsobj.RS))
        elif outlist[i] == "Zout Ohms":
            outvalueslist.append(zoutcalc(cascmatrix,termsobj.RS))
        elif outlist[i] == "Pout W":
            outvalueslist.append(poutcalc(cascmatrix,termsobj.RS,termsobj.RL,termsobj.VT))
        elif outlist[i] == "Zin Ohms":
            outvalueslist.append(zincalc(cascmatrix,termsobj.RL))
        elif outlist[i] == "Av":
            outvalueslist.append(avcalc(cascmatrix,termsobj.RS,termsobj.RL,termsobj.VT))
        elif outlist[i] == "Ai":
            outvalueslist.append(aicalc(cascmatrix,termsobj.VT,termsobj.RL,termsobj.RS))
    return outvalueslist
      
def reallist(outlist,termsobj,cascmatrix):
    
    #returns real components of the output values list
    
    rlist = []
    modlist = outputvaluelist(outlist,termsobj,cascmatrix)
    for i in range(len(modlist)):
        realn = modlist[i].real
        rlist.append(realn)
    return rlist

def imagelist(outlist,termsobj,cascmatrix):
    
    #returns imaginary components of the output values list
    
    ilist = []
    modlist = outputvaluelist(outlist,termsobj,cascmatrix)
    for i in range(len(modlist)):
        imagn = modlist[i].imag
        ilist.append(imagn)
    return ilist

def combinedlist(outlist,termsobj,cascmatrix):
    
    #returns the combined list of real and imaginary components of the output values
    
    newlist = []
    realist = reallist(outlist,termsobj,cascmatrix)
    imaglist = imagelist(outlist,termsobj,cascmatrix)
    for i in range(len(realist)):
        newlist.append("{:.3e}".format(realist[i]))
        newlist.append("{:.3e}".format(imaglist[i]))
    return newlist

def combinedliststring(outlist,termsobj,cascmatrix):
    
    #converts the combined list of output values to strings and formats them
    
    newlist = []
    for i in range(len(combinedlist(outlist,termsobj,cascmatrix))):
        for j in range(len(combinedlist(outlist,termsobj,cascmatrix)[i])):
            a = combinedlist(outlist,termsobj,cascmatrix)[i][j]
            a = str(a)
            a = a.rjust(11)
            newlist.append(a)
    return newlist
                    
def outputvaluelistwithfreq(freqlist,outlist,termsobj,cascmatrixlist):
    
    #creates a nested list of the output variables according to the different frequencies
    
    outvalfreqlist = []
    for i in range(len(freqlist)):
        a = combinedlist(outlist,termsobj,cascmatrixlist[i])
        outvalfreqlist.append(a)
    return outvalfreqlist
        
def blockformatter(listtobeformatted,length):
    
    #formats values in a list before they are written to a csv file
    
    for i in range(len(listtobeformatted)):
        holdvalue = listtobeformatted[i]
        holdvalue = str(holdvalue)
        holdvalue = holdvalue.rjust(length)
        listtobeformatted[i] = holdvalue
        
try:             

    pi = np.pi #creates pi variable

    if (len(sys.argv)<3): #if less than 3 files are provided the program does not run and outputs an error message
        
        print("\n\tError, program needs two arguments to run\n" )
        
        output_empty_csv(sys.argv[2])
        
        sys.exit(3)
    
    input_filename=sys.argv[1] #the 2nd file provided is the assigned to a variable


    try:
        
        fin=open( input_filename,'rt') #reads the input file
        
        
    except FileNotFoundError:
        
        print('File <%s> not found'%(input_filename))
        
        output_empty_csv(sys.argv[2])
        
        sys.exit(3)
        
     
    try:
        
        clean = commentremover(input_filename) #cleans the input file text of comments
        
        
    except:
        
        print("error cleaning input file")
        
        output_empty_csv(sys.argv[2])
        
        sys.exit(3)

    #checks if all the required blocks are present in the input file

    blockcheck("<CIRCUIT>",clean)
    
    blockcheck("<TERMS>",clean)
    
    blockcheck("<OUTPUT>",clean)
    
    
    try:
        
        #creates a list of components with the values provided in the input file
    
        componentlist = [] #initializes empty list
        for line in clean.split('\n'): #iterates through the string that contains the values needed for the objects
            A = line.find("n1=") #searches string for the desired value
            if A>-1:
                A3 = extractint("n1=",line)
                B3 = extractint("n2=",line)
                C3 = extract("R=",line)
                D3 = extract("C=",line)    
                E3 = extract("L=",line)
                F3 = extract("G=",line)
                #extracs values and assigns them to variables    
                componentlist.append(Component(A3,B3,C3,D3,E3,F3,float(0))) #creates objects with the extracted values and appends them to the previously empty list   
        
        #sorting order for the list of component objects
                        
        attributes_to_sort_by = ['n1', 'n2']
        
        #creates a new list of the sorted components
                    
        newcomplist = sorted(componentlist, key=lambda x: custom_sort_key(x, attributes_to_sort_by)) #creates a new list of components sorted by nodes with the start node taking priority over the end node
        
        
    except:
        
        print("error extracting components")
        
        output_empty_csv(sys.argv[2])
        
        sys.exit(3)
        
        
    try:    
        
        #creates an object of class terms
    
        VT1 = extract("VT=",clean)
        RS1 = extract("RS=",clean)
        RL1 = extract("RL=",clean)
        Fs1 = extract("Fstart=",clean)
        Fed1 = extract("Fend=",clean)
        Nfreq1 = extract("Nfreqs=",clean)
        
        Termsvalues = Terms(VT1,RS1,RL1,Fs1,Fed1,Nfreq1)
        
        
    except:
        
        print("error extracting terms")
        
        output_empty_csv(sys.argv[2])
        
        sys.exit(3)


    try:

        #creates a list of text/units found in the output block
        
        outputlist = outputscraper("<OUTPUT>", "</OUTPUT>",clean)
        
        #runs through the measurements found in the output block and creates a list of titles needed for the csv output file
        
        outlabelslist = []
        for i in range(len(outputlist)):
            if outputlist[i] == "Vin V":
                outlabelslist.append("Vin")
            elif outputlist[i] == "Vout V":
                outlabelslist.append("Vout")
            elif outputlist[i] == "Iin A":
                outlabelslist.append("Iin")
            elif outputlist[i] == "Iout A":
                outlabelslist.append("Iout")
            elif outputlist[i] == "Pin W":
                outlabelslist.append("Pin")
            elif outputlist[i] == "Zout Ohms":
                outlabelslist.append("Zout")
            elif outputlist[i] == "Pout W":
                outlabelslist.append("Pout")
            elif outputlist[i] == "Zin Ohms":
                outlabelslist.append("Zin")
            elif outputlist[i] == "Av":
                outlabelslist.append("Av")
            elif outputlist[i] == "Ai":
                outlabelslist.append("Ai")
            elif outputlist[i] == "Vin dBV":
                outlabelslist.append("Vin")
            elif outputlist[i] == "Vout dBV":
                outlabelslist.append("Vout")
            elif outputlist[i] == "Iin dBA":
                outlabelslist.append("Iin")
            elif outputlist[i] == "Pin dBW":
                outlabelslist.append("Pin")
            elif outputlist[i] == "Av dB":
                outlabelslist.append("Av")
        
        #adds real or imaginary tags to the list of titles
        
        outputlabel = []
        for i in range(len(outlabelslist)):
            outputlabel.append("Re("+outlabelslist[i]+")")
            outputlabel.append("Im("+outlabelslist[i]+")")
    
        blockformatter(outputlabel,11)
        
        #finds what units are necessary for the titles found in the output title list previously created
        
        outlabelunitlist = []
        for i in range(len(outputlist)):
            if outputlist[i] == "Vin V":
                outlabelunitlist.append("V")
            elif outputlist[i] == "Vout V":
                outlabelunitlist.append("V")
            elif outputlist[i] == "Iin A":
                outlabelunitlist.append("A")
            elif outputlist[i] == "Iout A":
                outlabelunitlist.append("A")
            elif outputlist[i] == "Pin W":
                outlabelunitlist.append("W")
            elif outputlist[i] == "Zout Ohms":
                outlabelunitlist.append("Ohms")
            elif outputlist[i] == "Pout W":
                outlabelunitlist.append("W")
            elif outputlist[i] == "Zin Ohms":
                outlabelunitlist.append("Ohms")
            elif outputlist[i] == "Av":
                outlabelunitlist.append("L")
            elif outputlist[i] == "Ai":
                outlabelunitlist.append("L")
            elif outputlist[i] == "Vin dBV":
                outlabelunitlist.append("dBV")
            elif outputlist[i] == "Vout dBV":
                outlabelunitlist.append("dBV")
            elif outputlist[i] == "Iin dBA":
                outlabelunitlist.append("dBA")
            elif outputlist[i] == "Pin dBW":
                outlabelunitlist.append("dBW")
            elif outputlist[i] == "Av dB":
                outlabelunitlist.append("dB")
                
        #duplicates the unit list to account for real and imaginary measurements
                
        outputlabelunitlistfinal = []
        for i in range(len(outputlist)):
            outputlabelunitlistfinal.append(outlabelunitlist[i])
            outputlabelunitlistfinal.append(outlabelunitlist[i])     
    
        blockformatter(outputlabelunitlistfinal,11)
        
        
    except:
        
        print("error extracting output file format")
        
        output_empty_csv(sys.argv[2])
        
        sys.exit(3)
            
            
    try:
        
        #creates a list of frequencies from the frequency attributes in the terms object
        
        freaklist = genfreaklist(Termsvalues.Fend,Termsvalues.Fstart,Termsvalues.Nfreqs)
        
        #formats the frequency values in scientific notation
        
        freaklist3e = formatlistnumerical(freaklist,"{:.3e}")
        
        blockformatter(freaklist3e,10)

        #generates a list of angular frequencies using the previously created frequency list
        
        anglefreaklist =  [i * 2 * pi for i in freaklist]
        
        
    except:
        
        print("error calculating frequencies")
        
        output_empty_csv(sys.argv[2])
        
        sys.exit(3)


    try:
        
        #creates a nested list of objects with the top list corresponding to the different frequency steps
        
        nestedlistofcomponents = []

        for angle_freq in anglefreaklist:
            temp_list1 = []  # Temporary list to hold components with impedances for the current angle frequency
            for comp in newcomplist:
                cloned_object = Component(comp.n1, comp.n2, comp.R, comp.C, comp.L, comp.G, comp.Z)  # Create a clone of the component
                cloned_object.Z = impedancecalc(comp, angle_freq)  # Calculate and assign impedance for the clone
                temp_list1.append(cloned_object)  # Append the clone to the temporary list
            nestedlistofcomponents.append(temp_list1)  # Append the temporary list to the storage list
    

    except:
        
        print("error creating list of list of components with impedance values")
        
        output_empty_csv(sys.argv[2])
        
        sys.exit(3)
    
    
    try:   
        
        #creates a nested list of abcd matrices calculated from the previously created nested components list
         
        matrixgrid = []
        for smallcomplist in nestedlistofcomponents:
            small_list = []
            for obj in smallcomplist:
                matrix = matrixcalc(obj)
                small_list.append(matrix)
            matrixgrid.append(small_list)

        #calculates the overall abcd matrix for the circuit and puts it in a list according to the frequencies
        
        cascadematrixlist = []
        for i in range(len(matrixgrid)):
            finmatrix = cascadematrix(matrixgrid[i])
            cascadematrixlist.append(finmatrix)
       
        #calculates the output values with real and imaginary parts using the previously calculated overall abcd matrices and places them in a nested list 
        
        combinedfinaloutputvalueslist = outputvaluelistwithfreq(freaklist,outputlist,Termsvalues,cascadematrixlist)
        
        #formats the lists for the output csv file
        
        for i in range(len(combinedfinaloutputvalueslist)):
        
            blockformatter(combinedfinaloutputvalueslist[i],11)
        
            (combinedfinaloutputvalueslist[i]).append("'")
    
    
    except:
        
        print("error creating list of overall cascade matrix values")
        
        output_empty_csv(sys.argv[2])
        
        sys.exit(3)
    

    try:
        
        #writes the ouputs to a csv file
                    
        with open(sys.argv[2], 'w', newline='') as file:
        
            writer = csv.writer(file)
        
            writer.writerow(["Freq".rjust(10)] + outputlabel)
        
            writer.writerow(["Hz".rjust(10)] + outputlabelunitlistfinal)
        
            for i in range(len(combinedfinaloutputvalueslist)):
            
                writer.writerow([freaklist3e[i]] + combinedfinaloutputvalueslist[i])
                
    except:
        print("error writing to output csv file")
        
        output_empty_csv(sys.argv[2])
        
        sys.exit(3)
    
    print(sys.argv[2] + " successfully created")
   
except Exception as e:
    
    print(f"An unexpected error occurred: {e}")
    
    output_empty_csv(sys.argv[2])
    
    sys.exit(3)
    
    
    