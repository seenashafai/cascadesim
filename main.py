######  Libraries  ######
import sys
import csv
import math

######  Functions  ######
def unitType(): #Determine output format: standard form or exponent prefixes.
    choice = "SF"
    while True: #This loop makes sure that only valid answers are accepted.
        if not(choice == "SF" or choice == "prefix"):
            choice = input("Input Error: Please type 'SF' or 'prefix' to continue.")
        else:
            break
    return choice #choice is returned for use in the prefix() function

def clean(text): #Standardises the code and removes unwanted characters.
    array = []
    for line in text:
        line = line.replace("\n", "") #Removes instances of "\n"
        #If statement removes commented lines.
        if "#" not in line:
            array.append(line) #if no hashtag, the whole line is appended to the output array
        else:
            array.append(line[:line.find("#")]) #if a hashtag is present, all data infront of the hashtag is appended to the output array.
    array = list(filter(None, array)) #empty lines are removed
    return array

def sortInput(text, header): #extracts data blocks from cleaned data
    i=0 #counter used to get line indexes of header and footer.
    for line in text:
        if str('<'+ header +'>') in line: #when header found, store the index in min.
            min = i
        elif str('</'+ header +'>') in line: #when footer found, store the index in max.
            max = i
        else:
            i=i+1  #if neither found, continue searching
    return text[min+1:max+1] #return the data between the header and footer.

def logOrLin(array): #determine whether a linear or logarithmic frequency sweep is required
    for index in array:
        if "LFend" in index and "LFstart" in index: #if logarithmic set boolean to True
            logBoo = True           
        else:   #else set it to False
            logBoo = False
    return logBoo #return for use in getFreq()

def extractTerms(array, word):  #extracts numeric data from terms data
    for index in array: #for each line...
        index = index.split(" ") #....split up values on the same line
        for sub in index:       #for each value....
            if word in sub:    
                mod = sub.find('=')
                return sub[mod+1:] #...extract the numeric data after the '='
    return "N/A" #if keyword not found, return 'N/A'

def extractOutput(array, word):  #extracts units from output data
    for index in array:
        if word in index and " " in index:
            indexClean = index[index.find(" ")+1:].replace(" ","") #removes excess whitespace found after the units.
            answer = [word, indexClean[indexClean.find(" ")+1:]] #stores the variable and its units in a 2 index array
            if answer[1] == '': #in the case of a none value for units...
                answer[1] = "L"  #...replace with 'L'
            return answer
        elif word in index:
            answer = [word, "L"] #if no units, set to 'L'
            return answer
    return "N/A"    #if keyword not found, return 'N/A'

def findData(data,searchArray,identifier):  #sets up the extractTerms and extractOutput functions.
    outputData = []
    if identifier == "terms": #if terms required...
        for index in searchArray:
            outputData.append(extractTerms(data, index)) #...call extractTerms for each index
    elif identifier == "output": #if output required...
        for index in searchArray:
            outputData.append(extractOutput(data, index)) #...call extractOutput for each index
    return outputData #return extracted data

def readprefix(string): #takes in a string and removes exponent prefixes.
    string = str(string).replace(" ","") #remove any spaces (i.e 100 k --> 100k)
    if not(string == 'N/A'): #checks if string is a numeric value
        if "G" in string:   #checks for exponent....
            num = float(string[:string.find("G")])*1000000000 #...if exponent present, take the data before it, and multiplies or divides by the respective magnitude
        elif "M" in string:                    #repeat for each exponent.
            num = float(string[:string.find("M")])*1000000
        elif "k" in string:
            num = float(string[:string.find("k")])*1000
        elif "m" in string:
            num = float(string[:string.find("m")])/1000
        elif "u" in string:
            num = float(string[:string.find("u")])/1000000
        elif "n" in string:
            num = float(string[:string.find("n")])/1000000000
        elif "p" in string:
            num = float(string[:string.find("p")])/1000000000000
        else: 
            num = float(string) #if no exponent found, convert string to float.
        return num #return the modified float.
    return string #if string = 'N/A', return unchanged.
        
def standardiseTerms(array):  #Convert Norton Current to Thevenin voltage and runs each index through readprefix()
   
    if not(array[3] == "N/A") and array[1] == "N/A": #if Norton source present and no Thevenin source
        try:
            array[1] = float(array[3])*float(array[2])  #Convert using V=IR, unless R not given...
        except:
            array[1] = float(array[3])*(1/float(array[4])) #.. in which case use V=I/G
    for index in array:
        array[array.index(index)] = readprefix(index) #convert prefixes to float.
    return array

def unitsPrefix(array):  #checks for exponent prefixes in the units in the output block.
    for index in array:
        if not(index == "N/A"): #if value is numeric...
            unit = index[1]
            if "Ohms" in unit:          #'Ohms' containes the letter 'm' which is an exponent prefix...
                unit = unit.replace("Ohms", "X")    #...so we replace 'Ohms' with an arbitary value 'X', so as not to cause a read error.
            if "G" in unit:     #for each exponent, the function checks the index.
                mod = 0.000000001   #if found, mod is set to its respective modifer.
            elif "M" in unit:
                mod = 0.000001
            elif "k" in unit:
                mod = 0.001
            elif "m" in unit:
                mod = 1000
            elif "u" in unit:
                mod = 1000000
            elif "n" in unit:
                mod = 1000000000
            elif "p" in unit:
                mod = 1000000000000
            else:
                mod = 1
            array[array.index(index)] = [index[0],index[1],mod] #expand the 2 index array, to now contain a third index, which is the exponent modifier.
    return array

def csvHeader(array): #generates the title and units line for the output file.
    #every line is up to 250 chars. 24 per variable (12 Re, 12 Im) and 10 for Freq. (Note: its 11 + a comma)
        #(I added Ap and T as potential outputs, meaning the maximum characters per line is not 274, but for auto tests it is 250)
    names = ['      Freq'] #frequency is always required and never changes so can be manually input.
    units = ['        Hz'] #same with units.
    for index in array:
        if not(index == "N/A"): #if output required
            nameReVal = "Re("+index[0]+")" #create real part header
            nameImVal = "Im("+index[0]+")" #create imaginary part header
            nameMod = 11-len(nameReVal)    #determine how much whitespace is needed to conform to the formatting requirements.
            extend = ""
            extend = extend.ljust(nameMod) #generate the whitespace
            names.append(extend+nameReVal) #add the whitespace and append to the output name array
            names.append(extend+nameImVal) #same again
            
            unitVal = index[1]             #retrive the units
            unitMod = 11-len(unitVal)      #calculate whitespace needed
            extend = ""
            extend  = extend.ljust(unitMod) #generate whitespace
            units.append(extend+unitVal)  #add whitespace and store units twice, once for real part and once for imaginary part
            units.append(extend+unitVal)
    return names,units #return completed output lines.

def getFreq(StartF,EndF,NF,logBoo): #generates a logarithmic or linear frequency sweep.
    if logBoo == True: #logarithmic sweep
        low = math.log(StartF,10) #find smallest power
        high = math.log(EndF,10)  #find largest power
        size = (high-low)/(NF-1)  #find size of power increment
        Freqs = [] #initialise output array
        while True:
            Freqs.append(StartF) #begin adding frequencies
            StartF = pow(10,(math.log(StartF,10)+size)) #generate next logarithmic frequency
            if StartF > EndF: #if current frequncy bigger than end frequency, then the function ends.
                return Freqs
    else:   #linear sweep
        size = (EndF-StartF)/(NF-1) #linear is a little simpler. This line create linear increments
        Freqs = [] #intialise output array
        while True:
            Freqs.append(StartF) #begin adding frequencies
            StartF = StartF+size #generate next linear frequency
            if StartF > EndF: #if current frequency bigger than end frequency then the functione ends.
                return Freqs

def sortResistors(data,freq): #calculate impedances and orders components
    i=1 #initialise counter and output arrays
    Zarray = []
    Yarray = []
    while True: #loop until endproccess remains True (i.e no more components found)
        holdValue = 0 #reset holdvalue to 0
        endproccess = True #reset endproccess to True
        for index in data: #begin looping through components
            n1 = index.find("n1=") #find node 1 value
            n2 = index.find("n2=") #find node 2 value
            G = index.find("G=") #find either G,R,C or L value if given.
            R = index.find("R=")
            C = index.find("C=")
            L = index.find("L=")
            if (int(index[n1+3:index.find(" ",n1)]) == i) and (int(index[n2+3]) == 0): #if current component is a shunt on current node
                endproccess = False #keep searching for resistors afterwards
                Zarray.append(0) #Z value for shunts is 0, we now calculate the Y value
                if "G=" in index:                               #All values are run through readprefix()
                    calcG = complex(readprefix(index[G+2:]),0) #for G we simply convert to complex. 
                    Yarray.append(calcG) #then append to Yarray
                elif "C=" in index:
                    calcC = complex(0,-1/(2*math.pi*float(readprefix(index[C+2:]))*freq)) #for C, addmitance has to be calculated, then it is converted to a complex
                    Yarray.append(1/calcC) #then append to Yarray
                elif "L=" in index:
                    calcL = complex(0,(2*math.pi*float(readprefix(index[L+2:]))*freq)) #L addmitance needs to be calculated, again turned into a complex
                    Yarray.append(1/calcL) #then append to Yarray
                else:
                    calcR = complex(float(readprefix(index[R+2:])),0) #R converted to G using R = 1/G
                    Yarray.append(1/calcR) #then append to Yarray
            elif (int(index[n1+3:index.find(" ",n1)]) == i):  #if current component is in series on current node (can only be 1)
                endproccess = False #keep searching for resistors afterwards
                if "R=" in index:                               #All values are run through readprefix()
                    holdValue = complex(float(readprefix(index[R+2:])),0) #R is stored straight away.
                elif "C=" in index:
                    holdValue = complex(0,-1/(2*math.pi*float(readprefix(index[C+2:]))*freq)) #C impedance calculated and stored in holdValue
                elif "L=" in index:
                    holdValue = complex(0,(2*math.pi*float(readprefix(index[L+2:]))*freq)) #L impedance calculated and stored
                else:
                    holdValue = complex(1/(float(readprefix(index[G+2:]))),0) # G= 1/R
        Zarray.append(holdValue) #if no more components on current node found, the series component is appended to Z and Y arrays
        Yarray.append(0) #Y value for series component is 0.
        i=i+1 #increment node counter by 1
        if endproccess == True: #only true if no resistors found on current node (i.e end of circuit)
            return Zarray, Yarray #return both arrays.

def simulator(Zarray, Yarray): #calculates transform matrix
    Q = 1 #set up initial identity matrix representation
    Zn = 0
    Yn = 0
    P = 1
    for Zindex,Yindex in zip(Zarray,Yarray): #sychornously use both Zarray and Yarray
        Qnew = Q + (Yindex*Zn) #have to use new variables, so as to not mess up the maths (i.e. Q would be reassgined before its used to calcule Zn)
        Znnew = Zn + (Zindex*Q)
        Ynnew = Yn + (Yindex*P)
        Pnew = P + (Zindex*Yn)
        Q = Qnew    #update ABCD values
        Zn = Znnew
        Yn = Ynnew
        P = Pnew
        #and repeat until all components simulated
    print(Q, Zn, Yn, P)
    return Q,Zn,Yn,P #return ABCD

def calculate(A,B,C,D,array): #calculates all possible output values
    Zin = (A*array[0]+B)/(C*array[0]+D)
    try:                        #If Rs given
        Zout = (D*array[2]+B)/(C*array[2]+A)
    except:                     #If Gs given
        Zout = (D*(1/array[4])+B)/(C*(1/array[4])+A)
    Av = 1/(A+B*(1/array[0]))
    Ai = 1/(C*array[0]+D)
    try:                         #If Rs given
        Vin = array[1]*(Zin/(Zin+array[2]))
    except:                     #If Gs given
        Vin = array[1]*(Zin/(Zin+(1/array[4])))
    Vout = Av*Vin
    Iin = Vin/Zin
    Iout = Vout/array[0]
    Pin = Vin*Iin.conjugate()
    Pout = Vout*Iout.conjugate()
    Ap = Av*Ai.conjugate()
    try:                         #If Rs given
        T=2/((A*array[0])+B+(C*array[0]*array[2])+D*array[2])
    except:                     #If Gs given
        T=2/((A*array[0])+B+(C*array[0]*(1/array[4]))+D*(1/array[4]))
    return [Vin,Vout,Iin,Iout,Pin,Zout,Pout,Zin,Av,Ai,Ap,T] #return values in a single array.

def splitComplex(array,arrayData,freq): #splits values into real and imaginary parts
    splitArray = [freq] #initialise output array with current frequency
    for index in array:
        if not(index == "N/A"): #check if value is needed in output
            num = arrayData[array.index(index)] #store value
            num = num*index[2]                  #multiple by exponent modifier from unitsPrefix()
            splitArray.append(num.real) #append real part
            splitArray.append(num.imag) #append imaginary part
    return splitArray

def prefix(value,choice):  #adds exponent prefixes to output values
    value = ("{:.3e}".format(value)) #round every value to standard form with 3.d.p
    if choice == "prefix": #if prefix selected
        factor = value.find("e")
        if value[factor+1] == '+' and int(value[factor+2:]) >= 9:  #checks order of magnitude for each possible exponent
            mod = int(value[factor+2:])-9                           #if within a certain bound
            value = str(round(float(value[:factor])*pow(10,mod),3)) + "G" #multiply by the correct value and add the expoent prefix
        elif value[factor+1] == '+' and int(value[factor+2:]) >= 6: #repeat for all exponents
            mod = int(value[factor+2:])-6                   
            value = str(round(float(value[:factor])*pow(10,mod),3)) + "M"
        elif value[factor+1] == '+' and int(value[factor+2:]) >= 3:
            mod = int(value[factor+2:])-3
            value = str(round(float(value[:factor])*pow(10,mod),3)) + "k"
        elif value[factor+1] == '+' and int(value[factor+2:]) >= 0:
            mod = int(value[factor+2:])
            value = str(round(float(value[:factor])*pow(10,mod),3))
        elif value[factor+1] == '-' and int(value[factor+2:]) >= 10:
            mod = abs(12-int(value[factor+2:]))
            value = str(round(float(value[:factor])*pow(10,mod),3)) + "p"
        elif value[factor+1] == '-' and int(value[factor+2:]) >= 7:
            mod = abs(9-int(value[factor+2:]))
            value = str(round(float(value[:factor])*pow(10,mod),3)) + "n"
        elif value[factor+1] == '-' and int(value[factor+2:]) >= 4:
            mod = abs(6-int(value[factor+2:]))
            value = str(round(float(value[:factor])*pow(10,mod),3)) + "u"
        elif value[factor+1] == '-' and int(value[factor+2:]) >= 1:
            mod = abs(3-int(value[factor+2:]))
            value = str(round(float(value[:factor])*pow(10,mod),3)) + "m"
        return value
    else:       #if SF selected...
        return value #...return the value without exponent conversion

def whiteSpace(array):  #formats output values with whitespace
    for index in array:
        if array.index(index) == 0:
            numSpace = 9 -len(index) #Frequency column has 1 less character than the rest.
        else:
            numSpace = 10 -len(index)
        whiteSpace1 = ''
        while numSpace > 0:
            whiteSpace1 = whiteSpace1 + " " #generate whitespace
            numSpace = numSpace-1
        if index[0] == '-':
            whiteSpace2 = "" #if there is a minus sign, less whitespace is given
        else:
            whiteSpace2 = " " #if positive, an extra space is given
        if array.index(index) == 0:
            array[array.index(index)] = whiteSpace1 + whiteSpace2 + str(index) #all whitespae is added
        else:
            array[array.index(index)] = whiteSpace1 + whiteSpace2 + str(index) #same again but for non frequency values.
    return array #return formatted output values

###### Main ######

try:

    prefixChoice = unitType() #determines output format

    inputFile = open(sys.argv[1], "rt") #opens doc from cmd
    inputLines = inputFile.readlines() #extracts lines
    cleanData = clean(inputLines) #standardises data

    circuitData=sortInput(cleanData,'CIRCUIT') #extract the three main data blocks
    termsData=sortInput(cleanData,'TERMS')
    outputData=sortInput(cleanData,'OUTPUT')


    logBoo = logOrLin(termsData) #determines which frequency sweep to use

    termsToFind = ['RL','VT','RS','IN','GS','Fstart','Fend','Nfreqs']
    termsValues = findData(termsData,termsToFind,"terms") #extracts numeric values for each term
    termsValues = standardiseTerms(termsValues) #read prefixes and converts Norton current to Thevenin voltage if required


    outputToFind = ['Vin','Vout','Iin','Iout','Pin','Zout','Pout','Zin','Av','Ai','Ap','T']
    outputValues = findData(outputData, outputToFind,"output") #extracts units from output values
    outputValues  = unitsPrefix(outputValues) #searches for prefixes and assigns modifer.

    names,units = csvHeader(outputValues) #generates first two lines of output doc ready for writing
    Frequencies = getFreq(termsValues[5],termsValues[6],termsValues[7],logBoo) #generates frequency sweep, taking logBoo to determine which type

    with open(sys.argv[2], 'w', newline='') as file: #opens output file for writing
        writer = csv.writer(file)
        writer.writerow(names) #writes first two lines to output file
        writer.writerow(units)
        
        for freq in Frequencies:                        #for each frequency...
            Z,Y = sortResistors(circuitData,freq)       #calculate impedances and order components
            An,Bn,Cn,Dn = simulator(Z,Y)                #generate transform matrix
            results = calculate(An,Bn,Cn,Dn,termsValues)    #calculate output values
            finalResults = splitComplex(outputValues,results,freq) #determine which outvalues are required
            for resultt in finalResults:
                finalResults[finalResults.index(resultt)] = prefix(resultt,prefixChoice) #add exponents to results if required
            finalResults = whiteSpace(finalResults) #add whitespace to format results
            writer.writerow(finalResults) #write results to output file.
        #if no more frequencies need simulating, the program ends.

except:     #in the event of an error, the program returns a blank output file, as specified by the brief.
    with open(sys.argv[2], 'w', newline='') as newfile:
        writer = csv.writer(newfile)


