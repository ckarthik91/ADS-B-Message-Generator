#!/usr/bin/env python

#Author : Karthik C

#//////////////////////////////////////////////////
#Code developed to genetare ADS-B message//////////
#//////////////////////////////////////////////////
import math 
import sys
import gmpy2 #multi precision module
from gmpy2 import mpz,mpq,mpfr,mpc
import numpy
from numpy import *

gmpy2.get_context().precision = 32#set gmpy precision to 128 bits

Dlat0 = 6 # Even message
Dlat1 = mpfr(6.101694915254237288135593220339) # Odd message
pi = mpfr(gmpy2.const_pi(),56)
evenMsg = mpz(1)
oddMsg = mpz(1)

#-----------For testing -------------------------
t_aid = ab7444
t_alt = 40000
t_lat = mpfr(49.925827)
t_lon = mpfr(89.193542)

#-------------------------------------------------


#to concatinate hex digits
def concat_hex(a,b):
    size_of_b = 0
    #print "a="+hex(a)
    #print "b="+hex(b)

    #get size of b in bits
    while((b >> size_of_b ) > 0):
        size_of_b += 1
        
    #print "size_of_b= "+str(size_of_b)

    #align answer to nearest 4 bits (hex digit)
    if ((size_of_b % 4) !=0):
        size_of_b += 4-(size_of_b % 4)

    return (a << size_of_b)|b 
#------------------------------------------------------------------------
#RTCA Special Committee 186 , Working Group 3 
def modulus(val,modval):
    #print "init val " + str(val)
    #print "init modval " + str(modval)
    if (val < 0 ):  #checking for negative values
        val = val + 360
        #print "if val " +str(val)
        result = int(math.floor((val/ modval)))
        #print "if Result : " + str(result) 
        #print "if  Modulus " + str(val - (result * modval))
        return val - (result * modval)
    
    else:
        result = int(math.floor(val/modval))
        #print "Else result " + str(result)
        #print "Else moduls " + str(val - (result * modval)) 
        return val - (result * modval)
        
#------------------------------------------------------------------------------

def genNlatValue(latitude):
    j = 2
    buff = numpy.empty(63,dtype = object)
    #print buff
    buff[:] = [mpfr('1.0') for _ in range(63)] #double array of length 63
    
    for i in xrange(2,60):
        Nlat = mpfr(0.0)
        sqrtHold = mpfr(0.0)
        Nlat = (180/pi)*(gmpy2.acos(math.sqrt((1-(math.cos(pi/30)))/(1-(gmpy2.cos((2*pi)/i))))));
        buff[j] = Nlat
        j = j+1
    
    for k in xrange(59,1,-1):
        if ( latitude > buff[k]):
            pass
        else:
            if (k<2):
                return 1
            return k
            
#--------------------------------------------------------------------------------
    

def calculateLatBitsOdd( lat,longitude):
    #calculate YZ which is to be put into message
    modHold = mpfr(modulus(lat,Dlat1))
    modHold = (modHold/Dlat1)
    modHold = modHold * math.pow(2,17)
    modHold = modHold + 0.5
    YZ = mpfr(math.floor(modHold))
    
    #calculate Rlatitude for Airborne
    Rlat1 = mpfr(1.0)
    floorHold = 1
    Rlat1 = mpfr(YZ) / math.pow(2,17)
    floorHold = (lat/mpfr(Dlat1))
    floorHold = math.floor(floorHold)
    Rlat1 = Rlat1 + mpfr(floorHold)
    Rlat1 = Rlat1 * Dlat1
    NlLat = 1
    NlLat = genNlatValue(Rlat1)
    
    #calculate Dlongitude
    Dlon1 = mpfr(1.0)
    if((NlLat-1) > 0):
        Dlon1 = (360/(mpfr(NlLat)-1));
    else:
        Dlon1 = 360
        
    #calculate XZ the decimal representation of longitude
    XZ = mpfr(1.0)
    modHold = 0
    modHold = modulus(longitude, Dlon1)
    modHold = (modHold/Dlon1)
    modHold = modHold * math.pow(2,17)
    modHold = modHold + 0.5
    XZ = math.floor(modHold)
    
    #Ensure this fits into our 17 bit space
    YZ1 = mpz(modulus(YZ,math.pow(2,17)))
    XZ1 = mpz(modulus(XZ,math.pow(2,17)))
    LatLon = mpz(1)
    YZ1 = YZ1 << 17
    LatLon = (YZ1 | XZ1)
    LatLon = (LatLon | 17179869184)
    
    #to be decided
    print "odd message : " + hex(LatLon)
    oddMsg = Latlon
    
    return 0

#---------------------------------------------------------------------------------


def calculateLatBitsEven():
    #lat = mpfr(raw_input(" Please enter required Latitude > "))
    #longitude =  mpfr(raw_input(" Please enter required longitude > "))
    lat = t_lat
    longitude = t_lon 
    
    calculateLatBitsOdd(lat,longitude)
    
    #Calculate YZ which will be what is put into our message
    modHold = mpfr(modulus(lat,Dlat0))
    modHold = (modHold/Dlat0)
    modHold = modHold * math.pow(2,17)
    modHold = modHold + 0.5
    YZ = mpfr(math.floor(modHold))
    
    #Calculate Rlatiude for airborne
    Rlat0 = mpfr(1.0)
    floorHold = 1
    Rlat0 = mpfr(YZ) / math.pow(2,17)
    floorHold = (lat/mpfr(Dlat0))
    floorHold = math.floor(floorHold)
    Rlat0 = Rlat0 + mpfr(floorHold)
    Rlat0 = Rlat0 * mpfr(Dlat0)
    NlLat = 1
    NlLat = genNlatValue(Rlat0)
    
    #Calcule DLongitude
    Dlon0 = mpfr(1.0)
    if((NlLat-1) > 0):
        Dlon0 = (360/(mpfr(NlLat)));
    else:
        Dlon0 = 360
        
    #calculate XZ the decimal representation of longitude
    XZ = mpfr(1.0)
    modHold = 0
    modHold = modulus(longitude, Dlon0)
    modHold = (modHold/Dlon0)
    modHold = modHold * math.pow(2,17)
    modHold = modHold + 0.5
    XZ = math.floor(modHold)
    #print mpfr(XZ)
    
    #Ensure this fits into our 17 bit space
    YZ1 = mpz(modulus(YZ,math.pow(2,17)))
    XZ1 = mpz(modulus(XZ,math.pow(2,17)))
    #print YZ1
    LatLon = mpz(1)
    YZ1 = YZ1 << 17
    LatLon = (YZ1 | XZ1)
    #LatLon = (LatLon | 17179869184)
    
    #to be decided
    print "even message: " + hex(LatLon)
    
    return 0
#--------------------------------------------------------------------------------



#---------------------------------------------------------------------------------

#get input Aircraft ID and append with DF type and CA

def generateAirID():
    airID = 0x8D # bit 1 to 5 DF format 17 and 6 to 8 CA type 3
    
    #ID = raw_input("Enter Aircraft ID (ex : AF4511) > ")
    ID = t_aid
    IDhex = int(ID,16)#convert input to hex
    
    print 'Aircrat ID: '+hex(IDhex)
    airID = airID << 24 # shift 24 bits to make space for CAO aircraft address
    airID = (airID | IDhex )
    return airID
    
#-----------------------------------------------------------------------------

#funcition to insert altitude with typecode 

def calAltitude():
    #altitude = int(raw_input("Enter the altitude (0 ft to 50000 ft) > "))
    altitude = t_alt
    altitude = (altitude + 1000 ) /25 # code in increments of 25ft
    hold = (altitude & 0x00F)
    altitude = (altitude & 0xFF0) << 1;
    altitude = (altitude | hold) #concatenate to the entire message
    altitude = (altitude | 0x010)
    altitude = (altitude | 0x58000) #Takes on the TC (0x58) field
        
    return altitude 
    
#--------------------------------------------------------------------------

    

    
#----------------------------------------------------------------------------
    
if __name__ == '__main__':
    print "This program generates ADS-B message "
    airid =  hex(generateAirID())
    alt = hex(calAltitude())
    calculateLatBitsEven()
    
    print "ADS-B message : " + str(airid) +" "+ str(alt)
    
    
    
    






