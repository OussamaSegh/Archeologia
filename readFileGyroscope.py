#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
Archeo-metrologia project
@author = Alban GOUGOUA & Oussama SEGHAIER

New treatment with the gyroscope (addition of angle)
"""

import csv
from math import floor
#from dataTreatmentsGyroscope import NUMBER_POINTS_IN_CIRCLE, INFINITY

#INITIAL_HEIGHT = 2 #: The device of measurement can go down in one vase until a height of 2 millimeter maximum to avoid a contact.

def readFile(file) :
    """
    @author = Alban GOUGOUA
    
    This function allows of reading a file .csv including the values of distances, associated angles and heights 
        measured by the device, and numbers of round. 
    It returns a list of lists of distances, associated angles and heights, and numbers of round.
    
    :param file : the file .csv including data.
    """
    
    listDistAngHeightsNum = [] #: The list of distances, associated angles and heights, and numbers of round.
    table = [] #: A table which will include values of file .csv.
    #correction = [] #: A list of values of heights to correcting.
    
    #: Opening of the file .csv and adding of values from file in table. Then we close this file .csv.
    with open(file, newline='') as csvfile :
        csvreader = csv.reader(csvfile) 
        for row in csvreader :
            table.append(row)
    csvfile.close()
    
    #: We take values from table in listDistAngHeightsNum in the form of floating numbers.
    for i in range(0, len(table)) :
        listDistAngHeightsNum.append([float(table[i][0]), float(table[i][1]), float(table[i][2]), float(table[i][4])])
    
    #: Sort of listDistAngHeightsNum in order ascending of the number.
    listDistAngHeightsNum = sorted(listDistAngHeightsNum, key=lambda d : d[3])
    
    """#: Correction of values less than INFINITY
    for l in listDistAngHeightsNum :
        if(l[2] < INFINITY) :
            correction.append((l[2], l[3])) #: We add each tuple of (height, number) in this list for height less than 10 mm.
    correction = sorted(correction, key=lambda d : d[1]) #: We sort this list in ascending order of number.
    counter = 0 #: We initialize a counter.
    initial_height = INITIAL_HEIGHT - 1 #: We initialize a variable to 1. This variable is incrementing of 1.
    while(counter <= int(floor(NUMBER_POINTS_IN_CIRCLE / 2))) : #: We compare counter to the whole part of the half of NUMBER_POINTS_IN_CIRCLE.
        if(len(correction) == 0) :
            counter = NUMBER_POINTS_IN_CIRCLE #: If the list correction is empty, we leave the loop while.
        else :
            value = correction[0] #: We use the first tuple of correction.
            for l in listDistAngHeightsNum :
                if(l[3] == value[1]) :
                    l[2] = initial_height + 1 #: We change the value of false height of listDistAngHeightsNum by initial_height + 1. We begin to 2 mm.
                    correction.remove(value) #: We remove the first tuple of correction.
                    counter += 1 #: We increment counter.
            if(counter == int(floor(NUMBER_POINTS_IN_CIRCLE / 2))) :
                const = NUMBER_POINTS_IN_CIRCLE % int(floor(NUMBER_POINTS_IN_CIRCLE / 2)) #: We figure out the rest of euclidian division of NUMBER_POINTS_IN_CIRCLE by the whole part of its half.
                if(const == 0) : #: If the rest is 0, we reinitialize counter to 0 and increment initial_height to 1.
                    counter = 0
                    initial_height += 1
                elif(const == 1) : #: If the rest is 1, we reinitialize counter to 0, increment initial_height to 1 and figure out the number of remaining points to make one round of measurement.
                    counter = 0
                    initial_height += 1
                    for i in range(NUMBER_POINTS_IN_CIRCLE - int(floor(NUMBER_POINTS_IN_CIRCLE / 2))) :
                        if(len(correction) == 0) :
                            break #: If the correction is empty, we leave the loop for.
                        else :
                            value = correction[0]
                            for l in listDistAngHeightsNum :
                                if(l[3] == value[1]) :
                                    l[2] = initial_height + 1
                                    correction.remove(value)
                    initial_height += 1

    #: We remove the column of number.
    listDistAngHeights1 = listDistAngHeightsNum[:]
    listDistAngHeightsNum = []
    for l in listDistAngHeights1 :
        listDistAngHeightsNum.append([l[0], l[1], l[2]])"""
    
    return listDistAngHeightsNum

