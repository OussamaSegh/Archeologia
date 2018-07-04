#! /usr/bin/env python3
# -*- coding: utf-8 -*-

"""
Archeo-metrologia project
@author = Alban GOUGOUA & Oussama SEGHAIER

New treatment with the gyroscope (addition of angle)
"""

from math import pi, fabs, cos, sin, sqrt #: Import a few of constant and functions of module math.
import numpy as np #: Import all functions of module numpy.
import random #: Import module random.


INFINITY = 11 #: All value less than 11 millemeters and greater than 2000 millimeters
                    #: aren't detected by the sensor.
                    #: We must complete with other values by approximation.
#SPEED_MOTOR = 2.09 #: The angular speed of motor in radians/second.
#FREQUENCY_SENSOR = 1/(238*(10**(-3))) #: The frequency of sensor in hertz (1/second).
#NUMBER_POINTS_IN_CIRCLE = int(FREQUENCY_SENSOR * ((2*pi) / SPEED_MOTOR)) #: The number of measured points when the screw doing one round on itself.

################################################### BEGINNING #################################################

def correction(listDistAngHeightsNum) :
    """
    @author = Alban GOUGOUA
    
    This function allows of averaging the values of the listDistAngHeightsNum which are infinites to values 
        near of their neighbours not infinite.
    It returns a couple composed of one list of lists (corlistDistAngHeightsNum) and another list (indexInfDistList1).
    
    :param listDistAngHeightsNum : The list of distances, associated angles and heights, and numbers of round from a file .csv.
        This list is a list of lists [distance, angle, height, number].
        Distance and height are in millemeter.
        Angle is in radian.
        Number is an integer.
    """
    
    #: Sort of listDistAngHeightsNum in order ascending of the number.
    listDistAngHeightsNum = sorted(listDistAngHeightsNum, key=lambda d : d[3])
    
    #: Share the listDistAngHeightsNum in three lists : first for the distances, second for the angles, third for the associated heights and last for the associated numbers of round.
    listDist = []
    listAng = []
    listHeight = []
    listNum = []
    for l in listDistAngHeightsNum :
        listDist.append(l[0])
        listAng.append(l[1])
        listHeight.append(l[2])
        listNum.append(l[3])
    
    #: Put the index of infinite distances in a list.
    indexInfDistList = []
    for i in range(len(listDist)) :
        if listDist[i] < INFINITY :
            indexInfDistList.append(i)
    
    #: Copy of the list of indexes of infinite distances in another list.
    indexInfDistList1 = indexInfDistList[:]
    
    #: First correction in putting all infinte distances near INFINITY.
    for i in indexInfDistList :
        listDist[i] = random.randint(10, 15) #: This value is abritrary.
    
    #: Give an averaged value to infinite values in counting the number of successive indexes.
    sucIndex = [] #: The list of successive indexes.
    nonsucIndex = [] #: The list of nonsuccessive indexes.
    i = 0 #: First index of the list. His value doesn't never change during the processing.
    while i < len(indexInfDistList) :
        #: Case 1 : the lenght of the list of indexes is greater than 1.
        if len(indexInfDistList) > 1 :
            #: We check if the list of successive indexes is empty.
            if len(sucIndex) == 0 :
                if indexInfDistList[i] == (indexInfDistList[i+1] - 1) : #: We check if the two first elements are successives.
                    #: We add the two first elements to the list of successive indexes.
                    sucIndex.append(indexInfDistList[i]) 
                    sucIndex.append(indexInfDistList[i+1])
                    #: We delete the two first elements to the list of indexes. When we remove the first element of the list, the second element locates at the index 0 now.
                    indexInfDistList.remove(indexInfDistList[i])
                    indexInfDistList.remove(indexInfDistList[i])
                else : #: If the two first elements aren't successives, do below.
                    nonsucIndex.append(indexInfDistList[i]) #: We add the first element in the list of nonsuccessive indexes.
                    indexInfDistList.remove(indexInfDistList[i]) #: And we remove it of the list of indexes.
            #: We check if the list of successive indexes isn't empty, that is to say, is greater than or equal to 2.
            elif len(sucIndex) >= 2 :
                sucIndex.sort() #: We sort this list in ascending order.
                #: We check if the first element of the list of indexes is equal to last element of the list of successive indexes.
                if indexInfDistList[i] == (sucIndex[-1] + 1) :
                    sucIndex.append(indexInfDistList[i]) #: We add it in the list of successive indexes.
                    indexInfDistList.remove(indexInfDistList[i]) #: And we remove it of the list of indexes.
                #: In the reverse case, do below.
                else :
                    minDist = listDist[sucIndex[0] - 1] #: We take the value of element locating before the first element of the list of successive indexes in the list of distances.
                    maxDist = listDist[sucIndex[-1] + 1] #: We take the value of element locating after the last element of the list of successive indexes in the list of distances.
                    newMin = min(minDist, maxDist) #: We take the smallest value.
                    #: We correct the value of infinite distances in finite values in doing a mean.
                    for s in sucIndex :
                        listDist[s] = newMin + (fabs(maxDist - minDist)/(len(sucIndex) + 1))
                        newMin = listDist[s]
                    sucIndex = [] #: We empty this list.
        #: Case 2 : the lenght of the list of indexes is equal to 1.
        else :
            #: We check if the list of successive indexes isn't empty, that is to say, is greater than or equal to 2.
            if len(sucIndex) >= 2 :
                sucIndex.sort() #: We sort this list in ascending order.
                #: We check if the first element of the list of indexes is equal to last element of the list of successive indexes.
                if indexInfDistList[i] == (sucIndex[-1] + 1) :
                    sucIndex.append(indexInfDistList[i]) #: We add it in the list of successive indexes.
                    indexInfDistList.remove(indexInfDistList[i]) #: And we remove it of the list of indexes.
                #: In the reverse case, do below.
                else :
                    minDist = listDist[sucIndex[0] - 1] #: We take the value of element locating before the first element of the list of successive indexes in the list of distances.
                    maxDist = listDist[sucIndex[-1] + 1] #: We take the value of element locating after the last element of the list of successive indexes in the list of distances.
                    newMin = min(minDist, maxDist) #: We take the smallest value.
                    #: We correct the value of infinite distances in finite values in doing a mean.
                    for s in sucIndex :
                        listDist[s] = newMin + (fabs(maxDist - minDist)/(len(sucIndex) + 1))
                        newMin = listDist[s]
                    sucIndex = [] #: We empty this list.
            #: We check if the list of successive indexes is empty.
            else :
                nonsucIndex.append(indexInfDistList[i]) #: We add the one element in the list of nonsuccessive indexes.
                indexInfDistList.remove(indexInfDistList[i]) #: And we remove it of the list of indexes.
    
    #: We correct the value of last infinite distances which are in the list of successive indexes.
    if len(sucIndex) > 0 :
        minDist = listDist[sucIndex[0] - 1]
        if (sucIndex[-1] + 1) > 0 :
             maxDist = listDist[0]
        else :
            maxDist = listDist[sucIndex[-1] + 1]
        newMin = min(minDist, maxDist)
        for s in sucIndex :
            listDist[s] = newMin + (fabs(maxDist - minDist)/(len(sucIndex) + 1))
            newMin = listDist[s]
    
    #: We correct the value of infinite distances which are in the list of nonsuccessive indexes.
    for i in nonsucIndex :
        if i+1 > len(listDist) :
            newMin = min(listDist[i-1], listDist[0])
            listDist[i] = newMin + (fabs(listDist[i-1] - listDist[0]) / 2)
        else :
            newMin = min(listDist[i-1], listDist[i+1])
            listDist[i] = newMin + (fabs(listDist[i-1] - listDist[i+1]) / 2)
    
    #: Put the list of distances, the list of angles, the list of associated heights and the list of numbers in a list.
    corlistDistAngHeightsNum = [] #: The corrected list of distances, angles, heights and numbers.
    for d in range(len(listDist)) :
        corlistDistAngHeightsNum.append([listDist[d], listAng[d], listHeight[d], listNum[d]])
    
    return (corlistDistAngHeightsNum, indexInfDistList1)

    ################################################################################################################################################################

def distancesToPoints(corlistDistAngHeightsNum, indexInfDistList) :
    """
    @author = Alban GOUGOUA
    
    This function allows to figure out coordinates of points from a distance.
    It does it for a list of distances and returns a list of coordinates (x,y,h,n) where x and y are calculated with distances 
        and angles ; h corresponds to the height associated to each distance ; n corresponds to the number of round associated to each distance.
    It returns also the list of indexes of infinite distances (indexInfDistList).
    
    :param corlistDistAngHeightsNum : The list of distances, associated angles and heights, and numbers of round from a file .csv which was corrected.
        This list is a list of lists [distance, angle, height, number].
        Distance and height are in millemeter.
        Angle is in radian.
        Number is an integer.

    :param indexInfDistList : The list of indexes of infinite distances.
        It is a list of integers.
    """
    
    #: Initialization of differents variables to 0.
    listCoordAndHeightsNum = []
    x = 0
    y = 0
    h = 0
    n = 0
    
    #: We figure out the coordinates (x, y, h) of each distance of the list of distances, angles and heights.
    for d in corlistDistAngHeightsNum :
        #: We figure out coordinates (x, y) from distance and associated angle. h receives the value of associated height.
        x = d[0] * cos(d[1])
        y = d[0] * sin(d[1])
        h = d[2]
        n = d[3]
        
        listCoordAndHeightsNum.append((x, y, h, n)) #: We add each calculated data in the list.
    
    #: We return the list of coordinates and associated heights and the list of indexes of infinite distances.
    return (listCoordAndHeightsNum, indexInfDistList)

    ################################################################################################################################################################
       
def pointsToCircles(listCoord, indexInfDistList) :
    """
    @author = Alban GOUGOUA & Oussama SEGHAIER
    
    This methods finds the least squares circle fitting a set of 2D points (x,y).
    It returns the list of the best circles using the least square method and the associated height ;
        it returns also, the list of gaps.
        
    INPUTS : 
      @ListCoord = [(x,y,h) .... ] : (x,y) the cartesien coordinates of a point , h : height of the point.

    source of methods : http://scipy-cookbook.readthedocs.io/items/Least_Squares_Circle.html
        We choose the method 2b : "leastsq with jacobian".
    """
    
        ################################################## Beginning of this method ##################################################
    # == METHOD 2b ==
    # Advanced usage, with jacobian
    
    def pointsToRadius(X, Y) :
        """
        @author = Alban GOUGOUA

        This function figures out from a numpy array list of coordinates the least squares circle.
        It returns the radius of this circle.

        :param X : a numpy array list including NUMBER_POINTS_IN_CIRCLE coordinates in x.

        :param Y : a numpy array list including NUMBER_POINTS_IN_CIRCLE coordinates in y.
        """

        from scipy import optimize as op


        def calc_R(xc, yc):
            """ calculate the distance of each 2D points from the center c=(xc, yc) """
            return np.sqrt((X-xc)**2 + (Y-yc)**2)
        
        def f_2b(c):
            """ calculate the algebraic distance between the 2D points and the mean circle centered at c=(xc, yc) """
            Ri = calc_R(*c)
            return Ri - Ri.mean()
        
        def Df_2b(c):
            """ Jacobian of f_2b
            The axis corresponding to derivatives must be coherent with the col_deriv option of leastsq"""
            xc, yc     = c
            df2b_dc    = np.empty((len(c), X.size))

            Ri = calc_R(xc, yc)
            df2b_dc[ 0] = (xc - X)/Ri                   # dR/dxc
            df2b_dc[ 1] = (yc - Y)/Ri                   # dR/dyc
            df2b_dc       = df2b_dc - df2b_dc.mean(axis=1)[:, np.newaxis]

            return df2b_dc


        # coordinates of the barycenter
        x_m = np.mean(X)
        y_m = np.mean(Y)
        
        #: We figure out the center of the circle.
        center_estimate = x_m, y_m
        center_2b = op.leastsq(f_2b, center_estimate, Dfun=Df_2b, col_deriv=True)
        
        #: We figure out the radius of the circle.
        xc_2b, yc_2b = center_2b[0]
        Ri_2b        = calc_R(xc_2b, yc_2b)
        R_2b         = Ri_2b.mean()

        return R_2b
    
        ################################################### End of this method #################################################
        
        ##########################################################################################        
        ##                           CoordsAndHeight(ListCoord)-->> X,Y,H                       ##   
        ##########################################################################################        
        
    def coordsAndHeight(listCoord) :
        """
        @author = Oussama SEGHAIER
        
        This function transforms a listCoord as explained above to two numpy arrays X and Y 
            constituted respectively of the x coordinates and the y coordinates.

        the funcion returns  
        @X : a numpy array of the x coordinates 
        @Y : -------------------- y -----------
        @H : a list of the h associated coordinates
        """
        
        #: We create the list of each coordinate x, y.
        X_r = []
        Y_r = []
        H = []
        
        #:We add each value of the list of coordinates and associated heights in lists above.
        for coord in listCoord :
            X_r.append(coord[0])
            Y_r.append(coord[1])
            H.append(coord[2])
        
        #: We transform these lists in numpy array lists.
        X = np.array(X_r)
        Y = np.array(Y_r)
        
        return X,Y,H
    
    
    #: We create a dictionnary of number of points in each round.
    listNumber = []
    for l in listCoord :
        listNumber.append(l[3])
    points = {}
    counter = 1
    value = 1
    while(counter != 0) :
        counter = listNumber.count(value)
        if(counter != 0) :
            points[value] = counter
        value += 1
    
    #: We figure out the radius of each circle to each NUMBER_POINTS_IN_CIRCLE points.
    pseudoListCoord = [] #: A list for storing the NUMBER_POINTS_IN_CIRCLE values of listCoord.
    listRadiAndHeights = [] #: A list which will return at the end of this function.
    for key in points.keys() :
        number = points[key]
        counter = 0 #: A counter for counting the number of points reaching in the listCoord.
        for coord in listCoord :
            counter += 1
            pseudoListCoord.append(coord) #: We add each value of listCoord.
            #: We figure out the radius when the counter is equal to NUMBER_POINTS_IN_CIRCLE.
            if(counter == number) :
                X = coordsAndHeight(pseudoListCoord)[0] #: We recover the numpy array list of coordinates x.
                Y = coordsAndHeight(pseudoListCoord)[1] #: We recover the numpy array list of coordinates y.
                H = np.mean(coordsAndHeight(pseudoListCoord)[2]) #: We recover the mean value of the list of heights.
                R = pointsToRadius(X, Y) #: We figure out the radius.
                listRadiAndHeights.append((R, H)) #: We add the value of radius and the associated height in the list.
                pseudoListCoord = [] #: We empty this list.
                #: We remove values from listCoord already used for the calculation of radius.
                while(counter > 0) :
                    for i in listCoord :
                        listCoord.remove(i)
                        break
                    counter -= 1
                break

    #: We sort the listRadiAndHeights in order descending of heights.
    listRadiAndHeights = sorted(listRadiAndHeights, key=lambda d : d[1], reverse=True)  
    
    #: We complete this list with the values of missed radius and heights. In fact, the device doesn't reach the bottom of the vase.
        #: We add the three last values of the listRadiAndHeights in sorting of the order ascending of the heights.
    corListRadiAndHeights = []
    for l in range(10) :
        corListRadiAndHeights.append(sorted(listRadiAndHeights, key=lambda d : d[1])[l])
    minHeight = corListRadiAndHeights[0][1] #: We recover the smallest value of heights of the corListRadiAndHeights.
    listRadius = [] #: The list of radius extracted of the corListRadiAndHeights.
    while(minHeight > 0) :
        #: We put all values of radius from corListRadiAndHeights in listRadius.
        for c in corListRadiAndHeights :
            listRadius.append(c[0])
        stdRadius = np.std(listRadius) #: We figure out the standard deviation of the list of radius.
        if((sum([listRadius[i] - listRadius[i+1] for i in range(len(listRadius)-1)])) < 0) : #: The form of the vase is convex.
            corListRadiAndHeights.append((listRadius[0] - stdRadius, minHeight - 1))
        else : #: The form of the vase is concave.
            corListRadiAndHeights.append((listRadius[0] + stdRadius, minHeight - 1))
        listRadiAndHeights.append(corListRadiAndHeights[-1])
        corListRadiAndHeights.remove(corListRadiAndHeights[-2]) #: We just remove the first element of the list.
        corListRadiAndHeights = sorted(corListRadiAndHeights, key=lambda d : d[1])
        listRadius = []
        minHeight -= 5


    return (listRadiAndHeights, indexInfDistList)

    ################################################################################################################################################################

def calculVolum(rayList) :
    """
    @author = Alban GOUGOUA & Oussama SEGHAIER

    This function figures out the volume in millimeter cube from a list of rays and heights.
    
    rayList = [(R,h).....] list of radius and associated heights that calculates the vase volume
    """
    
    def volume(r1, r2, h) :
        """calculates the volume of a cone truncated with parameters R = radius and H = height. 
        source : https://calculis.net/volume/cone-tronque
        """
        return (h*pi/3)*((r1**2)+(r2**2)+r1*r2)
    
    volum = 0
    #: We figure out the volume by sum of each cone's volume.
    for i in range(len(rayList)) :
        if(i != len(rayList) - 1) :
            h = rayList[i][1] - rayList[i+1][1]
            r1 = rayList[i][0]
            r2 = rayList[i+1][0]
            volum += volume(r1, r2, h)
    
    return round(volum, 3)

################################################### END #################################################
