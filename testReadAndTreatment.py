#!/usr/bin/env python3
# -*- coding: utf-8 -*-

from readFileGyroscope import readFile
from dataTreatmentsGyroscope import correction, distancesToPoints, pointsToCircles, calculVolum


file = 'VaseTest.CSV'

listDH = readFile(file)

a =  correction(listDH)

a = distancesToPoints(a[0], a[1])

a = pointsToCircles(a[0], a[1])

a = calculVolum(a[0])
print("\nVolume =", a, "mm3")

print("\nVOLUME =", a / (10**6), "L")

