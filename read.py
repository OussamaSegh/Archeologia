from readFileGyroscope import readFile
import sys 
from dataTreatmentsGyroscope import correction, distancesToPoints, pointsToCircles
file1 = sys.argv[1]

listDH = readFile(file1)

a =  correction(listDH)
a = distancesToPoints(a[0], a[1])
a = pointsToCircles(a[0], a[1])
for i in range(len(a[0])):
    print(a[0][i][0],',',a[0][i][1],'\n')
